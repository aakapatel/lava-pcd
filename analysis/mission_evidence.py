"""Mission-evidence extraction from the two day-2 long-flight bags (ROS 2, sqlite3).

Runs INSIDE the stage_cuda container (ROS 2 Humble). Reads, with a storage
filter, only the light topics needed to resolve the manuscript claims:

  - /shafterx2/dlio/odom_node/pose (geometry_msgs/PoseStamped): path length
    (0.1 m downsample, >2 m jump rejection), max Euclidean distance from the
    start pose ("forward reach"), first-to-last-pose closure error, z range,
    duration.
  - /shafterx2/bt_visualization/transitions (bt_msgs/BtTransition, custom, not
    installed): decoded manually from CDR. Layout verified against
    exploration_ws/src/bt_logger/bt_msgs/msg/BtTransition.msg:
      string visualization_name; builtin_interfaces/Time timestamp;
      string node_name; uint8 status; uint8 prev_status
      (SUCCESS=0 RUNNING=1 FAILURE=2 IDLE=3 SKIPPED=4)
  - /shafterx2/execution_complete (std_msgs/String): count + payloads.
  - /shafterx2/homing_path (nav_msgs/Path): timestamps (first = RTH/homing
    activation evidence) + waypoint counts.
  - /shafterx2/planner_status, /shafterx2/planner_mode, /shafterx2/takeoff
    (std_msgs/String): payloads.
  - /rosout (rcl_interfaces/Log): keyword scan (replan, homing, battery, ...).

PX4 topics (battery, vehicle_status) need px4_msgs, which the container lacks;
neither bag contains a battery topic at all (see metadata topic lists).

Run: python3 mission_evidence.py <out.json> <bag.db3> [<bag.db3> ...]
Pass the .db3 FILE path, not the bag directory: the day-2 bags were renamed
after recording (first_long_flight.db3 etc.), so their metadata.yaml still
points at the original rosbag2_* filename and directory-open fails.
"""
import json
import re
import struct
import sys
from collections import Counter

import numpy as np
import rosbag2_py
from rclpy.serialization import deserialize_message
from geometry_msgs.msg import PoseStamped
from nav_msgs.msg import Path
from rcl_interfaces.msg import Log
from std_msgs.msg import String

POSE = "/shafterx2/dlio/odom_node/pose"
BT = "/shafterx2/bt_visualization/transitions"
EXEC = "/shafterx2/execution_complete"
HOMING = "/shafterx2/homing_path"
STRINGS = ["/shafterx2/planner_status", "/shafterx2/planner_mode",
           "/shafterx2/takeoff"]
ROSOUT = "/rosout"
WANT = [POSE, BT, EXEC, HOMING, ROSOUT] + STRINGS

STATUS = {0: "SUCCESS", 1: "RUNNING", 2: "FAILURE", 3: "IDLE", 4: "SKIPPED"}
KEYWORDS = ["replan", "plan", "homing", "home", "return", "battery", "volt",
            "fail", "invalid", "collision", "land", "mission", "budget",
            "abort", "time"]


def parse_bt(data):
    """Manual CDR parse of bt_msgs/BtTransition (little-endian CDR)."""
    off = 4  # skip encapsulation header

    def align(n):
        nonlocal off
        pad = (n - (off - 4) % n) % n
        off += pad

    def read_string():
        nonlocal off
        align(4)
        (n,) = struct.unpack_from("<I", data, off)
        off += 4
        s = data[off:off + n - 1].decode(errors="replace")  # strip null
        off += n
        return s

    viz = read_string()
    align(4)
    sec, nsec = struct.unpack_from("<iI", data, off)
    off += 8
    node = read_string()
    status, prev = struct.unpack_from("<BB", data, off)
    return viz, sec + nsec * 1e-9, node, status, prev


def read_bag(bag_dir):
    reader = rosbag2_py.SequentialReader()
    reader.open(rosbag2_py.StorageOptions(uri=bag_dir, storage_id="sqlite3"),
                rosbag2_py.ConverterOptions("cdr", "cdr"))
    reader.set_filter(rosbag2_py.StorageFilter(topics=WANT))

    poses, pose_ts = [], []
    bt_count = 0
    bt_by_node = Counter()
    bt_transition_kinds = Counter()
    bt_first_ts = {}       # node -> first ts of any transition
    bt_first_running = {}  # node -> first ts entering RUNNING
    bt_failures = Counter()
    exec_msgs = []
    homing_ts, homing_npts = [], []
    small_strings = {}
    rosout_hits = {k: [] for k in KEYWORDS}
    rosout_count = 0
    t0 = t1 = None

    while reader.has_next():
        topic, data, ts = reader.read_next()
        t0 = ts if t0 is None else t0
        t1 = ts
        if topic == POSE:
            m = deserialize_message(data, PoseStamped)
            poses.append((m.pose.position.x, m.pose.position.y,
                          m.pose.position.z))
            pose_ts.append(ts)
        elif topic == BT:
            _, _, node, status, prev = parse_bt(data)
            bt_count += 1
            bt_by_node[node] += 1
            kind = f"{node}: {STATUS.get(prev, prev)}->{STATUS.get(status, status)}"
            bt_transition_kinds[kind] += 1
            bt_first_ts.setdefault(node, ts)
            if status == 1:
                bt_first_running.setdefault(node, ts)
            if status == 2:
                bt_failures[node] += 1
        elif topic == EXEC:
            m = deserialize_message(data, String)
            exec_msgs.append((ts, m.data))
        elif topic == HOMING:
            m = deserialize_message(data, Path)
            homing_ts.append(ts)
            homing_npts.append(len(m.poses))
        elif topic in STRINGS:
            m = deserialize_message(data, String)
            small_strings.setdefault(topic, []).append((ts, m.data))
        elif topic == ROSOUT:
            m = deserialize_message(data, Log)
            rosout_count += 1
            low = m.msg.lower()
            for k in KEYWORDS:
                if k in low and len(rosout_hits[k]) < 40:
                    rosout_hits[k].append((ts, m.name, m.msg[:200]))

    P = np.array(poses, dtype=float)
    T = np.array(pose_ts, dtype=float)
    rel = lambda ts: round((ts - T[0]) / 1e9, 1)

    # Path length: 0.1 m downsample, reject >2 m teleport jumps.
    keep = [P[0]]
    for p in P[1:]:
        if np.linalg.norm(p - keep[-1]) >= 0.1:
            keep.append(p)
    K = np.array(keep)
    steps = np.linalg.norm(np.diff(K, axis=0), axis=1)
    path_len = float(steps[steps <= 2.0].sum())

    d_from_start = np.linalg.norm(P - P[0], axis=1)
    i_max = int(d_from_start.argmax())
    closure = P[-1] - P[0]

    result = {
        "bag": bag_dir.rstrip("/").split("/")[-1],
        "read_span_s": round((t1 - t0) / 1e9, 1),
        "pose": {
            "n_msgs": len(P),
            "duration_s": round((T[-1] - T[0]) / 1e9, 1),
            "path_length_m": round(path_len, 1),
            "max_dist_from_start_m": round(float(d_from_start[i_max]), 2),
            "max_dist_time_rel_s": rel(T[i_max]),
            "closure_first_to_last": {
                "dx_m": round(float(closure[0]), 3),
                "dy_m": round(float(closure[1]), 3),
                "dz_m": round(float(closure[2]), 3),
                "norm_m": round(float(np.linalg.norm(closure)), 3),
            },
            "z_min_m": round(float(P[:, 2].min()), 2),
            "z_max_m": round(float(P[:, 2].max()), 2),
            "z_range_m": round(float(P[:, 2].ptp()), 2),
            "start_xyz": [round(float(v), 3) for v in P[0]],
            "end_xyz": [round(float(v), 3) for v in P[-1]],
        },
        "bt_transitions": {
            "n_total": bt_count,
            "by_node": dict(bt_by_node.most_common()),
            "top_transition_kinds": dict(bt_transition_kinds.most_common(40)),
            "failures_by_node": dict(bt_failures.most_common()),
            "first_seen_rel_s": {n: rel(t) for n, t in sorted(
                bt_first_ts.items(), key=lambda kv: kv[1])},
            "first_running_rel_s": {n: rel(t) for n, t in sorted(
                bt_first_running.items(), key=lambda kv: kv[1])},
        },
        "execution_complete": {
            "n": len(exec_msgs),
            "payloads": Counter(d for _, d in exec_msgs).most_common(),
            "times_rel_s": [rel(t) for t, _ in exec_msgs],
        },
        "homing_path": {
            "n": len(homing_ts),
            "first_rel_s": rel(homing_ts[0]) if homing_ts else None,
            "last_rel_s": rel(homing_ts[-1]) if homing_ts else None,
            "times_rel_s": [rel(t) for t in homing_ts],
            "n_waypoints": homing_npts,
        },
        "small_string_topics": {
            t: [(rel(ts), d) for ts, d in v]
            for t, v in small_strings.items()},
        "rosout": {
            "n_total": rosout_count,
            "keyword_hits": {
                k: [(rel(ts), name, msg) for ts, name, msg in v]
                for k, v in rosout_hits.items() if v},
        },
    }
    return result


def main():
    out = sys.argv[1]
    results = []
    for bag in sys.argv[2:]:
        print("reading", bag, flush=True)
        results.append(read_bag(bag))
        print("  done:", results[-1]["bag"],
              results[-1]["pose"]["path_length_m"], "m path", flush=True)
    json.dump(results, open(out, "w"), indent=1)
    print("wrote", out)


if __name__ == "__main__":
    main()
