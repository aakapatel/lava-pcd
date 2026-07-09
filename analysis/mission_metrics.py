"""T7: per-flight mission metrics from the ROS 2 bags. Runs INSIDE the stage_cuda
container (ROS 2 Humble + px4_msgs from the exploration workspace).

For each bag, using a storage filter so only the light topics are read:
  - duration (from message timestamps)
  - travelled path length + straight-line extent, from the in-flight DLIO pose
    (/shafterx2/dlio/odom_node/pose), the odometry the manuscript describes
  - autonomy / intervention check from PX4 /fmu/out/vehicle_status nav_state
    (OFFBOARD = 14 is autonomous; any MANUAL/POSCTL/ALTCTL would be a takeover)
  - event counts: BT transitions, execution_complete, planner replans

Writes JSON to the mounted output dir.
Run (in container): python3 mission_metrics.py <out.json> <bag_dir> [<bag_dir> ...]
"""
import json
import sys

import numpy as np
import rosbag2_py
from rclpy.serialization import deserialize_message
from geometry_msgs.msg import PoseStamped

POSE_TOPIC = "/shafterx2/dlio/odom_node/pose"
STATUS_TOPIC = "/fmu/out/vehicle_status"
MANUAL_TOPIC = "/fmu/out/manual_control_setpoint"
BT_TOPIC = "/shafterx2/bt_visualization/transitions"
EXEC_TOPIC = "/shafterx2/execution_complete"
# Only POSE is deserialized (geometry_msgs, always available). The PX4 topics
# use px4_msgs, which this container lacks, so they are only COUNTED here; the
# definitive nav_state intervention check needs px4_msgs built in the workspace.
WANT = [POSE_TOPIC, STATUS_TOPIC, MANUAL_TOPIC, BT_TOPIC, EXEC_TOPIC]


def read_bag(bag_dir):
    reader = rosbag2_py.SequentialReader()
    reader.open(rosbag2_py.StorageOptions(uri=bag_dir, storage_id="sqlite3"),
                rosbag2_py.ConverterOptions("cdr", "cdr"))
    reader.set_filter(rosbag2_py.StorageFilter(topics=WANT))
    poses = []
    counts = {STATUS_TOPIC: 0, MANUAL_TOPIC: 0, BT_TOPIC: 0, EXEC_TOPIC: 0}
    t0 = t1 = None
    while reader.has_next():
        topic, data, ts = reader.read_next()
        t0 = ts if t0 is None else t0
        t1 = ts
        if topic == POSE_TOPIC:
            m = deserialize_message(data, PoseStamped)
            poses.append((m.pose.position.x, m.pose.position.y, m.pose.position.z))
        elif topic in counts:
            counts[topic] += 1
    P = np.array(poses, dtype=float)
    if len(P) > 1:
        # Downsample the ~100 Hz pose stream to 0.1 m spacing before summing, so
        # per-sample jitter does not inflate the travelled path length; reject
        # outlier teleport jumps (> 2 m between kept samples).
        keep = [P[0]]
        for p in P[1:]:
            if np.linalg.norm(p - keep[-1]) >= 0.1:
                keep.append(p)
        K = np.array(keep)
        steps = np.linalg.norm(np.diff(K, axis=0), axis=1)
        path_len = float(steps[steps <= 2.0].sum())
        extent = float(np.linalg.norm(P - P[0], axis=1).max())
        z_desc = float(P[:, 2].max() - P[:, 2].min())
    else:
        path_len = extent = z_desc = 0.0
    return dict(
        bag=bag_dir.rstrip("/").split("/")[-1],
        duration_s=round((t1 - t0) / 1e9, 1) if t0 else None,
        n_pose=len(P),
        path_length_m=round(path_len, 1),
        straight_line_extent_m=round(extent, 1),
        z_range_m=round(z_desc, 1),
        n_vehicle_status=counts[STATUS_TOPIC],
        n_manual_control_setpoint=counts[MANUAL_TOPIC],
        n_bt_transitions=counts[BT_TOPIC],
        n_execution_complete=counts[EXEC_TOPIC],
    )


def main():
    out = sys.argv[1]
    results = []
    for bag in sys.argv[2:]:
        print("reading", bag, flush=True)
        try:
            results.append(read_bag(bag))
            print("  ", json.dumps(results[-1]), flush=True)
        except Exception as e:
            print("  ERROR", e, flush=True)
            results.append(dict(bag=bag, error=str(e)))
    json.dump(results, open(out, "w"), indent=2)
    print("wrote", out)


if __name__ == "__main__":
    main()
