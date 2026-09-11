"""Vector redraw of the mission behaviour tree (Extended Data Fig. 2b).

The tree is transcribed node for node from the Groot screenshot of the
flown tree (manuscript figures/stage_bt.png, 5883 x 4317 px): the same
node types, the same instance names and the same child order. Groot's
"Root" container is not a behaviour-tree node and is omitted; the tree
starts at RootSequence.

Drawing conventions (all vector, Nimbus Sans, 6 to 7 pt at final size):
    composite nodes  rounded box, type icon (arrow: Sequence, question
                     mark: Fallback) in colour, instance name in words
    decorator        rounded box, circular-arrow icon, "Retry until successful"
    action leaves    box with a lightning icon; consecutive leaf children
                     of one parent are stacked vertically in one column
                     (their order top to bottom is the execution order)
Children of a node execute left to right. Layout is a contour-based tidy
tree (Reingold-Tilford style with a 1 mm vertical raster) so that shallow
branches may sit above the wide, deep part of the tree, as in Groot.

Run standalone (repo root) to write stage_bt.{pdf,svg,png} into
analysis_out_v10/figures_nature/:
    env -u PYTHONPATH PYTHONPATH=src .venv/bin/python analysis/behaviour_tree.py
The same draw_tree() is called by make_nature_ed_figures.py for ED Fig. 2.
"""
from __future__ import annotations

import re
import sys
from dataclasses import dataclass, field
from pathlib import Path

import matplotlib
import numpy as np
from matplotlib.font_manager import FontProperties
from matplotlib.patches import FancyBboxPatch, PathPatch, Polygon
from matplotlib.path import Path as MPath
from matplotlib.textpath import TextToPath

MM = 1.0 / 25.4
FONT = "Nimbus Sans"

# Groot's palette, slightly darkened for print: composites purple,
# decorators blue, leaves black.
C_SEQ = "#7b1fa2"
C_FB = "#7b1fa2"
C_DEC = "#1565c0"
C_LEAF = "#000000"
C_BOX = "#f4f4f4"
C_EDGE = "#8a8a8a"
C_LINE = "#4a4a4a"

# ----------------------------------------------------------------------
# the tree (type, name, children); transcribed from the Groot screenshot
# ----------------------------------------------------------------------
SEQ, FB, DEC, ACT = "Sequence", "Fallback", "Decorator", "Action"


def S(name, *ch): return (SEQ, name, list(ch))
def F(name, *ch): return (FB, name, list(ch))
def R(*ch): return (DEC, "RetryUntilSuccessful", list(ch))
def A(name): return (ACT, name, [])


TREE = S("RootSequence",
         F("HealthCheckFallback",
           S("HealthSequence", A("CheckSensors"), A("CheckOdometry")),
           A("SetPlannerStatus")),
         F("PlannerStatusFallback",
           S("IfExecuting", A("CheckPlannerStatus"), A("DoNothing")),
           S("IfReady",
             A("CheckPlannerStatus"),
             F("PlannerModeSelector",
               S("WaypointNavSequence",
                 A("CheckPlannerMode"),
                 R(S("WaypointLocalPlanExec",
                     A("RunWaypointPlanner"), A("CheckPlanValidity"),
                     A("ExecutePath")))),
               S("AutonomousExplorationSequence",
                 A("CheckPlannerMode"),
                 F("LocalOrGlobalOrHomeFallback",
                   R(S("AutoLocalPlanExec",
                       A("RunLocalPlanner"), A("CheckPlanValidity"),
                       A("ExecutePath"))),
                   R(S("AutoGlobalPlanExec",
                       A("RunGlobalPlanner"), A("CheckPlanValidity"),
                       A("ExecutePath"))),
                   S("ReturnHomeSequence",
                     A("PlanHoming"), A("ExecutePath"),
                     A("SetPlannerStatus")))),
               S("ReturnHomeSequence",
                 A("CheckPlannerMode"), A("PlanHoming"), A("ExecutePath"),
                 A("SetPlannerStatus"))))))


def words(camel: str) -> str:
    """CamelCase node name -> words ("CheckPlanValidity" -> "Check plan validity")."""
    w = re.findall(r"[A-Z][a-z]*|[a-z]+", camel)
    out = [w[0]] + [x.lower() if x not in ("IMU",) else x for x in w[1:]]
    return " ".join(out)


# ----------------------------------------------------------------------
# text measurement and wrapping
# ----------------------------------------------------------------------
_T2P = TextToPath()


def text_w_mm(s: str, size: float, bold=False) -> float:
    fp = FontProperties(family=FONT, size=size, weight="bold" if bold else "normal")
    w, _, _ = _T2P.get_text_width_height_descent(s, fp, ismath=False)
    return w / 72.0 * 25.4


def wrap_mm(s: str, width_mm: float, size: float) -> list[str]:
    lines, cur = [], ""
    for w in s.split(" "):
        t = (cur + " " + w).strip()
        if cur and text_w_mm(t, size) > width_mm:
            lines.append(cur); cur = w
        else:
            cur = t
    if cur:
        lines.append(cur)
    return lines


# ----------------------------------------------------------------------
# layout
# ----------------------------------------------------------------------
@dataclass
class Box:
    kind: str
    name: str
    lines: list[str]
    w: float
    h: float
    x: float = 0.0          # centre x (mm), relative until placed
    y: float = 0.0          # top y (mm, downward)
    children: list = field(default_factory=list)
    stack: list = field(default_factory=list)   # stacked leaf boxes (kind == "stack")
    parent: "Box | None" = None


@dataclass
class Style:
    size: float = 6.0        # node text (pt)
    leaf_w_max: float = 14.0  # wrap width of leaf text (mm)
    pad_x: float = 1.3       # box padding (mm)
    icon_w: float = 2.9      # icon column inside a box (mm)
    line_h: float = 2.65     # text line pitch (mm) at 6 pt
    pad_y: float = 0.8
    gap_x: float = 2.0       # minimum horizontal gap between boxes (mm)
    gap_y: float = 3.8       # vertical gap parent bottom -> child top (mm)
    comp_w_shallow: float = 34.0   # wrap width of composite names, depth <= 3 (mm)
    comp_w_deep: float = 16.0      # wrap width of composite names, deeper rows (mm)
    stack_gap: float = 0.8   # gap between stacked leaves (mm)
    spine: float = 1.3       # offset of a stack's spine left of its boxes (mm)
    lw: float = 0.6


def build(node, st: Style, depth: int = 0) -> Box:
    kind, name, ch = node
    label = "Retry until successful" if kind == DEC else words(name)
    if kind == ACT:
        lines = wrap_mm(label, st.leaf_w_max, st.size)
    else:
        lines = wrap_mm(label, st.comp_w_shallow if depth <= 3 else st.comp_w_deep, st.size)
    tw = max(text_w_mm(l, st.size, bold=(kind == ACT)) for l in lines)
    w = st.icon_w + tw + 2 * st.pad_x
    h = len(lines) * st.line_h + 2 * st.pad_y
    b = Box(kind, name, lines, w, h)
    kids = [build(c, st, depth + 1) for c in ch]
    # group consecutive leaves into a stack column
    out, run = [], []

    def flush():
        if len(run) >= 2:
            w_s = max(k.w for k in run)
            for k in run:
                k.w = w_s
            h_s = sum(k.h for k in run) + st.stack_gap * (len(run) - 1)
            sb = Box("stack", "", [], w_s, h_s); sb.stack = list(run)
            out.append(sb)
        else:
            out.extend(run)
        run.clear()
    for k in kids:
        if k.kind == ACT:
            run.append(k)
        else:
            flush(); out.append(k)
    flush()
    b.children = out
    for k in out:
        k.parent = b
    return b


def layout(root: Box, st: Style):
    """Contour-based tidy layout. Returns (width, height) in mm and sets
    absolute x (centre) and y (top) on every box."""
    def place(b: Box, y: float):
        """Position subtree with b at top y; returns contour dict
        {row_mm: [xmin, xmax]} in coordinates relative to b.x = 0."""
        b.y = y
        cont = {}

        def add(c, x0, x1, y0, y1):
            for r in range(int(np.floor(y0)), int(np.ceil(y1))):
                if r in c:
                    c[r][0] = min(c[r][0], x0); c[r][1] = max(c[r][1], x1)
                else:
                    c[r] = [x0, x1]
        add(cont, -b.w / 2, b.w / 2, y, y + b.h)
        if b.kind == "stack":
            yy = y
            for k in b.stack:
                k.y = yy; k.x = 0.0; yy += k.h + st.stack_gap
            return cont
        if not b.children:
            return cont
        cy = y + b.h + st.gap_y
        subs = [place(k, cy) for k in b.children]
        # pack children left to right with contour separation
        offs = [0.0]
        merged = {r: list(v) for r, v in subs[0].items()}
        for sc in subs[1:]:
            need = -1e9
            for r, (x0, _) in sc.items():
                if r in merged:
                    need = max(need, merged[r][1] + st.gap_x - x0)
            off = need if need > -1e8 else offs[-1] + st.gap_x
            offs.append(off)
            for r, (x0, x1) in sc.items():
                if r in merged:
                    merged[r][0] = min(merged[r][0], x0 + off); merged[r][1] = max(merged[r][1], x1 + off)
                else:
                    merged[r] = [x0 + off, x1 + off]
        # centre parent over first and last child
        mid = (offs[0] + offs[-1]) / 2
        for k, off, sc in zip(b.children, offs, subs):
            k.x = off - mid            # relative to parent
        for r, (x0, x1) in merged.items():
            add(cont, x0 - mid, x1 - mid, r, r + 1)
        return cont

    cont = place(root, 0.0)
    xmin = min(v[0] for v in cont.values()); xmax = max(v[1] for v in cont.values())
    ymax = max(cont.keys()) + 1

    def absolute(b: Box, px: float):
        b.x = px + b.x
        for k in b.children:
            absolute(k, b.x)
        for k in b.stack:
            k.x = b.x
    root.x = -xmin
    for k in root.children:
        absolute(k, root.x)
    return xmax - xmin, ymax, cont


def all_boxes(b: Box):
    yield b
    for k in b.children:
        yield from all_boxes(k)
    for k in b.stack:
        yield k


# ----------------------------------------------------------------------
# drawing (axes in mm, y downward)
# ----------------------------------------------------------------------
def icon(ax, kind, cx, cy, s, st: Style):
    """Type icon centred at (cx, cy), nominal size s mm."""
    if kind == SEQ:
        ax.annotate("", xy=(cx + s / 2, cy), xytext=(cx - s / 2, cy),
                    arrowprops=dict(arrowstyle="-|>,head_width=0.22,head_length=0.32",
                                    color=C_SEQ, lw=1.1, shrinkA=0, shrinkB=0))
    elif kind == FB:
        ax.text(cx, cy, "?", ha="center", va="center", fontsize=st.size + 2.5,
                fontweight="bold", color=C_FB, family=FONT)
    elif kind == DEC:
        th = np.linspace(np.deg2rad(-40), np.deg2rad(250), 40)
        r = s * 0.36
        ax.plot(cx + r * np.cos(th), cy - r * np.sin(th), color=C_DEC, lw=1.0,
                solid_capstyle="round")
        # arrowhead at the arc start
        a0 = np.deg2rad(-40)
        px, py = cx + r * np.cos(a0), cy - r * np.sin(a0)
        tx, ty = -np.sin(a0), -np.cos(a0)       # tangent direction (y down)
        hd = s * 0.28
        ax.add_patch(Polygon([[px, py], [px - hd * tx + hd * 0.55 * ty, py - hd * ty - hd * 0.55 * tx],
                              [px - hd * tx - hd * 0.55 * ty, py - hd * ty + hd * 0.55 * tx]],
                             closed=True, color=C_DEC, lw=0))
    else:  # lightning bolt
        p = np.array([[0.15, -0.5], [-0.35, 0.08], [-0.02, 0.08], [-0.15, 0.5],
                      [0.35, -0.1], [0.02, -0.1]]) * s * 0.9
        ax.add_patch(Polygon(np.c_[cx + p[:, 0], cy + p[:, 1]], closed=True,
                             color=C_LEAF, lw=0))


def draw_box(ax, b: Box, st: Style):
    x0, y0 = b.x - b.w / 2, b.y
    ax.add_patch(FancyBboxPatch((x0, y0), b.w, b.h,
                                boxstyle="round,pad=0,rounding_size=0.8",
                                fc=C_BOX, ec=C_EDGE, lw=st.lw))
    icon(ax, b.kind, x0 + st.pad_x + st.icon_w * 0.42, y0 + b.h / 2, 2.6, st)
    col = {SEQ: C_SEQ, FB: C_FB, DEC: C_DEC}.get(b.kind, C_LEAF)
    bold = b.kind in (ACT, DEC)
    n = len(b.lines)
    for i, line in enumerate(b.lines):
        cy = y0 + b.h / 2 + (i - (n - 1) / 2) * st.line_h
        ax.text(x0 + st.pad_x + st.icon_w, cy, line, ha="left", va="center",
                fontsize=st.size, color=col if b.kind != ACT else C_LEAF,
                fontweight="bold" if bold else "normal", family=FONT)


def edge(ax, x0, y0, x1, y1, st: Style):
    """Groot-style S-curve from a parent's bottom port to a child's top port."""
    ym = (y0 + y1) / 2
    path = MPath([(x0, y0), (x0, ym), (x1, ym), (x1, y1)],
                 [MPath.MOVETO, MPath.CURVE4, MPath.CURVE4, MPath.CURVE4])
    ax.add_patch(PathPatch(path, fc="none", ec=C_LINE, lw=st.lw))
    ax.plot([x0, x1], [y0, y1], "o", ms=1.6, color=C_LINE, mec="none")


def draw_tree(ax, root: Box, st: Style):
    for b in all_boxes(root):
        if b.kind == "stack":
            # spine from the parent port down the stack's left side, with a
            # short connector into each leaf
            xs = b.x - b.w / 2 - st.spine
            top = b.stack[0].y
            for k in b.stack:
                ax.plot([xs, xs], [top if k is b.stack[0] else k.y - st.stack_gap, k.y + k.h / 2],
                        color=C_LINE, lw=st.lw, solid_capstyle="round")
                ax.plot([xs, b.x - b.w / 2], [k.y + k.h / 2, k.y + k.h / 2],
                        color=C_LINE, lw=st.lw, solid_capstyle="round")
            continue
        if b.kind != "stack" and b.name:
            draw_box(ax, b, st)
    for b in all_boxes(root):
        for k in b.children:
            if k.kind == "stack":
                xs = k.x - k.w / 2 - st.spine
                edge(ax, b.x, b.y + b.h, xs, k.stack[0].y, st)
            else:
                edge(ax, b.x, b.y + b.h, k.x, k.y, st)


def draw_key(ax, x, y, st: Style):
    """Small on-figure key: icon + type name, four rows."""
    rows = [(SEQ, "Sequence: runs children in order, fails on the first failure"),
            (FB, "Fallback: tries children in order until one succeeds"),
            (DEC, "Retry decorator: repeats its child until it succeeds"),
            (ACT, "Action: a check or a command executed on the robot")]
    for i, (kind, txt) in enumerate(rows):
        cy = y + i * 3.4
        icon(ax, kind, x + 1.4, cy, 2.6, st)
        ax.text(x + 3.8, cy, txt, ha="left", va="center", fontsize=st.size + 0.5,
                family=FONT, color="#222222")


def render(ax, w_mm: float, h_mm: float, st: Style | None = None, key=True):
    """Lay out and draw the tree into `ax`, whose data coordinates are set to
    millimetres (x right, y down) over w_mm x h_mm. Returns the natural
    (width, height) of the tree in mm."""
    st = st or Style()
    root = build(TREE, st)
    tw, th, cont = layout(root, st)
    # centre horizontally, top-align
    dx = (w_mm - tw) / 2 if tw < w_mm else 0.0
    for b in all_boxes(root):
        b.x += dx
    ax.set_xlim(0, w_mm); ax.set_ylim(h_mm, 0); ax.set_aspect("equal")
    ax.set_axis_off()
    draw_tree(ax, root, st)
    if key:
        draw_key(ax, w_mm - 74.0, 2.0, st)
    return tw, th


def main():
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    root = Path(__file__).resolve().parents[1]
    out = root / "analysis_out_v10" / "figures_nature"
    out.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "sans-serif", "font.sans-serif": [FONT],
                         "pdf.fonttype": 42, "svg.fonttype": "none"})
    st = Style()
    probe = build(TREE, st)
    tw, th, _ = layout(probe, st)
    W, H = max(180.0, tw), th + 1.0
    fig = plt.figure(figsize=(W * MM, H * MM))
    ax = fig.add_axes([0, 0, 1, 1])
    tw, th = render(ax, W, H, st)
    for ext in ("pdf", "svg"):
        fig.savefig(out / f"stage_bt.{ext}")
    fig.savefig(out / "stage_bt.png", dpi=300, facecolor="white")
    print(f"tree {tw:.1f} x {th:.1f} mm -> {out}/stage_bt.{{pdf,svg,png}}")


if __name__ == "__main__":
    main()
