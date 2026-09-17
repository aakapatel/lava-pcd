"""Web copy of the rim ground-truth panels (Extended Data Fig. 5g, h) as one
180 x 62 mm figure in the paper's style, for the project page.
Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/make_web_groundtruth_figure.py
Output: analysis_out_v10/figures_nature/web_rim_groundtruth.{pdf,png}
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import make_nature_ed_figures as ned  # noqa: E402  (sets fonts and palette)


def main() -> None:
    recs, _ = ned.load_rim_groundtruth()
    fig = ned.new_fig(62.0)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.55, 1])
    ax_a = fig.add_subplot(gs[0, 0]); ned.rim_thickness(ax_a, recs)
    ax_b = fig.add_subplot(gs[0, 1]); ned.rim_floor(ax_b, recs)
    ned.layout(fig)
    ned.freeze_layout(fig)
    ned.place_letters(fig, [(ax_a, "a"), (ax_b, "b")])
    ned.save(fig, "web_rim_groundtruth")


if __name__ == "__main__":
    main()
