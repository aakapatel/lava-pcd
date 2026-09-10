"""Extended Data Fig. 4 for the Nature build: rockfall inspection missions in
an operating underground mine with the same autonomy stack. The four panels
are the author's own figures 6.11, 6.13, 6.14 and 6.15 from the PhD thesis
(AP_PhD_Thesis.pdf in the project root), rendered at 250 dpi and cropped;
only vector panel letters are added. Output: analysis_out_v10/figures_nature/
ed4_deployments.{pdf,png} at 180 x <=170 mm in Nimbus Sans.

Run (from the analysis repo):
  env -u PYTHONPATH PYTHONPATH=src .venv/bin/python analysis/make_nature_ed4_deployments.py
"""
import os, subprocess, tempfile
from pathlib import Path
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
THESIS = ROOT.parent / "AP_PhD_Thesis.pdf"
OUT = ROOT / os.environ.get("NATURE_FIGS", "analysis_out_v10/figures_nature")
OUT.mkdir(parents=True, exist_ok=True)
MM = 1 / 25.4
plt.rcParams.update({"font.family": "Nimbus Sans", "font.size": 7, "pdf.fonttype": 3})

# (pdf page, fractional crop box x0,y0,x1,y1) for thesis figs 6.11, 6.13, 6.14, 6.15
PANELS = {
    "a": (161, (0.12, 0.12, 0.92, 0.83)),
    "b": (163, (0.14, 0.16, 0.91, 0.745)),
    "c": (164, (0.08, 0.11, 0.88, 0.40)),
    "d": (164, (0.08, 0.445, 0.88, 0.73)),
}


def render_page(page, dpi=250):
    tmp = Path(tempfile.mkdtemp())
    subprocess.run(["pdftoppm", "-r", str(dpi), "-f", str(page), "-l", str(page), "-png", str(THESIS), str(tmp / "p")], check=True)
    return Image.open(next(tmp.glob("p*.png"))).convert("RGB")


def crop(im, box):
    w, h = im.size
    c = im.crop((int(box[0] * w), int(box[1] * h), int(box[2] * w), int(box[3] * h)))
    a = np.asarray(c.convert("L")); m = a < 245; ys, xs = np.where(m)
    return c.crop((xs.min(), ys.min(), xs.max() + 1, ys.max() + 1))


pages = {}
imgs = {}
for k, (pg, box) in PANELS.items():
    pages.setdefault(pg, render_page(pg))
    imgs[k] = crop(pages[pg], box)

W, H = 180, 168
fig = plt.figure(figsize=(W * MM, H * MM))
# layout: top row a (left, tall progression) and b (right, third-person + graph);
# bottom row c and d (wide strips) side by side.
boxes = {"a": (0.0, 0.29, 0.41, 0.71), "b": (0.43, 0.29, 0.57, 0.71), "c": (0.0, 0.0, 0.49, 0.27), "d": (0.51, 0.0, 0.49, 0.27)}
for k, (x, y, w, h) in boxes.items():
    ax = fig.add_axes((x, y, w, h)); ax.axis("off")
    im = imgs[k]; iw, ih = im.size
    # fit inside the box preserving aspect, anchored top-left so letters sit on the image
    box_ar = (w * W) / (h * H); im_ar = iw / ih
    if im_ar > box_ar:      # wider than box: full width, image at the top
        hh = box_ar / im_ar; ext = (0, 1, 1 - hh, 1)
    else:                   # taller than box: full height, image at the left
        ww = im_ar / box_ar; ext = (0, ww, 0, 1)
    ax.imshow(im, extent=ext, aspect="auto", interpolation="lanczos")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    fig.text(x + 0.006, y + h - 0.008, k, fontsize=8, fontweight="bold", color="white", va="top", ha="left")
fig.savefig(OUT / "ed4_deployments.pdf", dpi=300)
fig.savefig(OUT / "ed4_deployments.png", dpi=300)
print("wrote", OUT / "ed4_deployments.pdf", {k: v.size for k, v in imgs.items()})
