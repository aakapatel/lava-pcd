"""Extended Data panel: independent validation of the surface model and the
orbital-quality transfer test (dem_validation_panel.pdf).

a) Difference map, photogrammetric DEM minus ArcticDEM v4.1 (2 m), constant
   datum removed: shows the absence of doming or tilt at the decimetre level.
b) Distribution of the difference with the measured bounds annotated.
c) The overburden profile recomputed with the surface degraded to HiRISE-like
   and LROC-NAC-like DTM quality (single realisations), against the baseline
   and its systematic band: the profile a planetary mission would recover.

Same style contract as make_figures.py.
Run:  PYTHONPATH=src ANALYSIS_OUT=analysis_out_v6 \
      .venv/bin/python analysis/make_dem_validation_figure.py
"""
from __future__ import annotations

import csv
import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from pyproj import Transformer

sys.path.insert(0, str(Path(__file__).parent))
import orbital_dem_test as odt

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v6")))
FIGS = Path(os.environ.get("FIGS_DIR", str(
    ROOT.parent / "_Nature__Autonomous_aerial_reconnaissance_of_a_basaltic_"
    "lava_tube_reveals_interior_morphology_and_roof_thickness_distribution_"
    "for_planetary_subsurface" / "figures")))   # FIGS_DIR: review copies (v9)
# Elevation-axis label follows the vertical datum recorded by the run
# (ANALYSIS_OUT/datum.json, written by seed_v9_ortho.py); legacy runs carry
# WGS84 ellipsoidal heights and say so.
_DATUM = (json.loads((OUT / "datum.json").read_text())
          if (OUT / "datum.json").exists() else {})
ELEV_LABEL = _DATUM.get("elevation_axis_label", "elevation (m, WGS84 ellipsoidal)")
ORIGIN = np.array([479158.0, 7089826.0])

C = dict(tube="#0072B2", surface="#E69F00", intact="#009E73",
         hirise="#CC79A7", lroc="#D55E00")
plt.rcParams.update({
    "font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8,
    "legend.fontsize": 7, "xtick.labelsize": 7, "ytick.labelsize": 7,
    "axes.spines.top": False, "axes.spines.right": False,
    "pdf.fonttype": 42,
})


def main() -> None:
    val = json.loads((OUT / "national_dem_check.json").read_text())["arcticdem"]
    z, xmin, ymin, res = odt.load_dem()

    # recompute the difference cloud (as in national_dem_check, best shift)
    step = 4
    ix, iy = np.meshgrid(np.arange(0, z.shape[0], step),
                         np.arange(0, z.shape[1], step), indexing="ij")
    zz = z[ix, iy]
    ok = np.isfinite(zz)
    xs = ORIGIN[0] + xmin + ix[ok] * res
    ys = ORIGIN[1] + ymin + iy[ok] * res
    tr = Transformer.from_crs(32627, 3413, always_xy=True)
    bdx, bdy = val["best_shift_m"]
    X, Y = tr.transform(xs + bdx, ys + bdy)
    with rasterio.open(ROOT / "maps/national_dem_site.tif") as r:
        v = np.array([q[0] for q in r.sample(np.c_[X, Y])], dtype=np.float64)
        if r.nodata is not None:
            v[v == r.nodata] = np.nan
    v[np.abs(v) > 1e4] = np.nan
    m = np.isfinite(v)
    d = zz[ok][m] - v[m]
    d0 = d - np.median(d)

    fig = plt.figure(figsize=(7.1, 6.6))
    gs = fig.add_gridspec(2, 2, hspace=0.45, wspace=0.40,   # was 0.30: colourbar label touched panel b's ylabel
                          height_ratios=[1, 0.9], width_ratios=[1.5, 1])

    ax = fig.add_subplot(gs[0, 0])   # (a) difference map
    sc = ax.scatter(xs[m] - ORIGIN[0], ys[m] - ORIGIN[1], c=d0, s=2.5,
                    cmap="RdBu_r", vmin=-1, vmax=1, lw=0, rasterized=True)
    ax.set_aspect("equal")
    ax.set_xlabel("x (m, local frame)")
    ax.set_ylabel("y (m)")
    ax.set_title("a  Photogrammetric surface minus ArcticDEM\n"
                 "(constant datum removed)", loc="left")
    plt.colorbar(sc, ax=ax, label=r"$\Delta z$ (m)", shrink=0.75, pad=0.02)

    ax = fig.add_subplot(gs[0, 1])   # (b) histogram
    ax.hist(d0, bins=np.arange(-1.5, 1.55, 0.1), color=C["tube"], alpha=0.9)
    ax.axvline(0, color="k", lw=0.8)
    ax.annotate(f"std {val['std_raw_m']:.2f} m\n"
                f"doming {val['doming_amplitude_m']:.2f} m\n"
                f"tilt {val['tilt_m_per_km']:.2f} m km$^{{-1}}$",
                (0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=7, color="0.25")
    ax.set_xlabel(r"$\Delta z$ (m)")
    ax.set_ylabel("samples")
    ax.set_title("b  Difference distribution", loc="left")

    # (c) profile under degraded surfaces
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    s, x, y, zceil, tau0 = g("s"), g("x"), g("y"), g("z_ceil"), g("tau")
    it = np.array([r["class"] for r in rows]) == "intact"
    ax = fig.add_subplot(gs[1, :])
    o = np.argsort(s[it])
    ax.fill_between(s[it][o], tau0[it][o] - 1.5, tau0[it][o] + 1.5,
                    color=C["intact"], alpha=0.15, lw=0,
                    label=r"baseline $\pm\sigma_\tau$")
    ax.plot(s[it][o], tau0[it][o], color=C["intact"], lw=1.3,
            label="drone surface model (baseline)")
    for name, grid, sig, col in (("HiRISE-like 1 m / 0.3 m", 1.0, 0.3,
                                  C["hirise"]),
                                 ("LROC-NAC-like 5 m / 1 m", 5.0, 1.0,
                                  C["lroc"])):
        rng = np.random.default_rng(104)
        zd = odt.degrade(z, res, grid, sig, 0.0, rng)
        zt = odt.sample(zd, xmin, ymin, res, x, y)
        taud = zt - zceil
        ax.plot(s[it][o], taud[it][o], lw=0.8, color=col, alpha=0.85,
                label=name)
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("distance along tube $s$ (m)")
    ax.set_ylabel(r"overburden $\tau$ (m)")
    ax.set_title("c  The same profile measured against surfaces of orbital "
                 "DTM quality (single realisations)", loc="left")
    ax.legend(frameon=False, ncols=2, fontsize=6.5, loc="upper left")

    FIGS.mkdir(exist_ok=True)
    fig.savefig(FIGS / "dem_validation_panel.pdf", bbox_inches="tight",
                dpi=300)
    plt.close(fig)
    print("dem_validation_panel.pdf")


if __name__ == "__main__":
    main()
