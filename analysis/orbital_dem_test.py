"""Orbital-quality DEM degradation test (internal; NOT a manuscript input).

Question (editor item E47): if the surface reference were an orbital stereo
DTM instead of a drone photogrammetric model, would the overburden profile
and its headline structure survive? We degrade our own 0.5 m DEM to the
published quality of planetary stereo products and recompute tau at the v6
stations, leaving the interior map and registration untouched (the interior
side of a planetary campaign is the vehicle's own lidar, unaffected by the
orbital DTM; the DTM enters through z_DEM and through the availability of
rim anchors).

Degradation model per configuration:
  - posting: block-average the 0.5 m grid to G m and sample nearest;
  - correlated vertical noise: white noise Gaussian-smoothed to a 25 m
    correlation length, scaled to sigma_z (stereo-matching error is spatially
    correlated, which is the pessimistic case for a profile);
  - optional doming: a radial quadratic warp of amplitude A (peak to trough)
    across the footprint, the classic uncontrolled-SfM/jitter systematic.

Configurations follow published product specs (parameter sources recorded in
maps/national_dem_site.json when the companion download runs, and in the
session report): HiRISE-like (1 m, 0.3 m), (2 m, 0.5 m); LROC-NAC-like
(5 m, 1 m), (5 m, 2 m); CTX-like (20 m, 3 m); doming 1 and 2 m variants.

Metrics over the 318 intact stations, 20 noise seeds per configuration:
median-tau shift, per-station RMS error, fraction beyond the 1.5 m budget,
false negative-roof fraction, and survival of the thick-thin contrast
(median tau beyond s=250 minus median before s=150). Plus the anchor
question: which skylight apertures remain resolvable at each posting
(3-post detection rule), against the catalogued planetary aperture sizes.

Output: analysis_out_v6/orbital_dem_test.json, fig_orbital_dem_test.png.
Run:  PYTHONPATH=src ANALYSIS_OUT=analysis_out_v6 \
      .venv/bin/python analysis/orbital_dem_test.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v6")))

CORR_LEN_M = 25.0
N_SEEDS = 20
SIGMA_TAU = 1.5
CONFIGS = [
    dict(name="control 0.5 m / 0 m", grid=0.5, sigma=0.0, dome=0.0),
    dict(name="HiRISE-like 1 m / 0.3 m", grid=1.0, sigma=0.3, dome=0.0),
    dict(name="HiRISE-like 2 m / 0.5 m", grid=2.0, sigma=0.5, dome=0.0),
    dict(name="LROC-NAC-like 5 m / 1 m", grid=5.0, sigma=1.0, dome=0.0),
    dict(name="LROC-NAC-like 5 m / 2 m", grid=5.0, sigma=2.0, dome=0.0),
    dict(name="CTX-like 20 m / 3 m", grid=20.0, sigma=3.0, dome=0.0),
    dict(name="doming 1 m (0.5 m grid)", grid=0.5, sigma=0.0, dome=1.0),
    dict(name="doming 2 m (0.5 m grid)", grid=0.5, sigma=0.0, dome=2.0),
    dict(name="LROC-NAC-like 5 m / 1 m + doming 2 m",
         grid=5.0, sigma=1.0, dome=2.0),
]
APERTURES = dict(S1=2.83, S2=4.56, S3=9.34)  # long axes, m (surface view)


def load_dem():
    d = np.load(OUT / "dem_grid_full_surface_and_subsurface_merged.npz")
    return (d["z"].astype(np.float64), float(d["xmin"]), float(d["ymin"]),
            float(d["res"]))


def load_stations():
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    return (g("s"), g("x"), g("y"), g("z_ceil"), g("tau"),
            np.array([r["class"] for r in rows]))


def degrade(z, res, grid, sigma, dome, rng):
    """Return a degraded copy of the DEM raster (NaN-aware)."""
    filled = np.isfinite(z)
    _, (ii, jj) = ndi.distance_transform_edt(~filled, return_indices=True)
    zf = z[ii, jj]
    f = max(int(round(grid / res)), 1)
    if f > 1:
        ny = (zf.shape[0] // f) * f
        nx = (zf.shape[1] // f) * f
        blocks = zf[:ny, :nx].reshape(ny // f, f, nx // f, f).mean(axis=(1, 3))
        zg = np.repeat(np.repeat(blocks, f, axis=0), f, axis=1)
        out = zf.copy()
        out[:ny, :nx] = zg
    else:
        out = zf.copy()
    if sigma > 0:
        k = CORR_LEN_M / res
        w = rng.standard_normal(out.shape)
        w = ndi.gaussian_filter(w, k)
        w /= max(w.std(), 1e-9)
        out = out + sigma * w
    if dome > 0:
        u = (np.arange(out.shape[0]) / out.shape[0]) * 2 - 1
        v = (np.arange(out.shape[1]) / out.shape[1]) * 2 - 1
        r2 = u[:, None] ** 2 + v[None, :] ** 2
        out = out + dome * (r2 / r2.max() - 0.5)
    out[~filled] = np.nan
    return out


def sample(z, xmin, ymin, res, x, y):
    # the DEM cache is x-major: z[ix, iy] (see run_roof.build_dem)
    i = ((x - xmin) / res).astype(int)
    j = ((y - ymin) / res).astype(int)
    ok = (i >= 0) & (i < z.shape[0]) & (j >= 0) & (j < z.shape[1])
    v = np.full(len(x), np.nan)
    v[ok] = z[i[ok], j[ok]]
    return v


def main() -> None:
    z0, xmin, ymin, res = load_dem()
    s, x, y, zceil, tau0, klass = load_stations()
    it = klass == "intact"
    thick0 = np.median(tau0[it & (s >= 250)]) - np.median(tau0[it & (s < 150)])

    results = []
    for cfg in CONFIGS:
        seeds = N_SEEDS if cfg["sigma"] > 0 else 1
        med_shift, rms, frac_big, frac_neg, contrast = [], [], [], [], []
        for k in range(seeds):
            rng = np.random.default_rng(100 + k)
            zd = degrade(z0, res, cfg["grid"], cfg["sigma"], cfg["dome"], rng)
            zt = sample(zd, xmin, ymin, res, x, y)
            taud = zt - zceil
            d = taud[it] - tau0[it]
            med_shift.append(np.nanmedian(taud[it]) - np.nanmedian(tau0[it]))
            rms.append(np.sqrt(np.nanmean(d ** 2)))
            frac_big.append(np.nanmean(np.abs(d) > SIGMA_TAU))
            frac_neg.append(np.nanmean(taud[it] <= 0))
            contrast.append(np.nanmedian(taud[it & (s >= 250)])
                            - np.nanmedian(taud[it & (s < 150)]))
        results.append(dict(
            config=cfg["name"], grid_m=cfg["grid"], sigma_m=cfg["sigma"],
            dome_m=cfg["dome"], n_seeds=seeds,
            median_shift_m=dict(mean=float(np.mean(med_shift)),
                                worst=float(np.max(np.abs(med_shift)))),
            per_station_rms_m=dict(mean=float(np.mean(rms)),
                                   worst=float(np.max(rms))),
            frac_error_gt_1p5=dict(mean=float(np.mean(frac_big)),
                                   worst=float(np.max(frac_big))),
            frac_false_negative_roof=float(np.max(frac_neg)),
            thick_thin_contrast_m=dict(baseline=float(thick0),
                                       mean=float(np.mean(contrast)),
                                       worst=float(np.min(contrast)))))

    anchors = {g: {k: bool(a >= 3 * g) for k, a in APERTURES.items()}
               for g in (1.0, 2.0, 5.0, 20.0)}

    out = dict(
        note="Internal E47 test: interior map and registration held fixed; "
             "only z_DEM degraded. Doming is applied to the whole footprint "
             "(pessimistic for the along-tube trend). Not a manuscript input.",
        corr_len_m=CORR_LEN_M, n_intact=int(it.sum()),
        baseline=dict(median=float(np.nanmedian(tau0[it])),
                      thick_thin_contrast_m=float(thick0)),
        configs=results,
        skylight_anchor_resolvable_3post=anchors,
        planetary_aperture_context=dict(
            mare_pit_inner_median_m=100.0, apc_outer_median_m=120.0,
            note="catalogued planetary anchors are 10-40x larger than the "
                 "terrestrial skylights used here, so anchor detectability "
                 "at orbital posting is easier at the targets than at the "
                 "analogue"))
    (OUT / "orbital_dem_test.json").write_text(json.dumps(out, indent=2))

    fig, axes = plt.subplots(1, 2, figsize=(10, 3.6))
    names = [r["config"] for r in results]
    axes[0].barh(names, [r["per_station_rms_m"]["mean"] for r in results],
                 color="#0072B2")
    axes[0].axvline(SIGMA_TAU, color="k", ls="--", lw=1,
                    label="stated budget 1.5 m")
    axes[0].set_xlabel("per-station RMS tau error (m)")
    axes[0].legend(frameon=False, fontsize=7)
    axes[1].barh(names, [r["thick_thin_contrast_m"]["mean"] for r in results],
                 color="#009E73")
    axes[1].axvline(thick0, color="k", ls="--", lw=1, label="baseline")
    axes[1].set_xlabel("thick-thin contrast (m)")
    axes[1].legend(frameon=False, fontsize=7)
    axes[1].set_yticklabels([])
    fig.tight_layout()
    fig.savefig(OUT / "fig_orbital_dem_test.png", dpi=140)
    print(json.dumps(out["configs"], indent=1))
    print("anchors:", json.dumps(anchors))


if __name__ == "__main__":
    main()
