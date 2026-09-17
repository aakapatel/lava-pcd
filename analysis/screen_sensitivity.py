"""Sensitivity of the overburden summary to the two exclusion screens.

Recomputes median, interquartile range and 5th/95th percentiles of tau over
the retained stations (a) with the aperture factor set to 1.2, 1.6 and 2.0
times the semi-major axis of each detected skylight (the pipeline uses 1.6)
and (b) with the four negative stations retained at their measured values.
Also the strength-demand tail (Equation 7, beta 0.5, rho 2,600) under the
tau > sigma_tau screen, a relaxed tau > 0.5 m screen and no screen.
Output: ANALYSIS_OUT/screen_sensitivity.json
Run:  ANALYSIS_OUT=analysis_out_v9 .venv/bin/python analysis/screen_sensitivity.py
"""
from __future__ import annotations
import csv, json, os
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
BETA, RHO, G = 0.5, 2600.0, 9.81


def main() -> None:
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    f = lambda k: np.array([float(r[k]) if r[k] != "nan" else np.nan for r in rows])
    x, y, tau, L = f("x"), f("y"), f("tau"), f("span")
    cls = np.array([r["class"] for r in rows])
    sigma_tau = float(rows[0]["sigma_tau"])
    holes = json.load(open(OUT / "aerial_holes.json"))["holes"][:3]

    def summary(k, keep_negative=False):
        ap = np.zeros(len(tau), bool)
        for h in holes:
            d = np.hypot(x - h["centroid"][0], y - h["centroid"][1])
            ap |= d <= k * h["semi_major"]
        m = np.isfinite(tau) & ~ap
        if not keep_negative:
            m &= tau > 0
        t = tau[m]
        return dict(n=int(m.sum()), median=float(np.median(t)), q25=float(np.percentile(t, 25)),
                    q75=float(np.percentile(t, 75)), p5=float(np.percentile(t, 5)),
                    p95=float(np.percentile(t, 95)))

    out = dict(aperture_factor={str(k): summary(k) for k in (1.2, 1.6, 2.0)},
               negatives_retained_factor_1p6=summary(1.6, True),
               pipeline=dict(n=int((cls == "intact").sum()), median=float(np.nanmedian(tau[cls == "intact"]))))
    it = cls == "intact"
    sig = BETA * RHO * G * L ** 2 / tau / 1e6
    out["strength_demand_MPa"] = {}
    for name, thr in (("tau_gt_sigma_tau", sigma_tau), ("tau_gt_0p5", 0.5), ("no_screen", 0.0)):
        m = it & (tau > thr)
        out["strength_demand_MPa"][name] = dict(threshold_m=thr, n=int(m.sum()), max=float(sig[m].max()),
                                                p95=float(np.percentile(sig[m], 95)),
                                                median=float(np.median(sig[m])))
    json.dump(out, open(OUT / "screen_sensitivity.json", "w"), indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
