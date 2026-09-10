"""Strength demand (clamped-strip model, Methods Eq. envelope) at the stations
approaching the road crossing, for the Results sentence that compares the
road reach with the survey median. Reads analysis_out_v9 (or ANALYSIS_OUT);
writes <out>/road_strength.json. Run:
  env -u PYTHONPATH PYTHONPATH=src .venv/bin/python analysis/road_strength_demand.py
"""
import json, os
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / os.environ.get("ANALYSIS_OUT", "analysis_out_v9")
RHO, G, BETA = 2600.0, 9.81, 0.5          # as in run_review_stats / Methods
S_FROM = float(os.environ.get("ROAD_S_FROM", "280"))

m = pd.read_csv(OUT / "morphometry.csv")
r = pd.read_csv(OUT / "roof_thickness.csv")
d = r.merge(m[["s", "width"]], on="s")
sig_tau = float(r["sigma_tau"].dropna().median())
ok = (d["class"] == "intact") & (d["tau"] > sig_tau)   # same rule as the paper's 248-station set
sig = BETA * RHO * G * d["width"] ** 2 / d["tau"] / 1e6
allv = sig[ok]
road = sig[ok & (d["s"] >= S_FROM)]
res = {
    "note": f"sigma_req = beta*rho*g*L^2/tau (MPa) with L = robust width, over intact stations with tau > sigma_tau; road reach = s >= {S_FROM:g} m",
    "rho": RHO, "g": G, "beta": BETA, "tau_threshold_m": sig_tau,
    "n_all": int(ok.sum()), "median_all_MPa": float(allv.median()), "p95_all_MPa": float(allv.quantile(0.95)),
    "n_road_reach": int(road.size), "min_road_reach_MPa": float(road.min()), "median_road_reach_MPa": float(road.median()), "max_road_reach_MPa": float(road.max()),
    "frac_all_below_road_median": float((allv < road.median()).mean()),
    "width_road_reach_m": [float(d.loc[ok & (d.s >= S_FROM), "width"].min()), float(d.loc[ok & (d.s >= S_FROM), "width"].max())],
}
(OUT / "road_strength.json").write_text(json.dumps(res, indent=2))
print(json.dumps(res, indent=2))
