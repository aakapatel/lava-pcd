"""Collect every manuscript number into analysis_out/paper_numbers.json.

The manuscript pulls its data numbers from this single file (plan/04 rule:
no data number is typed into the tex by hand). Re-run after any upstream
pipeline changes:  .venv/bin/python analysis/collect_paper_numbers.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np

OUT = Path(os.environ.get("ANALYSIS_OUT",
                          str(Path(__file__).resolve().parents[1] / "analysis_out")))


def main() -> None:
    reg = json.loads((OUT / "registration_report.json").read_text())
    bud = json.loads((OUT / "uncertainty_budget.json").read_text())
    mor = json.loads((OUT / "morphometry_summary.json").read_text())
    roof = json.loads((OUT / "roof_summary.json").read_text())
    pla = json.loads((OUT / "planetary_summary.json").read_text())
    vert = json.loads((OUT / "vertical_check.json").read_text())

    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    clean = [r for r in rows if r["class"] == "intact"]
    s_clean = [float(r["s"]) for r in clean]

    n = dict(
        # --- registration ---
        reg_mode=reg["landmark"]["mode"],
        reg_matches=reg["landmark"]["matches"],
        reg_rim_rms_m=reg["landmark"]["rms"],
        reg_vres_mm=[round(1000 * v) for v in
                     reg["landmark"]["vertical_residuals_m"]],
        reg_sigma_tau_m=bud["sigma_tau_m"],
        reg_sigma_reg_z_m=bud["sigma_reg_z_m"],
        vertical_q90_diffs_m=[e["q90_diff"] for e in vert["per_skylight"]],
        scale_ruler_diff_m=[p["diff_m"] for p in reg["scale_ruler"]],
        # --- survey scope (single long flight processed on the Mac) ---
        centreline_length_m=mor["total_centreline_m"],
        n_stations=mor["n_stations"],
        n_valid_sections=mor["n_valid_sections"],
        # --- morphometry ---
        area=mor["area_m2"], height=mor["height_m"], eta=mor["eta"],
        width=mor.get("width_m"), sinuosity=mor["sinuosity_50m"],
        # --- roof thickness (clean stations only; see multipass caveat) ---
        tau=roof["tau_m"],
        n_tau_clean=roof["n_intact"],
        n_multipass=roof["n_multipass"],
        n_skylight_stations=roof["n_skylight"],
        n_no_dem=roof["n_no_dem"],
        tau_s_range_m=[min(s_clean), max(s_clean)] if s_clean else None,
        # minimum tau/L retained for provenance only; the manuscript's
        # stability analysis is the per-station strength inversion in
        # review_stats.json (the single-station minimum was retired in v6)
        min_tau_over_L_retired=roof["kappa_env"],
        min_tau_over_L_station=roof["kappa_env_station"],
        shielding=roof["shielding_g_cm2"],
        # --- planetary ---
        span_scale=pla["span_scale"],
        mars=pla["mars"], moon=pla["moon"],
        # --- provenance ---
        note=(
            ("Tube map = offline FAST-LIO2 first_long_flight carried into the "
             "Llewellin slice-registered frame by map-to-map GICP "
             "(transform_flf_ed_chain.json); onboard DLIO map processed "
             "identically as cross-check (analysis_out_ed, basis_comparison)."
             if OUT.name.endswith("_v6") else
             "Tube map = DLIO first_long_flight in the Llewellin slice-based "
             "registration (transform_ed_slice.json); cross-checked against the "
             "FAST-LIO map chained into the same frame (ceilings agree to "
             "0.1-0.2 m, analysis/chain_flf_to_ed.py).")
            + " Surface DEM = full_surface_and_subsurface_merged.pcd. tau uses "
            "intact stations (tau>0, non-skylight, non-multipass, "
            "non-inconsistent)."
            if reg["landmark"]["mode"] == "slice_manual_ed" else
            "Tube map = offline FAST-LIO first_long_flight (flf_*_aerial), "
            "gravity-preserving 4dof skylight registration. Replaces the Mac "
            "DLIO map, which carried ~60 m of vertical drift that inflated the "
            "deep roof (see analysis/dlio_vs_fastlio_drift.py). Surface DEM = "
            "full_surface_and_subsurface_merged.pcd. tau uses intact stations "
            "(tau>0, non-skylight, non-multipass, non-inconsistent)."),
    )
    (OUT / "paper_numbers.json").write_text(json.dumps(n, indent=2))
    print(json.dumps(n, indent=2))


if __name__ == "__main__":
    main()
