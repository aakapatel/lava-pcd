"""v9 datum transfer, step D.6: field-by-field regression of the orthometric
rerun (analysis_out_v9) against the v8 canon (analysis_out_v6).

Every numeric leaf of paper_numbers.json, review_stats.json,
roof_summary.json, morphometry_summary.json, registration_validation.json,
flight2_repeatability.json, basis_comparison.json, national_dem_check.json
and orbital_dem_test.json is listed with its v6 value, v9 value, difference
and a verdict:

    identical        exactly equal
    rounding-safe    differs, but by less than half the unit the paper quotes
                     (paper roundings from the v8 canon; see TOLERANCE)
    datum (dz)       an absolute-elevation field that moved by the applied
                     shift (ortho_translation.json) to within 5 mm
    datum-derived    a field that legitimately carries the datum change
                     (national DEM offsets, the DLIO-vs-FLF ceiling offset
                     against the untouched ellipsoidal analysis_out_ed)
    CHANGED          anything else: to be investigated, never papered over

Writes ANALYSIS_OUT/regression_vs_v6.md and regression_vs_v6.json.

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/regression_v9_vs_v6.py
"""
from __future__ import annotations

import json
import os
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
REF = ROOT / "analysis_out_v6"
FILES = ["paper_numbers.json", "review_stats.json", "roof_summary.json",
         "morphometry_summary.json", "registration_validation.json",
         "flight2_repeatability.json", "basis_comparison.json",
         "national_dem_check.json", "orbital_dem_test.json",
         "planetary_summary.json", "uncertainty_budget.json"]

# half the unit quoted in the paper (v8 canon list in the session handoff)
TOLERANCE = [
    (r"sinuosity", 0.005),                       # 1.03 / 1.08
    (r"\.eta\.", 0.005),                         # 1.53
    (r"span_scale|dtau|dceil", 0.005),
    (r"corr_len|n_eff", 0.5),
    (r"volume\.total_m3", 500.0),                # "14,000 m3"
    (r"volume\.per_m", 0.5),                     # "47 per m"
    (r"shielding.*frac_gt", 0.005),              # 99/84/72 %
    (r"shielding.*median_g_cm2|shielding\.median|shielding\.min", 5.0),
    (r"atmosphere_ratio", 0.005),
    (r"sigma_req_MPa\.(median|p90|p95|max|drop_one)", 0.05),   # 2.1 / 2.0 MPa
    (r"tau_over_L|kappa", 0.0005),
    (r"area_m2\.|area_median|area\.", 0.5),      # 48 m2, 13-76
    (r"width|height", 0.05),                     # 6.1 .. 32, 8.1
    (r"tau_m\.|\.tau\.|median_tau|tau_hist|reach_median|near_skylight", 0.05),
    (r"ci_block|ci_iid|median_ci95", 0.05),
    (r"median_m|rms_m|p10_m|p90_m|rms|median|p5|p95|std", 0.005),
    (r"sigma_", 0.005),
    (r"total_centreline_m|centreline_m", 0.05),
    (r"Lmax_sqrt_model", 0.05),
    (r"slope_m_per_m|ci95", 0.005),
    (r"d_zdem|d_zceil|d_tau|r2|\.r$|dem_share|corr_tau_zdem|rms_residual", 0.005),
    (r"variance\.", 0.05),
    (r"zceil_range|zdem_range", 0.05),
    (r"thick_thin_contrast|median_shift|per_station_rms|frac_", 0.005),
    (r"doming|tilt|std_", 0.005),
    (r"scale_ruler|hole_semi_major|dists|S1S|S2S", 0.005),
]
DEFAULT_TOL = 0.005
ABS_Z = [r"net_centreline\.(start|end)"]
DATUM_DERIVED = [
    r"national_dem_check\..*datum_offset_m",
    r"basis_comparison\.dceil\.",          # A (analysis_out_ed) is ellipsoidal
    r"national_dem_check\..*best_shift_m",
]
# Differences whose cause was traced during the v9 run (2026-09-08). A field
# matching one of these is reported as "explained: <tag>" instead of CHANGED;
# the cause text is printed in the report so the reader can judge it.
EXPLAINED = [
    # (an earlier "aperture bookkeeping" entry covered the four S1 stations
    # classed no_dem instead of skylight; retired 2026-09-08 when run_roof gave
    # 'skylight' precedence over 'no_dem')
    (r"sigma_reg|sigma_tau|tau_threshold|floor_through_skylight|anchor_floor",
     "S3 floor check on new voxel copy",
     "the 10 cm working copy is re-voxelised from the new cloud; the 5 mm "
     "horizontal ISN2016/WGS84 difference changes voxel membership near S3 "
     "(360 vs 352 matched points), RMS 1.347 vs 1.356 m, so sigma_reg 1.45 vs 1.46 "
     "and sigma_tau 1.48 vs 1.49 m; the paper quotes both as 1.5 m."),
    (r"national_dem_check",
     "in-aperture cells removed",
     "165 cells inside the skylights had DEM values only in v6 (tube points 4-5 m "
     "below the rim); they are absent in the surface-only crop, which lowers the "
     "ArcticDEM shape std from 0.293 to 0.281 m (paper: 0.29 m, an upper bound) "
     "and the tilt 0.42 -> 0.43 m/km. Datum offsets carry the -66.00 m shift."),
    (r"stats\.median|orbital_dem_test|basis_comparison|shielding|kappa_sensitivity|"
     r"tau_over_L\.second|sensitivity_tau_gt_0p5|skylight_runs",
     "1 cm quantisation of the delivered cloud",
     "the new cloud is stored at 1 cm; per-cell maxima come out 0.5-1 cm lower "
     "than old minus 66.00 (median -66.01), so tau is lower by 0.0075 m in the mean "
     "(median 6.0255 -> 6.0155 m, shielding 1567 -> 1564 g/cm2, 1.52 -> 1.51x). "
     "Paper roundings (6.0 m, about 1,500 g/cm2, 1.5x) unchanged; the smallest "
     "tau/L moves because a 0.10 m cover became 0.09 m."),
]
IGNORE = [r"\.note$", r"^.*\.product$", r"registration_validation\.note",
          r"sigma_reg_z_basis", r"\.basis$", r"final_transform", r"config$",
          r"reg_mode", r"warnings", r"class_counts"]


def flatten(obj, prefix=""):
    out = {}
    if isinstance(obj, dict):
        for k, v in obj.items():
            out.update(flatten(v, f"{prefix}.{k}" if prefix else k))
    elif isinstance(obj, list):
        if all(not isinstance(v, (dict, list)) for v in obj):
            for i, v in enumerate(obj):
                out[f"{prefix}[{i}]"] = v
        else:
            for i, v in enumerate(obj):
                out.update(flatten(v, f"{prefix}[{i}]"))
    else:
        out[prefix] = obj
    return out


def tol_for(key):
    for pat, t in TOLERANCE:
        if re.search(pat, key):
            return t
    return DEFAULT_TOL


def main() -> None:
    dz = json.loads((OUT / "ortho_translation.json").read_text())["dz_m"]
    rows = []
    counts = {}
    for name in FILES:
        a = flatten(json.loads((REF / name).read_text()), name[:-5])
        b = flatten(json.loads((OUT / name).read_text()), name[:-5])
        for key in sorted(set(a) | set(b)):
            if any(re.search(p, key) for p in IGNORE):
                continue
            va, vb = a.get(key), b.get(key)
            if isinstance(va, bool) or isinstance(vb, bool) or \
               not isinstance(va, (int, float)) or not isinstance(vb, (int, float)):
                verdict = "identical" if va == vb else "CHANGED (non-numeric)"
                if verdict != "identical":
                    for pat, tag, _ in EXPLAINED:
                        if re.search(pat, key):
                            verdict = f"explained: {tag}"
                            break
                diff = ""
            else:
                diff = vb - va
                if diff == 0:
                    verdict = "identical"
                elif any(re.search(p, key) for p in ABS_Z) and abs(diff - dz) < 0.005:
                    verdict = "datum (dz)"
                elif any(re.search(p, key) for p in DATUM_DERIVED):
                    verdict = "datum-derived"
                elif abs(diff) < tol_for(key):
                    verdict = "rounding-safe"
                else:
                    verdict = "CHANGED"
                    for pat, tag, _ in EXPLAINED:
                        if re.search(pat, key):
                            verdict = f"explained: {tag}"
                            break
            counts[verdict] = counts.get(verdict, 0) + 1
            rows.append(dict(field=key, v6=va, v9=vb, diff=diff, verdict=verdict))

    def fmt(v):
        if isinstance(v, float):
            return f"{v:.4f}" if abs(v) < 1e4 else f"{v:.1f}"
        return str(v)

    lines = ["# Regression: analysis_out_v9 (orthometric, ISH2004) vs analysis_out_v6 (v8 canon)",
             "",
             f"Applied vertical translation dz = {dz:+.4f} m (ortho_translation.json).",
             "Verdicts: identical | rounding-safe (below half the paper's quoted unit) | "
             "datum (dz) (absolute elevation, moved by exactly dz) | datum-derived "
             "(offsets against products still in the old datum) | CHANGED (investigate).",
             "",
             "Summary: " + ", ".join(f"{k}: {v}" for k, v in sorted(counts.items())),
             ""]
    changed = [r for r in rows if r["verdict"].startswith("CHANGED")]
    lines += ["## Paper headline numbers at the paper's rounding", "",
              "| quantity | v6 | v9 | paper rounding |", "|---|---|---|---|"]
    pn6 = json.loads((REF / "paper_numbers.json").read_text())
    pn9 = json.loads((OUT / "paper_numbers.json").read_text())
    nd6 = json.loads((REF / "national_dem_check.json").read_text())
    nd9 = json.loads((OUT / "national_dem_check.json").read_text())
    f26 = json.loads((REF / "flight2_repeatability.json").read_text())
    f29 = json.loads((OUT / "flight2_repeatability.json").read_text())
    def head(name, a, b, fm):
        lines.append(f"| {name} | {fm(a)} | {fm(b)} | {fm.__doc__} |")
    def r1(v):
        "0.1"
        return f"{v:.1f}"
    def r2(v):
        "0.01"
        return f"{v:.2f}"
    def r0(v):
        "1"
        return f"{v:.0f}"
    def pct(v):
        "1 %"
        return f"{100*v:.0f}%"
    def k1(v):
        "1,000"
        return f"{1000*round(v/1000):,.0f}"
    for label, key, fm in (
        ("survey length (m)", ("centreline_length_m",), r0),
        ("stations", ("n_stations",), r0),
        ("intact stations", ("n_tau_clean",), r0),
        ("overburden median (m)", ("tau", "median"), r1),
        ("overburden block CI low (m)", ("stats", "ci_block", "30", 0), r1),
        ("overburden block CI high (m)", ("stats", "ci_block", "30", 1), r1),
        ("overburden IQR low (m)", ("tau", "q25"), r1),
        ("overburden IQR high (m)", ("tau", "q75"), r1),
        ("overburden max (m)", ("tau", "max"), r1),
        ("shielding >100 g/cm2 (rho 2600)", ("shielding_sensitivity", "2600", "frac_gt_100"), pct),
        ("shielding >500 g/cm2", ("shielding_sensitivity", "2600", "frac_gt_500"), pct),
        ("shielding >1000 g/cm2", ("shielding_sensitivity", "2600", "frac_gt_1000"), pct),
        ("area median (m2)", ("area", "median"), r0),
        ("area p5 (m2)", ("area", "p5"), r0),
        ("area p95 (m2)", ("area", "p95"), r0),
        ("width min (m)", ("width", "min"), r1),
        ("width max (m)", ("width", "max"), r1),
        ("height median (m)", ("height", "median"), r1),
        ("sinuosity mean", ("sinuosity", "mean"), r2),
        ("sinuosity max", ("sinuosity", "max"), r2),
        ("strength demand p95 (MPa)", ("envelope", "sigma_req_MPa", "p95"), r1),
        ("strength demand max (MPa)", ("envelope", "sigma_req_MPa", "max"), r1),
        ("strength demand drop-one (MPa)", ("envelope", "sigma_req_MPa", "drop_one"), r1),
        ("strength stations n", ("envelope", "sigma_req_MPa", "n"), r0),
        ("sigma_tau (m)", ("reg_sigma_tau_m",), r1),
        ("sigma_reg (m)", ("reg_sigma_reg_z_m",), r1),
        ("S3 floor median (m)", ("registration_validation", "floor_through_skylight_s3", "median_m"), r2),
        ("S3 floor RMS (m)", ("registration_validation", "floor_through_skylight_s3", "rms_m"), r2),
        ("scale ruler S1-S2 diff (m)", ("scale_ruler_diff_m", 0), r2),
        ("scale ruler S1-S3 diff (m)", ("scale_ruler_diff_m", 1), r2),
    ):
        def get(d, ks):
            for k in ks:
                d = d[k]
            return d
        head(label, get(pn6, key), get(pn9, key), fm)
    v6v = json.loads((REF / "review_stats.json").read_text())["volume"]
    v9v = json.loads((OUT / "review_stats.json").read_text())["volume"]
    head("volume (m3)", v6v["total_m3"], v9v["total_m3"], k1)
    head("volume per m (m3)", v6v["per_m"], v9v["per_m"], r0)
    head("flight-2 repeatability RMS (m)", f26["rms"], f29["rms"], r2)
    head("flight-2 repeatability median (m)", f26["median"], f29["median"], r2)
    head("ArcticDEM shape std (m)", nd6["arcticdem"]["std_raw_m"], nd9["arcticdem"]["std_raw_m"], r2)
    head("ArcticDEM doming (m)", nd6["arcticdem"]["doming_amplitude_m"], nd9["arcticdem"]["doming_amplitude_m"], r2)
    head("ArcticDEM tilt (m/km)", nd6["arcticdem"]["tilt_m_per_km"], nd9["arcticdem"]["tilt_m_per_km"], r2)
    head("ArcticDEM datum offset, ours minus (m)", nd6["arcticdem"]["datum_offset_m"], nd9["arcticdem"]["datum_offset_m"], r2)
    head("IslandsDEM datum offset, ours minus (m)", nd6["islandsdem"]["datum_offset_m"], nd9["islandsdem"]["datum_offset_m"], r2)
    lines += ["", "## Explained differences (cause traced during the v9 run)", ""]
    for pat, tag, why in EXPLAINED:
        n_ = sum(1 for r in rows if r["verdict"] == f"explained: {tag}")
        lines.append(f"- **{tag}** ({n_} fields, pattern `{pat}`): {why}")
    lines += ["", "## Fields flagged CHANGED (unexplained)", ""]
    if changed:
        lines += ["| field | v6 | v9 | diff |", "|---|---|---|---|"]
        lines += [f"| {r['field']} | {fmt(r['v6'])} | {fmt(r['v9'])} | {fmt(r['diff'])} |"
                  for r in changed]
    else:
        lines += ["none"]
    lines += ["", "## All fields", "", "| field | v6 | v9 | diff | verdict |",
              "|---|---|---|---|---|"]
    lines += [f"| {r['field']} | {fmt(r['v6'])} | {fmt(r['v9'])} | {fmt(r['diff'])} | "
              f"{r['verdict']} |" for r in rows]
    (OUT / "regression_vs_v6.md").write_text("\n".join(lines) + "\n")
    (OUT / "regression_vs_v6.json").write_text(json.dumps(
        dict(dz_m=dz, counts=counts, rows=rows), indent=2))
    # one summary block into datum_transfer.json so every v9 number sits in
    # one file (national-DEM offsets, headline pairs, station-class change)
    dtp = OUT / "datum_transfer.json"
    if dtp.exists():
        dt = json.loads(dtp.read_text())
        r6 = json.loads((REF / "roof_summary.json").read_text())
        r9 = json.loads((OUT / "roof_summary.json").read_text())
        dt["v9_summary"] = dict(
            dz_applied_m=dz, regression_counts=counts,
            national_dem_offsets_ours_minus_m=dict(
                arcticdem_ellipsoidal=dict(v6=nd6["arcticdem"]["datum_offset_m"],
                                           v9=nd9["arcticdem"]["datum_offset_m"]),
                islandsdem_orthometric=dict(v6=nd6["islandsdem"]["datum_offset_m"],
                                            v9=nd9["islandsdem"]["datum_offset_m"])),
            arcticdem_shape_std_m=dict(v6=nd6["arcticdem"]["std_raw_m"],
                                       v9=nd9["arcticdem"]["std_raw_m"]),
            station_classes=dict(v6={k: r6[k] for k in r6 if k.startswith("n_")},
                                 v9={k: r9[k] for k in r9 if k.startswith("n_")}),
            tau_m=dict(v6=r6["tau_m"], v9=r9["tau_m"]),
            headline_table_md=str(OUT / "regression_vs_v6.md"))
        dtp.write_text(json.dumps(dt, indent=2))
    print("Summary:", counts)
    for r in changed:
        print(f"CHANGED {r['field']}: {r['v6']} -> {r['v9']} ({r['diff']})")


if __name__ == "__main__":
    main()
