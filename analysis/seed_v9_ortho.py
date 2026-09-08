"""v9 datum transfer, step D.4: apply the measured vertical datum shift and
seed analysis_out_v9 for the orthometric pipeline run.

Reads ANALYSIS_OUT/datum_transfer.json (written by datum_transfer_check.py)
and applies the plan's rule: the old-to-new surface transform must be a pure
vertical translation (nearest-neighbour dz std < 0.05 m, |dx|, |dy| < 0.5 m).
If it is not, STOP.md is written and the script exits with status 2 without
touching any map.

If it is, dz = median nearest-neighbour (new minus old) is composed AFTER the
existing registration chain (Ed's slice matrix, the GICP chains): every
ed-frame tube map gets a copy translated by (0, 0, dz):

    maps/{flf,slf,tube}_{10,30}cm_ed.pcd  ->  maps/*_ed_ortho.pcd

(# LAVA_PCD_ORIGIN header preserved by lava_pcd.merge.apply_transform).
A 10 cm orthometric surface working copy over the v6 DEM extent is cut from
maps/aerial_isn16_ortho_10cm.pcd for the showcase render.

ANALYSIS_OUT is then seeded the way setup_ed_frame.py seeds an Ed-frame run,
but from analysis_out_v6 and with every absolute-z field shifted by dz:
    aerial_holes.json        centroid z (+dz)          (xy, ellipse unchanged)
    registration_report.json aerial_holes centroid z (+dz); vertical residuals
                             are differences, unchanged
    vertical_check.json      q*_aerial and q*_tube (+dz); q*_diff unchanged
    tube_holes.json, planetary_catalogue.csv, uncertainty_budget.json  copied
No DEM cache is copied: run_roof rebuilds it from the new surface.

Also written: ANALYSIS_OUT/ortho_translation.json (the 4x4 applied, dz, the
statistics it came from) and ANALYSIS_OUT/datum.json (the vertical-datum
label the figure scripts read for their elevation axes).

Run:  ANALYSIS_OUT=analysis_out_v9 env -u PYTHONPATH PYTHONPATH=src \
      .venv/bin/python analysis/seed_v9_ortho.py
"""
from __future__ import annotations

import json
import os
import shutil
import sys
from pathlib import Path

import numpy as np

from lava_pcd.crop import crop_pcd
from lava_pcd.merge import apply_transform

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_v9")))
REF = ROOT / "analysis_out_v6"
MAPS = ROOT / "maps"
STD_MAX = 0.05
SHIFT_MAX = 0.5
BBOX_DEM = (1200.7193603515625, 1520.7193603515625,
            500.17303466796875, 881.6730346679688)


def main() -> None:
    dt = json.loads((OUT / "datum_transfer.json").read_text())
    nn = dt["nearest_neighbour"]["dz_new_minus_old"]
    dz = float(nn["median"])
    std = float(nn["std"])
    shifts = dt["horizontal_shift"]
    dxdy = []
    for k, v in shifts.items():
        dxdy.append((k, "phase", v["dx_m"], v["dy_m"]))
        dxdy.append((k, "brute", *v["brute_force_best_shift_m"]))
    worst = max(max(abs(a), abs(b)) for _, _, a, b in dxdy)
    checks = dict(dz_median_m=dz, dz_std_m=std, std_limit_m=STD_MAX,
                  horizontal_shifts=[dict(grid=k, method=m, dx_m=a, dy_m=b)
                                     for k, m, a, b in dxdy],
                  worst_abs_horizontal_m=worst, shift_limit_m=SHIFT_MAX)
    passed = std < STD_MAX and worst < SHIFT_MAX
    if not passed:
        msg = ("# STOP: the new surface cloud is not a pure vertical translation "
               "of the old one within tolerance\n\n```\n"
               + json.dumps(checks, indent=2) + "\n```\n\nPer plan 11 section 5 "
               "the fallback is to keep the old cloud as the analysis basis and "
               "apply only the measured mean dz to report orthometric heights.\n")
        (OUT / "STOP.md").write_text(msg)
        print(msg)
        sys.exit(2)
    print("rule passed:", json.dumps(checks, indent=2))

    T = np.eye(4)
    T[2, 3] = dz
    written = []
    for base in ("flf", "slf", "tube"):
        for res in ("10cm", "30cm"):
            src = MAPS / f"{base}_{res}_ed.pcd"
            dst = MAPS / f"{base}_{res}_ed_ortho.pcd"
            apply_transform(src, dst, T, show_progress=False)
            written.append(dst.name)
            print("wrote", dst.name)
    surf = MAPS / "aerial_isn16_ortho_surface_10cm.pcd"
    r = crop_pcd(MAPS / "aerial_isn16_ortho_10cm.pcd", surf, BBOX_DEM,
                 show_progress=False)
    print(f"wrote {surf.name}: {r.point_count:,} of {r.source_count:,}")

    OUT.mkdir(exist_ok=True)
    for name in ("tube_holes.json", "planetary_catalogue.csv",
                 "uncertainty_budget.json"):
        shutil.copy2(REF / name, OUT / name)

    ah = json.loads((REF / "aerial_holes.json").read_text())
    for h in ah["holes"]:
        h["centroid"][2] = h["centroid"][2] + dz
    ah["source_path"] = str(MAPS / "aerial_isn16_ortho_crop.pcd")
    ah["note"] = (f"v6 detections (aerial_crop.pcd) with centroid z shifted by "
                  f"{dz:+.3f} m to the ISH2004 orthometric datum; xy and "
                  "ellipse parameters unchanged")
    (OUT / "aerial_holes.json").write_text(json.dumps(ah, indent=2))

    rep = json.loads((REF / "registration_report.json").read_text())
    for h in rep["aerial_holes"]:
        h["centroid"][2] = round(h["centroid"][2] + dz, 3)
    rep["inputs"]["aerial"] = "aerial_isn16_ortho_crop.pcd (v6 detections, z shifted)"
    rep["note"] = (rep.get("note", "") + f" v9: aerial hole z shifted by {dz:+.3f} m "
                   "(orthometric datum); vertical residuals are differences and unchanged.")
    (OUT / "registration_report.json").write_text(json.dumps(rep, indent=2))

    vc = json.loads((REF / "vertical_check.json").read_text())
    for e in vc["per_skylight"]:
        for q in (50, 75, 90, 95):
            e[f"q{q}_aerial"] = round(e[f"q{q}_aerial"] + dz, 3)
            e[f"q{q}_tube"] = round(e[f"q{q}_tube"] + dz, 3)
    vc["note"] = (vc.get("note", "") + f" v9: absolute quantiles shifted by {dz:+.3f} m "
                  "(orthometric datum); differences unchanged.")
    (OUT / "vertical_check.json").write_text(json.dumps(vc, indent=2))

    (OUT / "ortho_translation.json").write_text(json.dumps(dict(
        matrix=[list(map(float, r)) for r in T],
        dz_m=dz, dz_source="median of nearest-neighbour (new minus old) vertical "
        "differences over the tube footprint, datum_transfer.json",
        composition="applied AFTER the registration chain: "
        "T_ortho<-slam = Tz(dz) @ T_aerial<-slam (transform_ed_slice / "
        "transform_flf_ed_chain / transform_slf_ed_chain)",
        maps_written=written, checks=checks), indent=2))
    (OUT / "datum.json").write_text(json.dumps(dict(
        vertical_datum="ISH2004 geoid (orthometric), ISN2016 delivery of 2026-09-08",
        horizontal="EPSG:32627 minus (479158, 7089826), unchanged local frame",
        elevation_axis_label="elevation (m a.s.l., ISH2004)",
        legacy_frame="WGS84 ellipsoidal heights (v6/v8 canon)",
        dz_applied_m=dz), indent=2))
    print("seeded", OUT)


if __name__ == "__main__":
    main()
