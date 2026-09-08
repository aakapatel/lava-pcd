#!/usr/bin/env bash
# v9 datum-transfer rerun (plan 11, work package D): the whole pipeline on the
# orthometric (ISH2004) surface and the translated ed-frame maps, in the v8
# order. Prerequisites: convert_isn16_surface.py, datum_transfer_check.py,
# seed_v9_ortho.py have run (maps/*_ed_ortho.pcd, maps/aerial_isn16_ortho_*.pcd,
# analysis_out_v9 seeded). Run from the repo root:
#     bash analysis/run_v9_pipeline.sh 2>&1 | tee -a analysis_out_v9/pipeline_v9.log
set -euo pipefail
cd "$(dirname "$0")/.."
export ANALYSIS_OUT=analysis_out_v9
export ROOF_DEM_PCD=maps/aerial_isn16_ortho_dem.pcd     # also names the DEM cache
PY=".venv/bin/python"
run() { echo; echo "=== $(date '+%H:%M:%S') $*"; env -u PYTHONPATH PYTHONPATH=src "$@"; }

# SKIP_MORPHO=1 reruns from the roof stage (the centreline/morphometry cache is
# untouched by roof-class changes; used for the 2026-09-08 class-precedence rerun)
if [ "${SKIP_MORPHO:-0}" != "1" ]; then
run env MORPHO_SKEL_PCD=maps/flf_30cm_ed_ortho.pcd MORPHO_SECT_PCD=maps/flf_10cm_ed_ortho.pcd \
    $PY analysis/run_morphometry.py
fi
run env ROOF_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd $PY analysis/run_roof.py
run env VALID_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd VALID_AERIAL_PCD=maps/aerial_isn16_ortho_10cm.pcd \
    $PY analysis/validate_ed_registration.py
run env ROOF_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd $PY analysis/run_roof.py
run $PY analysis/run_planetary.py
# chain scripts: GICP between SLAM-frame maps is datum-free, so the committed
# chain transforms and ed maps are reused (translated by seed_v9_ortho.py);
# only the flight-2 repeatability is recomputed on the v9 centreline.
run env SLF_CHAIN_REUSE=1 SLF_F1_PCD=maps/flf_30cm_ed_ortho.pcd SLF_F2_PCD=maps/slf_30cm_ed_ortho.pcd \
    $PY analysis/chain_slf_to_ed.py
run $PY analysis/compare_bases.py
run $PY analysis/orbital_dem_test.py
run $PY analysis/national_dem_check.py
run $PY analysis/run_review_stats.py
run $PY analysis/collect_paper_numbers.py
run $PY analysis/geotiff_sensitivity.py
run $PY analysis/regression_v9_vs_v6.py
echo; echo "=== $(date '+%H:%M:%S') pipeline done"
