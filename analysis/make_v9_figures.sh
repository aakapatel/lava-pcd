#!/usr/bin/env bash
# v9 datum-transfer rerun: regenerate every data figure from analysis_out_v9
# on the orthometric maps, into analysis_out_v9/figures_for_manuscript/ (NOT
# the manuscript figures dir; the orchestrator copies after review).
# Run from the repo root after run_v9_pipeline.sh:
#     bash analysis/make_v9_figures.sh 2>&1 | tee -a analysis_out_v9/figures_v9.log
set -euo pipefail
cd "$(dirname "$0")/.."
export ANALYSIS_OUT=analysis_out_v9
export ROOF_DEM_PCD=maps/aerial_isn16_ortho_dem.pcd
export FIGS_DIR=analysis_out_v9/figures_for_manuscript
mkdir -p "$FIGS_DIR" analysis_out_v9/renders
PY=".venv/bin/python"
run() { echo; echo "=== $(date '+%H:%M:%S') $*"; env -u PYTHONPATH PYTHONPATH=src "$@"; }
render() { echo; echo "=== $(date '+%H:%M:%S') (render) $*"; env -u PYTHONPATH -u WAYLAND_DISPLAY DISPLAY=:1 PYTHONPATH=src "$@"; }

run $PY analysis/make_figures.py
run env EDFIG_SHELL_PCD=maps/flf_30cm_ed_ortho.pcd EDFIG_SECT_PCD=maps/flf_10cm_ed_ortho.pcd \
    EDFIG_AERIAL_PCD=maps/aerial_isn16_ortho_10cm.pcd $PY analysis/make_ed_figures.py
run $PY analysis/make_dem_validation_figure.py
render env SHOWCASE_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd \
    SHOWCASE_SURFACE_PCD=maps/aerial_isn16_ortho_surface_10cm.pcd $PY analysis/render_merged_showcase.py
run env SHOWCASE_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd \
    SHOWCASE_SURFACE_PCD=maps/aerial_isn16_ortho_surface_10cm.pcd $PY analysis/make_merged_showcase_figure.py
render env SKEL3D_TUBE_PCD=maps/flf_10cm_ed_ortho.pcd $PY analysis/render_skeleton_3d.py
echo; echo "=== $(date '+%H:%M:%S') figures done -> $FIGS_DIR"; ls -la "$FIGS_DIR"
