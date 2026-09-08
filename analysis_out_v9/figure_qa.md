# Figure QA, v9 orthometric regeneration (2026-09-08; round 3 after the class-precedence rerun)

Method (HANDOFF_VERIFICATION rule 0): every panel in
`analysis_out_v9/figures_for_manuscript/` rasterised at 110 dpi
(`pdftoppm -r 110`, PNGs read directly) into `analysis_out_v9/figure_qa/`,
inspected, fixed, regenerated with `analysis/make_v9_figures.sh`, and
re-inspected. The v8.1 manuscript copies were rasterised alongside for the
two panels where a defect might have been pre-existing. Nothing was copied
into the manuscript figures directory (timestamps there are still 14 Aug).

| panel | file | round 1 | fix | round 2 verdict |
|---|---|---|---|---|
| Fig 5 morphometry | morphometry_panel.pdf | panel g label "centreline elevation (m a.s.l., ISH2004)" clipped at the left edge and touching panel f's label | two-line label (`ELEV_LABEL_2L`) | PASS. Panel g axis 160-167.5 m a.s.l. Round 3: S1 band continuous again (s 29-36) after `skylight` took precedence over `no_dem` in run_roof |
| Fig 6 roof | roof_panel.pdf | no overlap, no clipping; same layout as v8.1 | none | PASS. Round 3: bands 29-36 / 54-56 / 89-99 as in v8.1; no grey ticks at S1 any more (the four stations are aperture-class) |
| Fig 7 planetary | planetary_panel.pdf | unchanged from v8.1 apart from orbital-test decimals | none | PASS |
| Fig 4 showcase | merged_showcase_panel.pdf/.png | (i) "roof, coloured by tau" clipped at the right axis edge (present in the shipped v8.1 as well); (ii) "mapped, excluded by screen" over the S1/S2 markers; (iii) panel d y-label touching panel c | (i) right-aligned inside the axes; (ii) anchor +10.5 m above the DEM as in v8.1; (iii) two-line labels; annotation anchors now computed from the profile (were hardcoded 229.5/242.5/248.5/238.5/246.5 m ellipsoidal) | PASS. Panels c and d axes 160-185 m a.s.l.; panel d tau = 14.1 m (CSV 14.140); S1 marker now visible (v8.1 anchored it on the in-aperture junk DEM value). Panel b: the three apertures render white (no surface DEM inside the openings) where v8.1 showed pits filled from the merged cloud's tube points; hillshade and road/pit geometry otherwise identical to v8.1 |
| ED Fig 5 registration | ed_registration_panel.pdf | none | none | PASS. Panel c axis 165-185 m a.s.l.; pink 4-DOF ceiling (shifted by dz) above the DEM in the skylight reach, blue slice-registered ceiling below everywhere; panel b median 0.87 / RMS 1.35 m |
| ED Fig 6 consistency | ed_consistency_panel.pdf | none | none | PASS on layout. Panel a axis 165-180 m a.s.l. CONTENT CHANGE for the caption: panel b shows ZERO stations above the 5% line (every plotted station sits at 0.000). v8.1's 'few points above it are skylight-class' were s 31-34 measured against the fake in-aperture DEM value of the merged cloud; under the surface-only DEM those four aperture stations carry no column statistic (ghost_frac = NaN, not plotted). Panel a: the DEM line has a 4 m gap at S1 (no surface cell inside the opening) |
| ED Fig 7 centreline | ed_centreline_panel.pdf | panel b label long but not overlapping | two-line label | PASS. Panel b axis 150-180 m a.s.l.; round 3 skylight stations 22 as v8.1 |
| ED Fig 11 DEM validation | dem_validation_panel.pdf | colourbar label "dz (m)" touching panel b's "samples" label (pre-existing in v8.1) | gridspec wspace 0.30 -> 0.40 | PASS. Road stripe runs N-S near x 1250-1300 as expected; panel b now reads std 0.28 m, tilt 0.43 m/km (was 0.29 / 0.42) |
| skeleton_3d | skeleton_3d.png | raster, no text | none | PASS (same view; a vertical translation does not change it) |

Elevation axes now read "elevation (m a.s.l., ISH2004)" from
`analysis_out_v9/datum.json` (legacy runs without that file label their axes
"elevation (m, WGS84 ellipsoidal)"); ranges are 150-185 m throughout.

Not verified here: the built PDF (.fls) embedding, because the manuscript
copy step is the orchestrator's after review.

Round 3 (2026-09-08, after run_roof gave `skylight` precedence over `no_dem`
and writes ghost_frac = NaN where there is no surface cell): all nine
panels regenerated with `make_v9_figures.sh`, rasterised and re-read; every
verdict above holds, and the S1 band/marker artefacts of rounds 1-2 are gone.
