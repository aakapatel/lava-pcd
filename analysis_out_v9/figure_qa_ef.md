# Figure QA, work packages E and F (2026-09-08)

Method (HANDOFF_VERIFICATION rule 0): each PDF rasterised at 110 dpi with
`pdftoppm -r 110` into `analysis_out_v9/figure_qa/`, read as an image, fixed,
regenerated, re-read; panels a and e of the ED panel were also inspected at
300 dpi crops. Scripts: `analysis/road_crossing.py`,
`analysis/survey1970_compare.py` (+ `--figure`).

## fig_road_crossing.pdf (task E; becomes a panel of ED Fig 12)

| round | defect | fix |
|---|---|---|
| 1 | DEM line stopped at s = 301 (only roof_thickness rows drawn); legend covered the floor dots; plan legend covered the ortho and hid the star | DEM joined from the cells beyond 301; legend moved into the empty 169-178 m band; plan legend outside the axes |
| 2 | legend star scaled by markerscale to a giant marker; ceiling line stopped at 280 | marker sizes set only for scatter handles; ceiling mask separated |
| 3 | third legend column ran into the "survey end" line | return counts removed from the legend (they live in road_crossing.json), labels shortened |
| 4 | after the ceiling-consistency flag, cross markers added for cells whose ceiling is not sampled | title text updated |

Verdict: PASS. No text over data or other text, nothing clipped, legends in
empty regions. Content check: DEM crest band 323-328 m; ceiling cells at
325 and 327.5 m carry 2838-4403 returns each; open circles absent (no cell
below 30 returns survived the trim); crosses at 317.5 (second flight,
pipeline ceiling 160.7 vs p99 167.5) and 320 m (all maps, shadowed cell).
Plan panel: the ortho road, the OSM line and the DEM crest coincide; the
far-range returns end at x = 1252 under the east edge of the crest.

## figures_for_manuscript/ed_survey1970_panel.pdf (task F)

| round | defect | fix |
|---|---|---|
| 1 | panel a title ran into panel b's; Fig. 3 raster cluttered by the section and cross-sections; panel c "tie points" label on the legend; panel e legend over the floor dots and "survey end" label on the legend | two-line titles; raster cropped to the plan area (Fig. 3 rows 470-960, cols 380-2470); label moved and boxed; panel e rebuilt as a copy of fig_road_crossing panel a |
| 2 | panel e legend top row on the DEM line (short panel) | y-limits 155-184, legend moved to the empty lower-left, two columns |
| 3 | panel e third column touched the "survey end" line | labels shortened ("ceiling, DLIO") |
| 4 | panel c: Fig. 4C curve labelled "SMCC 1970" was ambiguous against the map-sheet curve; "tie-point reach" label brushed the Munger curve | label "SMCC 1970 (Fig. 4C)"; label moved to s = 64 |

Alignment check (rule 1/2): at 300 dpi the warped Fig. 3 walls straddle the
lidar centreline over the whole shared reach, the "Roof Collapses" label sits
beside S1-S3 and the Tube 1 / Tube 2 loops lie west of the shell, so the
raster warp and the digitised medial (green dotted) agree with each other and
with the shell. Panel b: the map-sheet medial (tie-point fit) and the Fig. 4C
medial (shape fit) overlap the centreline; Munger 1954 and Cambridgeshire
1969 diverge as their RMS values say.

Verdict: PASS (all five panels).

Not done here: the built manuscript PDF (.fls) embedding, because the
manuscript is the orchestrator's; the figures are in `analysis_out_v9/`.
