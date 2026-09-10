# Figure QA, Nature main figures (plan/13 Step 4), 2026-09-10

Script: `analysis/make_nature_figures.py` (run from the repo root with
`env -u PYTHONPATH PYTHONPATH=src .venv/bin/python analysis/make_nature_figures.py`;
`... qa` alone re-runs the checks). Outputs in `analysis_out_v10/figures_nature/`
(PDF + 300 dpi PNG each); QA rasters (110 dpi, 300 dpi, deuteranopia
simulation of the 110 dpi image) in `analysis_out_v10/figure_qa_nature/`.
Data source for every number: `analysis_out_v9/` (paper_numbers.json,
review_stats.json, morphometry.csv, roof_thickness.csv, centreline.csv,
aerial_holes.json, planetary_catalogue.csv, orbital_dem_test.json,
osm_way_620673065.json, survey1970.json + survey1970_georef.npz).

## Global checks (pdfinfo, pdffonts, script output)

| Figure | Size (mm) | Fonts embedded | Verdict |
|---|---|---|---|
| fig1_anatomy_site | 180.0 x 110.0 | NimbusSans-Regular, -Bold, -Italic (Type 3) | OK |
| fig2_survey | 180.0 x 160.0 | NimbusSans-Regular, -Bold (Type 3) | OK |
| fig3_merged_model | 180.0 x 170.0 | NimbusSans-Regular, -Bold, -Italic (Type 3) | OK |
| fig4_geometry_cover | 180.0 x 170.0 | NimbusSans-Regular, -Bold, -Italic (Type 3) | OK |
| fig5_planetary | 180.0 x 118.0 | NimbusSans-Regular, -Bold, -Italic (Type 3) | OK |

No DejaVu or other font appears in any PDF. Type 3 rather than Type 42: the
installed Nimbus Sans is an OpenType CFF font, and matplotlib's fonttype 42
wraps CFF outlines in a TrueType wrapper that poppler flags ("Mismatch
between font type and embedded font file"); Type 3 embeds the subset
cleanly with a ToUnicode map (pdffonts: emb yes, uni yes). Text size 7 pt,
ticks 6 pt, panel letters bold 8 pt; labels start with a capital and carry a
space before the unit; the only tick labels above 999 (shielding axis of
Fig 4d) read 1,000 / 2,000 / 3,000. All embedded rasters are at 300 ppi at
final size (pdfimages -list), resampled from sources at 467 to 887 dpi.

## Colour-vision check

Deuteranopia (Machado 2009 matrix, severity 1) applied to each 110 dpi
raster and inspected:
- Fig 1: earth palette (browns, tan, greys) and the orange skylight rings on
  the ortho; all classes remain separable.
- Fig 2: rasters as supplied (height colouring baked into the onboard
  visualiser renders; see caveats).
- Fig 3: the overburden scale is a luminance ramp (bright thin, dark thick)
  and survives unchanged; the blue-grey "excluded" cells near S2 become
  close to the hillshade grey but are also keyed and cover only 4 stations.
- Fig 4: intact (green -> olive grey) against skylight (orange -> yellow),
  additionally separated by marker shape; width (blue) against height (green
  -> grey) separable.
- Fig 5: Earth (green -> grey), Mars (orange -> yellow), Moon (light blue)
  separable, and the ECDFs in a also differ by line style.
No red/green pair and no rainbow scale in any vector panel.

## Panel verdicts (110 dpi and 300 dpi rasters read; rule 0)

### Fig 1, anatomy and site
- a (long section + cross-section, own vector art after lt_fix.jpg): OK.
  Labels: Surface; Soil, tephra; Basalt roof (stacked flow units); Drained
  conduit; Snow cone; Breakdown blocks; Basalt; Skylight (the one feature
  seen from above and from within); What orbit sees: the aperture; What the
  robot sees: the void from inside; Collapse trench (entrance); Overburden
  tau; Span L; Highest ceiling. No label overlaps another or a leader.
- b (site map): OK. Landmark check against HANDOFF rules: the three dark
  collapse pits in the orthomosaic sit inside the S1, S2, S3 rings
  (aerial_holes.json centroids), the entrance collapse and visitor centre
  sit at the s = 0 station (the 1970 tie point "entrance line midpoint"),
  and the OSM Route 39 polyline runs along the road crest. Scale: S1 to S3
  straight-line separation measured on the 300 dpi raster is 56 m against
  57.0 m from the centroids and the 50 m bar. Ortho read through
  EPSG:32627 -> EPSG:8088 at 0.25 m; 50 m arc ticks 0 to 300; survey-end
  label "s = 301 m" from centreline.csv; 1970 map-sheet medial line drawn
  only beyond the survey end (dashed white). Iceland inset from the
  Natural Earth 50 m coastline (extract cached in
  analysis_out_v10/iceland_coastline_ne50m.json); site at 63.94 N, 21.40 W
  computed from the s = 0 station. Map labels carry a white halo because a
  map has no unshaded area to put them on.

### Fig 2, autonomous survey (rasters)
- a: top-left and bottom-left quadrants of mission_tpv_collage.png (chamber
  beneath a skylight; dark constriction). 467 dpi source, embedded at
  300 ppi. Vector captions below each photo. OK.
- b: quadrants A and D of exploration_progression.png with the baked
  letters blacked out (the map does not enter the masked corner in either
  stage). Vector captions "Early in the mission" / "End of the recorded
  mission". OK.
- c: bottom tier of skylights_frontiers.png (rows 2916 to 3998). The baked
  Times title "Global Path to Skylight Frontier" was masked (pure black
  background there) and replaced by a vector label "Route to the skylight
  frontier" in the same place; the baked arrow to the route end is kept
  and points at the red route. OK.
- Letters a, b, c present; bottom margin 2 mm.

### Fig 3, co-registered model
- a: showcase_oblique.png trimmed, 110 mm wide (887 dpi source). Tube runs
  parallel to the road as in v9. A single horizontal colourbar
  "Overburden (m)", 0 to 16, serves a, b and c. OK.
- b: panel_b of make_merged_showcase_figure re-used unchanged (DEM sampled
  x-major through msf.dem_at); vertical crop to the tube corridor; skylight
  pits (white, no-DEM cells) inside the S1 to S3 rings; arc labels and S
  labels given a white halo; the "mapped, excluded by screen" key has a
  white backing. OK.
- c: panel_c re-used; y label "Elevation (m a.s.l., ISH2004)"; the
  "vertical exaggeration 3x" note moved off the shaded void to below the
  axis. Max tau at the far end 14.1 m matches paper_numbers tau.max 14.14.
  OK.
- d: panel_d re-used; slice at s = 281 m, tau = 14.1 m; "surface
  (photogrammetry)" label moved off the surface points. OK.

### Fig 4, geometry and cover
Layout: along-tube panels in one aligned left column (a plan, b area with
width/height strip, c overburden, f elevation), d and e in the right
column (the brief said "two rows"; a shared s axis reads better and the
letter assignment of the brief is kept).
- a plan (principal-axis frame, 50 m ticks, S1 to S3, entrance): OK.
- b area, width and height (two strips, inline keys "width"/"height"): OK.
- c overburden profile: band = sigma_tau 1.48 m from roof_thickness.csv,
  asserted equal to paper_numbers reg_sigma_tau_m, shown as "±1.5 m";
  band broken at non-intact stations; 22 skylight stations as orange
  bands (s 29-36, 54-56, 89-99); 4 excluded stations as grey ticks; key in
  the top strip above the bands (bands stop at 84 % height so no text lies
  on a shaded area). OK.
- d distribution: 276 intact stations, median line at 6.02 m; top axis in
  g cm^-2 at rho_mid = 2,600 kg m^-3 from review_stats envelope, so the
  median maps to 1,564 g cm^-2 = paper_numbers shielding.median. OK.
- e overburden against span: intact and skylight classes, three
  clamped-strip curves (beta 0.5, rho 2,600 from review_stats) labelled
  inline; spot check 1 MPa at L = 20 m gives 5.10 m; key in the empty
  lower-right corner. OK.
- f centreline elevation, ISH2004 datum: OK.

### Fig 5, planetary
fig_planetary of make_figures.py re-used through a figure-size hook
(180 x 118 mm), then: left titles removed, letters added, labels
capitalised, the sqrt formula label replaced by "Maximum stable span
L_max (m)" (formula belongs in the legend; mathtext cannot set \sqrt in
Nimbus Sans), fonts clamped to 7 pt, the key of a moved above the axes
(inside the axes it hid the top of the impact-melt curve), log ticks with
plain numbers. n = 131 / 248 / 18 in the key come from the catalogue CSV;
the dashed 1.5 m line is asserted equal to round(reg_sigma_tau_m, 1).
- a, b, c: no overlaps, nothing clipped. OK.

## Numbers computed in the script (printed by the run)
- Entrance (s = 0 station) at 63.9402 N, 21.3972 W (pyproj 32627 -> 4326).
- Survey end s = 301 m; arc ticks 0 to 300.
- Fig 4: n intact 276, excluded 4, skylight 22, tau median 6.02 m,
  1 MPa curve at L = 20 m: 5.10 m; shielding factor 260 g cm^-2 per m.
- S1-S3 ring separation on the fig 3b raster: 56 m (57.0 m from data).

## Caveats to carry into the legends and the handoff
- Fig 2 b and c are onboard-visualiser renders with a height colour ramp
  (blue to red) baked into the rasters; Nature's "no rainbow" rule cannot
  be met without re-rendering the point clouds. Flag in the legend
  ("colour encodes height") or re-render.
- The old roof_panel.pdf (make_figures.fig_roof) scaled its shielding axis
  by 300 g cm^-2 per m (rho = 3,000), inconsistent with the 2,600 kg m^-3
  used everywhere else; the new Fig 4d uses rho_mid = 2,600 from
  review_stats.json.
- Fig 3b shows only the tube corridor (the 32 m hillshade pad of the v9
  panel is cropped to +15 / -24 m); the legend need not change.
- The Natural Earth coastline file (1.6 MB) was fetched once with curl;
  only the 30 kB Iceland extract is committed and the script uses it.
