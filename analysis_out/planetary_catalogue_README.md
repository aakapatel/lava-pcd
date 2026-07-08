# planetary_catalogue.csv — provenance and documentation

Stage S4.1 deliverable. Machine-readable catalogue of (a) Martian candidate cave
entrances / atypical pit craters and (b) lunar pits, for gravity-scaling analysis
of collapse apertures. Built 2026-07-08 from primary archives (USGS/PDS and
LROC/PDS); no values were estimated or extrapolated.

## Contents

- `planetary_catalogue.csv` — 1340 rows: 1062 Mars + 278 Moon.
- `source_data/` — verbatim copies of every source file used:
  - `mars_cave_catalog.zip` — full MGC3 PDS4 bundle as downloaded from USGS.
  - `Mars_Cave_Catalog_shapefile.csv` — the MGC3 table actually parsed (from
    `miscellaneous/` inside the bundle; correctly ordered header).
  - `Mars_Cave_Catalog_data.csv` — the `data/` collection copy of the same table.
    WARNING: in this file the header row is misaligned with the data (data rows
    begin with the ID/Label column while the header begins with `longitude`).
    Use the shapefile CSV, or the column order defined in the archive
    description PDF (ID, longitude, latitude, type code, priority, APC diameter,
    APC depth, comment).
  - `mars_cave_catalog_archive_description.pdf` — Cushing & Okubo, MGC3 archive
    description v1.0 (2016-09-22): column definitions, type codes, survey method.
  - `LUNAR_PIT_LOCATIONS.CSV` — LROC PDS RDR lunar pit table (updated
    2021-03-09), 278 pits. NOTE: file is MacRoman-encoded (e.g. "Catalán B 1",
    "Schlüter Pit" garble under UTF-8/cp1252).
  - `wagner_robinson_2021_lpsc2530.pdf` — LPSC 52 abstract #2530 describing the
    lunar pit catalog and its summary statistics.
  - `cushing_2012_jcks.pdf` — Cushing (2012), full paper.
  - `horvath_2022_grl.pdf` — Horvath et al. (2022), full paper.

## CSV columns

| column | meaning |
|---|---|
| `body` | `mars` or `moon` |
| `feature_id` | Mars: MGC3 ID (`APC###` = atypical pit crater, `CC####` = cave candidate). Moon: LROC pit name (e.g. `Marius Hills Pit`, `King 1a`) |
| `lat` | Latitude, decimal degrees. Mars: planetocentric (MGC3 C-3). Moon: from LROC table (`Latitude`) |
| `lon_east_0_360` | Longitude, decimal degrees east, 0–360, for both bodies (MGC3 C-2 native; Moon uses the `Longitude_360` column) |
| `aperture_long_axis_m` | Mars: MGC3 `APC_Diameter` — "Diameter of the APC rounded to the nearest 5 meters ... For the small number of pits with elliptical shapes, this value represents the pit's major axis" (archive description C-6). Populated only for APCs (131 of 132). Moon: LROC `Funnel Max Diam` — maximum diameter of the outer funnel (surface aperture) |
| `aperture_short_axis_m` | Mars: not provided by MGC3 (empty). Moon: LROC `Funnel Min Diam` |
| `inner_long_axis_m` | Moon only: LROC `Inner Max Diam` — maximum diameter of the inner, vertical-walled opening. Empty for Mars |
| `inner_short_axis_m` | Moon only: LROC `Inner Min Diam` |
| `depth_m` | Mars: MGC3 `APC_Depth`, shadow-derived, "D=((Shadow-length) / tan(incidence angle))" (C-7). Moon: LROC `Depth` (shadow or stereo derived), numeric part only |
| `depth_is_minimum` | `yes` when the source flags the depth as a lower bound. Moon: LROC depth string begins with `>`. Mars: comment contains "deep" — per MGC3 C-7, "an APC described as 'deep' in the comment section signifies that no sunlit floor can be seen, and that only a minimum depth can be calculated" |
| `type` | Mars: MGC3 3-letter type code (see below). Moon: `pit (impact melt)`, `pit (mare)`, or `pit (highland)` from the LROC `Terrain` column |
| `source` | Short citation key (`Cushing2015_MGC3v1`, `WagnerRobinson2021_LROCPitAtlas`) |
| `source_detail` | Archive/file/URL provenance |
| `notes` | Mars: targeting priority (0–3) + verbatim MGC3 comment. Moon: host feature, type qualifier (`Uncertain`, `Large Bowl`, `Fracture`), long-axis azimuth (deg), raw depth string |

Aperture semantics: for the Moon, both the funnel (outer, surface) and inner
(vertical-walled) openings are carried, because the atlas reports both and the
choice matters for scaling analyses. For Mars APCs the single MGC3 diameter is
measured at the surface; APC walls are near-vertical or overhanging (Cushing et
al. 2015), so it is closer in kind to the lunar *inner* diameter than to the
funnel diameter. Decide explicitly which comparison to use in the analysis.

MGC3 type codes (archive description C-4, verbatim definitions abridged):
`APC` atypical pit crater; `sky` skylight entrance into a lava tube; `srp`
small shallow circular rimless collapse pit; `crk` deep dark crack/fracture;
`pin` pinhole; `pit` generic amorphous collapse pit smaller than APCs;
`irr` irregular; `rim` on rim of pit/channel/trough; `lat` lateral entrance
(cliff walls); `end` distal end of channel/fracture/trough; `flr` floor of
channel/fracture/trough; `pol` polar ice region; `kst` karst-type terrain;
`con` contributed; `smrp` small rimless pit (1 entry).

Counts by type in this file: sky 354, srp 218, APC 132, crk 87, pin 60, flr 42,
pit 37, irr 36, rim 30, lat 27, end 21, pol 8, kst 7, con 2, smrp 1.

## Where each number came from

### Mars (1062 entries) — fully captured

Source: Mars Global Cave Candidate Catalog (MGC3) PDS4 archive bundle,
downloaded 2026-07-08 from USGS Astrogeology (Astropedia):
`https://astrogeology.usgs.gov/ckan/dataset/c2960f77-ee20-4cc6-ac85-abbecc0fc7f6/resource/ccd54fb2-9a50-46bf-8971-80c241bb2f40/download/mars_cave_catalog.zip`
(landing page: https://astrogeology.usgs.gov/search/map/mars_global_cave_candidate_catalog_v1_cushing ;
DOI 10.17189/1519222). All 1062 rows of the bundle's table were copied without
modification (values, comments, coordinates verbatim).

Version note: the archive description PDF (v1.0, 2016-09-22) states "The table
currently contains 1036 rows of data" and its Fig. 3 caption says "(1035
total)". The bundle README is dated 2017-06-21 and the shipped table contains
1062 rows — the data file was updated after the description was written. The
1062-row file is what USGS currently distributes and is what this catalogue uses.

Aperture/depth coverage: only APC entries carry dimensions. "Most candidates in
this catalog were identified during a survey for APC candidates; measuring
non-APC candidates was outside the scope of budgeted criteria for the
cave-detection survey" (archive description C-6). Computed from the file:

- APC diameters: n = 131, range 25–520 m, median 120 m.
- APC depths: n = 132, range 17–345 m, median 66 m.

Published aggregate for context, Cushing et al. (2015) abstract: "locations of
115 APCs", "surface diameters of ~50-350 m", "depth-to-diameter (d/D) ratios
that are usually greater than 0.3 ... and can exceed values of 1.8". The MGC3
bundle notes (item 4): "Most of the APC data in this catalog was previously
published as on-line supplementary information to [Cushing et al., 2015].
Updated APC measurements from subsequent HiRISE observations have been applied,
this release supersedes all previously published APC data." Therefore the MGC3
values, not the 2015 supplement, are used here.

Caveats from the archive description carried over verbatim: "NONE can be
verified as genuine caves"; "coordinates given in this catalog are approximate"
(CTX/HiRISE not geodetically controlled; points may sit up to a few hundred
meters from the feature).

### Moon (278 entries) — fully captured

Source: LROC PDS RDR "Lunar Pit Locations" product (Wagner & Robinson),
updated 2021-03-09, downloaded 2026-07-08 from
`https://pds.lroc.im-ldi.com/data/LRO-L-LROC-5-RDR-V1.0/LROLRC_2001/EXTRAS/SHAPEFILE/LUNAR_PIT_LOCATIONS/LUNAR_PIT_LOCATIONS.CSV`
(product page: https://data.lroc.im-ldi.com/lroc/view_rdr/SHAPEFILE_LUNAR_PIT_LOCATIONS ;
recommended citation per that page: Wagner & Robinson 2021, LPSC 52, #2530).

Cross-check: the live LROC "Pits Atlas" web table
(https://www.lroc.asu.edu/atlases/pits/list, served from lroc.im-ldi.com) was
scraped the same day and contains the same 278 pits with identical funnel/inner
diameters and depths (name-set difference is zero after fixing the PDS file's
MacRoman encoding). So as of 2026-07-08 the live atlas and the March 2021 PDS
release are in sync.

Terrain breakdown in the file: impact melt 257, mare 16, highland 5 (total 278).
Computed from the file:

- Funnel max diameter: n = 266, range 6–430 m, median 29 m.
- Inner max diameter: n = 266, range 4–385 m, median 16 m.
- Depth: reported for 277 pits (40 of them lower bounds, flagged
  `depth_is_minimum=yes`); `Jackson 2b` has no depth.

Published aggregates, Wagner & Robinson (2021) LPSC #2530, verbatim: "To date,
we have identified ~281 pits in melt deposits of impact craters, 15 pits in mare
basalts, and 5 pits in non-impact-melt highland terrain." — "the median impact
melt pit diameter is 15 m (range 5-385 m), the median mare pit is 100 m
(17-175 m), and the median highland pit is 45 m (16-65 m)." — "Impact melt pit
depth/diameter ratios range from 0.08 (Crookes 1) to 2.8 (King 8b) or >3.5
(Copernicus 7) (mean = 0.67), while mare and highland pits range from 0.17
(Mare Insularum) to 2.5 (Southwest Fecunditatis) (mean = 0.64)." — "Potential
pits with diameters < 5 m were excluded from the catalog".
(The abstract's ~281/15 counts differ slightly from the released file's 257/16;
the abstract predates the March 2021 release and gives approximate numbers. The
released file is authoritative for this catalogue. The abstract's median 15 m
matches the file's inner-diameter median of 16 m to within its precision.)

Earlier published aggregate, Wagner & Robinson (2014) abstract, verbatim:
"steep-walled pits in mare basalt (n= 8), impact melt deposits (n= 221), and
highland terrain (n= 2)" (total 231; superseded by the 2021 catalog).

Update, Wagner & Robinson (2022) abstract, verbatim: "Lunar pits are small
(~10–300 m wide) collapse features with vertical walls and sometimes overhangs.
We have identified almost 300 pits, mostly in ponds of cooled impact melt inside
large craters younger than ~1 billion years. ... We investigated the 21 known
pits outside of impact melt ponds ..." (no new machine-readable table located
on the open web beyond the 2021 PDS release; see below).

Single-feature checks: Marius Hills Pit row (funnel 92 x 79 m, inner
55 x 49 m, depth 40 m, stereo-derived) is consistent with the discovery paper,
Haruyama et al. (2009): hole "nearly circular, 65 m in diameter", depth
"80 to 88 m" from initial SELENE Terrain Camera imagery — the later LROC
NAC/DTM measurements in the atlas supersede these first estimates; both are
retained here (atlas values in the CSV, Haruyama values in this README) rather
than mixed.

## Full citations

- Cushing, G.E. (2012). Candidate cave entrances on Mars. *Journal of Cave and
  Karst Studies*, 74(1), 33–47. doi:10.4311/2010EX0167R.
  NOTE: this paper is in JCKS, not JGR Planets as assumed in the task brief;
  confirmed from the paper's own header (source_data/cushing_2012_jcks.pdf).
- Cushing, G.E. (2015). Mars Global Cave Candidate Catalog PDS4 Archive Bundle
  (MGC3 v1). PDS Cartography and Imaging Sciences Node. doi:10.17189/1519222.
  Bundle authored by G.E. Cushing and C.H. Okubo, USGS Astrogeology.
- Cushing, G.E., Okubo, C.H., & Titus, T.N. (2015). Atypical pit craters on
  Mars: New insights from THEMIS, CTX, and HiRISE observations. *Journal of
  Geophysical Research: Planets*, 120, 1023–1043. doi:10.1002/2014JE004735.
- Cushing, G.E., Titus, T.N., Wynne, J.J., & Christensen, P.R. (2007). THEMIS
  observes possible cave skylights on Mars. *Geophysical Research Letters*, 34,
  L17201. doi:10.1029/2007GL030709. (The task brief's "Cushing 2015 THEMIS
  candidate caves" conflates this 2007 GRL paper with the 2015 JGR paper above.)
- Wagner, R.V., & Robinson, M.S. (2014). Distribution, formation mechanisms,
  and significance of lunar pits. *Icarus*, 237, 52–60.
  doi:10.1016/j.icarus.2014.04.002.
- Wagner, R.V., & Robinson, M.S. (2021). Occurrence and Origin of Lunar Pits:
  Observations from a New Catalog. *52nd Lunar and Planetary Science
  Conference*, LPI Contrib. No. 2548, abstract #2530.
- Wagner, R.V., & Robinson, M.S. (2022). Lunar Pit Morphology: Implications for
  Exploration. *Journal of Geophysical Research: Planets*, 127, e2022JE007328.
  doi:10.1029/2022JE007328.
- Haruyama, J., et al. (2009). Possible lunar lava tube skylight observed by
  SELENE cameras. *Geophysical Research Letters*, 36, L21206.
  doi:10.1029/2009GL040635.
- Horvath, T., Hayne, P.O., & Paige, D.A. (2022). Thermal and Illumination
  Environments of Lunar Pits and Caves: Models and Observations From the
  Diviner Lunar Radiometer Experiment. *Geophysical Research Letters*, 49,
  e2022GL099710. doi:10.1029/2022GL099710. Key abstract numbers, verbatim:
  Mare Tranquillitatis and Mare Ingenii pits "exhibit elevated thermal emission
  during the night, ~100 K warmer than the surrounding surface"; shadowed pit
  interiors "would maintain a nearly constant temperature of ~290 K".

## Not obtained on the open web / fetch later

- Wagner & Robinson (2014) full text and its Table 1 (per-pit dimensions for
  the 2014 set): behind Elsevier paywall. Not blocking — every 2014 pit is in
  the 2021 PDS release captured here, with updated measurements.
- Wagner & Robinson (2022) full text, figures, and 3D-reconstruction
  supplementary data (six pit interiors): behind AGU/Wiley paywall (HTTP 402).
  Abstract captured. Fetch via institutional access if interior geometry
  (overhang extents) is needed.
- Cushing et al. (2015) online supplementary APC table: paywalled; explicitly
  superseded by MGC3 (bundle note 4), so not needed.
- LROC Pits Atlas PDF version (announced in LPSC #2530): not located; the CSV
  and live atlas pages carry the same data.
- Any MGC3 v2: none found on Astropedia, PDS IMG, or data.nasa.gov as of
  2026-07-08; v1 (2017 update, 1062 rows) is current.
- Per-pit uncertainty estimates: neither archive publishes formal errors.
  MGC3 diameters are rounded to 5 m; MGC3 note 3 warns coordinates are
  approximate. LROC lat/lon are stated (LPSC #2530) to be better than the
  ~30 m single-NAC uncertainty.
