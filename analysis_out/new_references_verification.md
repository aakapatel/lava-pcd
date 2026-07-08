# Verification log for new_references.bib

Date of verification session: 2026-07-08.
Method notation: "DOI resolved" = `https://doi.org/<doi>` returned a 302 redirect to the publisher landing page during this session; "CrossRef/DataCite record fetched" = full bibliographic metadata retrieved from `api.crossref.org` / `api.datacite.org` for that exact DOI; "page fetched" = the landing page / PDF itself was retrieved and read.

All 16 requested keys are **VERIFIED**. Nothing is listed under "UNVERIFIED".

---

## bib_l1skeleton — Huang et al. 2013
- **Verification:** DOI 10.1145/2461912.2461913 resolved (redirects to dl.acm.org). CrossRef record fetched: ACM Transactions on Graphics 32(4), Article 65 (pp. 1--8), 2013; authors Huang, Wu, Cohen-Or, Gong, Zhang, Li, Chen (exactly as expected). Presented at SIGGRAPH 2013.
- **Claim:** Introduces the L1-medial skeleton, a curve-skeleton representation computed directly from raw, unoriented, incomplete point clouds by applying L1-medians locally with a growing support radius — no strong requirements on point-cloud quality, normals, or watertightness.

## bib_yang_l1 — Yang et al. 2024
- **Verification:** DOI 10.1016/j.jag.2024.104062 resolved (redirects to Elsevier linkinghub, PII S1569843224004163). CrossRef record fetched: Int. J. Applied Earth Observation and Geoinformation, vol. 132, article 104062, 2024. Full author list from CrossRef: Yang, Kang, Hu, Cao, Ye, Liu, Hu, Shao. (ScienceDirect landing page itself returned 403 to automated fetch, but DOI resolution + CrossRef metadata confirm the record.)
- **Claim:** Uses the L1-medial skeleton of terrestrial lava-tube point clouds to extract a centreline, slice perpendicular cross-sections along it, and reconstruct a 3D tube model — i.e., the direct methodological precedent for skeleton-based lava-tube cross-sectioning.

## bib_lofthellir — Lee et al. 2019 (LPSC abstract)
- **Verification:** Abstract PDF fetched and read directly: https://www.hou.usra.edu/meetings/lpsc2019/pdf/3118.pdf (HTTP 200). Confirms title, abstract #3118, LPI Contribution No. 2132, and authors Pascal Lee (SETI Institute/Mars Institute/NASA Ames), Eirik Kommedal, Andrew Horchler (Astrobotic), Eric Amoroso, Kerry Snyder, Anton F. Birgisson. This is the peer-community LPSC citation preferred over the SETI/Astrobotic press release.
- **Claim:** Reports the first 3D mapping of a lava tube (the ice-rich Lofthellir cave, Iceland) by drone-borne LiDAR, documenting micro-glaciers and rockfalls, and discusses implications for exploring potentially ice-rich lava tubes on the Moon and Mars. Note: this is a conference abstract, not a peer-reviewed journal paper.

## bib_venus_tube — Carrer, Diana & Bruzzone 2026
- **Verification:** DOI 10.1038/s41467-026-68643-6 resolved to nature.com; article landing page fetched (metadata + abstract). Nature Communications 17, article 1147 (2026). The DOI guessed in the task brief is correct. Authors: Leonardo Carrer, Elena Diana, Lorenzo Bruzzone.
- **Claim:** By reanalysing 1990--1992 Magellan SAR images near candidate skylight collapses with a technique designed to detect subsurface conduits adjacent to skylights, the authors report the first radar-based observation of a lava tube on Venus (previously hypothesized but unconfirmed).

## bib_moon_radar — Ding et al. 2025
- **Verification:** URL from brief confirmed live; article page fetched (nature.com/articles/s44453-025-00013-w) and CrossRef record fetched for DOI 10.1038/s44453-025-00013-w. npj Space Exploration 1, article 11 (2025). Authors (obtained as requested): Chunyu Ding, Jiangwan Xu, Yuxiao Zhi, Ravi Sharma, Zihang Liang, Wentao Chen, Xiaohang Qiu, Changzhi Jiang.
- **Claim:** A perspective/argument piece: lunar exploration should turn to subsurface cavities (lava tubes, impact-formed voids) as ready-made radiation-shielded, thermally stable habitats, and proposes rover- and orbiter-deployed ground-penetrating radar for global-scale detection and mapping, highlighting Chinese mission (Chang'E) GPR heritage. Cite it as a proposal/perspective, not as a detection result.

## bib_horvath — Horvath, Hayne & Paige 2022
- **Verification:** DOI 10.1029/2022GL099710 resolved (redirects to agupubs.onlinelibrary.wiley.com). CrossRef record fetched, including abstract. Geophysical Research Letters 49(14), e2022GL099710, 2022.
- **Exact title (note — no "Tycho"):** "Thermal and Illumination Environments of Lunar Pits and Caves: Models and Observations From the Diviner Lunar Radiometer Experiment".
- **Claim:** Diviner measurements show the Mare Tranquillitatis and Mare Ingenii pit floors stay ~100 K warmer than the surrounding surface at night; thermophysical modelling characterizes pit/cave interiors, with cave interiors remaining at relatively stable, benign temperatures (~290 K class) compared with the surface extremes.

## bib_tabernacle — Grechi et al. 2024
- **Verification:** DOI 10.1007/s00603-024-03868-9 resolved (Springer). CrossRef record and Springer's own BibTeX export fetched; Springer landing-page metadata fetched. Rock Mechanics and Rock Engineering 57(9), 2024. Authors: Grechi, Moore, Jensen, McCreary, Czech, Festin.
- **Pagination caveat:** Springer's official metadata reports pages "1--10" (likely a metadata artifact); I omitted the page field from the .bib entry — volume/issue/DOI identify it unambiguously. Add pages only after checking the print PDF.
- **Claim:** Uses ambient-vibration recordings from a 43-station nodal geophone array plus 3D finite-element modal analysis to identify resonance modes of a partly collapsed lava-tube roof at Tabernacle Hill, Utah, showing that resonance-based methods can characterize the structure and mechanical state of tube roofs non-invasively.

## bib_carrer_moon_tube — Carrer et al. 2024
- **Verification:** DOI 10.1038/s41550-024-02302-y resolved to nature.com; article page fetched (metadata + abstract). CrossRef record fetched. Nature Astronomy 8, 1119--1126 (2024).
- **Exact title (differs from the "MiniRF" phrasing in the brief):** "Radar evidence of an accessible cave conduit on the Moon below the Mare Tranquillitatis pit".
- **Claim:** Reanalysis of 2010 LRO Mini-RF radar images shows that part of the radar reflections from the Mare Tranquillitatis pit are best explained by an accessible subsurface cave conduit tens of metres long extending from the pit floor — radar evidence, not direct imaging.

## bib_daedalus — Rossi et al. 2021 (DAEDALUS)
- **Verification:** DOI 10.25972/OPUS-22791 resolved (redirects to opus.bibliothek.uni-wuerzburg.de/22791). Full DataCite record fetched (all 21 creators, publisher Universität Würzburg, series "Würzburger Forschungsberichte in Robotik und Telematik; 21", type: Report). TU Delft research portal entry cross-checked (same title/authors/DOI).
- **Surprise:** This is NOT a Frontiers in Robotics and AI 2021 paper — no such CrossRef record exists. The citable DAEDALUS reference is a University of Würzburg technical report: Rossi, A. P., et al. (2021), "DAEDALUS -- Descent And Exploration in Deep Autonomy of Lava Underground Structures: Open Space Innovation Platform (OSIP) Lunar Caves-System Study". Cited as @techreport.
- **Claim:** Describes the ESA OSIP-studied DAEDALUS concept — a ~46 cm tether-deployed spherical robot with stereo cameras and LiDAR that descends through a lunar skylight and then rolls autonomously to 3D-map and characterize the entrance and first section of a lunar lava tube.

## bib_gpr_tube — Miyamoto et al. 2005
- **Verification:** DOI 10.1029/2005GL024159 resolved (Wiley/AGU). CrossRef record fetched: Geophysical Research Letters 32(21), L21316, 2005; 14 authors led by Hideaki Miyamoto.
- **Claim:** A shielded stepped-frequency GPR system surveyed the Komoriana (Fuji volcano, Japan) lava tube and successfully resolved the tube's roof structure, size, and depth from the surface — a terrestrial demonstration that GPR can map lava-tube roofs.

## bib_cushing2012 — Cushing 2012
- **Verification:** DOI 10.4311/2010EX0167R resolved (redirects to caves.org PDF of JCKS v74 no.1 p.33). CrossRef record fetched: Journal of Cave and Karst Studies 74(1), 33--47, 2012, single author Glen E. Cushing. (The caves.org PDF itself 404'd over HTTPS during the session, but the DOI resolution and CrossRef record fully confirm the citation.)
- **Claim:** Reviews and catalogues candidate cave entrances on Mars identified in orbital imagery (atypical pit craters, skylights, fissures), discussing their likely volcanic/tectonic origins and exploration significance.

## bib_cushing2007 — Cushing et al. 2007
- **Verification:** DOI 10.1029/2007GL030709 resolved (Wiley/AGU). CrossRef record fetched: Geophysical Research Letters 34, L17201 (issue 17), 2007; authors Cushing, Titus, Wynne, Christensen.
- **Claim:** THEMIS visible and thermal-IR observations reveal seven dark, nearly circular features ("seven sisters") on the flanks of Arsia Mons whose thermal behaviour is consistent with deep subsurface voids — the first reported possible cave skylights on Mars.

## bib_cushing2015 — Cushing, Okubo & Titus 2015
- **Verification:** DOI 10.1002/2014JE004735 resolved (Wiley/AGU). CrossRef record fetched. Journal of Geophysical Research: Planets 120(6), 1023--1043, 2015.
- **Exact title:** "Atypical pit craters on Mars: New insights from THEMIS, CTX, and HiRISE observations".
- **Claim:** Characterizes Martian atypical pit craters (APCs) using THEMIS, CTX, and HiRISE, constraining their morphologies, depths, and formation mechanisms, and their potential as entrances to subsurface cavities.

## bib_mgc3_data — MGC3 dataset
- **Verification:** DOI 10.17189/1519222 resolved (redirects to the NASA PDS bundle viewer for urn:nasa:pds:mars_mro.odyssey_multi_cavecatalog_cushing_2016, v1.0). DataCite record fetched: title "Mars Global Cave Candidate Catalog Archive Bundle", creator Glen E. Cushing, publisher NASA Planetary Data System, publication year 2016.
- **Note:** DataCite lists the publisher as NASA PDS (not USGS; Cushing is a USGS scientist) and year 2016.
- **Claim (dataset):** The Mars Global Cave Candidate Catalog (MGC3) — over 1,000 candidate cave entrances/subsurface-access features across Mars compiled from MRO and Odyssey imagery, archived as a PDS4 bundle.

## bib_lroc_pits — LROC Lunar Pit Locations shapefile + Wagner & Robinson 2021
- **Verification:** Product page fetched (HTTP 200): https://data.lroc.im-ldi.com/lroc/view_rdr/SHAPEFILE_LUNAR_PIT_LOCATIONS — "Lunar Pit Locations Shapefile", and the page itself instructs: "When citing this product, use the following reference: Wagner, R. V. and Robinson, M. S. (2021)... 52nd LPSC, Abstract #2530." The LPSC abstract PDF was also fetched and read (https://www.hou.usra.edu/meetings/lpsc2021/pdf/2530.pdf, HTTP 200; LPI Contribution No. 2548).
- **Claim:** A catalogue/shapefile of lunar pit locations, dimensions, and morphologic descriptors (roughly 281 impact-melt pits, 15 mare pits, 5 highland pits), with depths from shadow measurements; the abstract discusses pit occurrence and origin. No DOI exists for this product; the stable product URL above is the citation anchor.

---

## UNVERIFIED — do not cite
(none — all 16 requested keys verified)

---

## Surprises / deviations from the brief
1. **bib_daedalus** is a Universität Würzburg technical report (DOI 10.25972/OPUS-22791, Rossi et al. 2021, 21 authors), not a Frontiers in Robotics and AI paper.
2. **bib_horvath** exact title contains no "Tycho": it is about lunar pits/caves generally, anchored on Mare Tranquillitatis and Mare Ingenii Diviner data (GRL 49, e2022GL099710).
3. **bib_carrer_moon_tube** exact title is "Radar evidence of an accessible cave conduit on the Moon below the Mare Tranquillitatis pit" (Nature Astronomy 8, 1119--1126, 2024) — no "MiniRF" in the title.
4. **bib_venus_tube** DOI guess in the brief was correct: Nat. Commun. 17, 1147 (2026); author list is short (Carrer, Diana, Bruzzone).
5. **bib_moon_radar** is a perspective/proposal (8 authors, Ding et al.), not an observational detection paper — word citations accordingly.
6. **bib_mgc3_data** publisher per DataCite is NASA PDS, year 2016 (brief called it "USGS ... data product"; creator is USGS's G. Cushing).
7. **bib_tabernacle** page numbers are ambiguous in publisher metadata ("1--10"); pages omitted from the .bib entry.
8. **bib_lofthellir**: a proper LPSC abstract exists (#3118, 2019), so the press release is not needed.
