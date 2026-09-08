# Regression: analysis_out_v9 (orthometric, ISH2004) vs analysis_out_v6 (v8 canon)

Applied vertical translation dz = -66.0000 m (ortho_translation.json).
Verdicts: identical | rounding-safe (below half the paper's quoted unit) | datum (dz) (absolute elevation, moved by exactly dz) | datum-derived (offsets against products still in the old datum) | CHANGED (investigate).

Summary: datum (dz): 4, datum-derived: 4, explained: 1 cm quantisation of the delivered cloud: 34, explained: S3 floor check on new voxel copy: 30, explained: aperture bookkeeping: 52, explained: in-aperture cells removed: 10, identical: 427, rounding-safe: 163

## Paper headline numbers at the paper's rounding

| quantity | v6 | v9 | paper rounding |
|---|---|---|---|
| survey length (m) | 301 | 301 | 1 |
| stations | 302 | 302 | 1 |
| intact stations | 276 | 276 | 1 |
| overburden median (m) | 6.0 | 6.0 | 0.1 |
| overburden block CI low (m) | 3.1 | 3.1 | 0.1 |
| overburden block CI high (m) | 8.7 | 8.7 | 0.1 |
| overburden IQR low (m) | 3.2 | 3.2 | 0.1 |
| overburden IQR high (m) | 9.0 | 9.0 | 0.1 |
| overburden max (m) | 14.2 | 14.1 | 0.1 |
| shielding >100 g/cm2 (rho 2600) | 99% | 99% | 1 % |
| shielding >500 g/cm2 | 84% | 84% | 1 % |
| shielding >1000 g/cm2 | 72% | 72% | 1 % |
| area median (m2) | 48 | 48 | 1 |
| area p5 (m2) | 13 | 13 | 1 |
| area p95 (m2) | 76 | 76 | 1 |
| width min (m) | 6.1 | 6.1 | 0.1 |
| width max (m) | 31.6 | 31.6 | 0.1 |
| height median (m) | 8.1 | 8.1 | 0.1 |
| sinuosity mean | 1.03 | 1.03 | 0.01 |
| sinuosity max | 1.08 | 1.08 | 0.01 |
| strength demand p95 (MPa) | 1.6 | 1.6 | 0.1 |
| strength demand max (MPa) | 2.0 | 2.0 | 0.1 |
| strength demand drop-one (MPa) | 2.0 | 2.0 | 0.1 |
| strength stations n | 248 | 248 | 1 |
| sigma_tau (m) | 1.5 | 1.5 | 0.1 |
| sigma_reg (m) | 1.5 | 1.4 | 0.1 |
| S3 floor median (m) | 0.89 | 0.87 | 0.01 |
| S3 floor RMS (m) | 1.36 | 1.35 | 0.01 |
| scale ruler S1-S2 diff (m) | -0.18 | -0.18 | 0.01 |
| scale ruler S1-S3 diff (m) | 2.18 | 2.18 | 0.01 |
| volume (m3) | 14,000 | 14,000 | 1,000 |
| volume per m (m3) | 47 | 47 | 1 |
| flight-2 repeatability RMS (m) | 0.37 | 0.37 | 0.01 |
| flight-2 repeatability median (m) | 0.08 | 0.08 | 0.01 |
| ArcticDEM shape std (m) | 0.29 | 0.28 | 0.01 |
| ArcticDEM doming (m) | 0.13 | 0.13 | 0.01 |
| ArcticDEM tilt (m/km) | 0.42 | 0.43 | 0.01 |
| ArcticDEM datum offset, ours minus (m) | 0.61 | -65.40 | 0.01 |
| IslandsDEM datum offset, ours minus (m) | 66.63 | 0.63 | 0.01 |

## Explained differences (cause traced during the v9 run)

- **aperture bookkeeping** (52 fields, pattern `n_skylight|n_no_dem|skylight_runs|near_skylight_median|ceiling_thinning`): stations s=31..34 lie inside the S1 opening (2.3 x 1.6 m); the surface-only orthometric crop has no return there, so run_roof classes them 'no_dem' instead of 'skylight' (the v6 DEM had values in those cells only because the merged July cloud carried tube points from the superseded 4-DOF registration inside the hole). The geometric aperture count is unchanged (18 + 4 = 22), the 276 intact stations are unchanged, and the run-based terrain diagnostics regroup (S1 splits into two runs).
- **S3 floor check on new voxel copy** (30 fields, pattern `sigma_reg|sigma_tau|tau_threshold|floor_through_skylight|anchor_floor`): the 10 cm working copy is re-voxelised from the new cloud; the 5 mm horizontal ISN2016/WGS84 difference changes voxel membership near S3 (360 vs 352 matched points), RMS 1.347 vs 1.356 m, so sigma_reg 1.45 vs 1.46 and sigma_tau 1.48 vs 1.49 m; the paper quotes both as 1.5 m.
- **in-aperture cells removed** (10 fields, pattern `national_dem_check`): 165 cells inside the skylights had DEM values only in v6 (tube points 4-5 m below the rim); they are absent in the surface-only crop, which lowers the ArcticDEM shape std from 0.293 to 0.281 m (paper: 0.29 m, an upper bound) and the tilt 0.42 -> 0.43 m/km. Datum offsets carry the -66.00 m shift.
- **1 cm quantisation of the delivered cloud** (34 fields, pattern `stats\.median|orbital_dem_test|basis_comparison|shielding|kappa_sensitivity|tau_over_L\.second|sensitivity_tau_gt_0p5`): the new cloud is stored at 1 cm; per-cell maxima come out 0.5-1 cm lower than old minus 66.00 (median -66.01), so tau is lower by 0.0075 m in the mean (median 6.0255 -> 6.0155 m, shielding 1567 -> 1564 g/cm2, 1.52 -> 1.51x). Paper roundings (6.0 m, about 1,500 g/cm2, 1.5x) unchanged; the smallest tau/L moves because a 0.10 m cover became 0.09 m.

## Fields flagged CHANGED (unexplained)

none

## All fields

| field | v6 | v9 | diff | verdict |
|---|---|---|---|---|
| paper_numbers.area.max | 81.4000 | 81.4000 | 0.0000 | identical |
| paper_numbers.area.mean_boot_ci95[0] | 45.2000 | 45.2000 | 0.0000 | identical |
| paper_numbers.area.mean_boot_ci95[1] | 49.2000 | 49.2000 | 0.0000 | identical |
| paper_numbers.area.median | 48.0000 | 48.0000 | 0.0000 | identical |
| paper_numbers.area.min | 5.1000 | 5.1000 | 0.0000 | identical |
| paper_numbers.area.p5 | 13.4000 | 13.4000 | 0.0000 | identical |
| paper_numbers.area.p95 | 75.8000 | 75.8000 | 0.0000 | identical |
| paper_numbers.area.q25 | 36.0000 | 36.0000 | 0.0000 | identical |
| paper_numbers.area.q75 | 59.2000 | 59.2000 | 0.0000 | identical |
| paper_numbers.centreline_length_m | 301.0000 | 301.0000 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.10MPa.at_tau_14m | 104.7750 | 104.7750 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.10MPa.at_tau_6m | 69.7252 | 69.7252 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.1MPa.at_tau_14m | 33.1328 | 33.1328 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.1MPa.at_tau_6m | 22.0490 | 22.0490 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.5MPa.at_tau_14m | 74.0871 | 74.0871 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.earth.5MPa.at_tau_6m | 49.3031 | 49.3031 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.10MPa.at_tau_14m | 170.3748 | 170.3748 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.10MPa.at_tau_6m | 113.3802 | 113.3802 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.1MPa.at_tau_14m | 53.8772 | 53.8772 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.1MPa.at_tau_6m | 35.8540 | 35.8540 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.5MPa.at_tau_14m | 120.4732 | 120.4732 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.mars.5MPa.at_tau_6m | 80.1719 | 80.1719 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.10MPa.at_tau_14m | 257.8308 | 257.8308 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.10MPa.at_tau_6m | 171.5800 | 171.5800 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.1MPa.at_tau_14m | 81.5333 | 81.5333 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.1MPa.at_tau_6m | 54.2584 | 54.2584 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.5MPa.at_tau_14m | 182.3139 | 182.3139 | 0.0000 | identical |
| paper_numbers.envelope.Lmax_sqrt_model.moon.5MPa.at_tau_6m | 121.3254 | 121.3254 | 0.0000 | identical |
| paper_numbers.envelope.beta | 0.5000 | 0.5000 | 0.0000 | identical |
| paper_numbers.envelope.kappa_sensitivity.drop_one | 0.0080 | 0.0072 | -0.0008 | explained: 1 cm quantisation of the delivered cloud |
| paper_numbers.envelope.kappa_sensitivity.median | 0.5455 | 0.5451 | -0.0004 | rounding-safe |
| paper_numbers.envelope.kappa_sensitivity.p5 | 0.0881 | 0.0884 | 0.0003 | rounding-safe |
| paper_numbers.envelope.kappa_sensitivity.single_min | 0.0034 | 0.0034 | 0.0000 | identical |
| paper_numbers.envelope.rho_mid | 2600.0000 | 2600.0000 | 0.0000 | identical |
| paper_numbers.envelope.sigma_req_MPa.drop_one | 1.9539 | 1.9669 | 0.0130 | rounding-safe |
| paper_numbers.envelope.sigma_req_MPa.max | 2.0462 | 2.0462 | 0.0000 | identical |
| paper_numbers.envelope.sigma_req_MPa.median | 0.2140 | 0.2142 | 0.0002 | rounding-safe |
| paper_numbers.envelope.sigma_req_MPa.n | 248 | 248 | 0 | identical |
| paper_numbers.envelope.sigma_req_MPa.p90 | 0.7710 | 0.7722 | 0.0013 | rounding-safe |
| paper_numbers.envelope.sigma_req_MPa.p95 | 1.5647 | 1.5647 | 0.0000 | identical |
| paper_numbers.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.drop_one | 3.0571 | 3.0866 | 0.0295 | explained: 1 cm quantisation of the delivered cloud |
| paper_numbers.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.max | 3.7041 | 3.7041 | 0.0000 | identical |
| paper_numbers.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.n | 270 | 270 | 0 | identical |
| paper_numbers.envelope.span_scale.mars | 1.6261 | 1.6261 | 0.0000 | identical |
| paper_numbers.envelope.span_scale.moon | 2.4608 | 2.4608 | 0.0000 | identical |
| paper_numbers.envelope.tau_over_L.median | 0.5455 | 0.5451 | -0.0004 | rounding-safe |
| paper_numbers.envelope.tau_over_L.min | 0.0034 | 0.0034 | 0.0000 | identical |
| paper_numbers.envelope.tau_over_L.p25 | 0.2974 | 0.2970 | -0.0004 | rounding-safe |
| paper_numbers.envelope.tau_over_L.p5 | 0.0881 | 0.0884 | 0.0003 | rounding-safe |
| paper_numbers.envelope.tau_over_L.second | 0.0080 | 0.0072 | -0.0008 | explained: 1 cm quantisation of the delivered cloud |
| paper_numbers.envelope.tau_threshold_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
| paper_numbers.eta.median | 1.5300 | 1.5300 | 0.0000 | identical |
| paper_numbers.eta.q25 | 1.3100 | 1.3100 | 0.0000 | identical |
| paper_numbers.eta.q75 | 1.8400 | 1.8400 | 0.0000 | identical |
| paper_numbers.height.max | 11.4100 | 11.4100 | 0.0000 | identical |
| paper_numbers.height.median | 8.0900 | 8.0900 | 0.0000 | identical |
| paper_numbers.height.min | 3.5300 | 3.5300 | 0.0000 | identical |
| paper_numbers.mars.aperture_median_m | 120.0000 | 120.0000 | 0.0000 | identical |
| paper_numbers.mars.counts.ge_10m | 131 | 131 | 0 | identical |
| paper_numbers.mars.counts.ge_20m | 131 | 131 | 0 | identical |
| paper_numbers.mars.counts.ge_50m | 120 | 120 | 0 | identical |
| paper_numbers.mars.counts.ge_5m | 131 | 131 | 0 | identical |
| paper_numbers.mars.counts.ge_narrowest | 131 | 131 | 0 | identical |
| paper_numbers.mars.n_dimensioned | 131 | 131 | 0 | identical |
| paper_numbers.mars.n_total | 1062 | 1062 | 0 | identical |
| paper_numbers.min_tau_over_L_retired | 0.0034 | 0.0034 | 0.0000 | identical |
| paper_numbers.min_tau_over_L_station.s | 51.0000 | 51.0000 | 0.0000 | identical |
| paper_numbers.min_tau_over_L_station.span | 12.2100 | 12.2100 | 0.0000 | identical |
| paper_numbers.min_tau_over_L_station.tau | 0.0400 | 0.0400 | 0.0000 | identical |
| paper_numbers.moon.aperture_median_m | 29.0000 | 29.0000 | 0.0000 | identical |
| paper_numbers.moon.counts.ge_10m | 258 | 258 | 0 | identical |
| paper_numbers.moon.counts.ge_20m | 190 | 190 | 0 | identical |
| paper_numbers.moon.counts.ge_50m | 63 | 63 | 0 | identical |
| paper_numbers.moon.counts.ge_5m | 266 | 266 | 0 | identical |
| paper_numbers.moon.counts.ge_narrowest | 265 | 265 | 0 | identical |
| paper_numbers.moon.n_dimensioned | 266 | 266 | 0 | identical |
| paper_numbers.moon.n_total | 278 | 278 | 0 | identical |
| paper_numbers.n_multipass | 0 | 0 | 0 | identical |
| paper_numbers.n_no_dem | 0 | 4 | 4 | explained: aperture bookkeeping |
| paper_numbers.n_skylight_stations | 22 | 18 | -4 | explained: aperture bookkeeping |
| paper_numbers.n_stations | 302 | 302 | 0 | identical |
| paper_numbers.n_tau_clean | 276 | 276 | 0 | identical |
| paper_numbers.n_valid_sections | 296 | 296 | 0 | identical |
| paper_numbers.reg_matches | 3 | 3 | 0 | identical |
| paper_numbers.reg_rim_rms_m | 1.7331 | 1.7331 | 0.0000 | identical |
| paper_numbers.reg_sigma_reg_z_m | 1.4600 | 1.4500 | -0.0100 | explained: S3 floor check on new voxel copy |
| paper_numbers.reg_sigma_tau_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
| paper_numbers.reg_vres_mm[0] | -5240 | -5240 | 0 | identical |
| paper_numbers.reg_vres_mm[1] | -5864 | -5864 | 0 | identical |
| paper_numbers.reg_vres_mm[2] | -3582 | -3582 | 0 | identical |
| paper_numbers.registration_validation.floor_through_skylight_s3.annulus_outside_cone.median_m | 0.9072 | 0.8983 | -0.0088 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.annulus_outside_cone.n | 321 | 327 | 6 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.annulus_outside_cone.p10_m | -0.0687 | -0.0789 | -0.0102 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.annulus_outside_cone.p90_m | 1.3719 | 1.3727 | 0.0009 | rounding-safe |
| paper_numbers.registration_validation.floor_through_skylight_s3.annulus_outside_cone.rms_m | 1.3603 | 1.3516 | -0.0087 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.median_m | 0.8869 | 0.8738 | -0.0131 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.n | 352 | 360 | 8 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.p10_m | -1.1877 | -1.2324 | -0.0447 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.p90_m | 1.3423 | 1.3372 | -0.0051 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.floor_through_skylight_s3.rms_m | 1.3561 | 1.3469 | -0.0092 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.n_negative_tau_stations | 4 | 4 | 0 | identical |
| paper_numbers.registration_validation.rim_horizontal_rms_m | 1.7331 | 1.7331 | 0.0000 | identical |
| paper_numbers.registration_validation.sigma_reg_components.anchor_floor_rms_m | 1.3561 | 1.3469 | -0.0092 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.sigma_reg_components.drift_rss_m | 0.5318 | 0.5318 | -0.0000 | rounding-safe |
| paper_numbers.registration_validation.sigma_reg_z_m | 1.4600 | 1.4500 | -0.0100 | explained: S3 floor check on new voxel copy |
| paper_numbers.registration_validation.two_slam_zones.s0-80.mean | -0.2324 | -0.2324 | 0.0000 | identical |
| paper_numbers.registration_validation.two_slam_zones.s0-80.n | 80 | 80 | 0 | identical |
| paper_numbers.registration_validation.two_slam_zones.s0-80.std | 0.3220 | 0.3220 | 0.0000 | identical |
| paper_numbers.registration_validation.two_slam_zones.s160-240.mean | -0.1664 | -0.1664 | -0.0000 | rounding-safe |
| paper_numbers.registration_validation.two_slam_zones.s160-240.n | 80 | 80 | 0 | identical |
| paper_numbers.registration_validation.two_slam_zones.s160-240.std | 0.4783 | 0.4783 | -0.0000 | rounding-safe |
| paper_numbers.registration_validation.two_slam_zones.s240-340.mean | -0.1012 | -0.1012 | -0.0000 | rounding-safe |
| paper_numbers.registration_validation.two_slam_zones.s240-340.n | 99 | 99 | 0 | identical |
| paper_numbers.registration_validation.two_slam_zones.s240-340.std | 0.1277 | 0.1277 | 0.0000 | rounding-safe |
| paper_numbers.registration_validation.two_slam_zones.s80-160.mean | -0.1352 | -0.1352 | -0.0000 | rounding-safe |
| paper_numbers.registration_validation.two_slam_zones.s80-160.n | 80 | 80 | 0 | identical |
| paper_numbers.registration_validation.two_slam_zones.s80-160.std | 0.1281 | 0.1281 | 0.0000 | rounding-safe |
| paper_numbers.scale_ruler_diff_m[0] | -0.1810 | -0.1810 | 0.0000 | identical |
| paper_numbers.scale_ruler_diff_m[1] | 2.1790 | 2.1790 | 0.0000 | identical |
| paper_numbers.scale_ruler_diff_m[2] | 2.4250 | 2.4250 | 0.0000 | identical |
| paper_numbers.shielding.atmosphere_ratio_median | 1.5200 | 1.5100 | -0.0100 | explained: 1 cm quantisation of the delivered cloud |
| paper_numbers.shielding.median | 1567.0000 | 1564.0000 | -3.0000 | rounding-safe |
| paper_numbers.shielding.min | 11.0000 | 11.0000 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2200.frac_gt_100 | 0.9819 | 0.9783 | -0.0036 | rounding-safe |
| paper_numbers.shielding_sensitivity.2200.frac_gt_1000 | 0.6884 | 0.6884 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2200.frac_gt_500 | 0.8043 | 0.8043 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2200.median_g_cm2 | 1325.6100 | 1323.4100 | -2.2000 | rounding-safe |
| paper_numbers.shielding_sensitivity.2600.frac_gt_100 | 0.9855 | 0.9855 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2600.frac_gt_1000 | 0.7210 | 0.7210 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2600.frac_gt_500 | 0.8370 | 0.8370 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.2600.median_g_cm2 | 1566.6300 | 1564.0300 | -2.6000 | rounding-safe |
| paper_numbers.shielding_sensitivity.3000.frac_gt_100 | 0.9855 | 0.9855 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.3000.frac_gt_1000 | 0.7428 | 0.7428 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.3000.frac_gt_500 | 0.8659 | 0.8659 | 0.0000 | identical |
| paper_numbers.shielding_sensitivity.3000.median_g_cm2 | 1807.6500 | 1804.6500 | -3.0000 | rounding-safe |
| paper_numbers.sinuosity.max | 1.0790 | 1.0790 | 0.0000 | identical |
| paper_numbers.sinuosity.mean | 1.0320 | 1.0320 | 0.0000 | identical |
| paper_numbers.span_scale.mars | 1.6300 | 1.6300 | 0.0000 | identical |
| paper_numbers.span_scale.moon | 2.4600 | 2.4600 | 0.0000 | identical |
| paper_numbers.stats.ci_block.20[0] | 3.5100 | 3.5000 | -0.0100 | rounding-safe |
| paper_numbers.stats.ci_block.20[1] | 8.3260 | 8.3160 | -0.0100 | rounding-safe |
| paper_numbers.stats.ci_block.30[0] | 3.1400 | 3.1400 | 0.0000 | identical |
| paper_numbers.stats.ci_block.30[1] | 8.6555 | 8.6505 | -0.0050 | rounding-safe |
| paper_numbers.stats.ci_block.50[0] | 2.4390 | 2.4290 | -0.0100 | rounding-safe |
| paper_numbers.stats.ci_block.50[1] | 8.9805 | 8.9755 | -0.0050 | rounding-safe |
| paper_numbers.stats.ci_iid[0] | 5.7305 | 5.7255 | -0.0050 | rounding-safe |
| paper_numbers.stats.ci_iid[1] | 6.4015 | 6.3920 | -0.0095 | rounding-safe |
| paper_numbers.stats.corr_len_1e_m | 50 | 50 | 0 | identical |
| paper_numbers.stats.corr_len_zero_m | 102 | 102 | 0 | identical |
| paper_numbers.stats.median | 6.0255 | 6.0155 | -0.0100 | explained: 1 cm quantisation of the delivered cloud |
| paper_numbers.stats.n_eff | 3.3673 | 3.3678 | 0.0005 | rounding-safe |
| paper_numbers.stats.sigma_tau_systematic_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
| paper_numbers.tau.max | 14.1500 | 14.1400 | -0.0100 | rounding-safe |
| paper_numbers.tau.median | 6.0300 | 6.0200 | -0.0100 | rounding-safe |
| paper_numbers.tau.median_ci95[0] | 5.7400 | 5.7300 | -0.0100 | rounding-safe |
| paper_numbers.tau.median_ci95[1] | 6.4000 | 6.3900 | -0.0100 | rounding-safe |
| paper_numbers.tau.min | 0.0400 | 0.0400 | 0.0000 | identical |
| paper_numbers.tau.q25 | 3.2400 | 3.2400 | 0.0000 | identical |
| paper_numbers.tau.q75 | 9.0100 | 9.0000 | -0.0100 | rounding-safe |
| paper_numbers.tau.sigma | 1.4900 | 1.4800 | -0.0100 | rounding-safe |
| paper_numbers.tau_s_range_m[0] | 0.0000 | 0.0000 | 0.0000 | identical |
| paper_numbers.tau_s_range_m[1] | 301.0000 | 301.0000 | 0.0000 | identical |
| paper_numbers.terrain.ceiling_thinning.ci95[0] | -0.1699 | -0.1490 | 0.0209 | explained: aperture bookkeeping |
| paper_numbers.terrain.ceiling_thinning.ci95[1] | -0.1093 | -0.0933 | 0.0159 | explained: aperture bookkeeping |
| paper_numbers.terrain.ceiling_thinning.n | 108 | 108 | 0 | identical |
| paper_numbers.terrain.ceiling_thinning.slope_m_per_m | -0.1345 | -0.1164 | 0.0181 | explained: aperture bookkeeping |
| paper_numbers.terrain.const_ceiling_control.r | 0.9054 | 0.9053 | -0.0001 | rounding-safe |
| paper_numbers.terrain.const_ceiling_control.r2 | 0.8155 | 0.8153 | -0.0002 | rounding-safe |
| paper_numbers.terrain.const_ceiling_control.rms_residual_m | 1.6953 | 1.6953 | -0.0000 | rounding-safe |
| paper_numbers.terrain.corr_tau_zdem | 0.9054 | 0.9053 | -0.0001 | rounding-safe |
| paper_numbers.terrain.dem_share | 0.7047 | 0.7045 | -0.0002 | rounding-safe |
| paper_numbers.terrain.near_skylight_median.10to20 | 2.2400 | 2.6655 | 0.4255 | explained: aperture bookkeeping |
| paper_numbers.terrain.near_skylight_median.20to30 | 4.6300 | 4.1080 | -0.5220 | explained: aperture bookkeeping |
| paper_numbers.terrain.near_skylight_median.lt10 | 1.3560 | 1.3560 | 0.0000 | identical |
| paper_numbers.terrain.net_centreline.end | 234.2190 | 168.2190 | -66.0000 | datum (dz) |
| paper_numbers.terrain.net_centreline.start | 232.5970 | 166.5970 | -66.0000 | datum (dz) |
| paper_numbers.terrain.reach_median.0-50 | 1.8355 | 1.8305 | -0.0050 | rounding-safe |
| paper_numbers.terrain.reach_median.100-150 | 4.6860 | 4.6710 | -0.0150 | rounding-safe |
| paper_numbers.terrain.reach_median.150-200 | 6.4750 | 6.4700 | -0.0050 | rounding-safe |
| paper_numbers.terrain.reach_median.200-250 | 7.9910 | 7.9810 | -0.0100 | rounding-safe |
| paper_numbers.terrain.reach_median.250-300 | 12.9425 | 12.9325 | -0.0100 | rounding-safe |
| paper_numbers.terrain.reach_median.50-100 | 1.5440 | 1.5390 | -0.0050 | rounding-safe |
| paper_numbers.terrain.skylight_runs[0].d_tau | -5.6004 | -5.5566 | 0.0438 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[0].d_zceil | 0.3957 | 0.2732 | -0.1224 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[0].d_zdem | -5.2048 | -5.2834 | -0.0786 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[0].run_s[0] | 29.0000 | 29.0000 | 0.0000 | identical |
| paper_numbers.terrain.skylight_runs[0].run_s[1] | 36.0000 | 30.0000 | -6.0000 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[1].d_tau | -7.2389 | -5.9394 | 1.2996 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[1].d_zceil | 2.2604 | 0.8009 | -1.4594 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[1].d_zdem | -4.9786 | -5.1385 | -0.1599 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[1].run_s[0] | 54.0000 | 35.0000 | -19.0000 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[1].run_s[1] | 56.0000 | 36.0000 | -20.0000 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[2].d_tau | -5.8345 | -7.2355 | -1.4010 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[2].d_zceil | 2.2244 | 2.2604 | 0.0360 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[2].d_zdem | -3.6101 | -4.9751 | -1.3650 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[2].run_s[0] | 89.0000 | 54.0000 | -35.0000 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[2].run_s[1] | 99.0000 | 56.0000 | -43.0000 | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[3].d_tau | None | -5.8329 |  | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[3].d_zceil | None | 2.2244 |  | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[3].d_zdem | None | -3.6085 |  | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[3].run_s[0] | None | 89.0000 |  | explained: aperture bookkeeping |
| paper_numbers.terrain.skylight_runs[3].run_s[1] | None | 99.0000 |  | explained: aperture bookkeeping |
| paper_numbers.terrain.tau_hist_2m.counts[0] | 47 | 47 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[1] | 30 | 30 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[2] | 59 | 59 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[3] | 44 | 44 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[4] | 34 | 34 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[5] | 22 | 22 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[6] | 38 | 38 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.counts[7] | 2 | 2 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[0] | 0 | 0 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[1] | 2 | 2 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[2] | 4 | 4 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[3] | 6 | 6 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[4] | 8 | 8 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[5] | 10 | 10 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[6] | 12 | 12 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[7] | 14 | 14 | 0 | identical |
| paper_numbers.terrain.tau_hist_2m.edges[8] | 16 | 16 | 0 | identical |
| paper_numbers.terrain.variance.tau | 15.5763 | 15.5609 | -0.0155 | rounding-safe |
| paper_numbers.terrain.variance.z_ceil | 2.8742 | 2.8742 | 0.0000 | identical |
| paper_numbers.terrain.variance.z_dem | 10.9771 | 10.9631 | -0.0140 | rounding-safe |
| paper_numbers.terrain.zceil_range | 10.5590 | 10.5590 | 0.0000 | identical |
| paper_numbers.terrain.zdem_range | 12.6600 | 12.6500 | -0.0100 | rounding-safe |
| paper_numbers.vertical_q90_diffs_m[0] | -0.0430 | -0.0430 | 0.0000 | identical |
| paper_numbers.vertical_q90_diffs_m[1] | -1.3690 | -1.3690 | 0.0000 | identical |
| paper_numbers.vertical_q90_diffs_m[2] | 0.4710 | 0.4710 | 0.0000 | identical |
| paper_numbers.width.max | 31.6000 | 31.6000 | 0.0000 | identical |
| paper_numbers.width.median | 10.5300 | 10.5300 | 0.0000 | identical |
| paper_numbers.width.min | 6.0700 | 6.0700 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.10MPa.at_tau_14m | 104.7750 | 104.7750 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.10MPa.at_tau_6m | 69.7252 | 69.7252 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.1MPa.at_tau_14m | 33.1328 | 33.1328 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.1MPa.at_tau_6m | 22.0490 | 22.0490 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.5MPa.at_tau_14m | 74.0871 | 74.0871 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.earth.5MPa.at_tau_6m | 49.3031 | 49.3031 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.10MPa.at_tau_14m | 170.3748 | 170.3748 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.10MPa.at_tau_6m | 113.3802 | 113.3802 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.1MPa.at_tau_14m | 53.8772 | 53.8772 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.1MPa.at_tau_6m | 35.8540 | 35.8540 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.5MPa.at_tau_14m | 120.4732 | 120.4732 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.mars.5MPa.at_tau_6m | 80.1719 | 80.1719 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.10MPa.at_tau_14m | 257.8308 | 257.8308 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.10MPa.at_tau_6m | 171.5800 | 171.5800 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.1MPa.at_tau_14m | 81.5333 | 81.5333 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.1MPa.at_tau_6m | 54.2584 | 54.2584 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.5MPa.at_tau_14m | 182.3139 | 182.3139 | 0.0000 | identical |
| review_stats.envelope.Lmax_sqrt_model.moon.5MPa.at_tau_6m | 121.3254 | 121.3254 | 0.0000 | identical |
| review_stats.envelope.beta | 0.5000 | 0.5000 | 0.0000 | identical |
| review_stats.envelope.kappa_sensitivity.drop_one | 0.0080 | 0.0072 | -0.0008 | explained: 1 cm quantisation of the delivered cloud |
| review_stats.envelope.kappa_sensitivity.median | 0.5455 | 0.5451 | -0.0004 | rounding-safe |
| review_stats.envelope.kappa_sensitivity.p5 | 0.0881 | 0.0884 | 0.0003 | rounding-safe |
| review_stats.envelope.kappa_sensitivity.single_min | 0.0034 | 0.0034 | 0.0000 | identical |
| review_stats.envelope.rho_mid | 2600.0000 | 2600.0000 | 0.0000 | identical |
| review_stats.envelope.sigma_req_MPa.drop_one | 1.9539 | 1.9669 | 0.0130 | rounding-safe |
| review_stats.envelope.sigma_req_MPa.max | 2.0462 | 2.0462 | 0.0000 | identical |
| review_stats.envelope.sigma_req_MPa.median | 0.2140 | 0.2142 | 0.0002 | rounding-safe |
| review_stats.envelope.sigma_req_MPa.n | 248 | 248 | 0 | identical |
| review_stats.envelope.sigma_req_MPa.p90 | 0.7710 | 0.7722 | 0.0013 | rounding-safe |
| review_stats.envelope.sigma_req_MPa.p95 | 1.5647 | 1.5647 | 0.0000 | identical |
| review_stats.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.drop_one | 3.0571 | 3.0866 | 0.0295 | explained: 1 cm quantisation of the delivered cloud |
| review_stats.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.max | 3.7041 | 3.7041 | 0.0000 | identical |
| review_stats.envelope.sigma_req_MPa.sensitivity_tau_gt_0p5.n | 270 | 270 | 0 | identical |
| review_stats.envelope.span_scale.mars | 1.6261 | 1.6261 | 0.0000 | identical |
| review_stats.envelope.span_scale.moon | 2.4608 | 2.4608 | 0.0000 | identical |
| review_stats.envelope.tau_over_L.median | 0.5455 | 0.5451 | -0.0004 | rounding-safe |
| review_stats.envelope.tau_over_L.min | 0.0034 | 0.0034 | 0.0000 | identical |
| review_stats.envelope.tau_over_L.p25 | 0.2974 | 0.2970 | -0.0004 | rounding-safe |
| review_stats.envelope.tau_over_L.p5 | 0.0881 | 0.0884 | 0.0003 | rounding-safe |
| review_stats.envelope.tau_over_L.second | 0.0080 | 0.0072 | -0.0008 | explained: 1 cm quantisation of the delivered cloud |
| review_stats.envelope.tau_threshold_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
| review_stats.mars_pits.apc_median | 120.0000 | 120.0000 | 0.0000 | identical |
| review_stats.mars_pits.n_apc | 132 | 132 | 0 | identical |
| review_stats.mars_pits.n_apc_dimensioned | 131 | 131 | 0 | identical |
| review_stats.mars_pits.n_sky_type | 354 | 354 | 0 | identical |
| review_stats.mars_pits.n_total | 1062 | 1062 | 0 | identical |
| review_stats.moon_pits.all_inner.frac_lt_10 | 0.2068 | 0.2068 | 0.0000 | identical |
| review_stats.moon_pits.all_inner.frac_lt_narrowest | 0.0902 | 0.0902 | 0.0000 | identical |
| review_stats.moon_pits.all_inner.median | 16.0000 | 16.0000 | 0.0000 | identical |
| review_stats.moon_pits.all_inner.n | 266 | 266 | 0 | identical |
| review_stats.moon_pits.by_type.pit (highland).depth_median | 27.0000 | 27.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (highland).inner_ge_narrowest | 5 | 5 | 0 | identical |
| review_stats.moon_pits.by_type.pit (highland).inner_max | 45.0000 | 45.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (highland).inner_median | 37.0000 | 37.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (highland).inner_min | 19.0000 | 19.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (highland).n | 5 | 5 | 0 | identical |
| review_stats.moon_pits.by_type.pit (highland).n_inner | 5 | 5 | 0 | identical |
| review_stats.moon_pits.by_type.pit (highland).n_outer | 5 | 5 | 0 | identical |
| review_stats.moon_pits.by_type.pit (highland).outer_median | 63.0000 | 63.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).depth_median | 9.0000 | 9.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).inner_ge_narrowest | 224 | 224 | 0 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).inner_max | 385.0000 | 385.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).inner_median | 15.5000 | 15.5000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).inner_min | 4.0000 | 4.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).n | 257 | 257 | 0 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).n_inner | 248 | 248 | 0 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).n_outer | 246 | 246 | 0 | identical |
| review_stats.moon_pits.by_type.pit (impact melt).outer_median | 26.5000 | 26.5000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (mare).depth_median | 40.0000 | 40.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (mare).inner_ge_narrowest | 13 | 13 | 0 | identical |
| review_stats.moon_pits.by_type.pit (mare).inner_max | 175.0000 | 175.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (mare).inner_median | 100.0000 | 100.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (mare).inner_min | 14.0000 | 14.0000 | 0.0000 | identical |
| review_stats.moon_pits.by_type.pit (mare).n | 16 | 16 | 0 | identical |
| review_stats.moon_pits.by_type.pit (mare).n_inner | 13 | 13 | 0 | identical |
| review_stats.moon_pits.by_type.pit (mare).n_outer | 15 | 15 | 0 | identical |
| review_stats.moon_pits.by_type.pit (mare).outer_median | 175.0000 | 175.0000 | 0.0000 | identical |
| review_stats.moon_pits.narrowest_section_m | 6.0700 | 6.0700 | 0.0000 | identical |
| review_stats.scale_ruler.aerial.S1S2 | 22.4430 | 22.4430 | 0.0000 | identical |
| review_stats.scale_ruler.aerial.S1S3 | 57.0158 | 57.0158 | 0.0000 | identical |
| review_stats.scale_ruler.aerial.S2S3 | 34.7254 | 34.7254 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_map.S1S2 | 22.2622 | 22.2622 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_map.S1S3 | 59.1951 | 59.1951 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_map.S2S3 | 37.1500 | 37.1500 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_aerial.S1S2 | -0.1808 | -0.1808 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_aerial.S1S3 | 2.1793 | 2.1793 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_aerial.S2S3 | 2.4246 | 2.4246 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_flf.S1S2 | -0.3044 | -0.3044 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_flf.S1S3 | 1.9356 | 1.9356 | 0.0000 | identical |
| review_stats.scale_ruler.dlio_minus_flf.S2S3 | 2.1814 | 2.1814 | 0.0000 | identical |
| review_stats.scale_ruler.flf_map.S1S2 | 22.5666 | 22.5666 | 0.0000 | identical |
| review_stats.scale_ruler.flf_map.S1S3 | 57.2595 | 57.2595 | 0.0000 | identical |
| review_stats.scale_ruler.flf_map.S2S3 | 34.9686 | 34.9686 | 0.0000 | identical |
| review_stats.scale_ruler.flf_minus_aerial.S1S2 | 0.1236 | 0.1236 | 0.0000 | identical |
| review_stats.scale_ruler.flf_minus_aerial.S1S3 | 0.2437 | 0.2437 | 0.0000 | identical |
| review_stats.scale_ruler.flf_minus_aerial.S2S3 | 0.2432 | 0.2432 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.aerial[0] | 2.2791 | 2.2791 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.aerial[1] | 1.4142 | 1.4142 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.aerial[2] | 4.6683 | 4.6683 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.dlio[0] | 2.5626 | 2.5626 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.dlio[1] | 2.4495 | 2.4495 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.dlio[2] | 7.6375 | 7.6375 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.flf[0] | 2.1839 | 2.1839 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.flf[1] | 2.0976 | 2.0976 | 0.0000 | identical |
| review_stats.scale_ruler.hole_semi_major.flf[2] | 5.1679 | 5.1679 | 0.0000 | identical |
| review_stats.shielding.2200.frac_gt_100 | 0.9819 | 0.9783 | -0.0036 | rounding-safe |
| review_stats.shielding.2200.frac_gt_1000 | 0.6884 | 0.6884 | 0.0000 | identical |
| review_stats.shielding.2200.frac_gt_500 | 0.8043 | 0.8043 | 0.0000 | identical |
| review_stats.shielding.2200.median_g_cm2 | 1325.6100 | 1323.4100 | -2.2000 | rounding-safe |
| review_stats.shielding.2600.frac_gt_100 | 0.9855 | 0.9855 | 0.0000 | identical |
| review_stats.shielding.2600.frac_gt_1000 | 0.7210 | 0.7210 | 0.0000 | identical |
| review_stats.shielding.2600.frac_gt_500 | 0.8370 | 0.8370 | 0.0000 | identical |
| review_stats.shielding.2600.median_g_cm2 | 1566.6300 | 1564.0300 | -2.6000 | rounding-safe |
| review_stats.shielding.3000.frac_gt_100 | 0.9855 | 0.9855 | 0.0000 | identical |
| review_stats.shielding.3000.frac_gt_1000 | 0.7428 | 0.7428 | 0.0000 | identical |
| review_stats.shielding.3000.frac_gt_500 | 0.8659 | 0.8659 | 0.0000 | identical |
| review_stats.shielding.3000.median_g_cm2 | 1807.6500 | 1804.6500 | -3.0000 | rounding-safe |
| review_stats.stats.ci_block.20[0] | 3.5100 | 3.5000 | -0.0100 | rounding-safe |
| review_stats.stats.ci_block.20[1] | 8.3260 | 8.3160 | -0.0100 | rounding-safe |
| review_stats.stats.ci_block.30[0] | 3.1400 | 3.1400 | 0.0000 | identical |
| review_stats.stats.ci_block.30[1] | 8.6555 | 8.6505 | -0.0050 | rounding-safe |
| review_stats.stats.ci_block.50[0] | 2.4390 | 2.4290 | -0.0100 | rounding-safe |
| review_stats.stats.ci_block.50[1] | 8.9805 | 8.9755 | -0.0050 | rounding-safe |
| review_stats.stats.ci_iid[0] | 5.7305 | 5.7255 | -0.0050 | rounding-safe |
| review_stats.stats.ci_iid[1] | 6.4015 | 6.3920 | -0.0095 | rounding-safe |
| review_stats.stats.corr_len_1e_m | 50 | 50 | 0 | identical |
| review_stats.stats.corr_len_zero_m | 102 | 102 | 0 | identical |
| review_stats.stats.median | 6.0255 | 6.0155 | -0.0100 | explained: 1 cm quantisation of the delivered cloud |
| review_stats.stats.n_eff | 3.3673 | 3.3678 | 0.0005 | rounding-safe |
| review_stats.stats.sigma_tau_systematic_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
| review_stats.terrain.ceiling_thinning.ci95[0] | -0.1699 | -0.1490 | 0.0209 | explained: aperture bookkeeping |
| review_stats.terrain.ceiling_thinning.ci95[1] | -0.1093 | -0.0933 | 0.0159 | explained: aperture bookkeeping |
| review_stats.terrain.ceiling_thinning.n | 108 | 108 | 0 | identical |
| review_stats.terrain.ceiling_thinning.slope_m_per_m | -0.1345 | -0.1164 | 0.0181 | explained: aperture bookkeeping |
| review_stats.terrain.const_ceiling_control.r | 0.9054 | 0.9053 | -0.0001 | rounding-safe |
| review_stats.terrain.const_ceiling_control.r2 | 0.8155 | 0.8153 | -0.0002 | rounding-safe |
| review_stats.terrain.const_ceiling_control.rms_residual_m | 1.6953 | 1.6953 | -0.0000 | rounding-safe |
| review_stats.terrain.corr_tau_zdem | 0.9054 | 0.9053 | -0.0001 | rounding-safe |
| review_stats.terrain.dem_share | 0.7047 | 0.7045 | -0.0002 | rounding-safe |
| review_stats.terrain.near_skylight_median.10to20 | 2.2400 | 2.6655 | 0.4255 | explained: aperture bookkeeping |
| review_stats.terrain.near_skylight_median.20to30 | 4.6300 | 4.1080 | -0.5220 | explained: aperture bookkeeping |
| review_stats.terrain.near_skylight_median.lt10 | 1.3560 | 1.3560 | 0.0000 | identical |
| review_stats.terrain.net_centreline.end | 234.2190 | 168.2190 | -66.0000 | datum (dz) |
| review_stats.terrain.net_centreline.start | 232.5970 | 166.5970 | -66.0000 | datum (dz) |
| review_stats.terrain.reach_median.0-50 | 1.8355 | 1.8305 | -0.0050 | rounding-safe |
| review_stats.terrain.reach_median.100-150 | 4.6860 | 4.6710 | -0.0150 | rounding-safe |
| review_stats.terrain.reach_median.150-200 | 6.4750 | 6.4700 | -0.0050 | rounding-safe |
| review_stats.terrain.reach_median.200-250 | 7.9910 | 7.9810 | -0.0100 | rounding-safe |
| review_stats.terrain.reach_median.250-300 | 12.9425 | 12.9325 | -0.0100 | rounding-safe |
| review_stats.terrain.reach_median.50-100 | 1.5440 | 1.5390 | -0.0050 | rounding-safe |
| review_stats.terrain.skylight_runs[0].d_tau | -5.6004 | -5.5566 | 0.0438 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[0].d_zceil | 0.3957 | 0.2732 | -0.1224 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[0].d_zdem | -5.2048 | -5.2834 | -0.0786 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[0].run_s[0] | 29.0000 | 29.0000 | 0.0000 | identical |
| review_stats.terrain.skylight_runs[0].run_s[1] | 36.0000 | 30.0000 | -6.0000 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[1].d_tau | -7.2389 | -5.9394 | 1.2996 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[1].d_zceil | 2.2604 | 0.8009 | -1.4594 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[1].d_zdem | -4.9786 | -5.1385 | -0.1599 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[1].run_s[0] | 54.0000 | 35.0000 | -19.0000 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[1].run_s[1] | 56.0000 | 36.0000 | -20.0000 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[2].d_tau | -5.8345 | -7.2355 | -1.4010 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[2].d_zceil | 2.2244 | 2.2604 | 0.0360 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[2].d_zdem | -3.6101 | -4.9751 | -1.3650 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[2].run_s[0] | 89.0000 | 54.0000 | -35.0000 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[2].run_s[1] | 99.0000 | 56.0000 | -43.0000 | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[3].d_tau | None | -5.8329 |  | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[3].d_zceil | None | 2.2244 |  | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[3].d_zdem | None | -3.6085 |  | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[3].run_s[0] | None | 89.0000 |  | explained: aperture bookkeeping |
| review_stats.terrain.skylight_runs[3].run_s[1] | None | 99.0000 |  | explained: aperture bookkeeping |
| review_stats.terrain.tau_hist_2m.counts[0] | 47 | 47 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[1] | 30 | 30 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[2] | 59 | 59 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[3] | 44 | 44 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[4] | 34 | 34 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[5] | 22 | 22 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[6] | 38 | 38 | 0 | identical |
| review_stats.terrain.tau_hist_2m.counts[7] | 2 | 2 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[0] | 0 | 0 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[1] | 2 | 2 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[2] | 4 | 4 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[3] | 6 | 6 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[4] | 8 | 8 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[5] | 10 | 10 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[6] | 12 | 12 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[7] | 14 | 14 | 0 | identical |
| review_stats.terrain.tau_hist_2m.edges[8] | 16 | 16 | 0 | identical |
| review_stats.terrain.variance.tau | 15.5763 | 15.5609 | -0.0155 | rounding-safe |
| review_stats.terrain.variance.z_ceil | 2.8742 | 2.8742 | 0.0000 | identical |
| review_stats.terrain.variance.z_dem | 10.9771 | 10.9631 | -0.0140 | rounding-safe |
| review_stats.terrain.zceil_range | 10.5590 | 10.5590 | 0.0000 | identical |
| review_stats.terrain.zdem_range | 12.6600 | 12.6500 | -0.0100 | rounding-safe |
| review_stats.volume.n_sections | 302 | 302 | 0 | identical |
| review_stats.volume.per_m | 46.5780 | 46.5780 | 0.0000 | identical |
| review_stats.volume.total_m3 | 14066.6 | 14066.6 | 0.0000 | identical |
| roof_summary.kappa_env | 0.0034 | 0.0034 | 0.0000 | identical |
| roof_summary.kappa_env_station.s | 51.0000 | 51.0000 | 0.0000 | identical |
| roof_summary.kappa_env_station.span | 12.2100 | 12.2100 | 0.0000 | identical |
| roof_summary.kappa_env_station.tau | 0.0400 | 0.0400 | 0.0000 | identical |
| roof_summary.n_inconsistent | 4 | 4 | 0 | identical |
| roof_summary.n_intact | 276 | 276 | 0 | identical |
| roof_summary.n_low_coverage | 0 | 0 | 0 | identical |
| roof_summary.n_multipass | 0 | 0 | 0 | identical |
| roof_summary.n_no_dem | 0 | 4 | 4 | explained: aperture bookkeeping |
| roof_summary.n_skylight | 22 | 18 | -4 | explained: aperture bookkeeping |
| roof_summary.n_stations | 302 | 302 | 0 | identical |
| roof_summary.rho_basalt | 2600.0000 | 2600.0000 | 0.0000 | identical |
| roof_summary.shielding_g_cm2.atmosphere_ratio_median | 1.5200 | 1.5100 | -0.0100 | explained: 1 cm quantisation of the delivered cloud |
| roof_summary.shielding_g_cm2.median | 1567.0000 | 1564.0000 | -3.0000 | explained: 1 cm quantisation of the delivered cloud |
| roof_summary.shielding_g_cm2.min | 11.0000 | 11.0000 | 0.0000 | identical |
| roof_summary.tau_m.max | 14.1500 | 14.1400 | -0.0100 | rounding-safe |
| roof_summary.tau_m.median | 6.0300 | 6.0200 | -0.0100 | rounding-safe |
| roof_summary.tau_m.median_ci95[0] | 5.7400 | 5.7300 | -0.0100 | rounding-safe |
| roof_summary.tau_m.median_ci95[1] | 6.4000 | 6.3900 | -0.0100 | rounding-safe |
| roof_summary.tau_m.min | 0.0400 | 0.0400 | 0.0000 | identical |
| roof_summary.tau_m.q25 | 3.2400 | 3.2400 | 0.0000 | identical |
| roof_summary.tau_m.q75 | 9.0100 | 9.0000 | -0.0100 | rounding-safe |
| roof_summary.tau_m.sigma | 1.4900 | 1.4800 | -0.0100 | rounding-safe |
| morphometry_summary.area_m2.max | 81.4000 | 81.4000 | 0.0000 | identical |
| morphometry_summary.area_m2.mean_boot_ci95[0] | 45.2000 | 45.2000 | 0.0000 | identical |
| morphometry_summary.area_m2.mean_boot_ci95[1] | 49.2000 | 49.2000 | 0.0000 | identical |
| morphometry_summary.area_m2.median | 48.0000 | 48.0000 | 0.0000 | identical |
| morphometry_summary.area_m2.min | 5.1000 | 5.1000 | 0.0000 | identical |
| morphometry_summary.area_m2.p5 | 13.4000 | 13.4000 | 0.0000 | identical |
| morphometry_summary.area_m2.p95 | 75.8000 | 75.8000 | 0.0000 | identical |
| morphometry_summary.area_m2.q25 | 36.0000 | 36.0000 | 0.0000 | identical |
| morphometry_summary.area_m2.q75 | 59.2000 | 59.2000 | 0.0000 | identical |
| morphometry_summary.eta.median | 1.5300 | 1.5300 | 0.0000 | identical |
| morphometry_summary.eta.q25 | 1.3100 | 1.3100 | 0.0000 | identical |
| morphometry_summary.eta.q75 | 1.8400 | 1.8400 | 0.0000 | identical |
| morphometry_summary.height_m.max | 11.4100 | 11.4100 | 0.0000 | identical |
| morphometry_summary.height_m.median | 8.0900 | 8.0900 | 0.0000 | identical |
| morphometry_summary.height_m.min | 3.5300 | 3.5300 | 0.0000 | identical |
| morphometry_summary.min_coverage | 0.5000 | 0.5000 | 0.0000 | identical |
| morphometry_summary.n_stations | 302 | 302 | 0 | identical |
| morphometry_summary.n_valid_sections | 296 | 296 | 0 | identical |
| morphometry_summary.sinuosity_50m.max | 1.0790 | 1.0790 | 0.0000 | identical |
| morphometry_summary.sinuosity_50m.mean | 1.0320 | 1.0320 | 0.0000 | identical |
| morphometry_summary.total_centreline_m | 301.0000 | 301.0000 | 0.0000 | identical |
| morphometry_summary.width_m.max | 31.6000 | 31.6000 | 0.0000 | identical |
| morphometry_summary.width_m.median | 10.5300 | 10.5300 | 0.0000 | identical |
| morphometry_summary.width_m.min | 6.0700 | 6.0700 | 0.0000 | identical |
| morphometry_summary.width_min_m | 6.0700 | 6.0700 | 0.0000 | identical |
| registration_validation.floor_through_skylight_s3.annulus_outside_cone.median_m | 0.9072 | 0.8983 | -0.0088 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.annulus_outside_cone.n | 321 | 327 | 6 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.annulus_outside_cone.p10_m | -0.0687 | -0.0789 | -0.0102 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.annulus_outside_cone.p90_m | 1.3719 | 1.3727 | 0.0009 | rounding-safe |
| registration_validation.floor_through_skylight_s3.annulus_outside_cone.rms_m | 1.3603 | 1.3516 | -0.0087 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.median_m | 0.8869 | 0.8738 | -0.0131 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.n | 352 | 360 | 8 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.p10_m | -1.1877 | -1.2324 | -0.0447 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.p90_m | 1.3423 | 1.3372 | -0.0051 | explained: S3 floor check on new voxel copy |
| registration_validation.floor_through_skylight_s3.rms_m | 1.3561 | 1.3469 | -0.0092 | explained: S3 floor check on new voxel copy |
| registration_validation.n_negative_tau_stations | 4 | 4 | 0 | identical |
| registration_validation.rim_horizontal_rms_m | 1.7331 | 1.7331 | 0.0000 | identical |
| registration_validation.sigma_reg_components.anchor_floor_rms_m | 1.3561 | 1.3469 | -0.0092 | explained: S3 floor check on new voxel copy |
| registration_validation.sigma_reg_components.drift_rss_m | 0.5318 | 0.5318 | -0.0000 | rounding-safe |
| registration_validation.sigma_reg_z_m | 1.4600 | 1.4500 | -0.0100 | explained: S3 floor check on new voxel copy |
| registration_validation.two_slam_zones.s0-80.mean | -0.2324 | -0.2324 | 0.0000 | identical |
| registration_validation.two_slam_zones.s0-80.n | 80 | 80 | 0 | identical |
| registration_validation.two_slam_zones.s0-80.std | 0.3220 | 0.3220 | 0.0000 | identical |
| registration_validation.two_slam_zones.s160-240.mean | -0.1664 | -0.1664 | -0.0000 | rounding-safe |
| registration_validation.two_slam_zones.s160-240.n | 80 | 80 | 0 | identical |
| registration_validation.two_slam_zones.s160-240.std | 0.4783 | 0.4783 | -0.0000 | rounding-safe |
| registration_validation.two_slam_zones.s240-340.mean | -0.1012 | -0.1012 | -0.0000 | rounding-safe |
| registration_validation.two_slam_zones.s240-340.n | 99 | 99 | 0 | identical |
| registration_validation.two_slam_zones.s240-340.std | 0.1277 | 0.1277 | 0.0000 | rounding-safe |
| registration_validation.two_slam_zones.s80-160.mean | -0.1352 | -0.1352 | -0.0000 | rounding-safe |
| registration_validation.two_slam_zones.s80-160.n | 80 | 80 | 0 | identical |
| registration_validation.two_slam_zones.s80-160.std | 0.1281 | 0.1281 | 0.0000 | rounding-safe |
| flight2_repeatability.median | 0.0809 | 0.0809 | 0.0000 | rounding-safe |
| flight2_repeatability.n_stations | 302 | 302 | 0 | identical |
| flight2_repeatability.p5 | -0.2060 | -0.2059 | 0.0000 | rounding-safe |
| flight2_repeatability.p95 | 0.7125 | 0.7125 | -0.0000 | rounding-safe |
| flight2_repeatability.rms | 0.3707 | 0.3707 | -0.0000 | rounding-safe |
| flight2_repeatability.s_range[0] | 0.0000 | 0.0000 | 0.0000 | identical |
| flight2_repeatability.s_range[1] | 301.0000 | 301.0000 | 0.0000 | identical |
| basis_comparison.dceil.median | 0.2050 | -65.7950 | -66.0000 | datum-derived |
| basis_comparison.dceil.rms | 0.8652 | 65.5968 | 64.7315 | datum-derived |
| basis_comparison.dtau.median | -0.2320 | -0.2400 | -0.0080 | explained: 1 cm quantisation of the delivered cloud |
| basis_comparison.dtau.p5 | -1.4166 | -1.4220 | -0.0054 | explained: 1 cm quantisation of the delivered cloud |
| basis_comparison.dtau.p95 | 0.1734 | 0.1734 | 0.0000 | identical |
| basis_comparison.dtau.rms | 0.8547 | 0.8585 | 0.0038 | rounding-safe |
| basis_comparison.dtau_zones.s0-80.mean | -0.3382 | -0.3427 | -0.0045 | rounding-safe |
| basis_comparison.dtau_zones.s0-80.n | 63 | 63 | 0 | identical |
| basis_comparison.dtau_zones.s0-80.std | 0.4873 | 0.4876 | 0.0003 | rounding-safe |
| basis_comparison.dtau_zones.s160-240.mean | -0.5547 | -0.5627 | -0.0080 | explained: 1 cm quantisation of the delivered cloud |
| basis_comparison.dtau_zones.s160-240.n | 80 | 80 | 0 | identical |
| basis_comparison.dtau_zones.s160-240.std | 0.8085 | 0.8084 | -0.0001 | rounding-safe |
| basis_comparison.dtau_zones.s240-360.mean | -0.3961 | -0.4067 | -0.0106 | explained: 1 cm quantisation of the delivered cloud |
| basis_comparison.dtau_zones.s240-360.n | 62 | 62 | 0 | identical |
| basis_comparison.dtau_zones.s240-360.std | 0.7376 | 0.7372 | -0.0004 | rounding-safe |
| basis_comparison.dtau_zones.s80-160.mean | -0.3421 | -0.3492 | -0.0071 | explained: 1 cm quantisation of the delivered cloud |
| basis_comparison.dtau_zones.s80-160.n | 69 | 69 | 0 | identical |
| basis_comparison.dtau_zones.s80-160.std | 0.8475 | 0.8478 | 0.0003 | rounding-safe |
| basis_comparison.median_horiz_offset_m | 0.4723 | 0.4723 | 0.0000 | identical |
| basis_comparison.morphometry.area_median[0] | 48.9000 | 48.9000 | 0.0000 | identical |
| basis_comparison.morphometry.area_median[1] | 48.0000 | 48.0000 | 0.0000 | identical |
| basis_comparison.morphometry.centreline_m[0] | 356.0000 | 356.0000 | 0.0000 | identical |
| basis_comparison.morphometry.centreline_m[1] | 301.0000 | 301.0000 | 0.0000 | identical |
| basis_comparison.morphometry.height_median[0] | 8.0500 | 8.0500 | 0.0000 | identical |
| basis_comparison.morphometry.height_median[1] | 8.0900 | 8.0900 | 0.0000 | identical |
| basis_comparison.morphometry.width_median[0] | 10.4600 | 10.4600 | 0.0000 | identical |
| basis_comparison.morphometry.width_median[1] | 10.5300 | 10.5300 | 0.0000 | identical |
| basis_comparison.n_matched | 274 | 274 | 0 | identical |
| national_dem_check.arcticdem.best_shift_m[0] | 2.0000 | 2.0000 | 0.0000 | identical |
| national_dem_check.arcticdem.best_shift_m[1] | -2.0000 | -2.0000 | 0.0000 | identical |
| national_dem_check.arcticdem.datum_offset_m | 0.6100 | -65.3966 | -66.0066 | datum-derived |
| national_dem_check.arcticdem.doming_amplitude_m | 0.1280 | 0.1291 | 0.0011 | rounding-safe |
| national_dem_check.arcticdem.n | 30559 | 30550 | -9 | explained: in-aperture cells removed |
| national_dem_check.arcticdem.std_after_plane_m | 0.2895 | 0.2780 | -0.0115 | explained: in-aperture cells removed |
| national_dem_check.arcticdem.std_after_quadratic_m | 0.2877 | 0.2760 | -0.0117 | explained: in-aperture cells removed |
| national_dem_check.arcticdem.std_raw_m | 0.2926 | 0.2814 | -0.0112 | explained: in-aperture cells removed |
| national_dem_check.arcticdem.tilt_m_per_km | 0.4214 | 0.4294 | 0.0080 | explained: in-aperture cells removed |
| national_dem_check.islandsdem.best_shift_m[0] | 0.0000 | 0.0000 | 0.0000 | identical |
| national_dem_check.islandsdem.best_shift_m[1] | 0.0000 | 0.0000 | 0.0000 | identical |
| national_dem_check.islandsdem.datum_offset_m | 66.6325 | 0.6269 | -66.0056 | datum-derived |
| national_dem_check.islandsdem.doming_amplitude_m | 0.1332 | 0.1337 | 0.0005 | rounding-safe |
| national_dem_check.islandsdem.n | 30559 | 30550 | -9 | explained: in-aperture cells removed |
| national_dem_check.islandsdem.std_after_plane_m | 0.3585 | 0.3498 | -0.0087 | explained: in-aperture cells removed |
| national_dem_check.islandsdem.std_after_quadratic_m | 0.3569 | 0.3482 | -0.0087 | explained: in-aperture cells removed |
| national_dem_check.islandsdem.std_raw_m | 0.3635 | 0.3555 | -0.0080 | explained: in-aperture cells removed |
| national_dem_check.islandsdem.tilt_m_per_km | 0.5487 | 0.5752 | 0.0265 | explained: in-aperture cells removed |
| orbital_dem_test.baseline.median | 6.0255 | 6.0155 | -0.0100 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.baseline.thick_thin_contrast_m | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[0].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[0].frac_error_gt_1p5.mean | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[0].frac_error_gt_1p5.worst | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[0].frac_false_negative_roof | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[0].grid_m | 0.5000 | 0.5000 | 0.0000 | identical |
| orbital_dem_test.configs[0].median_shift_m.mean | -0.0000 | -0.0000 | -0.0000 | rounding-safe |
| orbital_dem_test.configs[0].median_shift_m.worst | 0.0000 | 0.0000 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[0].n_seeds | 1 | 1 | 0 | identical |
| orbital_dem_test.configs[0].per_station_rms_m.mean | 0.0002 | 0.0002 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[0].per_station_rms_m.worst | 0.0002 | 0.0002 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[0].sigma_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[0].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[0].thick_thin_contrast_m.mean | 9.9515 | 9.9465 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[0].thick_thin_contrast_m.worst | 9.9515 | 9.9465 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[1].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[1].frac_error_gt_1p5.mean | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[1].frac_error_gt_1p5.worst | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[1].frac_false_negative_roof | 0.0217 | 0.0217 | 0.0000 | identical |
| orbital_dem_test.configs[1].grid_m | 1.0000 | 1.0000 | 0.0000 | identical |
| orbital_dem_test.configs[1].median_shift_m.mean | 0.0879 | 0.0906 | 0.0026 | rounding-safe |
| orbital_dem_test.configs[1].median_shift_m.worst | 0.4289 | 0.4289 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[1].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[1].per_station_rms_m.mean | 0.2884 | 0.2883 | -0.0000 | rounding-safe |
| orbital_dem_test.configs[1].per_station_rms_m.worst | 0.4183 | 0.4186 | 0.0003 | rounding-safe |
| orbital_dem_test.configs[1].sigma_m | 0.3000 | 0.3000 | 0.0000 | identical |
| orbital_dem_test.configs[1].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[1].thick_thin_contrast_m.mean | 9.9676 | 9.9643 | -0.0033 | rounding-safe |
| orbital_dem_test.configs[1].thick_thin_contrast_m.worst | 9.2341 | 9.2307 | -0.0035 | rounding-safe |
| orbital_dem_test.configs[2].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[2].frac_error_gt_1p5.mean | 0.0016 | 0.0013 | -0.0004 | rounding-safe |
| orbital_dem_test.configs[2].frac_error_gt_1p5.worst | 0.0145 | 0.0145 | 0.0000 | identical |
| orbital_dem_test.configs[2].frac_false_negative_roof | 0.0326 | 0.0326 | 0.0000 | identical |
| orbital_dem_test.configs[2].grid_m | 2.0000 | 2.0000 | 0.0000 | identical |
| orbital_dem_test.configs[2].median_shift_m.mean | 0.1360 | 0.1389 | 0.0029 | rounding-safe |
| orbital_dem_test.configs[2].median_shift_m.worst | 0.6900 | 0.6928 | 0.0028 | rounding-safe |
| orbital_dem_test.configs[2].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[2].per_station_rms_m.mean | 0.4835 | 0.4826 | -0.0009 | rounding-safe |
| orbital_dem_test.configs[2].per_station_rms_m.worst | 0.6982 | 0.6982 | -0.0000 | rounding-safe |
| orbital_dem_test.configs[2].sigma_m | 0.5000 | 0.5000 | 0.0000 | identical |
| orbital_dem_test.configs[2].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[2].thick_thin_contrast_m.mean | 10.0242 | 10.0200 | -0.0042 | rounding-safe |
| orbital_dem_test.configs[2].thick_thin_contrast_m.worst | 8.8375 | 8.8363 | -0.0012 | rounding-safe |
| orbital_dem_test.configs[3].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[3].frac_error_gt_1p5.mean | 0.1395 | 0.1386 | -0.0009 | rounding-safe |
| orbital_dem_test.configs[3].frac_error_gt_1p5.worst | 0.3333 | 0.3333 | 0.0000 | identical |
| orbital_dem_test.configs[3].frac_false_negative_roof | 0.1486 | 0.1449 | -0.0036 | rounding-safe |
| orbital_dem_test.configs[3].grid_m | 5.0000 | 5.0000 | 0.0000 | identical |
| orbital_dem_test.configs[3].median_shift_m.mean | 0.2879 | 0.2910 | 0.0030 | rounding-safe |
| orbital_dem_test.configs[3].median_shift_m.worst | 1.5767 | 1.5782 | 0.0015 | rounding-safe |
| orbital_dem_test.configs[3].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[3].per_station_rms_m.mean | 0.9693 | 0.9662 | -0.0031 | rounding-safe |
| orbital_dem_test.configs[3].per_station_rms_m.worst | 1.4033 | 1.4016 | -0.0017 | rounding-safe |
| orbital_dem_test.configs[3].sigma_m | 1.0000 | 1.0000 | 0.0000 | identical |
| orbital_dem_test.configs[3].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[3].thick_thin_contrast_m.mean | 10.0219 | 10.0093 | -0.0126 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[3].thick_thin_contrast_m.worst | 7.4294 | 7.4221 | -0.0073 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[4].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[4].frac_error_gt_1p5.mean | 0.4589 | 0.4578 | -0.0011 | rounding-safe |
| orbital_dem_test.configs[4].frac_error_gt_1p5.worst | 0.7899 | 0.7899 | 0.0000 | identical |
| orbital_dem_test.configs[4].frac_false_negative_roof | 0.2500 | 0.2500 | 0.0000 | identical |
| orbital_dem_test.configs[4].grid_m | 5.0000 | 5.0000 | 0.0000 | identical |
| orbital_dem_test.configs[4].median_shift_m.mean | 0.4878 | 0.4905 | 0.0027 | rounding-safe |
| orbital_dem_test.configs[4].median_shift_m.worst | 2.7585 | 2.7612 | 0.0027 | rounding-safe |
| orbital_dem_test.configs[4].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[4].per_station_rms_m.mean | 1.9172 | 1.9153 | -0.0018 | rounding-safe |
| orbital_dem_test.configs[4].per_station_rms_m.worst | 2.7937 | 2.7929 | -0.0009 | rounding-safe |
| orbital_dem_test.configs[4].sigma_m | 2.0000 | 2.0000 | 0.0000 | identical |
| orbital_dem_test.configs[4].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[4].thick_thin_contrast_m.mean | 9.6739 | 9.6591 | -0.0148 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[4].thick_thin_contrast_m.worst | 4.8412 | 4.8376 | -0.0036 | rounding-safe |
| orbital_dem_test.configs[5].dome_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[5].frac_error_gt_1p5.mean | 0.5996 | 0.5989 | -0.0007 | rounding-safe |
| orbital_dem_test.configs[5].frac_error_gt_1p5.worst | 0.8841 | 0.8841 | 0.0000 | identical |
| orbital_dem_test.configs[5].frac_false_negative_roof | 0.2826 | 0.2826 | 0.0000 | identical |
| orbital_dem_test.configs[5].grid_m | 20.0000 | 20.0000 | 0.0000 | identical |
| orbital_dem_test.configs[5].median_shift_m.mean | 0.5898 | 0.5980 | 0.0083 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[5].median_shift_m.worst | 3.9250 | 3.9276 | 0.0025 | rounding-safe |
| orbital_dem_test.configs[5].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[5].per_station_rms_m.mean | 2.9027 | 2.9025 | -0.0003 | rounding-safe |
| orbital_dem_test.configs[5].per_station_rms_m.worst | 4.2004 | 4.2004 | -0.0000 | rounding-safe |
| orbital_dem_test.configs[5].sigma_m | 3.0000 | 3.0000 | 0.0000 | identical |
| orbital_dem_test.configs[5].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[5].thick_thin_contrast_m.mean | 9.1844 | 9.1581 | -0.0263 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[5].thick_thin_contrast_m.worst | 1.9619 | 1.9281 | -0.0337 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[6].dome_m | 1.0000 | 1.0000 | 0.0000 | identical |
| orbital_dem_test.configs[6].frac_error_gt_1p5.mean | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[6].frac_error_gt_1p5.worst | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[6].frac_false_negative_roof | 0.0217 | 0.0217 | 0.0000 | identical |
| orbital_dem_test.configs[6].grid_m | 0.5000 | 0.5000 | 0.0000 | identical |
| orbital_dem_test.configs[6].median_shift_m.mean | -0.4806 | -0.4806 | -0.0000 | rounding-safe |
| orbital_dem_test.configs[6].median_shift_m.worst | 0.4806 | 0.4806 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[6].n_seeds | 1 | 1 | 0 | identical |
| orbital_dem_test.configs[6].per_station_rms_m.mean | 0.4016 | 0.4016 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[6].per_station_rms_m.worst | 0.4016 | 0.4016 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[6].sigma_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[6].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[6].thick_thin_contrast_m.mean | 10.2267 | 10.2217 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[6].thick_thin_contrast_m.worst | 10.2267 | 10.2217 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[7].dome_m | 2.0000 | 2.0000 | 0.0000 | identical |
| orbital_dem_test.configs[7].frac_error_gt_1p5.mean | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[7].frac_error_gt_1p5.worst | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[7].frac_false_negative_roof | 0.0326 | 0.0326 | 0.0000 | identical |
| orbital_dem_test.configs[7].grid_m | 0.5000 | 0.5000 | 0.0000 | identical |
| orbital_dem_test.configs[7].median_shift_m.mean | -0.8819 | -0.8769 | 0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[7].median_shift_m.worst | 0.8819 | 0.8769 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[7].n_seeds | 1 | 1 | 0 | identical |
| orbital_dem_test.configs[7].per_station_rms_m.mean | 0.8032 | 0.8032 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[7].per_station_rms_m.worst | 0.8032 | 0.8032 | 0.0000 | rounding-safe |
| orbital_dem_test.configs[7].sigma_m | 0.0000 | 0.0000 | 0.0000 | identical |
| orbital_dem_test.configs[7].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[7].thick_thin_contrast_m.mean | 10.4529 | 10.4479 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[7].thick_thin_contrast_m.worst | 10.4529 | 10.4479 | -0.0050 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[8].dome_m | 2.0000 | 2.0000 | 0.0000 | identical |
| orbital_dem_test.configs[8].frac_error_gt_1p5.mean | 0.2375 | 0.2371 | -0.0004 | rounding-safe |
| orbital_dem_test.configs[8].frac_error_gt_1p5.worst | 0.6014 | 0.6051 | 0.0036 | rounding-safe |
| orbital_dem_test.configs[8].frac_false_negative_roof | 0.2029 | 0.2029 | 0.0000 | identical |
| orbital_dem_test.configs[8].grid_m | 5.0000 | 5.0000 | 0.0000 | identical |
| orbital_dem_test.configs[8].median_shift_m.mean | -0.6150 | -0.6120 | 0.0030 | rounding-safe |
| orbital_dem_test.configs[8].median_shift_m.worst | 1.4177 | 1.4144 | -0.0033 | rounding-safe |
| orbital_dem_test.configs[8].n_seeds | 20 | 20 | 0 | identical |
| orbital_dem_test.configs[8].per_station_rms_m.mean | 1.2139 | 1.2085 | -0.0054 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[8].per_station_rms_m.worst | 1.9377 | 1.9312 | -0.0066 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[8].sigma_m | 1.0000 | 1.0000 | 0.0000 | identical |
| orbital_dem_test.configs[8].thick_thin_contrast_m.baseline | 9.9515 | 9.9465 | -0.0050 | rounding-safe |
| orbital_dem_test.configs[8].thick_thin_contrast_m.mean | 10.4625 | 10.4556 | -0.0069 | explained: 1 cm quantisation of the delivered cloud |
| orbital_dem_test.configs[8].thick_thin_contrast_m.worst | 7.8999 | 7.8955 | -0.0043 | rounding-safe |
| orbital_dem_test.corr_len_m | 25.0000 | 25.0000 | 0.0000 | identical |
| orbital_dem_test.n_intact | 276 | 276 | 0 | identical |
| orbital_dem_test.planetary_aperture_context.apc_outer_median_m | 120.0000 | 120.0000 | 0.0000 | identical |
| orbital_dem_test.planetary_aperture_context.mare_pit_inner_median_m | 100.0000 | 100.0000 | 0.0000 | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.1.0.S1 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.1.0.S2 | True | True |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.1.0.S3 | True | True |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.2.0.S1 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.2.0.S2 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.2.0.S3 | True | True |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.20.0.S1 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.20.0.S2 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.20.0.S3 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.5.0.S1 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.5.0.S2 | False | False |  | identical |
| orbital_dem_test.skylight_anchor_resolvable_3post.5.0.S3 | False | False |  | identical |
| planetary_summary.demonstrated_min_span_m | 6.0700 | 6.0700 | 0.0000 | identical |
| planetary_summary.kappa_env_earth | 0.0034 | 0.0034 | 0.0000 | identical |
| planetary_summary.mars.aperture_median_m | 120.0000 | 120.0000 | 0.0000 | identical |
| planetary_summary.mars.counts.ge_10m | 131 | 131 | 0 | identical |
| planetary_summary.mars.counts.ge_20m | 131 | 131 | 0 | identical |
| planetary_summary.mars.counts.ge_50m | 120 | 120 | 0 | identical |
| planetary_summary.mars.counts.ge_5m | 131 | 131 | 0 | identical |
| planetary_summary.mars.counts.ge_narrowest | 131 | 131 | 0 | identical |
| planetary_summary.mars.n_dimensioned | 131 | 131 | 0 | identical |
| planetary_summary.mars.n_total | 1062 | 1062 | 0 | identical |
| planetary_summary.moon.aperture_median_m | 29.0000 | 29.0000 | 0.0000 | identical |
| planetary_summary.moon.counts.ge_10m | 258 | 258 | 0 | identical |
| planetary_summary.moon.counts.ge_20m | 190 | 190 | 0 | identical |
| planetary_summary.moon.counts.ge_50m | 63 | 63 | 0 | identical |
| planetary_summary.moon.counts.ge_5m | 266 | 266 | 0 | identical |
| planetary_summary.moon.counts.ge_narrowest | 265 | 265 | 0 | identical |
| planetary_summary.moon.n_dimensioned | 266 | 266 | 0 | identical |
| planetary_summary.moon.n_total | 278 | 278 | 0 | identical |
| planetary_summary.span_scale.mars | 1.6300 | 1.6300 | 0.0000 | identical |
| planetary_summary.span_scale.moon | 2.4600 | 2.4600 | 0.0000 | identical |
| uncertainty_budget.horizontal_rim_rms_m | 1.7331 | 1.7331 | 0.0000 | identical |
| uncertainty_budget.sigma_dem_basis | MEASURED upper bound: std of (photogrammetric DEM minus ArcticDEM v4.1 2 m) over the survey footprint, 0.29 m, including ArcticDEM's own error and snow/temporal change (analysis/national_dem_check.py). Supersedes the assumed 2x-GSD value; doming bounded at 0.13 m, tilt 0.42 m/km. | MEASURED upper bound: std of (photogrammetric DEM minus ArcticDEM v4.1 2 m) over the survey footprint, 0.29 m, including ArcticDEM's own error and snow/temporal change (analysis/national_dem_check.py). Supersedes the assumed 2x-GSD value; doming bounded at 0.13 m, tilt 0.42 m/km. |  | identical |
| uncertainty_budget.sigma_dem_z_m | 0.2900 | 0.2900 | 0.0000 | identical |
| uncertainty_budget.sigma_lidar_basis | Ouster range precision (~3 cm) plus aggregated-map surface thickness (~10 cm at the 10 cm voxel working copy). | Ouster range precision (~3 cm) plus aggregated-map surface thickness (~10 cm at the 10 cm voxel working copy). |  | identical |
| uncertainty_budget.sigma_lidar_z_m | 0.1000 | 0.1000 | 0.0000 | identical |
| uncertainty_budget.sigma_reg_z_m | 1.4600 | 1.4500 | -0.0100 | explained: S3 floor check on new voxel copy |
| uncertainty_budget.sigma_tau_m | 1.4900 | 1.4800 | -0.0100 | explained: S3 floor check on new voxel copy |
