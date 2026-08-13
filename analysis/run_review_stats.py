"""Phase A review-response statistics (tri-role panel, 2026-08-13).

Implements the quantitative fixes the review panel called for, all from the
committed Ed-frame outputs, and writes analysis_out_ed/review_stats.json:

  A2  terrain decomposition: variance split between surface DEM and interior
      ceiling, constant-ceiling control, near-skylight thinning on the
      ceiling term alone, distribution shape (bimodality, reach-conditioned
      medians).
  A3  spatial statistics: tau autocorrelation length, Bartlett effective
      sample size, moving-block bootstrap CIs on the median; shielding as
      exceedance fractions at a defensible density range instead of a single
      atmosphere ratio.
  A1  stability: per-station required tensile strength
      sigma_req = beta rho g L^2 / tau (clamped-strip beta = 0.5), tau/L
      percentiles, sensitivity of the old single-point kappa, and the
      sqrt-model L_max(tau) curves used by the rebuilt planetary panel.
  A4  pit catalogue split by feature type with inner apertures.
  A5  scale-ruler forensics: inter-skylight distances in the surface model
      and in BOTH independent SLAM maps (rigid-invariant), plus per-view
      hole sizes, to separate map scale error from rim-centroid definition.
  +   enclosed volume per unit length from the section areas.

Run:  PYTHONPATH=src ANALYSIS_OUT=analysis_out_ed \
      .venv/bin/python analysis/run_review_stats.py
"""
from __future__ import annotations

import csv
import json
import os
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = Path(os.environ.get("ANALYSIS_OUT", str(ROOT / "analysis_out_ed")))
OUT_BASE = ROOT / "analysis_out"

G = dict(earth=9.81, mars=3.71, moon=1.62)
BETA = 0.5                  # clamped strip: sigma = 0.5 rho g L^2 / tau
RHOS = (2200.0, 2600.0, 3000.0)   # kg/m3: vesicular crust+soil .. dense basalt
RHO_MID = 2600.0
SIGMA_T_MPA = (1.0, 5.0, 10.0)    # rock-mass .. intact basalt tensile strength
RNG = np.random.default_rng(11)


def load_roof():
    rows = list(csv.DictReader(open(OUT / "roof_thickness.csv")))
    g = lambda k: np.array([float(r[k]) if r[k] else np.nan for r in rows])
    return (g("s"), g("x"), g("y"), g("z_dem"), g("z_ceil"), g("tau"),
            g("span"), np.array([r["class"] for r in rows]))


def sky_runs(s, klass):
    runs, start = [], None
    for i in range(len(s)):
        if klass[i] == "skylight" and start is None:
            start = s[i]
        elif klass[i] != "skylight" and start is not None:
            runs.append((start, s[i - 1])); start = None
    if start is not None:
        runs.append((start, s[-1]))
    return runs


def moving_block_ci(x, block, n=10_000):
    """95% CI on the median by moving-block bootstrap along the station order."""
    m = len(x)
    nblocks = int(np.ceil(m / block))
    starts = np.arange(0, m - block + 1)
    meds = np.empty(n)
    for i in range(n):
        idx = np.concatenate([np.arange(st, st + block)
                              for st in RNG.choice(starts, nblocks)])[:m]
        meds[i] = np.median(x[idx])
    return [float(np.percentile(meds, 2.5)), float(np.percentile(meds, 97.5))]


def main() -> None:
    s, x, y, zdem, zceil, tau, span, klass = load_roof()
    it = klass == "intact"
    si, taui, zdi, zci, Li = s[it], tau[it], zdem[it], zceil[it], span[it]
    out = {}

    # ---------------- A2: terrain decomposition ----------------
    var = dict(tau=float(np.var(taui)), z_dem=float(np.var(zdi)),
               z_ceil=float(np.var(zci)))
    tau_const_ceil = zdi - zci.mean()
    ss_res = float(np.sum((taui - tau_const_ceil) ** 2))
    ss_tot = float(np.sum((taui - taui.mean()) ** 2))
    const_ceiling = dict(
        r=float(np.corrcoef(taui, tau_const_ceil)[0, 1]),
        r2=1.0 - ss_res / ss_tot,
        rms_residual_m=float(np.sqrt(np.mean((taui - tau_const_ceil) ** 2))))
    runs = sky_runs(s, klass)
    d_sky = np.min(np.abs(si[:, None] -
                          np.array([(a + b) / 2 for a, b in runs])[None, :]),
                   axis=1)
    near, far = d_sky < 25, d_sky > 60
    per_run = []
    for a, b in runs:
        dn = np.abs(si - (a + b) / 2)
        nr, fr = dn < 25, dn > 60
        per_run.append(dict(run_s=[float(a), float(b)],
                            d_zdem=float(np.mean(zdi[nr]) - np.mean(zdi[fr])),
                            d_zceil=float(np.mean(zci[nr]) - np.mean(zci[fr])),
                            d_tau=float(np.mean(taui[nr]) - np.mean(taui[fr]))))
    # ceiling-only thinning: linear fit of z_ceil against distance to nearest
    # skylight within 40 m, with a station-bootstrap CI on the slope
    m40 = d_sky < 40
    if m40.sum() > 10:
        A = np.c_[d_sky[m40], np.ones(m40.sum())]
        sl = float(np.linalg.lstsq(A, zci[m40], rcond=None)[0][0])
        boots = []
        pos = np.where(m40)[0]
        for _ in range(2000):
            j = RNG.choice(pos, len(pos))
            boots.append(float(np.linalg.lstsq(
                np.c_[d_sky[j], np.ones(len(j))], zci[j], rcond=None)[0][0]))
        slope_ci = [float(np.percentile(boots, 2.5)),
                    float(np.percentile(boots, 97.5))]
        ceil_thinning = dict(slope_m_per_m=sl, ci95=slope_ci,
                             n=int(m40.sum()), note="z_ceil vs distance to "
                             "nearest skylight, stations within 40 m")
    else:
        ceil_thinning = None
    hist, edges = np.histogram(taui, bins=np.arange(0, 17, 2))
    reach_med = {f"{lo}-{lo+50}": float(np.median(taui[(si >= lo) & (si < lo + 50)]))
                 for lo in range(0, 350, 50)
                 if ((si >= lo) & (si < lo + 50)).sum() > 5}
    near_sky_med = {"lt10": float(np.median(taui[d_sky < 10])),
                    "10to20": float(np.median(taui[(d_sky >= 10) & (d_sky < 20)])),
                    "20to30": float(np.median(taui[(d_sky >= 20) & (d_sky < 30)]))}
    out["terrain"] = dict(variance=var,
                          dem_share=var["z_dem"] / var["tau"],
                          corr_tau_zdem=float(np.corrcoef(taui, zdi)[0, 1]),
                          const_ceiling_control=const_ceiling,
                          skylight_runs=per_run,
                          ceiling_thinning=ceil_thinning,
                          tau_hist_2m=dict(counts=hist.tolist(),
                                           edges=edges.tolist()),
                          reach_median=reach_med,
                          near_skylight_median=near_sky_med,
                          zceil_range=float(np.ptp(zci)),
                          zdem_range=float(np.ptp(zdi)),
                          net_centreline=dict(start=float(zci[0]), end=float(zci[-1])))

    # ---------------- A3: spatial statistics ----------------
    t0 = taui - taui.mean()
    ac = np.correlate(t0, t0, "full")[len(t0) - 1:]
    ac = ac / ac[0]
    e_len = int(np.argmax(ac < 1 / np.e))
    zero_len = int(np.argmax(ac < 0))
    n_eff = len(taui) / (1 + 2 * float(np.sum(ac[1:zero_len])))
    iid = [float(v) for v in np.percentile(
        [np.median(RNG.choice(taui, len(taui))) for _ in range(10_000)],
        [2.5, 97.5])]
    blocks = {b: moving_block_ci(taui, b) for b in (20, 30, 50)}
    out["stats"] = dict(corr_len_1e_m=e_len, corr_len_zero_m=zero_len,
                        n_eff=float(n_eff), median=float(np.median(taui)),
                        ci_iid=iid, ci_block={str(k): v for k, v in blocks.items()},
                        sigma_tau_systematic_m=1.49,
                        note="sigma_tau is dominated by the common registration "
                             "datum term and does not average down")

    # shielding exceedance at a density range
    shield = {}
    for rho in RHOS:
        gcm2 = taui * rho / 10.0
        shield[str(int(rho))] = dict(
            median_g_cm2=float(np.median(gcm2)),
            frac_gt_100=float(np.mean(gcm2 > 100)),
            frac_gt_500=float(np.mean(gcm2 > 500)),
            frac_gt_1000=float(np.mean(gcm2 > 1000)))
    out["shielding"] = shield

    # ---------------- A1: strength inversion + envelope ----------------
    ok = taui > 0.5   # exclude sub-uncertainty roofs from the strength stats
    sig_req = BETA * RHO_MID * 9.81 * Li[ok] ** 2 / taui[ok]  # Pa
    ratio = taui / Li
    kappa_sorted = np.sort(ratio)
    out["envelope"] = dict(
        beta=BETA, rho_mid=RHO_MID,
        sigma_req_MPa=dict(median=float(np.median(sig_req) / 1e6),
                           p95=float(np.percentile(sig_req, 95) / 1e6),
                           max=float(sig_req.max() / 1e6),
                           n=int(ok.sum()),
                           note="required tensile strength for stability of "
                                "each intact station under the clamped-strip "
                                "model; the max is the lower bound the intact "
                                "roof places on rock-mass strength"),
        tau_over_L=dict(min=float(kappa_sorted[0]),
                        second=float(kappa_sorted[1]),
                        p5=float(np.percentile(ratio, 5)),
                        p25=float(np.percentile(ratio, 25)),
                        median=float(np.median(ratio))),
        kappa_sensitivity=dict(
            single_min=float(kappa_sorted[0]),
            drop_one=float(kappa_sorted[1]),
            p5=float(np.percentile(ratio, 5)),
            median=float(np.median(ratio))),
        Lmax_sqrt_model={
            body: {f"{st:.0f}MPa": dict(
                at_tau_6m=float(np.sqrt(st * 1e6 * 6.2 / (BETA * RHO_MID * g))),
                at_tau_14m=float(np.sqrt(st * 1e6 * 14.0 / (BETA * RHO_MID * g))))
                for st in SIGMA_T_MPA}
            for body, g in G.items()},
        span_scale=dict(mars=float(np.sqrt(G["earth"] / G["mars"])),
                        moon=float(np.sqrt(G["earth"] / G["moon"]))))

    # ---------------- A4: pit catalogue by type + inner apertures ----------------
    rows = list(csv.DictReader(open(OUT / "planetary_catalogue.csv")))
    narrowest = 5.87
    moon = [r for r in rows if r["body"] == "moon"]
    def f(r, k):
        v = (r.get(k) or "").strip()
        return float(v) if v else None
    by_type = {}
    for t in ("pit (mare)", "pit (impact melt)", "pit (highland)"):
        sub = [r for r in moon if r["type"] == t]
        outer = [f(r, "aperture_long_axis_m") for r in sub]
        inner = [f(r, "inner_long_axis_m") for r in sub]
        depth = [f(r, "depth_m") for r in sub]
        outer = [v for v in outer if v]; inner = [v for v in inner if v]
        depth = [v for v in depth if v]
        by_type[t] = dict(
            n=len(sub), n_outer=len(outer), n_inner=len(inner),
            outer_median=float(np.median(outer)) if outer else None,
            inner_median=float(np.median(inner)) if inner else None,
            inner_min=float(np.min(inner)) if inner else None,
            inner_max=float(np.max(inner)) if inner else None,
            depth_median=float(np.median(depth)) if depth else None,
            inner_ge_narrowest=int(np.sum(np.array(inner) >= narrowest)) if inner else None)
    all_inner = [f(r, "inner_long_axis_m") for r in moon]
    all_inner = np.array([v for v in all_inner if v])
    out["moon_pits"] = dict(
        by_type=by_type,
        all_inner=dict(n=len(all_inner), median=float(np.median(all_inner)),
                       frac_lt_narrowest=float(np.mean(all_inner < narrowest)),
                       frac_lt_10=float(np.mean(all_inner < 10.0))),
        narrowest_section_m=narrowest)
    mars = [r for r in rows if r["body"] == "mars"]
    apc = [r for r in mars if r["type"] == "APC"]
    apc_dim = [f(r, "aperture_long_axis_m") for r in apc]
    apc_dim = [v for v in apc_dim if v]
    out["mars_pits"] = dict(n_total=len(mars), n_apc=len(apc),
                            n_sky_type=sum(1 for r in mars if r["type"] == "sky"),
                            n_apc_dimensioned=len(apc_dim),
                            apc_median=float(np.median(apc_dim)),
                            note="'sky' candidates (n=354) carry no published "
                                 "dimensions; the dimensioned set is APCs only")

    # ---------------- A5: scale-ruler forensics ----------------
    aer = json.loads((OUT / "aerial_holes.json").read_text())["holes"][:3]
    dlio = json.loads((OUT_BASE / "tube_holes.json").read_text())["holes"]
    flf = json.loads((OUT_BASE / "flf_tube_holes.json").read_text())["holes"]
    # correspondence: aerial ids 0,1,2 = S1,S2,S3; dlio ids 0,1,2 match 0,1,2;
    # flf ids map (1,2,0) -> aerial (0,1,2) per flf_transform_landmark inliers
    A = np.array([h["centroid"] for h in aer])
    D = np.array([dlio[i]["centroid"] for i in (0, 1, 2)])
    F = np.array([flf[i]["centroid"] for i in (1, 2, 0)])
    def dists(P):
        return dict(S1S2=float(np.linalg.norm(P[0] - P[1])),
                    S1S3=float(np.linalg.norm(P[0] - P[2])),
                    S2S3=float(np.linalg.norm(P[1] - P[2])))
    da, dd, df = dists(A), dists(D), dists(F)
    out["scale_ruler"] = dict(
        aerial=da, dlio_map=dd, flf_map=df,
        dlio_minus_aerial={k: dd[k] - da[k] for k in da},
        flf_minus_aerial={k: df[k] - da[k] for k in da},
        dlio_minus_flf={k: dd[k] - df[k] for k in da},
        hole_semi_major=dict(
            aerial=[h["semi_major"] for h in aer],
            dlio=[dlio[i]["semi_major"] for i in (0, 1, 2)],
            flf=[flf[i]["semi_major"] for i in (1, 2, 0)]),
        note="inter-hole distances are rigid-transform invariant; agreement "
             "between the two independent SLAM maps against a common offset "
             "to the aerial model implicates the rim-centroid definition "
             "(S3 is a 9-15 m elongated opening seen from opposite sides), "
             "not map scale")

    # ---------------- volume per metre ----------------
    mrows = list(csv.DictReader(open(OUT / "morphometry.csv")))
    areas = np.array([float(r["area"]) for r in mrows if r["area"]])
    out["volume"] = dict(total_m3=float(np.sum(areas)),   # 1 m stations
                         per_m=float(np.mean(areas)),
                         n_sections=len(areas),
                         note="enclosed cross-section area integrated at 1 m "
                              "stations; floor bounded by top of breakdown, "
                              "so a lower bound on the void volume")

    (OUT / "review_stats.json").write_text(json.dumps(out, indent=2))
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
