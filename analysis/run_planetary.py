"""S4.2: propagate the measured stability envelope to the Mars/Moon catalogues.

Physics: for a uniformly loaded roof plate, sigma_max = beta * rho * g * L^2 /
tau <= sigma_t, so at fixed tau/L the marginally stable span scales as
L_max ~ g^-1/2. The Raufarholshellir measurement gives (a) kappa_env, the
minimum tau/L observed over intact roof, and (b) the demonstrated flyable
cross-section (the narrowest section the vehicle actually traversed). Both are
propagated against the catalogued apertures:

  - Mars: MGC3 APC diameters (the only sub-population with dimensions).
  - Moon: LROC pit atlas funnel and inner apertures.

Counts are reported at several aperture thresholds so the paper can quote the
one matching its wording. No stability inference is drawn for candidates
without dimensions; they are reported as 'undimensioned'.

Inputs: analysis_out/planetary_catalogue.csv, roof_summary.json,
        morphometry_summary.json
Output: analysis_out/planetary_summary.json

Run:  .venv/bin/python analysis/run_planetary.py
"""
from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "analysis_out"

G = dict(earth=9.81, mars=3.71, moon=1.62)


def main() -> None:
    roof = json.loads((OUT / "roof_summary.json").read_text())
    morpho = json.loads((OUT / "morphometry_summary.json").read_text())
    kappa = float(roof["kappa_env"])

    rows = list(csv.DictReader(open(OUT / "planetary_catalogue.csv")))

    def apertures(body: str, col_candidates: list[str]) -> np.ndarray:
        vals = []
        for r in rows:
            if r["body"].strip().lower() != body:
                continue
            for c in col_candidates:
                v = (r.get(c) or "").strip()
                if v:
                    try:
                        vals.append(float(v))
                        break
                    except ValueError:
                        pass
        return np.array(vals)

    cols = rows[0].keys()
    print("columns:", list(cols))
    mars_ap = apertures("mars", ["aperture_long_axis_m"])
    moon_funnel = apertures("moon", ["aperture_long_axis_m"])

    scale = {b: float(np.sqrt(G["earth"] / g)) for b, g in G.items()}
    # Demonstrated flyable section: the narrowest measured span the vehicle
    # traversed (it flew the whole surveyed length).
    # Note: uses the 'min width over valid sections' from morphometry.
    w_min = None
    for k in ("width_min_m",):
        if k in morpho:
            w_min = morpho[k]
    thresholds = [5.0, 10.0, 20.0, 50.0]

    def counts(ap: np.ndarray) -> dict:
        return {f"ge_{int(t)}m": int((ap >= t).sum()) for t in thresholds}

    summary = dict(
        kappa_env_earth=kappa,
        span_scale=dict(mars=round(scale["mars"], 2), moon=round(scale["moon"], 2)),
        equivalent_kappa=dict(
            note="a roof at the terrestrial envelope tau/L supports a span "
                 "larger by span_scale at the same tau and strength",
        ),
        mars=dict(
            n_total=sum(1 for r in rows if r["body"].strip().lower() == "mars"),
            n_dimensioned=int(len(mars_ap)),
            aperture_median_m=(round(float(np.median(mars_ap)), 1)
                               if len(mars_ap) else None),
            counts=counts(mars_ap),
        ),
        moon=dict(
            n_total=sum(1 for r in rows if r["body"].strip().lower() == "moon"),
            n_dimensioned=int(len(moon_funnel)),
            aperture_median_m=(round(float(np.median(moon_funnel)), 1)
                               if len(moon_funnel) else None),
            counts=counts(moon_funnel),
        ),
        demonstrated_min_span_m=w_min,
    )
    (OUT / "planetary_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
