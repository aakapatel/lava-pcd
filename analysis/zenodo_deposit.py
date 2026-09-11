"""Create (or update) the Zenodo deposition for the data release and reserve
its DOI. The deposition is left UNPUBLISHED so the DOI can go into the
manuscript now and the record is published at acceptance (or earlier).

Usage (repo root of lava-pcd, or anywhere):
    ZENODO_TOKEN=... python analysis/zenodo_deposit.py \
        --bundle "../data_release" [--sandbox] [--deposition-id N] [--skip-upload]

Token: zenodo.org > account > Applications > Personal access tokens, scopes
"deposit:write" and "deposit:actions". Use --sandbox (sandbox.zenodo.org
token) for a dry run.

What it does:
  1. creates a deposition (or reuses --deposition-id) with the metadata below
     and prereserve_doi, and prints the reserved DOI;
  2. uploads every file of the bundle directory that is not yet in the
     deposition, through the bucket API (files up to 50 GB each; the record
     limit is 50 GB in total unless Zenodo raises it on request);
  3. never publishes. Publish from the Zenodo web page when ready.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from pathlib import Path

import requests

AUTHORS = [
    {"name": "Patel, Akash", "affiliation": "Lulea University of Technology", "orcid": "0000-0002-0020-6020"},
    {"name": "Stathoulopoulos, Nikolaos", "affiliation": "Lulea University of Technology", "orcid": "0000-0002-0108-6286"},
    {"name": "Llewellin, Edward W.", "affiliation": "Durham University", "orcid": "0000-0003-2165-7426"},
    {"name": "Oskarsson, Birgir V.", "affiliation": "Natural Science Institute of Iceland", "orcid": "0000-0001-5653-670X"},
    {"name": "Nikolakopoulos, George", "affiliation": "Lulea University of Technology", "orcid": "0000-0003-0126-1897"},
]

DESCRIPTION = """<p>Point clouds, surface model, derived products and an analysis-code
snapshot for the manuscript "Autonomous interior survey of a lava tube for
planetary habitat assessment" (under review, 2026). An autonomous aerial robot
mapped 300 m of the closed section of Raufarholshellir lava tube, Iceland;
the interior map is registered to an independent photogrammetric surface
model through the three skylights, giving a continuous roof-thickness
profile along the tube.</p>
<p>Contents: registered interior lidar maps (first and second long flights,
onboard and offline reconstructions) in ISN2016 / ISH2004 orthometric
coordinates; the photogrammetric surface model crops (dense cloud, ground
surface, DEM); national DEM crops used for validation; every derived number
(paper_numbers.json), centreline, cross-sections, overburden profile,
registration and validation reports; the figures; and a snapshot of the
analysis code (https://github.com/aakapatel/lava-pcd). See DATA_README.md for
coordinate frames, file contents and the reassembly of split files.</p>
<p>Project page: https://aakapatel.github.io/lavatube-survey/</p>"""


def metadata(sandbox: bool) -> dict:
    return {"metadata": {
        "title": "Autonomous interior survey of a lava tube for planetary habitat assessment: maps, surface model, derived products and analysis code",
        "upload_type": "dataset",
        "description": DESCRIPTION,
        "creators": AUTHORS,
        "access_right": "open",
        "license": "cc-by-4.0",
        "keywords": ["lava tube", "planetary caves", "autonomous exploration", "lidar", "roof thickness",
                     "Raufarholshellir", "Iceland", "photogrammetry", "aerial robot"],
        "related_identifiers": [
            {"identifier": "https://github.com/aakapatel/lava-pcd", "relation": "isSupplementTo", "resource_type": "software"},
            {"identifier": "https://aakapatel.github.io/lavatube-survey/", "relation": "isDescribedBy"},
        ],
        "prereserve_doi": True,
        "notes": "Raw lidar and IMU logs (about 450 GB) are available from the corresponding author on request.",
    }}


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bundle", required=True, help="directory with the release files")
    ap.add_argument("--sandbox", action="store_true")
    ap.add_argument("--deposition-id", type=int, default=None)
    ap.add_argument("--skip-upload", action="store_true")
    a = ap.parse_args()
    tok = os.environ.get("ZENODO_TOKEN")
    if not tok:
        sys.exit("set ZENODO_TOKEN (personal access token with deposit:write, deposit:actions)")
    base = "https://sandbox.zenodo.org/api" if a.sandbox else "https://zenodo.org/api"
    s = requests.Session(); s.headers["Authorization"] = f"Bearer {tok}"

    if a.deposition_id:
        r = s.get(f"{base}/deposit/depositions/{a.deposition_id}"); r.raise_for_status(); dep = r.json()
        r = s.put(dep["links"]["self"], json=metadata(a.sandbox)); r.raise_for_status(); dep = r.json()
    else:
        r = s.post(f"{base}/deposit/depositions", json=metadata(a.sandbox)); r.raise_for_status(); dep = r.json()
    doi = dep["metadata"].get("prereserve_doi", {}).get("doi") or dep.get("doi")
    print(f"deposition {dep['id']}: {dep['links']['html']}")
    print(f"reserved DOI: {doi}  ->  https://doi.org/{doi}")

    if a.skip_upload:
        return
    bucket = dep["links"]["bucket"]
    have = {f["filename"] for f in dep.get("files", [])}
    files = sorted(p for p in Path(a.bundle).iterdir() if p.is_file() and p.name != "upload.log")
    for p in files:
        if p.name in have:
            print(f"  skip (present) {p.name}"); continue
        size = p.stat().st_size / 1e6
        print(f"  upload {p.name} ({size:.0f} MB) ...", end="", flush=True)
        with open(p, "rb") as fh:
            r = s.put(f"{bucket}/{p.name}", data=fh)
        r.raise_for_status()
        ok = r.json()["checksum"].split(":")[-1] == hashlib.md5(open(p, "rb").read()).hexdigest() if size < 500 else True
        print(" done" if ok else " CHECKSUM MISMATCH")
    print("all files uploaded; the record is still a draft. Publish it on the Zenodo page when ready.")
    json.dump({"deposition_id": dep["id"], "doi": doi, "html": dep["links"]["html"]},
              open(Path(a.bundle) / "zenodo_deposition.json", "w"), indent=2)


if __name__ == "__main__":
    main()
