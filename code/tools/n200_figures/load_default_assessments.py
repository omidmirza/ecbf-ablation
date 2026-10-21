#!/usr/bin/env python3
"""Per-trajectory assessments of the N200 default setting.

Default (alpha0, alpha, delta, K) = (0.15, 0.10, 0.05, 50), b = 2.11, phi = 1, L_mu = 0.0005, rho_buf = 0;
nine methods (four variants + five baselines at alpha0 = 0.15), Corridor and Maze, seeds 970001-970200.
Read only. For every (configuration, seed) the row of N200_INDEX.sqlite (folder ECBF_N200_DATA) gives the source
and the assessment SHA-256; the saved assessment is read from
  n100_prefix_sealed  -> sealed N100 catalog observations.assessment_json (same raw_sha256 and seed)
  suffix_mac_original / suffix_mac_overflow -> <source>/results/<raw stem>.json
  suffix_windows      -> suffix_windows/staged_source/results/<raw stem>.json
  suffix_historical   -> suffix_historical/records.jsonl
and its canonical SHA-256 must equal the index. The per-cell means of every exported metric are then recomputed
and must equal N200_numbers/final/common/metrics.json (folder ECBF_N200_NUMBERS) exactly. Writes one JSON cache.
The N200 index, the raw assessments, the N = 100 run catalog and the exports are not part of this repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sqlite3
import statistics
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_figures/<this file>)
N200 = Path(os.environ.get("ECBF_N200_DATA", REPO / "N200_data"))          # N200_INDEX.sqlite and raw assessments
INDEX = N200 / "N200_INDEX.sqlite"
N100_CATALOG = Path(os.environ.get("ECBF_N100_CATALOG", REPO / "N100_data" / "runs.sqlite"))   # N = 100 run catalog
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))     # N200 summary exports
METRICS = NUM / "final/common/metrics.json"
OUT = Path(__file__).resolve().parent / "default_assessments.json"
SEEDS = list(range(970001, 970201))
CERTIFIED = ("V-U", "V-D", "V-Cd", "V-C")
BASELINES = ("Zero margin", "Conservative", "ACI", "QT", "Unfiltered")
KEYS = ("T_mis", "m_c", "m_Y", "q_CBF", "q_unsafe", "q_use", "q_nocert", "q_notrans", "p_max", "h_min", "tau_use", "N")


def recorded_path(p):
    """A path as recorded in manifests and caches: '<default folder>/<rest>' under an input folder named by an ECBF_*
    variable, else relative to the repository root, else the file name (no absolute local paths)."""
    p = Path(p).resolve()
    for var, name in (("ECBF_N200_NUMBERS", "N200_numbers"), ("ECBF_N200_DATA", "N200_data"),
                      ("ECBF_FIGURE_STYLE", "figure_style"), ("ECBF_DISPLAYS", "data")):
        if os.environ.get(var):
            try:
                return f"{name}/{p.relative_to(Path(os.environ[var]).resolve()).as_posix()}"
            except ValueError:
                pass
    try:
        return p.relative_to(REPO.resolve()).as_posix()
    except ValueError:
        return p.name


def canonical(v):
    return json.dumps(v, sort_keys=True, separators=(",", ":"))


def sha(text):
    return hashlib.sha256(text.encode()).hexdigest()


def config_of(scenario, method):
    if method in CERTIFIED:
        return {"K": 50, "L_mu": 0.0005, "alpha": 0.1, "alpha0": 0.15, "benchmark": scenario, "bound_scale": 1.0,
                "budget_normalizer": 2.11, "delta": 0.05, "method": method, "rho_buf": 0.0}
    return {"L_mu": 0.0005, "alpha0": 0.15, "benchmark": scenario, "method": method, "rho_buf": 0.0}


def main():
    idx = sqlite3.connect(f"file:{INDEX}?mode=ro&immutable=1", uri=True)
    n100 = sqlite3.connect(f"file:{N100_CATALOG}?mode=ro&immutable=1", uri=True)
    hist = None
    out = {}
    wanted = {canonical(config_of(s, m)): (s, m) for s in ("corridor", "maze") for m in CERTIFIED + BASELINES}
    index_rows = {k: [] for k in wanted}
    for row in idx.execute("SELECT configuration_json, seed, source, location, raw_sha256, assessment_sha256 FROM runs"):
        if row[0] in index_rows:
            index_rows[row[0]].append(row[1:])
    n100_cid = {}
    for cid, cj in n100.execute("SELECT configuration_id, configuration_json FROM configurations"):
        key = canonical({k: v for k, v in json.loads(cj).items() if v is not None})
        if key in wanted:
            n100_cid[key] = cid
    print("index scanned", flush=True)
    for scenario in ("corridor", "maze"):
        for method in CERTIFIED + BASELINES:
            cfg = config_of(scenario, method)
            rows = sorted(index_rows[canonical(cfg)])
            if [r[0] for r in rows] != SEEDS:
                raise ValueError(f"index rows for {cfg}: {len(rows)}")
            per_seed = []
            for seed, source, location, raw, asha in rows:
                if source == "n100_prefix_sealed":
                    found = n100.execute("SELECT assessment_json, raw_sha256 FROM observations WHERE configuration_id=? AND seed=?",
                                         (n100_cid[canonical(cfg)], seed)).fetchall()
                    cands = [json.loads(t) for (t, r) in found if r == raw]
                elif source in ("suffix_mac_original", "suffix_mac_overflow", "suffix_windows"):
                    stem = Path(location).stem
                    base = N200 / (location.split("/raw/")[0])
                    rec = json.loads((base / "results" / f"{stem}.json").read_text())
                    if rec["raw_sha256"] != raw:
                        raise ValueError(f"raw hash {location}")
                    cands = [rec["assessment"]]
                elif source == "suffix_historical":
                    if hist is None:
                        hist = {}
                        for line in (N200 / "suffix_historical/records.jsonl").open():
                            r = json.loads(line)
                            hist.setdefault(sha(canonical(r["assessment"])), r["assessment"])
                    cands = [hist[asha]] if asha in hist else []
                else:
                    raise ValueError(source)
                match = [a for a in cands if sha(canonical(a)) == asha]
                if len(match) < 1 or match[0]["seed"] != seed:
                    raise ValueError(f"assessment not found or hash differs: {cfg} seed {seed} ({source})")
                per_seed.append({k: match[0][k] for k in KEYS})
            out[f"{scenario}|{method}"] = dict(configuration=cfg, seeds=SEEDS, values=per_seed)
            print(scenario, method, "ok", flush=True)
    # exact parity with the verified export
    data = json.loads(METRICS.read_text())
    cells = {}
    for v in data["views"]:
        for c in v["cells"]:
            cells[canonical({k: x for k, x in c["configuration"].items() if x is not None})] = c
    checked = 0
    for key, entry in out.items():
        cell = cells[canonical(entry["configuration"])]
        for k in KEYS:
            if k == "N":
                continue
            vals = [r[k] for r in entry["values"] if r[k] is not None]
            st = cell["statistics"][k]
            if not st["applicable"]:
                continue
            mean = statistics.mean(vals) if vals else None
            if st["available_mean"] != mean:
                raise ValueError(f"parity {key} {k}: {mean} vs {st['available_mean']}")
            checked += 1
        entry["view_cell"] = dict(configuration_id=cell["configuration_id"], population_sha256=cell["population_sha256"],
                                  deployment_decision_fraction=dict(confirmed_mean=cell["deployment_decision_fraction"]["confirmed"]["mean"]))
    OUT.write_text(json.dumps(dict(generator=recorded_path(Path(__file__)), index=recorded_path(INDEX), metrics=recorded_path(METRICS),
                                   parity_means_checked=checked, cells=out), indent=0) + "\n")
    print("parity means checked:", checked, "->", OUT)


if __name__ == "__main__":
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    main()
