#!/usr/bin/env python3
"""Checks 50 CSV values against the original metrics.json (exact float equality). Writes VERIFICATION.md.
The exports are read from ECBF_N200_NUMBERS (default: N200_numbers/ at the repository root; not included).
Command (inside data/): python3 tools/verify_sample.py"""
import csv, json, os, random, statistics, sys
from pathlib import Path
sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", OUT.parent / "N200_numbers"))   # N200 summary exports
rng = random.Random(20261007)
cells = {}
for src in ("final", "new_anchors"):
    d = json.load(open(NUM / src / "common/metrics.json"))
    for v in d["views"]:
        cells[(src, v["view_id"])] = v["cells"][0]
    del d
def rows(p):
    with open(OUT / p) as f:
        return list(csv.DictReader(f))
def num(x):
    return None if x == "" else float(x)
checks = []
def cmp(kind, where, csv_val, ref):
    ok = (csv_val is None and ref is None) or (csv_val is not None and ref is not None and csv_val == float(ref))
    checks.append((kind, where, csv_val, ref, ok))
# 15 summary values
sm = rows("summary/common_metrics_final.csv") + rows("summary/common_metrics_new_anchors.csv")
cols = [c for c in sm[0] if c.endswith(".mean") or c.endswith(".sd") or c.endswith("possible") or c.endswith("upper95_possible")]
for _ in range(15):
    r = rng.choice(sm); c = rng.choice([c for c in cols if r.get(c, "") != ""])
    ref = cells[(r["source"], r["view_id"])]
    for part in c.split("."):
        ref = ref[part]
    cmp("summary", f'{r["source"]}/{r["view_id"][:10]} {c}', num(r[c]), ref)
# 15 sweep values
sw = rows("sweeps/sweep_series.csv")
for _ in range(15):
    r = rng.choice(sw)
    ref = cells[(r["source"], r["view_id"])]
    if r["kind"] == "statistic":
        col, key = ("mean", "mean") if rng.random() < .5 else ("sd", "sd")
        cmp("sweep", f'{r["setting"]} {r["scenario"]} {r["swept_parameter"]}={r["x_value"]} {r["variant"]} {r["metric"]} {col}', num(r[col]), ref["statistics"][r["metric"]][key])
    else:
        f = {"q_cert": "qcert", "q_dep": "qdep"}[r["metric"]]
        cmp("sweep", f'{r["setting"]} {r["scenario"]} {r["swept_parameter"]}={r["x_value"]} {r["variant"]} {r["metric"]} possible', num(r["count"]), ref["trajectory_events"][f]["possible"])
# 15 heatmap values
hm = rows("heatmaps/heatmap_grid.csv")
for _ in range(15):
    r = rng.choice(hm); ref = cells[("final", r["view_id"])]["statistics"][r["metric"]]
    col = rng.choice(["mean", "sd"])
    cmp("heatmap", f'{r["scenario"]} {r["metric"]} a0={r["alpha0"]} a={r["alpha"]} d={r["delta"]} K={r["K"]} {r["variant"]} {col}', num(r[col]), ref[col])
# 5 per-trajectory means vs metrics.json (cell mean = mean of the 200 per-trajectory values)
pt = rows("per_trajectory/per_trajectory_values.csv")
groups = {}
for r in pt:
    groups.setdefault((r["setting"], r["scenario"], r["method_key"]), []).append(r)
want = {"Default": (0.15, 0.10, 0.05, 50), "T1": (0.30, 0.25, 0.20, 50), "T2": (0.30, 0.25, 0.15, 50), "T3": (0.30, 0.20, 0.20, 50), "T4": (0.30, 0.25, 0.025, 100), "T5": (0.30, 0.25, 0.10, 10)}
byconf = {}
for (src, vid), c in cells.items():
    cf = c["configuration"]
    if src == "final" and cf.get("L_mu") == 0.0005:
        byconf.setdefault((cf["benchmark"], cf["method"], cf["alpha0"], cf.get("alpha"), cf.get("delta"), str(cf.get("K")), cf.get("bound_scale"), cf.get("budget_normalizer")), []).append(c)
for _ in range(5):
    key = rng.choice(list(groups)); g = groups[key]; sett, scen, m = key
    a0, a, d, K = want[sett]
    if m in ("V-U", "V-D", "V-Cd", "V-C"):
        ck = (scen, m, a0, a, d, str(K), 1.0, 2.11)
    else:
        ck = None
    if ck is None:
        ref = [c for (src, vid), c in cells.items() if src == "final" and c["configuration"].get("benchmark") == scen and c["configuration"]["method"] == m and c["configuration"]["alpha0"] == a0 and c["configuration"]["L_mu"] == 0.0005 and "K" not in c["configuration"]]
    else:
        ref = byconf[ck]
    metric = rng.choice(["T_mis", "m_c", "q_CBF", "q_unsafe", "h_min"])
    vals = [float(r[metric]) for r in g if r[metric] != ""]
    cmp("per-trajectory", f"{sett} {scen} {m} {metric}: mean of 200 CSV values", statistics.mean(vals), ref[0]["statistics"][metric]["available_mean"])
n_ok = sum(c[4] for c in checks)
lines = ["# Verification of a sample of 50 CSV values against metrics.json", "",
         "Script `tools/verify_sample.py` (random seed 20261007); exact float equality with the original `final` and `new_anchors` `common/metrics.json` (per-trajectory rows: mean of the 200 CSV values equals the cell mean).", "",
         f"Result: {n_ok} of {len(checks)} equal.", "", "| # | file | value checked | CSV | metrics.json | equal |", "|---|---|---|---|---|---|"]
for i, (k, w, a, b, ok) in enumerate(checks, 1):
    lines.append(f"| {i} | {k} | {w} | {a!r} | {b!r} | {'yes' if ok else 'NO'} |")
(OUT / "VERIFICATION.md").write_text("\n".join(lines) + "\n")
print(n_ok, "of", len(checks))
