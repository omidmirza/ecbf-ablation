#!/usr/bin/env python3
"""Independent check of the tables of build_tables_v10_original.py: re-reads metrics.json and asserts every printed
string.

For every cell in tables_values.json it
  * finds the view by view_id in metrics.json and checks the recorded configuration,
  * recomputes the printed string (mean +- SD, event count/N (upper %), N-A, not-simulated dash),
  * recomputes the bold marks (best displayed mean among the four V-variants of one setting, ties all bold,
    arrow direction; never for m_cert, conditional metrics or SDs),
  * checks that the strings occur in the fragment file in the recorded order,
  * checks that every used view is an N200 view (family *_200, n = 200) whose definition.n100_view_id equals the
    recorded n100_view_id (final/ export: N200 view IDs, linked to N100 by n100_view_id),
  * checks that every used view has exactly n trajectories in population.jsonl (seeds 970001-970200), keyed by
    n200_view_id, each row carrying the same n100_view_id.
Rounding rule (stated in build_tables_v10_original.py): half-up after rounding the double to d+9 decimals; upper bounds ceil to 0.1 %.
Usage: python3 check_tables_v10_original.py [--values tables_values.json] [--metrics PATH] [--population PATH]
Inputs: ECBF_N200_NUMBERS (default <repository>/N200_numbers; not included); a relative recorded metrics path
"N200_numbers/<rest>" resolves against it.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
from collections import defaultdict
from decimal import ROUND_CEILING, ROUND_HALF_UP, Decimal
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_tables/<this file>)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports (not included)
DEF_POP = NUM / "final/common/population.jsonl"
PCT = {"q_CBF", "q_Y", "q_use", "q_unsafe", "q_nocert", "q_notrans", "p_max"}
DEC = {"m_c": 2, "m_Y": 2, "m_int": 2, "n_certificates": 2, "q_CBF": 2, "q_Y": 2, "q_use": 1, "q_unsafe": 2,
       "q_nocert": 1, "q_notrans": 1, "h_min": 2, "clearance_min": 2, "tau_use": 2, "p_max": 1, "T_mis": 2}
UP = {"q_use", "h_min", "clearance_min"}
NOBOLD = {"n_certificates", "tau_use", "p_max"}


def rnd(v, d):
    s = "{:.%df}" % (d + 9)
    q = Decimal(s.format(v)).quantize(Decimal(10) ** -d, rounding=ROUND_HALF_UP)
    t = format(q, f".{d}f")
    return t[1:] if t.startswith("-") and q.is_zero() else t


def ceil_pct(u):
    return format((Decimal(repr(u)) * 100).quantize(Decimal("0.1"), rounding=ROUND_CEILING), ".1f")


def cp_upper(x, n):
    """Exact one-sided 95% Clopper-Pearson bound via exact binomial tails (fractions) and bisection."""
    if x >= n:
        return 1.0
    lo, hi = Fraction(0), Fraction(1)
    for _ in range(60):
        p = (lo + hi) / 2
        tail = sum(math.comb(n, k) * p ** k * (1 - p) ** (n - k) for k in range(x + 1))
        if tail > Fraction(1, 20):
            lo = p
        else:
            hi = p
        lo, hi = Fraction(float(lo)), Fraction(float(hi))
    return float(hi)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--values", type=Path, default=HERE / "tables_values.json")
    ap.add_argument("--metrics", type=Path)
    ap.add_argument("--population", type=Path, default=DEF_POP)
    ap.add_argument("--skip-population", action="store_true")
    a = ap.parse_args()
    V = json.loads(a.values.read_text())
    rec = Path(V["metrics"])   # relative "N200_numbers/<rest>" resolves against NUM; absolute as is
    mpath = a.metrics or (rec if rec.is_absolute() else NUM.joinpath(*rec.parts[1:]))
    raw = mpath.read_bytes()
    assert hashlib.sha256(raw).hexdigest() == V["metrics_sha256"], "metrics.json changed since the build"
    M = json.loads(raw)
    del raw
    views = {v["view_id"]: v for v in M["views"]}
    base_present = set()
    for v in M["views"]:
        c = v["cells"][0]["configuration"]
        if v["definition"]["family"] in ("matched_baselines_200", "matched_ACI_200"):
            base_present.add((c["benchmark"], round(c["alpha0"], 6), round(c["L_mu"], 6), c["method"]))
    cells = V["cells"]
    nerr = 0

    def fail(msg, e):
        nonlocal nerr
        nerr += 1
        if nerr <= 20:
            print("FAIL", msg, {k: e.get(k) for k in ("fragment", "table", "method", "metric", "printed", "view_id")})

    expected_bold = {}
    groups = defaultdict(list)
    for i, e in enumerate(cells):
        kind = e["kind"]
        if kind == "not_simulated":
            r = e["row"]
            if (r["scenario"], round(r["alpha0"], 6), round(r["L_mu"], 6), e["method"]) in base_present:
                fail("dash although simulated", e)
            if e["printed"] != "--":
                fail("dash string", e)
            continue
        v = views.get(e["view_id"])
        if v is None or len(v["cells"]) != 1:
            fail("view missing", e)
            continue
        c = v["cells"][0]
        if c["configuration"] != e["configuration"] or c["configuration"]["method"] != e["method"]:
            fail("configuration mismatch", e)
        df = v["definition"]
        if (not df["family"].endswith("_200") or c["n"] != 200 or not df.get("n100_view_id")
                or df["n100_view_id"] != e.get("n100_view_id") or df["n100_view_id"] == e["view_id"]):
            fail("N200 view / n100_view_id link", e)
        grp, sub = e["source"].split(".")
        if grp == "statistics":
            s = c["statistics"][sub]
            if not s["applicable"] or s["observed_n"] == 0:
                if kind != "NA" or e["printed"] != r"\na":
                    fail("should be N-A", e)
                continue
            if kind != "meanSD":
                fail("should be mean+-SD", e)
                continue
            sc = 100 if sub in PCT else 1
            m, sd = rnd(s["mean"] * sc, DEC[sub]), rnd(s["sd"] * sc, DEC[sub])
            if e["method"] in ("V-U", "V-D", "V-Cd", "V-C") and e["row"].get("table") != "topfive" or \
                    (e["row"].get("table") == "topfive" and e["row"].get("anchor") == "current_default"):
                key = (e["fragment"], e["table"], json.dumps({k: x for k, x in e["row"].items()}, sort_keys=True), sub)
                groups[key].append((i, m, s["conditional"]))
            expected_bold[i] = False
            e["_m"], e["_sd"] = m, sd
        else:
            x = c[grp][sub]
            if grp == "trajectory_events":
                n, lo, hi, u = x["applicable_n"], x["confirmed"], x["possible"], x["upper95_possible"]
            elif grp == "realized_transition_events":
                n, lo, hi, u = x["assessed_n"], x["count"], x["count"], x["upper95"]
            else:
                n, lo, hi = x["assessed_trajectory"], x["any_excess"], x["any_excess"]
                u = cp_upper(hi, n) if n else None
            if not n:
                if e["printed"] != r"\na":
                    fail("zero assessed must be N-A", e)
                continue
            want = (f"{lo}" if lo == hi else f"{lo}--{hi}") + f"/{n} ({ceil_pct(u)})"
            if e["printed"] != want:
                fail(f"event string, want {want}", e)
    for key, lst in groups.items():
        if len(lst) != 4 or any(cond for _, _, cond in lst) or key[3] in NOBOLD:
            continue
        vals = [Decimal(m) for _, m, _ in lst]
        best = max(vals) if key[3] in UP else min(vals)
        for (i, m, _), val in zip(lst, vals):
            expected_bold[i] = val == best
    for i, b in expected_bold.items():
        e = cells[i]
        mm = rf"\mathbf{{{e['_m']}}}" if b else e["_m"]
        want = f"${mm}\\vpm {e['_sd']}$"
        if e["printed"] != want:
            fail(f"mean/SD string, want {want}", e)
    # order of appearance in the fragments
    byfrag = defaultdict(list)
    for e in cells:
        byfrag[e["fragment"]].append(e["printed"])
    for f, strs in byfrag.items():
        txt = (a.values.parent / "fragments" / f).read_text()
        body = txt[txt.find(r"\begin{tabular}"):]
        pos = 0
        for s in strs:
            j = body.find(s, pos)
            if j < 0:
                nerr += 1
                print("FAIL order", f, s)
                break
            pos = j + len(s)
    used = {e["view_id"] for e in cells if e.get("view_id")}
    print(f"{len(cells)} cells, {len(used)} views, {sum(expected_bold.values())} bold marks checked")
    if not a.skip_population:
        cnt, seeds, link = defaultdict(int), defaultdict(set), defaultdict(set)
        tag = '"n200_view_id":"'
        with open(a.population) as fh:
            for line in fh:
                k = line.find(tag)
                if k < 0:
                    raise SystemExit("population.jsonl has no n200_view_id field (old common_02 layout?)")
                vid = line[k + len(tag):k + len(tag) + 64]
                if vid in used:
                    r = json.loads(line)
                    assert r["n200_view_id"] == vid
                    cnt[vid] += 1
                    seeds[vid].add(r["seed"])
                    link[vid].add(r["n100_view_id"])
        for vid in used:
            n = views[vid]["cells"][0]["n"]
            if cnt[vid] != n or seeds[vid] != set(range(970001, 970001 + n)):
                nerr += 1
                print("FAIL population", vid, cnt[vid], n)
            if link[vid] != {views[vid]["definition"]["n100_view_id"]}:
                nerr += 1
                print("FAIL population n100 link", vid, link[vid])
        print(f"population rows checked for {len(used)} views")
    print("ERRORS:", nerr)
    raise SystemExit(1 if nerr else 0)


if __name__ == "__main__":
    main()
