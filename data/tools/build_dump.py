#!/usr/bin/env python3
"""Data dump behind every N200 display (the CSV and JSON files of data/). Never simulates; reads only saved exports.

Inputs (read only):
  N200_numbers/{final,new_anchors}/common/metrics.json          (folder ECBF_N200_NUMBERS, not included)
  N200_numbers/{final,new_anchors}/temporal/temporal_metrics.json
  <SRC>/tools/n200_figures/setting_assessments.json   (per-trajectory saved assessments of Default and T1-T5,
      written by load_setting_assessments.py; hash-checked there against N200_INDEX.sqlite, means equal to metrics.json)
  <SRC>/tools/n200_tables/tables_values.json.gz        (every printed table cell with provenance)
  <SRC>/tools/n200_figures/sweepbox_ext_values.json    (per-trajectory values of the sweep box figures, written and
      checked by build_sweeps_box_ext.py)
  <map builder folder>/*manifest.json                  (map manifests; folder of ECBF_MAPS_BUILDER, not included)
<SRC> is code/ of this repository (the tools/ and catalog/ folders of the paper source tree).
Outputs (in the folder that contains tools/): tables/, catalog/, summary/, per_trajectory/, sweeps/, heatmaps/, maps/.
sweeps/sweepbox_ext_values.csv.gz is the gzip-compressed form of the CSV that build_sweeps_box_ext.py writes
(sweeps/sweepbox_ext_values.csv, same rows); sweeps/sweepbox_T3_values.csv is written by build_sweeps_box.py.
Command (inside data/):  python3 tools/build_dump.py
"""
from __future__ import annotations

import csv
import gc
import gzip
import json
import os
import shutil
import statistics
import sys
from pathlib import Path

import numpy as np

sys.dont_write_bytecode = True
OUT = Path(__file__).resolve().parents[1]
REPO = Path(__file__).resolve().parents[2]   # repository root (data/tools/build_dump.py)
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))   # N200 summary exports
SRC = REPO / "code"                           # tools/ and catalog/ of the paper source tree
ASSESS = OUT / "per_trajectory/source/setting_assessments.json"   # copy of SRC/tools/n200_figures/setting_assessments.json
TABLES = SRC / "tools/n200_tables/tables_values.json.gz"
MAPS = Path(os.environ.get("ECBF_MAPS_BUILDER", REPO / "maps_top5" / "build_maps_top5.py")).parent   # map manifests
HEATVALS = SRC / "tools/n200_heatmaps"
SWEEPBOX_EXT = SRC / "tools/n200_figures/sweepbox_ext_values.json"

# display names: internal method key -> name printed in the documents
DISPLAY = {"Zero margin": "c=0", "Conservative": "c_max", "ACI": "ACI", "QT": "QT", "Unfiltered": "Unfiltered",
           "V-U": "V-U", "V-D": "V-Delta", "V-Cd": "V-L", "V-C": "V-C"}
CERT = ("V-U", "V-D", "V-Cd", "V-C")
BASE = ("Zero margin", "Conservative", "ACI", "QT", "Unfiltered")
# setting label -> (alpha0, alpha, delta, K, paired variant); Default, New1..New5 = T1..T5
SETTINGS = {"Default": (0.15, 0.10, 0.05, 50, "V-C"), "New1": (0.30, 0.25, 0.20, 50, "V-C"),
            "New2": (0.30, 0.25, 0.15, 50, "V-C"), "New3": (0.30, 0.20, 0.20, 50, "V-Cd"),
            "New4": (0.30, 0.25, 0.025, 100, "V-D"), "New5": (0.30, 0.25, 0.10, 10, "V-C")}
TNAME = {"Default": "Default", "New1": "T1", "New2": "T2", "New3": "T3", "New4": "T4", "New5": "T5"}
STAT_KEYS = ("T_mis", "m_c", "m_Y", "m_int", "q_CBF", "q_Y", "q_use", "q_unsafe", "q_nocert", "q_notrans",
             "tau_use", "p_max", "runtime_ms", "h_min", "clearance_min", "n_certificates")
SWEEP_STATS = ("m_c", "q_CBF", "q_unsafe", "q_use", "T_mis", "p_max", "h_min")
KVALS = ["10", "50", "100", "200", "infinity"]
SKIP = {"observed_seeds", "seeds", "assessed_seeds"}


def write_csv(path, rows, fields=None):
    if fields is None:
        fields = []
        seen = set()
        for r in rows:
            for k in r:
                if k not in seen:
                    seen.add(k)
                    fields.append(k)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: ("" if r.get(k) is None else (json.dumps(r[k], sort_keys=True) if isinstance(r[k], (dict, list)) else
                                                         repr(r[k]) if isinstance(r[k], float) else r[k])) for k in fields})
    print("wrote", path, len(rows), "rows,", len(fields), "columns", flush=True)
    return fields


def flatten(o, prefix="", out=None):
    out = {} if out is None else out
    for k, v in o.items():
        if k in SKIP:
            continue
        key = f"{prefix}.{k}" if prefix else k
        if isinstance(v, dict) and k not in ("missing_reasons", "termination_counts", "source_counts"):
            flatten(v, key, out)
        else:
            out[key] = v
    return out


def cell_rows_common(source, data):
    rows = []
    for v in data["views"]:
        d = v["definition"]
        for c in v["cells"]:
            r = dict(source=source, view_id=v["view_id"], family=d["family"], view_label=d.get("label"),
                     anchor=d.get("anchor"), anchor_id=d.get("anchor_id"), role=d.get("role"),
                     anchor_alpha0=(d.get("anchor_setting") or {}).get("alpha0"),
                     anchor_alpha=(d.get("anchor_setting") or {}).get("alpha"),
                     anchor_delta=(d.get("anchor_setting") or {}).get("delta"),
                     anchor_K=(d.get("anchor_setting") or {}).get("K"), n100_view_id=d.get("n100_view_id"),
                     seed_first=min(c["seeds"]), seed_last=max(c["seeds"]), n_seeds=len(c["seeds"]))
            cf = c["configuration"]
            r.update(scenario=cf.get("benchmark"), method=cf.get("method"), method_display=DISPLAY.get(cf.get("method")),
                     alpha0=cf.get("alpha0"), alpha=cf.get("alpha"), delta=cf.get("delta"), K=cf.get("K"),
                     L_mu=cf.get("L_mu"), phi_bound_scale=cf.get("bound_scale"), b_budget_normalizer=cf.get("budget_normalizer"),
                     rho_buf=cf.get("rho_buf"), configuration_id=c.get("configuration_id"),
                     population_sha256=c.get("population_sha256"), n=c["n"])
            body = {k: x for k, x in c.items() if k not in ("configuration", "seeds", "n", "configuration_id", "population_sha256")}
            r.update(flatten(body))
            rows.append(r)
    return rows


def cell_rows_temporal(source, data):
    rows = []
    for v in data["views"]:
        d = v["definition"]
        for c in v["cells"]:
            cf = c["configuration"]
            r = dict(source=source, view_id=v["view_id"], family=d["family"], view_label=d.get("label"),
                     anchor=d.get("anchor"), anchor_id=d.get("anchor_id"),
                     anchor_alpha0=(d.get("anchor_setting") or {}).get("alpha0"),
                     anchor_alpha=(d.get("anchor_setting") or {}).get("alpha"),
                     anchor_delta=(d.get("anchor_setting") or {}).get("delta"),
                     anchor_K=(d.get("anchor_setting") or {}).get("K"),
                     scenario=cf.get("benchmark"), method=cf.get("method"), method_display=DISPLAY.get(cf.get("method")),
                     alpha0=cf.get("alpha0"), alpha=cf.get("alpha"), delta=cf.get("delta"), K=cf.get("K"),
                     L_mu=cf.get("L_mu"), phi_bound_scale=cf.get("bound_scale"), b_budget_normalizer=cf.get("budget_normalizer"),
                     rho_buf=cf.get("rho_buf"), configuration_id=c.get("configuration_id"), n=c["n"],
                     comparison_eligible=c.get("comparison_eligible"))
            r.update(flatten(c["statistics"], "statistics"))
            rows.append(r)
    return rows


def eq(have, want, key):
    if have is None:
        return False
    if key == "K" or isinstance(want, str):
        return str(have) == str(want)
    return abs(float(have) - float(want)) < 1e-12


class Reduced:
    """Compact cells for sweeps / heatmaps (same fields as extract_sweep_cells.py + sd, observed_n)."""

    def __init__(self):
        self.cells = []

    def add(self, source, data):
        for v in data["views"]:
            d = v["definition"]
            for c in v["cells"]:
                st = c["statistics"]
                ev = c.get("trajectory_events") or {}
                self.cells.append(dict(
                    source=source, view_id=v["view_id"], family=d["family"], anchor=d.get("anchor"),
                    configuration=c["configuration"], n=c["n"],
                    stats={k: dict(mean=st[k].get("mean"), sd=st[k].get("sd"), observed_n=st[k].get("observed_n"),
                                   applicable=st[k].get("applicable"), conditional=st[k].get("conditional"))
                           for k in STAT_KEYS if k in st},
                    events={f: {x: ev[f].get(x) for x in ("confirmed", "possible", "applicable_n", "upper95_possible",
                                                          "lower95_confirmed")} for f in ("qcert", "qdep") if ev.get(f)}))

    def pick(self, fam, source, scen, **want):
        out = []
        for c in self.cells:
            if c["family"] != fam or (source and c["source"] != source):
                continue
            cf = c["configuration"]
            if cf.get("benchmark") != scen:
                continue
            if all(eq(cf.get(k), v, k) for k, v in want.items()):
                out.append(c)
        return out

    def one(self, fam, source, scen, **want):
        r = self.pick(fam, source, scen, **want)
        if len(r) != 1:
            raise ValueError(f"{fam} {source} {scen} {want}: {len(r)} cells")
        return r[0]

    def sweep(self, label, scen, col, method):
        a0, a, d, K, _ = SETTINGS[label]
        base = dict(alpha0=a0, alpha=a, delta=d, K=str(K), method=method, L_mu=0.0005)
        src = "final" if label == "Default" else "new_anchors"
        if col == "alpha":
            return [(x, x, self.one("main_grid_200", "final", scen, **dict(base, alpha=x)))
                    for x in (0.05, 0.10, 0.15, 0.20, 0.25) if x < a0 - 1e-12]
        if col == "delta":
            return [(x, x, self.one("main_grid_200", "final", scen, **dict(base, delta=x))) for x in (0.025, 0.05, 0.10, 0.15, 0.20)]
        if col == "K":
            return [(k, i, self.one("main_grid_200", "final", scen, **dict(base, K=k))) for i, k in enumerate(KVALS)]
        if col == "b":
            return [(x, x, self.one("budget_normalizer_200", src, scen, **dict(base, budget_normalizer=x))) for x in (2.11, 2.5, 3.0, 4.0)]
        if col == "phi":
            return {lm: [(x, x, self.one("stress_drift_200", src, scen, **dict(base, L_mu=lm, bound_scale=x, budget_normalizer=2.11)))
                         for x in (0.0, 0.25, 0.5, 0.75, 1.0)] for lm in (0.0005, 0.005)}
        raise ValueError(col)


def sweep_rows(red):
    rows = []

    def emit(label, scen, col, x, xpos, method, lm, cell, is_setting):
        a0, a, d, K, paired = SETTINGS[label]
        common = dict(setting=TNAME[label], asset_label=label, alpha0=a0, alpha=a, delta=d, K=K, paired_variant=DISPLAY[paired],
                      scenario=scen, swept_parameter=col, x_value=x, x_position=xpos, variant=DISPLAY[method], method_key=method,
                      L_mu=lm, is_setting_value=is_setting, source=cell["source"], view_id=cell["view_id"])
        for k in SWEEP_STATS:
            s = cell["stats"].get(k)
            if not s:
                continue
            pv = s["mean"] * 100 if (k in ("q_CBF", "q_unsafe") and s["mean"] is not None) else s["mean"]
            rows.append(dict(common, metric=k, kind="statistic", mean=s["mean"], sd=s["sd"], count=None, count_confirmed=None,
                             upper95=None, n_assessed=s["observed_n"], n_total=cell["n"], applicable=s["applicable"],
                             plot_value=pv if s["applicable"] else None,
                             plot_unit="percent" if k in ("q_CBF", "q_unsafe") else ("s" if k == "T_mis" else "")))
        for f, name in (("qcert", "q_cert"), ("qdep", "q_dep")):
            e = cell["events"].get(f)
            if not e:
                continue
            ok = e.get("applicable_n") not in (None, 0)
            rows.append(dict(common, metric=name, kind="event_count", mean=None, sd=None, count=e["possible"],
                             count_confirmed=e["confirmed"], upper95=e["upper95_possible"], n_assessed=e["applicable_n"],
                             n_total=cell["n"], applicable=ok, plot_value=float(e["possible"]) if ok else None, plot_unit="count"))

    for label in SETTINGS:
        a0, a, d, K, _ = SETTINGS[label]
        for scen in ("corridor", "maze"):
            for col in ("alpha", "delta", "K", "b", "phi"):
                for m in CERT:
                    sw = red.sweep(label, scen, col, m)
                    series = sw.items() if isinstance(sw, dict) else [(0.0005, sw)]
                    for lm, pts in series:
                        for x, xpos, cell in pts:
                            sx = {"alpha": a, "delta": d, "K": str(K), "b": 2.11, "phi": 1.0}[col]
                            emit(label, scen, col, x, xpos, m, lm, cell, eq(x, sx, col) and (col != "phi" or lm == 0.0005 or True))
            for m in BASE:
                for lm in (0.0005, 0.005):
                    cell = red.one("matched_baselines_200", "final", scen, method=m, alpha0=a0, L_mu=lm)
                    emit(label, scen, "baseline", None, None, m, lm, cell, None)
    return rows


def heatmap_rows(red):
    rows = []
    for c in red.cells:
        if c["family"] != "main_grid_200":
            continue
        cf = c["configuration"]
        for k in STAT_KEYS:
            s = c["stats"].get(k)
            if not s:
                continue
            rows.append(dict(scenario=cf["benchmark"], metric=k, alpha0=cf["alpha0"], alpha=cf["alpha"], delta=cf["delta"],
                             K=cf["K"], variant=DISPLAY[cf["method"]], method_key=cf["method"], mean=s["mean"], sd=s["sd"],
                             n_assessed=s["observed_n"], n_total=c["n"], applicable=s["applicable"],
                             conditional=s["conditional"], view_id=c["view_id"]))
    return rows


def per_trajectory():
    cache = json.loads(ASSESS.read_text())
    rows, box = [], []
    keys = ("T_mis", "m_c", "m_Y", "q_CBF", "q_unsafe", "q_use", "q_nocert", "q_notrans", "p_max", "h_min", "tau_use", "N")
    for key, entry in cache["cells"].items():
        label, scen, method = key.split("|")
        a0, a, d, K, paired = SETTINGS[label]
        cf = entry["configuration"]
        vals = entry["values"]
        for seed, v in zip(entry["seeds"], vals):
            r = dict(setting=TNAME[label], asset_label=label, scenario=scen, method=DISPLAY[method], method_key=method,
                     is_baseline=method in BASE, setting_alpha0=a0, setting_alpha=a, setting_delta=d, setting_K=K,
                     setting_paired_variant=DISPLAY[paired], L_mu=cf["L_mu"], seed=seed)
            for k in keys:
                r[k] = v[k]
            r["q_CBF_pct"] = None if v["q_CBF"] is None else v["q_CBF"] * 100
            r["q_unsafe_pct"] = None if v["q_unsafe"] is None else v["q_unsafe"] * 100
            rows.append(r)
        if method != "Unfiltered" or True:
            for k in ("T_mis", "m_c", "q_CBF", "q_unsafe", "p_max", "h_min"):
                x = np.array([v[k] for v in vals if v[k] is not None], dtype=float)
                scale = 100.0 if k in ("q_CBF", "q_unsafe") else 1.0
                if len(x) == 0:
                    box.append(dict(setting=TNAME[label], asset_label=label, scenario=scen, method=DISPLAY[method], metric=k,
                                    n=0, plot_unit="percent" if scale == 100 else ""))
                    continue
                x = x * scale
                q1, med, q3 = np.percentile(x, [25, 50, 75])
                iqr = q3 - q1
                lo = x[x >= q1 - 1.5 * iqr].min()
                hi = x[x <= q3 + 1.5 * iqr].max()
                box.append(dict(setting=TNAME[label], asset_label=label, scenario=scen, method=DISPLAY[method], metric=k,
                                plot_unit="percent" if scale == 100 else "", n=int(len(x)), mean=float(x.mean()), sd=float(x.std(ddof=1)) if len(x) > 1 else None,
                                min=float(x.min()), q1=float(q1), median=float(med), q3=float(q3), max=float(x.max()),
                                whisker_low=float(lo), whisker_high=float(hi)))
    write_csv(OUT / "per_trajectory/per_trajectory_values.csv", rows)
    write_csv(OUT / "per_trajectory/boxplot_statistics.csv", box)
    return rows


def tables():
    shutil.copy2(TABLES, OUT / "tables/tables_values.json.gz")
    (OUT / "catalog").mkdir(exist_ok=True)
    shutil.copy2(SRC / "tools/n200_tables/catalog_values.json.gz", OUT / "catalog/catalog_values.json.gz")
    shutil.copy2(SRC / "tools/n200_tables/starred_setting.json", OUT / "tables/starred_setting.json")
    shutil.copy2(SRC / "catalog/README.md", OUT / "catalog/CATALOG_README.md")
    d = json.load(gzip.open(TABLES))
    rows = []
    for c in d["cells"]:
        r = {k: v for k, v in c.items() if k != "row"}
        r["row"] = c.get("row")
        r["row_label"] = " | ".join(f"{k}={v}" for k, v in (c.get("row") or {}).items())
        rows.append(r)
    write_csv(OUT / "tables/tables_values_flat.csv", rows)
    (OUT / "tables/tables_values_meta.json").write_text(json.dumps({k: v for k, v in d.items() if k != "cells"}, indent=1) + "\n")


def sweepbox_ext():
    """sweeps/sweepbox_ext_values.csv.gz: the rows of write_csv() in build_sweeps_box_ext.py, gzip-compressed (mtime 0)."""
    import io
    d = json.loads(SWEEPBOX_EXT.read_text())
    buf = io.StringIO(newline="")
    w = csv.writer(buf)
    keys = ["alpha0", "alpha", "delta", "K", "budget_normalizer", "bound_scale", "L_mu", "benchmark", "method"]
    w.writerow(["configuration_id"] + keys + ["seed", "index_source", "assessment_sha256", "m_c_over_cmax",
                                              "q_CBF_fraction", "q_unsafe_fraction", "T_mis_s",
                                              "qcert_applicable", "qcert_confirmed", "qcert_possible",
                                              "qdep_applicable", "qdep_confirmed", "qdep_possible"])
    n = 0
    for cid, c in sorted(d["configs"].items()):
        v, cf = c["values"], c["configuration"]
        for i, seed in enumerate(c["seeds"]):
            w.writerow([cid] + [cf.get(k) for k in keys] + [seed, c["index_source"][i],
                       c["assessment_sha256"][i], v["m_c"][i], v["q_CBF"][i], v["q_unsafe"][i], v["T_mis"][i],
                       v["q_cert_applicable"][i], v["q_cert_confirmed"][i], v["q_cert_possible"][i],
                       v["q_dep_applicable"][i], v["q_dep_confirmed"][i], v["q_dep_possible"][i]])
            n += 1
    dst = OUT / "sweeps/sweepbox_ext_values.csv.gz"
    with open(dst, "wb") as f, gzip.GzipFile(filename="sweepbox_ext_values.csv", mode="wb", fileobj=f, mtime=0,
                                             compresslevel=9) as g:
        g.write(buf.getvalue().encode())
    print("wrote", dst, n, "rows", flush=True)


def main():
    gc.disable()
    tables()
    sweepbox_ext()
    per_trajectory()
    red = Reduced()
    for src in ("final", "new_anchors"):
        data = json.loads((NUM / src / "common/metrics.json").read_text())
        rows = cell_rows_common(src, data)
        write_csv(OUT / f"summary/common_metrics_{src}.csv", rows)
        red.add(src, data)
        del data, rows
        t = json.loads((NUM / src / "temporal/temporal_metrics.json").read_text())
        write_csv(OUT / f"summary/temporal_metrics_{src}.csv", cell_rows_temporal(src, t))
        del t
    write_csv(OUT / "sweeps/sweep_series.csv", sweep_rows(red))
    write_csv(OUT / "heatmaps/heatmap_grid.csv", heatmap_rows(red))
    for p in sorted(HEATVALS.glob("out_B4_*/*/heatmap_values.json")):
        dst = OUT / "heatmaps/printed_values" / f"{p.parent.parent.name}__{p.parent.name}__heatmap_values.json"
        dst.parent.mkdir(exist_ok=True)
        shutil.copy2(p, dst)
    for p in sorted(MAPS.glob("*manifest.json")):
        shutil.copy2(p, OUT / "maps" / p.name)


if __name__ == "__main__":
    main()
