#!/usr/bin/env python3
"""Box-plot sweep figure at T3 for the arXiv edition.

sweepbox_T3_compact.pdf: setting T3 = (alpha0, alpha, delta, K) = (0.30, 0.20, 0.20, 50), b = 2.11,
phi = 1, L_mu = 0.0005. Columns: the swept parameter alpha, delta, K (10, 50, 100, 200, inf), b; the column header
states what varies and the fixed values of the other three parameters and alpha0. Two blocks, Corridor then Maze
(scenario named once per block), four rows each: m_c,k/c_max, q+_CBF,k (%), T_mis,k (s) as box plots over the
N_MC = 200 per-trajectory values, and q_dep as the count of trajectories with a deployment event (confirmed or
possible, as in the tables) with a bar to its one-sided 95 % Clopper-Pearson upper bound (in trajectories).
Shared x per column, shared y per row; column widths proportional to the number of swept values.
Boxes: V-L (#009E73, left) and V-C (#E69F00, right); box = quartiles, black median line, whiskers = 5th and 95th
percentiles, no outlier points, mean = the variant's marker (V-L pentagon, V-C triangle). No line connects swept
values. The T3 value: light grey column band and bold tick label.
Scales: each box row's y range spans the 5th-95th percentile whiskers of V-L and V-C over all columns of the row
plus MARGIN of that span on each side (never widened for a baseline or a label). Baselines QT and ACI at
alpha0 = 0.30 (matched_baselines_200 means): a horizontal reference line when inside the row's range; otherwise an
edge annotation inside a free corner of one panel of the row ("ACI 0.22, QT 0.30" with an arrow towards the value).
c_max is not drawn (its values are known: m_c/c_max = 1, q+_CBF = 0, q_dep = 0; checked here).

Style of this figure: text 11 pt (ticks, labels, headers, legend, annotations), scenario names 12 pt, mean markers
5.5 pt with a different shape per variant (so that the annotations are distinguishable), box, whisker, median and
q_dep bar lines 1.5 pt; one legend row above all panels, outside the axes, centred, so that no legend affects any
scale. The rest follows the figure helper through n200_style / ieee_figure.py: vector PDF, usetex
Computer Modern, no text under 10 pt (sub- and superscripts at 10 pt via \\DeclareMathSizes), reference lines
1.8 pt, strict helper save (zero margin, registered width 505.89 TeX pt = arXiv \\textwidth), text-overlap check.

Data (read only, no simulation): the sweep cells of build_sweeps.Data (sweep_cells.json: main_grid_200 cells of
N200_numbers/final for alpha, delta, K; budget_normalizer_200 cells of N200_numbers/new_anchors for b);
baselines from final matched_baselines_200. Per-trajectory values from the stored assessments, read with the method
of load_setting_assessments.py: the row of N200_INDEX.sqlite per (configuration, seed) gives
the source and the assessment SHA-256 (the canonical SHA-256 of the saved assessment must equal it); the
new_anchor_sweeps source (results/<raw stem>.json) is added. Checks: seeds 970001-970200 per cell; the (seed,
assessment SHA-256) pairs equal those of {final,new_anchors}/common/population.jsonl for the cell's configuration
id; n = 200; the means of the per-trajectory m_c, q_CBF and T_mis equal the metrics.json means (available_mean) to
1e-9 relative; the per-trajectory q_dep counts (possible, confirmed) equal metrics.json trajectory_events.qdep and the
Clopper-Pearson bound equals its upper95_possible to 1e-9.
Outputs: sweepbox_T3_values.json (cache next to this script, with the checks), the CSV
sweeps/sweepbox_T3_values.csv of the display data (ECBF_DISPLAYS, default data/; one row per box and seed), the PDF and a
manifest with every box statistic, q_dep point, reference line and annotation as drawn.

Usage: python3 build_sweeps_box.py [--reload] [--load-only] [--out out_paper] [--omit-zero-rows] [--name STEM]
(arXiv: --omit-zero-rows, name sweepbox_T3_compact; full: --name sweepbox_T3_full with the q_dep rows.)
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
import sqlite3
import statistics
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs
import build_sweeps as B
import load_setting_assessments as L

import matplotlib.colors as mcolors  # noqa: E402
import matplotlib.patheffects as pe  # noqa: E402
from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402
from matplotlib.transforms import ScaledTranslation  # noqa: E402

HERE = Path(__file__).resolve().parent
LABEL = "New3"
VARIANTS = ("V-Cd", "V-C")
SCENS = ("corridor", "maze")
COLS = ("alpha", "delta", "K", "b")
METRICS = ("m_c", "q_CBF", "T_mis")            # box rows
ROWS = ("m_c", "q_CBF", "T_mis", "q_dep")
QDEP_KEYS = ("qdep_applicable", "qdep_confirmed", "qdep_possible")
REF_BASE = ("QT", "ACI")
NAME = "sweepbox_T3_compact"
CACHE = HERE / "sweepbox_T3_values.json"
DISPLAYS = Path(os.environ.get("ECBF_DISPLAYS", L.REPO / "data"))   # display data folder (data/ of this repository)
DUMP_CSV = DISPLAYS / "sweeps/sweepbox_T3_values.csv"
EXPORTS = {"final": L.NUM / "final/common",
           "new_anchors": L.NUM / "new_anchors/common"}
SEEDS = L.SEEDS
N_MC = 200
SCALE = {"m_c": 1.0, "q_CBF": 100.0, "T_mis": 1.0}

# figure style (see header)
TXT = 11.0
TXT_SCEN = 12.0
LW_BOX = 1.5
MS_MEAN = 5.5
HEIGHT_PT = 600.0
HEIGHT_PER_UNIT = 600.0 * 0.80 / (2 * (3 * 1.0 + 0.72) + 0.34)   # pt per row-height unit (about 80 % of the page is plot area)
HEIGHT_CAP = 640.0                 # arXiv \textheight 696 pt minus a four-line caption
R_BOX, R_DEP, R_GAP = 1.0, 0.72, 0.34   # row heights: box rows, q_dep row, gap between the blocks
MARGIN = 0.07                      # y margin per side, fraction of the whisker span of the row
QDEP_YLIM = (-0.9, 6.6)
OFFSET, WIDTH = 0.2, 0.34          # box centre offset and box width in category units
BAND = "#E4E4E4"                   # T3 column band
WSPACE, HSPACE = 0.04, 0.10
COL_RATIO = {"alpha": 5.0, "delta": 5.6, "K": 5.0, "b": 4.0}   # per swept value; delta wider for "0.025" "0.05"
YLAB = {"m_c": r"$m_{c,k}/c_{\max}$", "q_CBF": r"$q^{+}_{\mathrm{CBF},k}$ (\%)", "T_mis": r"$T_{\mathrm{mis},k}$ (s)",
        "q_dep": r"$\widehat q_{\mathrm{dep}}$ (count)"}
PNAME = {"alpha0": r"\alpha_0", "alpha": r"\alpha", "delta": r"\delta", "K": "K", "b": "b"}
QT_EFFECT = [pe.Stroke(linewidth=S.LW_REF + 1.0, foreground="#7A6A00"), pe.Normal()]
REF_LS = {"QT": S.STYLE["QT"][2], "ACI": S.STYLE["ACI"][2]}


# ----------------------------------------------------------------------------------------------- data
def x_label(col, x):
    if col == "K":
        return B.KVALS[int(x)]
    return f"{x:g}"


def sweep_boxes(data):
    """One record per (scenario, column, swept value, variant): the sweep cell of build_sweeps.Data."""
    out = []
    for scen in SCENS:
        for col in COLS:
            for m in VARIANTS:
                for pos, (x, cell) in enumerate(data.sweep(LABEL, scen, col, m)):
                    out.append(dict(scen=scen, col=col, pos=pos, x=x, x_label=x_label(col, x), method=m,
                                    is_setting=abs(float(x) - float(B.setting_x(LABEL, col))) < 1e-12,
                                    cell=cell))
    return out


def cp_upper(k, n, level=0.95):
    """Exact one-sided Clopper-Pearson upper bound of a binomial proportion (bisection on the binomial CDF)."""
    if k >= n:
        return 1.0

    def cdf(p):
        return sum(math.comb(n, i) * p ** i * (1 - p) ** (n - i) for i in range(k + 1))
    lo, hi = 0.0, 1.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if cdf(mid) > 1 - level:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def read_assessments(cfgs):
    """{canonical config: [(seed, index source, assessment_sha256, assessment dict)]} for every config, read as in
    load_setting_assessments.main (hash-checked against N200_INDEX.sqlite), plus the new_anchor_sweeps source."""
    idx = sqlite3.connect(f"file:{L.INDEX}?mode=ro&immutable=1", uri=True)
    n100 = sqlite3.connect(f"file:{L.N100_CATALOG}?mode=ro&immutable=1", uri=True)
    n100_cid = {}
    for cid, cj in n100.execute("SELECT configuration_id, configuration_json FROM configurations"):
        key = L.canonical({k: v for k, v in json.loads(cj).items() if v is not None})
        if key in cfgs:
            n100_cid[key] = cid
    hist = None
    out = {}
    for key in sorted(cfgs):
        rows = sorted(idx.execute("SELECT seed, source, location, raw_sha256, assessment_sha256 FROM runs "
                                  "WHERE configuration_key=? OR configuration_json=?", (key, key)).fetchall())
        if [r[0] for r in rows] != SEEDS:
            raise ValueError(f"index rows for {key}: {len(rows)}")
        per = []
        for seed, source, location, raw, asha in rows:
            if source == "n100_prefix_sealed":
                found = n100.execute("SELECT assessment_json, raw_sha256 FROM observations WHERE configuration_id=? "
                                     "AND seed=?", (n100_cid[key], seed)).fetchall()
                cands = [json.loads(t) for (t, r) in found if r == raw]
            elif source in ("suffix_mac_original", "suffix_mac_overflow", "suffix_windows", "new_anchor_sweeps"):
                stem = Path(location).stem
                base = L.N200 / (location.split("/raw/")[0])
                rec = json.loads((base / "results" / f"{stem}.json").read_text())
                if rec["raw_sha256"] != raw:
                    raise ValueError(f"raw hash {location}")
                cands = [rec["assessment"]]
            elif source == "suffix_historical":
                if hist is None:
                    hist = {}
                    for line in (L.N200 / "suffix_historical/records.jsonl").open():
                        r = json.loads(line)
                        hist.setdefault(L.sha(L.canonical(r["assessment"])), r["assessment"])
                cands = [hist[asha]] if asha in hist else []
            else:
                raise ValueError(source)
            match = [a for a in cands if L.sha(L.canonical(a)) == asha]
            if len(match) < 1 or match[0]["seed"] != seed:
                raise ValueError(f"assessment not found or hash differs: {key} seed {seed} ({source})")
            per.append((seed, source, asha, match[0]))
        out[key] = per
        print("read", json.loads(key)["benchmark"], json.loads(key)["method"], key[:90], flush=True)
    return out


def population_pairs(cids):
    """{configuration_id: {seed: set(assessment_sha256)}} and {configuration_id: set(view ids)} from both
    population.jsonl files."""
    pairs = {c: {} for c in cids}
    views = {c: set() for c in cids}
    for name, folder in EXPORTS.items():
        with open(folder / "population.jsonl") as f:
            for line in f:
                i = line.find('"configuration_id":"')
                cid = line[i + 20:i + 84]
                if cid not in pairs:
                    continue
                r = json.loads(line)
                pairs[cid].setdefault(r["seed"], set()).add(r["assessment_sha256"])
                views[cid].add(f"{name}:{r.get('n200_view_id') or r.get('view_id')}")
    return pairs, views


def metrics_cells(view_ids):
    """{(source, view_id, canonical config): metrics.json cell} for the wanted views."""
    out = {}
    for name, folder in EXPORTS.items():
        d = json.loads((folder / "metrics.json").read_text())
        for v in d["views"]:
            if (name, v.get("view_id")) not in view_ids:
                continue
            for c in v["cells"]:
                cfg = L.canonical({k: x for k, x in c["configuration"].items() if x is not None})
                out[(name, v["view_id"], cfg)] = dict(n=c["n"], configuration_id=c.get("configuration_id"),
                                                      stats={k: c["statistics"][k] for k in METRICS},
                                                      qdep=c["trajectory_events"]["qdep"])
        del d
        print("metrics read:", name, flush=True)
    return out


def _mean_of(st):
    return st.get("available_mean", st.get("mean"))


def load(reload=False):
    if CACHE.exists() and not reload:
        return json.loads(CACHE.read_text())
    data = B.Data()
    boxes = sweep_boxes(data)
    cfgs = {}
    for b in boxes:
        cfg = L.canonical({k: v for k, v in b["cell"]["configuration"].items() if v is not None})
        b["cfg"] = cfg
        b["configuration_id"] = L.sha(cfg)
        cfgs[cfg] = b["configuration_id"]
    per = read_assessments(set(cfgs))
    pairs, views = population_pairs(set(cfgs.values()))
    mcells = metrics_cells({(b["cell"]["source"], b["cell"]["view_id"]) for b in boxes})
    checks = []
    for cfg, cid in sorted(cfgs.items()):
        idx_pairs = {s: a for s, _, a, _ in per[cfg]}
        pop = pairs[cid]
        ok = sorted(pop) == SEEDS and all(pop[s] == {idx_pairs[s]} for s in SEEDS)
        if not ok:
            raise ValueError(f"population pairs differ for {cfg}")
        checks.append(dict(kind="population", configuration_id=cid, seeds=len(pop), views=sorted(views[cid])))
    records = []
    max_rel, max_cp = 0.0, 0.0
    for b in boxes:
        cell = b["cell"]
        mc = mcells.get((cell["source"], cell["view_id"], b["cfg"]))
        if mc is None:
            raise ValueError(f"metrics cell missing: {cell['source']} {cell['view_id']}")
        if mc["n"] != N_MC or (mc["configuration_id"] and mc["configuration_id"] != b["configuration_id"]):
            raise ValueError(f"n or configuration id differs: {cell['view_id']}")
        assess = [a for _, _, _, a in per[b["cfg"]]]
        vals = {}
        for k in METRICS:
            v = [a[k] for a in assess]
            if any(x is None for x in v):
                raise ValueError(f"{k} missing for some seed in {cell['view_id']}")
            want = _mean_of(mc["stats"][k])
            rel = abs(statistics.mean(v) - want) / max(abs(want), 1e-300)
            if rel > 1e-9:
                raise ValueError(f"mean parity {cell['view_id']} {k}: {statistics.mean(v)} vs {want}")
            max_rel = max(max_rel, rel)
            vals[k] = v
        for k in QDEP_KEYS:
            vals[k] = [bool(a[k]) for a in assess]
        if not all(vals["qdep_applicable"]):
            raise ValueError(f"q_dep not applicable for some seed in {cell['view_id']}")
        if any(c and not p for c, p in zip(vals["qdep_confirmed"], vals["qdep_possible"])):
            raise ValueError(f"confirmed but not possible q_dep event in {cell['view_id']}")
        qd = dict(possible=sum(vals["qdep_possible"]), confirmed=sum(vals["qdep_confirmed"]), applicable_n=N_MC)
        ev = mc["qdep"]
        if (qd["possible"], qd["confirmed"], qd["applicable_n"]) != (ev["possible"], ev["confirmed"],
                                                                     ev["applicable_n"]):
            raise ValueError(f"q_dep counts differ in {cell['view_id']}: {qd} vs {ev}")
        qd["upper95"] = cp_upper(qd["possible"], N_MC)
        if abs(qd["upper95"] - ev["upper95_possible"]) > 1e-9:
            raise ValueError(f"q_dep bound differs in {cell['view_id']}: {qd['upper95']} vs {ev['upper95_possible']}")
        max_cp = max(max_cp, abs(qd["upper95"] - ev["upper95_possible"]))
        records.append(dict(scen=b["scen"], col=b["col"], pos=b["pos"], x=b["x"], x_label=b["x_label"],
                            method=b["method"], is_setting=b["is_setting"], source=cell["source"],
                            family=cell["family"], view_id=cell["view_id"], configuration_id=b["configuration_id"],
                            configuration=json.loads(b["cfg"]), n=mc["n"],
                            metrics_mean={k: _mean_of(mc["stats"][k]) for k in METRICS}, qdep=qd, qdep_export=ev,
                            seeds=[s for s, _, _, _ in per[b["cfg"]]],
                            index_source=[src for _, src, _, _ in per[b["cfg"]]],
                            assessment_sha256=[a for _, _, a, _ in per[b["cfg"]]], values=vals))
    base = {}
    for scen in SCENS:
        for m in REF_BASE + ("Conservative",):
            c = data.baseline(scen, S.SETTINGS[LABEL][0], m)
            ev = c["events"]["qdep"]
            base.setdefault(scen, {})[m] = dict(view_id=c["view_id"], q_dep=ev["possible"],
                                                q_dep_confirmed=ev["confirmed"], q_dep_n=ev["applicable_n"],
                                                **{k: B.value(c, k) for k in METRICS})
        cm = base[scen]["Conservative"]
        if not (cm["m_c"] == 1.0 and cm["q_CBF"] == 0.0 and cm["q_dep"] == 0):
            raise ValueError(f"c_max values differ from m_c = 1, q_CBF = 0, q_dep = 0: {cm}")
    out = dict(generator=S.recorded_path(Path(__file__)), index=S.recorded_path(L.INDEX), cells_file=S.recorded_path(B.CELLS),
               metrics={k: S.recorded_path(v / "metrics.json") for k, v in EXPORTS.items()},
               checks=dict(boxes=len(records), distinct_configurations=len(cfgs), seeds_per_box=N_MC,
                           population_pairs_equal=len(checks), mean_parity_metrics=list(METRICS),
                           mean_parity_max_relative=max_rel, mean_parity_tolerance=1e-9,
                           qdep_counts_equal=len(records), qdep_upper95_max_abs_difference=max_cp,
                           cmax_known_values="m_c/c_max = 1, q_CBF = 0, q_dep = 0 (both scenarios)"),
               population=checks, baselines=base, boxes=records)
    CACHE.write_text(json.dumps(out) + "\n")
    write_csv(out)
    print("cache:", CACHE, "boxes", len(records), "configs", len(cfgs), "max rel", max_rel, "max cp", max_cp,
          flush=True)
    return out


def write_csv(d):
    DUMP_CSV.parent.mkdir(parents=True, exist_ok=True)
    a0, a, dl, K, paired = S.SETTINGS[LABEL]
    with open(DUMP_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["setting", "asset_label", "alpha0", "alpha", "delta", "K", "scenario", "swept_parameter",
                    "x_value", "x_position", "is_setting_value", "variant", "method_key", "source", "family",
                    "view_id", "configuration_id", "seed", "index_source", "assessment_sha256", "m_c_over_cmax",
                    "q_CBF_fraction", "q_CBF_percent", "T_mis_s", "qdep_applicable", "qdep_confirmed",
                    "qdep_possible", "cell_mean_m_c", "cell_mean_q_CBF_fraction", "cell_mean_T_mis_s",
                    "cell_qdep_possible_count", "cell_qdep_confirmed_count", "cell_qdep_upper95"])
        for r in d["boxes"]:
            v = r["values"]
            for i, seed in enumerate(r["seeds"]):
                w.writerow([S.TLABEL[LABEL], LABEL, a0, a, dl, K, r["scen"], r["col"], r["x_label"], r["pos"],
                            r["is_setting"], S.LABEL[r["method"]], r["method"], r["source"], r["family"],
                            r["view_id"], r["configuration_id"], seed, r["index_source"][i],
                            r["assessment_sha256"][i], repr(v["m_c"][i]), repr(v["q_CBF"][i]),
                            repr(100.0 * v["q_CBF"][i]), repr(v["T_mis"][i]), v["qdep_applicable"][i],
                            v["qdep_confirmed"][i], v["qdep_possible"][i], repr(r["metrics_mean"]["m_c"]),
                            repr(r["metrics_mean"]["q_CBF"]), repr(r["metrics_mean"]["T_mis"]),
                            r["qdep"]["possible"], r["qdep"]["confirmed"], repr(r["qdep"]["upper95"])])
    print("csv:", DUMP_CSV, flush=True)


# ----------------------------------------------------------------------------------------------- figure
def light(c, f=0.30):
    r, g, b = mcolors.to_rgb(c)
    return (1 - f * (1 - r), 1 - f * (1 - g), 1 - f * (1 - b))


def setup_text():
    width_in, width_tex, _ = S.setup("paper")
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize",
                "legend.fontsize", "legend.title_fontsize"):
        plt.rcParams[key] = TXT
    # 11 pt text; sub- and superscripts at 10 pt (never below the 10 pt minimum)
    plt.rcParams["text.latex.preamble"] = (r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}"
                                           rf"\DeclareMathSizes{{{TXT:g}}}{{{TXT:g}}}{{10}}{{10}}"
                                           rf"\DeclareMathSizes{{{TXT_SCEN:g}}}{{{TXT_SCEN:g}}}{{10}}{{10}}")
    return width_in, width_tex


def tick_labels(labels, setting_pos):
    out = []
    for i, t in enumerate(labels):
        t = r"\infty" if t == "infinity" else t
        if t == "2.11":   # the default budget normalizer b = kappa (2.11) prints as kappa
            out.append(r"$\bm{\kappa}$" if i == setting_pos else r"$\kappa$")
            continue
        out.append(rf"$\mathbf{{{t}}}$" if i == setting_pos else rf"${t}$")
    return out


def header(col):
    """Column header: what varies; the fixed values of alpha0 and the other three parameters (three lines)."""
    a0, a, d, K, _ = S.SETTINGS[LABEL]
    fixed = dict(alpha0=f"{a0:.2f}", alpha=f"{a:.2f}", delta=f"{d:g}" if d < 0.05 else f"{d:.2f}", K=f"{K}",
                 b=r"\kappa")   # b = kappa by default; delta = 0.025 keeps its third decimal (as the _ext script)
    others = [k for k in ("alpha0", "alpha", "delta", "K", "b") if k != col]
    items = [rf"${PNAME[k]}{{=}}{fixed[k]}$" for k in others]
    return rf"${PNAME[col]}$ varies" + "\n" + ", ".join(items[:2]) + "\n" + ", ".join(items[2:])


def box_stats(v):
    v = np.asarray(v, dtype=float)
    p5, q1, med, q3, p95 = np.percentile(v, [5, 25, 50, 75, 95])
    return dict(n=int(v.size), mean=float(v.mean()), p5=float(p5), q1=float(q1), median=float(med), q3=float(q3),
                p95=float(p95), min=float(v.min()), max=float(v.max()))


def draw_boxes(ax, recs, metric, info):
    for k, m in enumerate(VARIANTS):
        c, mk, _ = S.STYLE[m]
        rs = sorted(recs[m], key=lambda r: r["pos"])
        vals = [SCALE[metric] * np.asarray(r["values"][metric], dtype=float) for r in rs]
        pos = [r["pos"] + (2 * k - 1) * OFFSET for r in rs]
        ax.boxplot(vals, positions=pos, widths=WIDTH, whis=(5, 95), showfliers=False, patch_artist=True,
                   showmeans=True, manage_ticks=False, zorder=3 + k * 0.1,
                   boxprops=dict(facecolor=light(c), edgecolor=c, lw=LW_BOX),
                   whiskerprops=dict(color=c, lw=LW_BOX, solid_capstyle="butt"),
                   capprops=dict(color=c, lw=LW_BOX, solid_capstyle="butt"),
                   medianprops=dict(color="black", lw=LW_BOX, solid_capstyle="butt"),
                   meanprops=dict(marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.7, ls=""))
        for r, v, x in zip(rs, vals, pos):
            st = box_stats(v)
            if abs(st["mean"] - SCALE[metric] * r["metrics_mean"][metric]) > 1e-9 * max(1.0, abs(st["mean"])):
                raise ValueError("drawn mean differs from the export mean")
            info["boxes"].append(dict(scen=r["scen"], col=r["col"], x=r["x_label"], variant=S.LABEL[m],
                                      metric=metric, unit={"m_c": "ratio", "q_CBF": "percent", "T_mis": "s"}[metric],
                                      position=x, view_id=r["view_id"], is_setting=r["is_setting"], **st))


def draw_qdep(ax, recs, info):
    for k, m in enumerate(VARIANTS):
        c, mk, _ = S.STYLE[m]
        for r in sorted(recs[m], key=lambda r: r["pos"]):
            x = r["pos"] + (2 * k - 1) * OFFSET
            y, ub = r["qdep"]["possible"], N_MC * r["qdep"]["upper95"]
            ax.plot([x, x], [y, ub], color=c, lw=LW_BOX, solid_capstyle="butt", zorder=3)
            ax.plot([x - 0.1, x + 0.1], [ub, ub], color=c, lw=LW_BOX, solid_capstyle="butt", zorder=3)
            ax.plot([x], [y], ls="", marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.7, zorder=4)
            info["qdep"].append(dict(scen=r["scen"], col=r["col"], x=r["x_label"], variant=S.LABEL[m], position=x,
                                     count=y, confirmed=r["qdep"]["confirmed"], n=N_MC, upper95=r["qdep"]["upper95"],
                                     upper95_count=ub, delta=r["configuration"]["delta"], view_id=r["view_id"]))


def row_limits(d, scen, metric, h_pt=None):
    """y range of a row: the 5th-95th percentile whiskers of V-L and V-C over all columns plus MARGIN of that span on
    each side, widened so that whisker caps keep CLEAR_PT and mean markers CLEAR_PT from the
    frame at the row's height h_pt (pt; None in the first layout pass). Ticks: none in the top tenth (it meets the
    row above)."""
    if metric == "q_dep":
        top = max(N_MC * r["qdep"]["upper95"] for r in d["boxes"] if r["scen"] == scen)
        if top > QDEP_YLIM[1] - 1.5:
            raise ValueError("q_dep bound near the top of the q_dep axis")
        return QDEP_YLIM, [0, 3, 6]
    st = [box_stats(SCALE[metric] * np.asarray(r["values"][metric], dtype=float))
          for r in d["boxes"] if r["scen"] == scen]
    lo, hi = min(s["p5"] for s in st), max(s["p95"] for s in st)
    pad = MARGIN * (hi - lo)
    if h_pt is None:
        lim = (lo - pad, hi + pad)
    else:
        mlo, mhi = min(s["mean"] for s in st), max(s["mean"] for s in st)
        refs = [d["baselines"][scen][m][metric] for m in REF_BASE if d["baselines"][scen][m][metric] is not None]
        lim = limits_with_refs([(hi + pad, 0.0), (hi, line_pad(LW_BOX)), (mhi, marker_pad(MS_MEAN, 0.7))],
                               [(lo - pad, 0.0), (lo, line_pad(LW_BOX)), (mlo, marker_pad(MS_MEAN, 0.7))], refs,
                               h_pt, S.LW_REF)
    span = lim[1] - lim[0]
    ticks = [float(np.round(t, 6)) for t in MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(*lim)
             if lim[0] <= t <= lim[1] - 0.1 * span]
    return lim, ticks


CLEAR_PT = 3.0           # data keep this distance from every frame
FRAME_HALF_PT = 0.6      # half of the 1.2 pt frame


def pt_limits(tops, bottoms, h_pt, iters=60):
    """Axis limits for an axes h_pt tall such that every (value, pad_pt) of `tops` lies at least pad_pt below the top
    frame and every one of `bottoms` at least pad_pt above the bottom frame (pads in pt; a pad of 0 with a value
    beyond the data gives a margin in data units). Fixed point of span = top - bottom."""
    hi = max(v for v, _ in tops)
    lo = min(v for v, _ in bottoms)
    span = max(hi - lo, 1e-12)
    for _ in range(iters):
        top = max(v + p * span / h_pt for v, p in tops)
        bot = min(v - p * span / h_pt for v, p in bottoms)
        new = top - bot
        if abs(new - span) <= 1e-12 * max(1.0, abs(new)):
            break
        span = new
    return bot, top


def limits_with_refs(tops, bottoms, refs, h_pt, ref_lw):
    """pt_limits, then every reference value (QT, ACI) that would lie within line_pad(ref_lw) of a frame, inside the
    range, is kept that far from the frame by a small widening (a reference line never runs along a frame); values
    beyond the range stay outside (right-hand labels)."""
    tops, bottoms = list(tops), list(bottoms)
    lim = pt_limits(tops, bottoms, h_pt)
    for _ in range(6):
        band = line_pad(ref_lw) * (lim[1] - lim[0]) / h_pt
        add_t = [(v, line_pad(ref_lw)) for v in refs if lim[1] - band < v <= lim[1] and (v, line_pad(ref_lw)) not in tops]
        add_b = [(v, line_pad(ref_lw)) for v in refs if lim[0] <= v < lim[0] + band
                 and (v, line_pad(ref_lw)) not in bottoms]
        if not add_t and not add_b:
            break
        tops += add_t
        bottoms += add_b
        lim = pt_limits(tops, bottoms, h_pt)
    return lim


def marker_pad(ms, mew=0.6):
    """pt between a marker's centre and the frame line so that the marker keeps CLEAR_PT from the frame."""
    return ms / 2 + mew / 2 + CLEAR_PT + FRAME_HALF_PT


def line_pad(lw):
    return lw / 2 + CLEAR_PT + FRAME_HALF_PT


def fmt_base(metric, v):
    """Baseline value in the unit of the row's y label (the label carries the unit)."""
    return {"m_c": f"{v:.2f}", "q_CBF": f"{v:.1f}", "T_mis": f"{v:.1f}", "q_dep": f"{int(round(v))}"}[metric]


def fmt_label(metric, m, v, where, lines=2):
    """Baseline label: name (with an arrow towards the value when it is off the axis) and value."""
    arrow = {"line": "", "up": r" $\blacktriangle$", "down": r" $\blacktriangledown$"}[where]
    return f"{S.LABEL[m]}{arrow}" + ("\n" if lines == 2 else " ") + fmt_base(metric, v)


def free_box(ax, bb, pad_px, skip=()):
    """True when the display box bb lies inside the axes and covers no line point, no box patch and no text."""
    fig = ax.figure
    r = fig.canvas.get_renderer()
    abb = ax.get_window_extent(r)
    if bb.x0 < abb.x0 + 1 or bb.x1 > abb.x1 - 1 or bb.y0 < abb.y0 + 1 or bb.y1 > abb.y1 - 1:
        return False
    for ln in ax.get_lines():
        if not ln.get_visible():
            continue
        p = S._display_points(ax, ln)
        ms_px = ln.get_markersize() * fig.dpi / 72.0 / 2 if ln.get_marker() not in (None, "None", "", " ") else 0
        m = pad_px + ms_px + ln.get_linewidth() * fig.dpi / 72.0 / 2
        if ((p[:, 0] > bb.x0 - m) & (p[:, 0] < bb.x1 + m) & (p[:, 1] > bb.y0 - m) & (p[:, 1] < bb.y1 + m)).any():
            return False
    for p in ax.patches:
        if p.get_visible() and p.get_zorder() > 0 and p.get_window_extent(r).padded(pad_px).overlaps(bb):
            return False
    for t in ax.texts:
        if t not in skip and t.get_visible() and t.get_text().strip() and t.get_window_extent(r).overlaps(bb):
            return False
    return True


def place_annotation(row_axes, text, up):
    """In-panel edge annotation (q_dep rows, whose upper part is free): top edge (up) or bottom edge, first free
    corner of the last, then earlier panels; checked against every line point, marker and box."""
    fig = row_axes[0].figure
    pad_px = 2.0 * fig.dpi / 72.0
    y, va = (0.96, "top") if up else (0.04, "bottom")
    for ax in (row_axes[3], row_axes[2], row_axes[1], row_axes[0]):
        for x, ha in ((0.98, "right"), (0.02, "left"), (0.5, "center")):
            t = ax.text(x, y, text, transform=ax.transAxes, ha=ha, va=va, fontsize=TXT, zorder=6)
            fig.canvas.draw()
            if free_box(ax, t.get_window_extent(fig.canvas.get_renderer()), pad_px, skip=(t,)):
                return ax, t, ha
            t.remove()
    raise ValueError(f"no free edge for the annotation {text!r}")


def place_right_labels(fig, axes, labels, gap_pt=5.0, sep_pt=4.0, edge_pt=4.0):
    """Baseline labels right of the last panel of each row, in the order: values above the range (arrow up, highest
    first, from the top of the row), values inside it (at the line's height), values below it (arrow down, from the
    bottom of the row); sep_pt apart, never closer than edge_pt to the row's top or bottom edge (so that labels of
    adjacent rows do not meet at the row boundary), gap_pt right of
    the frame (all in pt). Two passes: top-down (each label at most at its target), then bottom-up."""
    r = fig.canvas.get_renderer()
    Wp, Hp = fig.get_size_inches() * fig.dpi
    px = fig.dpi / 72.0
    sep, edge = sep_pt * px, edge_pt * px
    for i, items in labels.items():
        ax = axes[i, -1]
        bb = ax.get_window_extent(r)
        x = bb.x1 + gap_pt * px
        hs = [t.get_window_extent(r).height for t, _, _ in items]
        top_b, bot_b = bb.y1 - edge, bb.y0 + edge
        target, key = {}, {}
        for k, ((t, where, v), h) in enumerate(zip(items, hs)):
            if where == "line":
                target[k] = ax.transData.transform((0, v))[1]
                key[k] = (1, -target[k])
            elif where == "up":
                target[k] = top_b - h / 2
                key[k] = (0, -v)
            else:
                target[k] = bot_b + h / 2
                key[k] = (2, -v)
        order = sorted(range(len(items)), key=lambda k: key[k])
        ys = {}
        prev = None
        for k in order:                         # top-down: at most at the target, below the previous label
            y = min(target[k], top_b - hs[k] / 2)
            if prev is not None:
                y = min(y, ys[prev] - (hs[prev] + hs[k]) / 2 - sep)
            ys[k] = y
            prev = k
        nxt = None
        for k in reversed(order):               # bottom-up: inside the row, above the next label
            y = max(ys[k], bot_b + hs[k] / 2)
            if nxt is not None:
                y = max(y, ys[nxt] + (hs[nxt] + hs[k]) / 2 + sep)
            ys[k] = y
            nxt = k
        if ys[order[0]] > top_b - hs[order[0]] / 2 + 0.5:
            raise ValueError(f"row {i}: baseline labels do not fit")
        for k, (t, _, _) in enumerate(items):
            t.set_position((x / Wp, ys[k] / Hp))


def separate_tick_labels(fig, bottom_axes, gap_pt=3.0):
    """Where the first two x tick labels of a column are closer than gap_pt (delta: "0.025" next to "0.05"), move
    the first label left by the missing amount (it has the free half category and the panel gap on its left)."""
    r = fig.canvas.get_renderer()
    fig.canvas.draw()
    done = []
    for ax in bottom_axes:
        labs = [t for t in ax.get_xticklabels() if t.get_visible() and t.get_text()]
        if len(labs) < 2:
            continue
        a, b = labs[0].get_window_extent(r), labs[1].get_window_extent(r)
        gap = (b.x0 - a.x1) * 72.0 / fig.dpi
        if gap < gap_pt:
            dx = gap_pt - gap       # 30 % to the left (the previous column's last label is near), 70 % to the right
            labs[0].set_transform(labs[0].get_transform() + ScaledTranslation(-0.3 * dx / 72.0, 0, fig.dpi_scale_trans))
            labs[1].set_transform(labs[1].get_transform() + ScaledTranslation(0.7 * dx / 72.0, 0, fig.dpi_scale_trans))
            done.append((labs[0].get_text(), round(-0.3 * dx, 2)))
            done.append((labs[1].get_text(), round(0.7 * dx, 2)))
    fig.canvas.draw()
    for ax in bottom_axes:      # verify after the redraw (the shift must survive it)
        labs = [t for t in ax.get_xticklabels() if t.get_visible() and t.get_text()]
        for a, b in zip(labs[:-1], labs[1:]):
            g = (b.get_window_extent(r).x0 - a.get_window_extent(r).x1) * 72.0 / fig.dpi
            if g < gap_pt - 0.05:
                raise ValueError(f"tick labels {a.get_text()} and {b.get_text()} only {g:.2f} pt apart")
    return done


def check_cross_gaps(fig, bottom_axes, gap_pt=3.0):
    """The last tick label of one column and the first of the next are at least gap_pt apart (0.25 | 0.025)."""
    r = fig.canvas.get_renderer()
    fig.canvas.draw()
    for ax1, ax2 in zip(bottom_axes[:-1], bottom_axes[1:]):
        a = [t for t in ax1.get_xticklabels() if t.get_visible() and t.get_text()][-1]
        b = [t for t in ax2.get_xticklabels() if t.get_visible() and t.get_text()][0]
        g = (b.get_window_extent(r).x0 - a.get_window_extent(r).x1) * 72.0 / fig.dpi
        print("cross gap", a.get_text(), b.get_text(), round(g, 2), flush=True)
        if g < gap_pt - 0.05:
            raise ValueError(f"tick labels {a.get_text()} and {b.get_text()} of neighbouring columns only {g:.2f} pt apart")


def legend_handles():
    def box_handle(m):
        c, mk, _ = S.STYLE[m]
        return (Patch(facecolor=light(c), edgecolor=c, lw=LW_BOX),
                Line2D([], [], color="black", lw=LW_BOX, solid_capstyle="butt"),
                Line2D([], [], ls="", marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.7))
    handles = [box_handle(m) for m in VARIANTS]
    labels = [S.LABEL[m] for m in VARIANTS]
    for m in REF_BASE:
        h = Line2D([], [], color=S.STYLE[m][0], ls=REF_LS[m], lw=S.LW_REF)
        if m == "QT":
            h.set_path_effects(QT_EFFECT)
        handles.append(h)
        labels.append(S.LABEL[m] + rf" ($\alpha_0{{=}}{S.SETTINGS[LABEL][0]:.2f}$)")   # alpha0 of the setting
    handles.append(Patch(facecolor=BAND, edgecolor="none"))
    labels.append(r"value at S3 (bold)")
    return handles, labels


def build(d, out: Path):
    """Two layout passes: the first measures each row's height in pt, the second sets the y ranges from it."""
    heights = _build(d, out, None)
    return _build(d, out, heights)


def _build(d, out: Path, row_pt):
    width_in, width_tex = setup_text()
    ncat = {c: max(r["pos"] for r in d["boxes"] if r["col"] == c) + 1 for c in COLS}
    rows = [(s, m) for s in SCENS for m in ROWS]
    blk_hr = [R_DEP if m == "q_dep" else R_BOX for m in ROWS]
    hr = blk_hr + [R_GAP] + blk_hr
    grid_rows = list(range(len(blk_hr))) + [len(blk_hr) + 1 + k for k in range(len(blk_hr))]
    full_hr = [R_BOX] * 3 + [R_DEP]
    height_pt = HEIGHT_PT - (sum(full_hr) - sum(blk_hr)) * 2 * HEIGHT_PER_UNIT   # same box-row height without q_dep rows
    fig = plt.figure(figsize=(width_in, height_pt / 72.27))
    gs = fig.add_gridspec(len(hr), len(COLS), height_ratios=hr, width_ratios=[COL_RATIO[c] for c in COLS],
                          hspace=HSPACE, wspace=WSPACE)
    axes = np.empty((len(rows), len(COLS)), dtype=object)
    for i, gi in enumerate(grid_rows):
        for j in range(len(COLS)):
            axes[i, j] = fig.add_subplot(gs[gi, j], sharex=axes[0, j] if i else None,
                                         sharey=axes[i, 0] if j else None)
    info = dict(figure=NAME, version=2, setting=LABEL, tuple=S.SETTINGS[LABEL][:4],
                variants=[S.LABEL[m] for m in VARIANTS], rows=[f"{s}:{m}" for s, m in rows], cols=list(COLS),
                override=dict(text_pt=TXT, scenario_pt=TXT_SCEN, box_lw=LW_BOX, mean_marker_pt=MS_MEAN,
                              legend="one row above the panels, outside the axes"),
                whiskers="5th and 95th percentiles", boxes=[], qdep=[], ranges={}, reference_lines=[],
                annotations=[], baselines=d["baselines"])
    right = {}
    for i, (scen, metric) in enumerate(rows):
        lim, ticks = row_limits(d, scen, metric, None if row_pt is None else row_pt[i])
        info["ranges"][f"{scen}:{metric}"] = dict(ylim=list(lim), ticks=ticks)
        for j, col in enumerate(COLS):
            ax = axes[i, j]
            recs = {m: [r for r in d["boxes"] if r["scen"] == scen and r["col"] == col and r["method"] == m]
                    for m in VARIANTS}
            n = len(recs[VARIANTS[0]])
            sp = [r["pos"] for r in recs[VARIANTS[0]] if r["is_setting"]]
            for p in sp:
                ax.axvspan(p - 0.5, p + 0.5, color=BAND, lw=0, zorder=0)
            for m in REF_BASE:
                v = d["baselines"][scen][m][metric]
                if lim[0] <= v <= lim[1]:
                    ln = ax.axhline(v, color=S.STYLE[m][0], ls=REF_LS[m], lw=S.LW_REF, zorder=2)
                    if m == "QT":
                        ln.set_path_effects(QT_EFFECT)
                    if j == 0:
                        info["reference_lines"].append(dict(row=f"{scen}:{metric}", baseline=m, value=v))
                        right.setdefault(i, []).append((m, v, "line"))
                elif j == 0:
                    right.setdefault(i, []).append((m, v, "up" if v > lim[1] else "down"))
            if metric == "q_dep":
                draw_qdep(ax, recs, info)
            else:
                draw_boxes(ax, recs, metric, info)
            ax.set_xlim(-0.5, n - 0.5)
            labels = [r["x_label"] for r in sorted(recs[VARIANTS[0]], key=lambda r: r["pos"])]
            ax.set_xticks(range(n), tick_labels(labels, sp[0] if sp else -1))
            ax.tick_params(axis="x", length=0, pad=3)
            ax.grid(False)
            ax.grid(True, which="major", axis="y")
            ax.set_xlabel(B.XLAB[col])
            ax.set_ylabel(YLAB[metric])
            if j > 0:
                ax.yaxis.label.set_visible(False)
                ax.tick_params(labelleft=False)
            if i != len(rows) - 1:
                ax.xaxis.label.set_visible(False)
            # the swept values under the last row of each block (the q_dep row, or T_mis when the all-zero q_dep
            # rows are omitted)
            ax.tick_params(labelbottom=(metric == ROWS[-1]))
            if i == 0:
                ax.set_title(header(col), fontsize=TXT, pad=5, linespacing=1.15)
        axes[i, 0].set_ylim(*lim)
        axes[i, 0].set_yticks(ticks, [f"{t:g}" for t in ticks])
    # q_dep bars must stay below every delta of the figure (the caption says so)
    dmin = min(r["configuration"]["delta"] for r in d["boxes"])
    if max(r["qdep"]["upper95"] for r in d["boxes"]) > dmin:
        raise ValueError("a q_dep upper bound exceeds the smallest delta")
    for blk in range(len(SCENS)):
        fig.align_ylabels(axes[blk * len(ROWS):(blk + 1) * len(ROWS), 0])
    scen_txt = [fig.text(0, 0.5, rf"\textbf{{{S.SCEN[s]}}}", rotation=90, ha="left", va="center",
                         fontsize=TXT_SCEN) for s in SCENS]
    rlabels, inpanel = {}, {}
    for i, items in right.items():
        metric = rows[i][1]
        if metric == "q_dep":           # the upper part of the q_dep panels is free: in-panel annotation
            inpanel[i] = items
            continue
        rlabels[i] = [(fig.text(0, 0, fmt_label(metric, m, v, where), ha="left", va="center", fontsize=TXT,
                                linespacing=1.0), where, v)
                      for m, v, where in sorted(items, key=lambda z: -z[1])]
        for (t, where, v), (m, _, _) in zip(rlabels[i], sorted(items, key=lambda z: -z[1])):
            info["annotations"].append(dict(row=info["rows"][i], baseline=m, value=v, kind=where, text=t.get_text()))
    handles, labels = legend_handles()
    leg = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=len(labels),
                     frameon=False, fontsize=TXT, handlelength=1.6, handleheight=0.9, handletextpad=0.35,
                     columnspacing=1.0, borderaxespad=0.0, borderpad=0.1,
                     handler_map={tuple: HandlerTuple(ndivide=1, pad=0.0)})
    if os.environ.get("SWEEP_DEBUG"):
        fig.subplots_adjust(left=0.12, right=0.99, bottom=0.05, top=0.90)
        fig.savefig(out / (NAME + "_debug.png"), dpi=150)
        plt.close(fig)
        return None
    layout(fig, axes, leg, scen_txt, rlabels)
    if row_pt is None:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        heights = [axes[i, 0].get_window_extent(r).height * 72.0 / fig.dpi for i in range(len(rows))]
        plt.close(fig)
        return heights
    for i, items in inpanel.items():
        metric = rows[i][1]
        if any(w == "line" for _, _, w in items) or len({w for _, _, w in items}) > 1:
            raise ValueError("q_dep baselines expected on one side of the axis")
        up = items[0][2] == "up"
        txt = ", ".join(fmt_label(metric, m, v, "line", lines=1) for m, v, _ in sorted(items)) + \
            (r" $\blacktriangle$" if up else r" $\blacktriangledown$")
        ax, t, ha = place_annotation(axes[i], txt, up)
        info["annotations"].append(dict(row=info["rows"][i], baselines={m: v for m, v, _ in items},
                                        kind="up" if up else "down", text=txt, panel=COLS[list(axes[i]).index(ax)],
                                        ha=ha))
    shifts = separate_tick_labels(fig, axes[len(ROWS) - 1, :]) + separate_tick_labels(fig, axes[-1, :])
    check_cross_gaps(fig, axes[len(ROWS) - 1, :]); check_cross_gaps(fig, axes[-1, :])
    bad = S.text_overlaps(fig)
    if bad:
        raise ValueError(f"overlapping texts {bad[:6]}")
    pdf = fs.save(fig, out / NAME, png=False, strict=True)
    plt.close(fig)
    w, h = S.check_width(pdf, width_tex)
    if h > HEIGHT_CAP + 0.5:
        raise ValueError(f"height {h} pt above {HEIGHT_CAP}")
    info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=h, text_pt=TXT,
                tick_label_shifts_pt=shifts)
    print(NAME, f"{w} x {h} TeX pt", "labels:", [(a["row"], a["text"]) for a in info["annotations"]], flush=True)
    return info


def layout(fig, axes, leg, scen_txt, rlabels, gap_pt=4.0, iters=12):
    """Subplot box such that the scenario names touch the left edge, the bottom labels the bottom edge, the
    baseline labels the right edge, the legend the top edge, and the column headers sit gap_pt below the legend."""
    fig.subplots_adjust(left=0.12, right=0.99, bottom=0.05, top=0.90)
    W, H = fig.get_size_inches()
    for _ in range(iters):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        Wp, Hp = W * fig.dpi, H * fig.dpi
        for blk, t in enumerate(scen_txt):
            rows = axes[blk * len(ROWS):(blk + 1) * len(ROWS), 0]
            x0 = min(ax.yaxis.label.get_window_extent(r).x0 for ax in rows)
            ytop, ybot = rows[0].get_window_extent(r).y1, rows[-1].get_window_extent(r).y0
            tw = t.get_window_extent(r).width
            t.set_position(((x0 - gap_pt * fig.dpi / 72 - tw) / Wp, 0.5 * (ytop + ybot) / Hp))
        place_right_labels(fig, axes, rlabels)
        fig.canvas.draw()
        bb = fig.get_tightbbox(r)
        sp = fig.subplotpars
        dl, dr, db = -bb.x0 / W, (W - bb.x1) / W, -bb.y0 / H
        title_top = max(ax.title.get_window_extent(r).y1 for ax in axes[0]) / Hp
        leg_bottom = leg.get_window_extent(r).y0 / Hp
        dt = (leg_bottom - gap_pt / 72.0 / H) - title_top
        if max(abs(dl), abs(dr), abs(db), abs(dt)) < 1e-4:
            break
        fig.subplots_adjust(left=sp.left + dl, right=sp.right + dr, bottom=sp.bottom + db, top=sp.top + dt)


def main():
    global ROWS, NAME
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=HERE / "out_paper")
    ap.add_argument("--reload", action="store_true", help="re-read the assessments and redo the checks")
    ap.add_argument("--load-only", action="store_true")
    ap.add_argument("--omit-zero-rows", action="store_true",
                    help="omit rows whose plotted variants are all zero (the q_dep rows; arXiv figure)")
    ap.add_argument("--name", default=NAME, help="output file stem (default sweepbox_T3_compact)")
    a = ap.parse_args()
    d = load(a.reload)
    NAME = a.name
    if a.omit_zero_rows:
        zero = [m for m in ROWS if m == "q_dep" and all(
            r["qdep"]["possible"] == 0 and r["qdep"]["confirmed"] == 0 for r in d["boxes"] if r["method"] in VARIANTS)]
        # box rows are never omitted here (their values are not all zero); only a verified all-zero q_dep row is
        ROWS = tuple(m for m in ROWS if m not in zero)
        print("omitted all-zero rows:", zero, flush=True)
    if a.load_only:
        return
    a.out.mkdir(parents=True, exist_ok=True)
    info = build(d, a.out)
    if info is None:
        return
    info["checks"] = d["checks"]
    (a.out / ("sweeps_manifest_" + ("paper_compact_box" if NAME == "sweepbox_T3_compact" else NAME) + ".json")).write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(B.CELLS): S.sha_file(B.CELLS),
                                                                 S.recorded_path(CACHE): S.sha_file(CACHE),
                                                                 S.recorded_path(DUMP_CSV): S.sha_file(DUMP_CSV)},
        figures=[info]), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
