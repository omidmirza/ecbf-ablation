#!/usr/bin/env python3
"""Box-plot sweep figures with all four variants, extended edition and report, in one uniform layout.

sweepbox_{scen}_{label}.pdf (label Default, New1..New5 = T1..T5), one figure per scenario and setting, in the visual
language of the arXiv box figure (build_sweeps_box.py): box plots over the N_MC = 200 per-trajectory
values, no connecting lines, per-row zoom from the variants' whiskers, fixed parameters in the column headers, ACI/QT
reference lines inside the range or right-hand edge labels with an arrow, one legend strip above all panels.
Paper style (extended edition) and report style: Default and T1-T5. Width 505.89 TeX pt (both).

Layout: two bands of columns (one band of five or six columns would leave about 3 pt per box):
  band 1: alpha, delta, K (10, 50, 100, 200, inf);
  band 2: b, phi at L_mu = 0.0005, phi at L_mu = 0.005 (separate column; phi < 1: Assumption 2 may fail).
Rows per band: m_{c,k}/c_max, q+_{CBF,k} (%), q+_{unsafe,k} (%) [only if some value is non-zero], T_{mis,k} (s),
q_cert and q_dep counts [each only if some count is non-zero]. Boxes: V-U, V-Delta, V-L, V-C side by side in the
paper colours; box = quartiles, black median, whiskers = 5th and 95th percentiles, no outliers, mean = the variant's
marker. Count rows: marker at the conservative count (confirmed + numerically possible events) with a bar to the
one-sided 95 % Clopper-Pearson upper bound (in trajectories); the level N_MC*delta (finite K) or N_MC*delta_1,
delta_1 = delta/(2 b ln^2 2) (K = inf), as a short dashed segment per value where it lies inside the row's range.
The setting's own value: grey column band and bold tick label. Baselines QT and ACI at the setting's alpha0
(matched_baselines_200). c_max is not drawn (m_c/c_max = 1, q_CBF = 0, q_unsafe = 0, q_dep = 0; checked).

Style (as in build_sweeps_box.py): text 11 pt (paper; sub- and superscripts 10 pt) or
11 pt throughout (report body 11 pt), mean markers 5 pt with one shape per variant, box lines 1.3 pt,
legend strip above the panels outside the axes. Everything else: vector PDF, usetex Computer Modern, strict helper
save (zero margin, registered width), text-overlap check.

Data (read only, no simulation): sweep cells of build_sweeps.Data (sweep_cells.json); per-trajectory assessments
read through N200_INDEX.sqlite and hash-checked (canonical SHA-256 = index), as in build_sweeps_box.read_assessments.
Checks: seeds 970001-970200; (seed, assessment SHA-256) pairs equal population.jsonl; means of m_c, q_CBF, q_unsafe,
T_mis equal metrics.json available_mean to 1e-9 relative; q_cert and q_dep counts (possible, confirmed) equal
trajectory_events and the Clopper-Pearson bound equals upper95_possible to 1e-9.
Cache: sweepbox_ext_values.json (next to this script; it holds the same per-trajectory values). CSV, written whenever
the cache is rebuilt (no cache or --reload): sweeps/sweepbox_ext_values.csv of the display data (ECBF_DISPLAYS, default
data/; one row per configuration and seed; this repository ships it gzip-compressed as
data/sweeps/sweepbox_ext_values.csv.gz).

Usage: python3 build_sweeps_box_ext.py [--reload] [--load-only] --kind paper|report --out DIR
       [--labels Default New3] [--scens corridor maze]
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import sqlite3
import statistics
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs
import build_sweeps as B
import build_sweeps_box as X
import load_setting_assessments as L

from matplotlib.legend_handler import HandlerTuple  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.ticker import MaxNLocator  # noqa: E402

HERE = Path(__file__).resolve().parent
VARIANTS = S.CERTIFIED                       # V-U, V-D, V-Cd (V-L), V-C
SCENS = ("corridor", "maze")
LABELS = tuple(S.SETTINGS)
COLS = ("alpha", "delta", "K", "b", "phi", "phi5")
BANDS = (("alpha", "delta", "K"), ("b", "phi", "phi5"))
BOXM = ("m_c", "q_CBF", "q_unsafe", "T_mis")
CNTM = ("q_cert", "q_dep")
ALL_ROWS = ("m_c", "q_CBF", "q_unsafe", "T_mis", "q_cert", "q_dep")
EVK = {"q_cert": "qcert", "q_dep": "qdep"}
SCALE = {"m_c": 1.0, "q_CBF": 100.0, "q_unsafe": 100.0, "T_mis": 1.0}
CACHE = HERE / "sweepbox_ext_values.json"
DUMP_CSV = X.DISPLAYS / "sweeps/sweepbox_ext_values.csv"
EXPORTS = X.EXPORTS
SEEDS = L.SEEDS
N_MC = 200
LEVEL_COLOUR = "#D55E00"   # theorem-check level colour

TXT_SCEN = 12.0
LW_BOX = 1.3
MS_MEAN = 5.0
R_BOX, R_CNT, R_GAP = 1.0, 0.78, 0.62
HEIGHT_CAP = 640.0
MARGIN = 0.07
SLOT, WIDTH = 0.205, 0.17                  # box centre spacing and box width in category units (four boxes)
BAND_COLOUR = X.BAND
WSPACE, HSPACE = 0.05, 0.10
YLAB = {"m_c": r"$m_{c,k}/c_{\max}$", "q_CBF": "$q^{+}_{\\mathrm{CBF},k}$\n(\\%)",
        "q_unsafe": "$q^{+}_{\\mathrm{unsafe},k}$\n(\\%)", "T_mis": "$T_{\\mathrm{mis},k}$\n(s)",
        "q_cert": "$\\widehat q_{\\mathrm{cert}}$\n(count)", "q_dep": "$\\widehat q_{\\mathrm{dep}}$\n(count)"}
XLAB = {"alpha": r"$\alpha$", "delta": r"$\delta$", "K": r"$K$", "b": r"$b$", "phi": r"$\phi$", "phi5": r"$\phi$"}
PNAME = {"alpha0": r"\alpha_0", "alpha": r"\alpha", "delta": r"\delta", "K": "K", "b": "b", "phi": r"\phi",
         "L_mu": r"L_\mu"}


# ----------------------------------------------------------------------------------------------- data
def sweep_boxes(data, labels):
    out = []
    for label in labels:
        for scen in SCENS:
            for col in COLS:
                for m in VARIANTS:
                    if col in ("phi", "phi5"):
                        pts = data.sweep(label, scen, "phi", m)[0.0005 if col == "phi" else 0.005]
                        sx = 1.0 if col == "phi" else None
                    else:
                        pts = data.sweep(label, scen, col, m)
                        sx = B.setting_x(label, col)
                    for pos, (x, cell) in enumerate(pts):
                        out.append(dict(label=label, scen=scen, col=col, pos=pos, x=x,
                                        x_label=X.x_label("K" if col == "K" else col, x), method=m,
                                        is_setting=sx is not None and abs(float(x) - float(sx)) < 1e-12, cell=cell))
    return out


def read_assessments(cfgs):
    """{canonical config: [(seed, source, assessment_sha256, assessment)]}: one scan of the index, then the saved
    assessment per seed, hash-checked (method of load_setting_assessments.main)."""
    idx = sqlite3.connect(f"file:{L.INDEX}?mode=ro&immutable=1", uri=True)
    n100 = sqlite3.connect(f"file:{L.N100_CATALOG}?mode=ro&immutable=1", uri=True)
    rows = {k: [] for k in cfgs}
    for r in idx.execute("SELECT configuration_json, seed, source, location, raw_sha256, assessment_sha256 FROM runs"):
        if r[0] in rows:
            rows[r[0]].append(r[1:])
    print("index scanned", flush=True)
    n100_cid = {}
    for cid, cj in n100.execute("SELECT configuration_id, configuration_json FROM configurations"):
        key = L.canonical({k: v for k, v in json.loads(cj).items() if v is not None})
        if key in cfgs:
            n100_cid[key] = cid
    hist = None
    out = {}
    for i, key in enumerate(sorted(cfgs)):
        rs = sorted(rows[key])
        if [r[0] for r in rs] != SEEDS:
            raise ValueError(f"index rows for {key}: {len(rs)}")
        n100_rows = None
        per = []
        for seed, source, location, raw, asha in rs:
            if source == "n100_prefix_sealed":
                if n100_rows is None:
                    n100_rows = {}
                    for s, t, rr in n100.execute("SELECT seed, assessment_json, raw_sha256 FROM observations WHERE "
                                                 "configuration_id=?", (n100_cid[key],)):
                        n100_rows.setdefault(s, []).append((t, rr))
                cands = [json.loads(t) for (t, rr) in n100_rows.get(seed, []) if rr == raw]
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
        if i % 50 == 0:
            print(f"read {i + 1}/{len(cfgs)}", flush=True)
    return out


def metrics_cells(view_ids):
    out = {}
    for name, folder in EXPORTS.items():
        d = json.loads((folder / "metrics.json").read_text())
        for v in d["views"]:
            if (name, v.get("view_id")) not in view_ids:
                continue
            for c in v["cells"]:
                cfg = L.canonical({k: x for k, x in c["configuration"].items() if x is not None})
                out[(name, v["view_id"], cfg)] = dict(n=c["n"], configuration_id=c.get("configuration_id"),
                                                      stats={k: c["statistics"][k] for k in BOXM},
                                                      events={k: c["trajectory_events"][EVK[k]] for k in CNTM})
        del d
        print("metrics read:", name, flush=True)
    return out


def load(labels, reload=False):
    if CACHE.exists() and not reload:
        d = json.loads(CACHE.read_text())
        if set(labels) <= set(d["labels"]):
            return d
    labels = LABELS        # the cache always holds every setting
    data = B.Data()
    boxes = sweep_boxes(data, labels)
    cfgs = {}
    for b in boxes:
        cfg = L.canonical({k: v for k, v in b["cell"]["configuration"].items() if v is not None})
        b["cfg"] = cfg
        b["configuration_id"] = L.sha(cfg)
        cfgs[cfg] = b["configuration_id"]
    print("boxes", len(boxes), "configurations", len(cfgs), flush=True)
    per = read_assessments(set(cfgs))
    pairs, views = X.population_pairs(set(cfgs.values()))
    for cfg, cid in cfgs.items():
        idx_pairs = {s: a for s, _, a, _ in per[cfg]}
        pop = pairs[cid]
        if not (sorted(pop) == SEEDS and all(pop[s] == {idx_pairs[s]} for s in SEEDS)):
            raise ValueError(f"population pairs differ for {cfg}")
    mcells = metrics_cells({(b["cell"]["source"], b["cell"]["view_id"]) for b in boxes})
    configs = {}
    max_rel, max_cp = 0.0, 0.0
    for cfg, cid in cfgs.items():
        assess = [a for _, _, _, a in per[cfg]]
        vals = {k: [a[k] for a in assess] for k in BOXM}
        for k in CNTM:
            ek = EVK[k]
            vals[k + "_applicable"] = [bool(a[ek + "_applicable"]) for a in assess]
            vals[k + "_confirmed"] = [bool(a[ek + "_confirmed"]) for a in assess]
            vals[k + "_possible"] = [bool(a[ek + "_possible"]) for a in assess]
        configs[cid] = dict(configuration=json.loads(cfg), seeds=[s for s, _, _, _ in per[cfg]],
                            index_source=[s for _, s, _, _ in per[cfg]],
                            assessment_sha256=[a for _, _, a, _ in per[cfg]], values=vals)
    records = []
    for b in boxes:
        cell = b["cell"]
        mc = mcells.get((cell["source"], cell["view_id"], b["cfg"]))
        if mc is None:
            raise ValueError(f"metrics cell missing: {cell['source']} {cell['view_id']}")
        if mc["n"] != N_MC or (mc["configuration_id"] and mc["configuration_id"] != b["configuration_id"]):
            raise ValueError(f"n or configuration id differs: {cell['view_id']}")
        vals = configs[b["configuration_id"]]["values"]
        means = {}
        for k in BOXM:
            st = mc["stats"][k]
            v = [x for x in vals[k] if x is not None]
            if not st.get("applicable", True):
                means[k] = None
                continue
            want = X._mean_of(st)
            have = statistics.mean(v) if v else None
            if have is None or want is None:
                if have != want:
                    raise ValueError(f"mean parity {cell['view_id']} {k}: {have} vs {want}")
                means[k] = None
                continue
            rel = abs(have - want) / max(abs(want), 1e-300) if want else abs(have - want)
            if rel > 1e-9:
                raise ValueError(f"mean parity {cell['view_id']} {k}: {have} vs {want}")
            max_rel = max(max_rel, rel)
            means[k] = want
        counts = {}
        for k in CNTM:
            ev = mc["events"][k]
            app = sum(vals[k + "_applicable"])
            if any(c and not p for c, p in zip(vals[k + "_confirmed"], vals[k + "_possible"])):
                raise ValueError(f"confirmed but not possible {k} event in {cell['view_id']}")
            qd = dict(possible=sum(p for p, a in zip(vals[k + "_possible"], vals[k + "_applicable"]) if a),
                      confirmed=sum(c for c, a in zip(vals[k + "_confirmed"], vals[k + "_applicable"]) if a),
                      applicable_n=app)
            if (qd["possible"], qd["confirmed"], qd["applicable_n"]) != (ev.get("possible"), ev.get("confirmed"),
                                                                         ev.get("applicable_n")):
                raise ValueError(f"{k} counts differ in {cell['view_id']}: {qd} vs {ev}")
            if app:
                qd["upper95"] = X.cp_upper(qd["possible"], app)
                if abs(qd["upper95"] - ev["upper95_possible"]) > 1e-9:
                    raise ValueError(f"{k} bound differs in {cell['view_id']}")
                max_cp = max(max_cp, abs(qd["upper95"] - ev["upper95_possible"]))
            else:
                qd["upper95"] = None
            counts[k] = qd
        records.append(dict(label=b["label"], scen=b["scen"], col=b["col"], pos=b["pos"], x=b["x"],
                            x_label=b["x_label"], method=b["method"], is_setting=b["is_setting"],
                            source=cell["source"], family=cell["family"], view_id=cell["view_id"],
                            configuration_id=b["configuration_id"], configuration=json.loads(b["cfg"]),
                            metrics_mean=means, counts=counts))
    base = {}
    for label in labels:
        a0 = S.SETTINGS[label][0]
        for scen in SCENS:
            for m in X.REF_BASE + ("Conservative",):
                c = data.baseline(scen, a0, m)
                ev = c["events"]["qdep"]
                base.setdefault(f"{a0:g}", {}).setdefault(scen, {})[m] = dict(
                    view_id=c["view_id"], q_dep=ev["possible"], q_dep_confirmed=ev["confirmed"],
                    q_dep_n=ev["applicable_n"], **{k: B.value(c, k) for k in BOXM})
            cm = base[f"{a0:g}"][scen]["Conservative"]
            if not (cm["m_c"] == 1.0 and cm["q_CBF"] == 0.0 and cm["q_unsafe"] == 0.0 and cm["q_dep"] == 0):
                raise ValueError(f"c_max values differ from the known values: {cm}")
    out = dict(generator=S.recorded_path(Path(__file__)), index=S.recorded_path(L.INDEX), cells_file=S.recorded_path(B.CELLS),
               labels=list(labels),
               checks=dict(boxes=len(records), distinct_configurations=len(cfgs), seeds_per_configuration=N_MC,
                           population_pairs_equal=len(cfgs), mean_parity_metrics=list(BOXM),
                           mean_parity_max_relative=max_rel, mean_parity_tolerance=1e-9,
                           count_parity=list(CNTM), upper95_max_abs_difference=max_cp,
                           cmax_known_values="m_c/c_max = 1, q_CBF = 0, q_unsafe = 0, q_dep = 0"),
               baselines=base, configs=configs, boxes=records)
    CACHE.write_text(json.dumps(out) + "\n")
    write_csv(out)
    print("cache:", CACHE, "boxes", len(records), "configs", len(cfgs), "max rel", max_rel, "max cp", max_cp,
          flush=True)
    return out


def write_csv(d):
    DUMP_CSV.parent.mkdir(parents=True, exist_ok=True)
    with open(DUMP_CSV, "w", newline="") as f:
        w = csv.writer(f)
        keys = ["alpha0", "alpha", "delta", "K", "budget_normalizer", "bound_scale", "L_mu", "benchmark", "method"]
        w.writerow(["configuration_id"] + keys + ["seed", "index_source", "assessment_sha256", "m_c_over_cmax",
                                                  "q_CBF_fraction", "q_unsafe_fraction", "T_mis_s",
                                                  "qcert_applicable", "qcert_confirmed", "qcert_possible",
                                                  "qdep_applicable", "qdep_confirmed", "qdep_possible"])
        for cid, c in sorted(d["configs"].items()):
            v, cf = c["values"], c["configuration"]
            for i, seed in enumerate(c["seeds"]):
                w.writerow([cid] + [cf.get(k) for k in keys] + [seed, c["index_source"][i],
                           c["assessment_sha256"][i], v["m_c"][i], v["q_CBF"][i], v["q_unsafe"][i], v["T_mis"][i],
                           v["q_cert_applicable"][i], v["q_cert_confirmed"][i], v["q_cert_possible"][i],
                           v["q_dep_applicable"][i], v["q_dep_confirmed"][i], v["q_dep_possible"][i]])
    print("csv:", DUMP_CSV, flush=True)


# ----------------------------------------------------------------------------------------------- figure
def setup_text(kind):
    width_in, width_tex, body = S.setup(kind)
    txt = 11.0
    for key in ("font.size", "axes.labelsize", "axes.titlesize", "xtick.labelsize", "ytick.labelsize",
                "legend.fontsize", "legend.title_fontsize"):
        plt.rcParams[key] = txt
    sub = 10 if body <= 10 else 11
    plt.rcParams["text.latex.preamble"] = (r"\usepackage{amsmath}\usepackage{amssymb}\usepackage{bm}"
                                           rf"\DeclareMathSizes{{{txt:g}}}{{{txt:g}}}{{{sub}}}{{{sub}}}"
                                           rf"\DeclareMathSizes{{{TXT_SCEN:g}}}{{{TXT_SCEN:g}}}{{{sub}}}{{{sub}}}")
    return width_in, width_tex, txt


def header(label, col):
    a0, a, d, K, _ = S.SETTINGS[label]
    fixed = dict(alpha0=f"{a0:.2f}", alpha=f"{a:.2f}", delta=f"{d:g}" if d < 0.05 else f"{d:.2f}", K=f"{K}",
                 b=r"\kappa")   # the default budget normalizer b = kappa
    c = "phi" if col == "phi5" else col
    others = [k for k in ("alpha0", "alpha", "delta", "K", "b") if k != c]
    items = [rf"${PNAME[k]}{{=}}{fixed[k]}$" for k in others]
    first = rf"${PNAME[c]}$ varies"
    if c == "phi":
        first += rf", $L_\mu{{=}}{'0.0005' if col == 'phi' else '0.005'}$"
        return first + "\n" + ", ".join(items[:3]) + "\n" + ", ".join(items[3:])
    return first + "\n" + ", ".join(items[:2]) + "\n" + ", ".join(items[2:])


def recs_of(d, label, scen, col, m):
    return sorted([r for r in d["boxes"] if r["label"] == label and r["scen"] == scen and r["col"] == col
                   and r["method"] == m], key=lambda r: r["pos"])


def values_of(d, r, metric):
    v = d["configs"][r["configuration_id"]]["values"][metric]
    return SCALE[metric] * np.asarray([x for x in v if x is not None], dtype=float)


def offsets():
    n = len(VARIANTS)
    return [(k - (n - 1) / 2) * SLOT for k in range(n)]


def level_of(r):
    cf = r["configuration"]
    if str(cf["K"]) == "infinity":
        return cf["delta"] / (2 * cf["budget_normalizer"] * math.log(2) ** 2)
    return cf["delta"]


def draw_boxes(ax, d, label, scen, col, metric, info):
    for k, m in enumerate(VARIANTS):
        c, mk, _ = S.STYLE[m]
        rs = recs_of(d, label, scen, col, m)
        vals = [values_of(d, r, metric) for r in rs]
        pos = [r["pos"] + offsets()[k] for r in rs]
        ax.boxplot(vals, positions=pos, widths=WIDTH, whis=(5, 95), showfliers=False, patch_artist=True,
                   showmeans=True, manage_ticks=False, zorder=3 + k * 0.1,
                   boxprops=dict(facecolor=X.light(c), edgecolor=c, lw=LW_BOX),
                   whiskerprops=dict(color=c, lw=LW_BOX, solid_capstyle="butt"),
                   capprops=dict(color=c, lw=LW_BOX, solid_capstyle="butt"),
                   medianprops=dict(color="black", lw=LW_BOX, solid_capstyle="butt"),
                   meanprops=dict(marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.6, ls=""))
        for r, v, x in zip(rs, vals, pos):
            st = X.box_stats(v)
            if abs(st["mean"] - SCALE[metric] * r["metrics_mean"][metric]) > 1e-9 * max(1.0, abs(st["mean"])):
                raise ValueError("drawn mean differs from the export mean")
            info["boxes"].append(dict(col=col, x=r["x_label"], variant=S.LABEL[m], metric=metric, position=x,
                                      view_id=r["view_id"], is_setting=r["is_setting"], **st))


def draw_counts(ax, d, label, scen, col, metric, ylim, info, level_cut):
    drawn_level = False
    for k, m in enumerate(VARIANTS):
        c, mk, _ = S.STYLE[m]
        for r in recs_of(d, label, scen, col, m):
            q = r["counts"][metric]
            if not q["applicable_n"]:
                continue
            x = r["pos"] + offsets()[k]
            n = q["applicable_n"]
            y, ub = q["possible"], n * q["upper95"]
            ax.plot([x, x], [y, ub], color=c, lw=LW_BOX, solid_capstyle="butt", zorder=3)
            ax.plot([x - 0.06, x + 0.06], [ub, ub], color=c, lw=LW_BOX, solid_capstyle="butt", zorder=3)
            ax.plot([x], [y], ls="", marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.6, zorder=4)
            lev = level_of(r)
            info["counts"].append(dict(col=col, x=r["x_label"], variant=S.LABEL[m], metric=metric, position=x,
                                       count=y, confirmed=q["confirmed"], n=n, upper95=q["upper95"],
                                       upper95_count=ub, level=lev, level_count=n * lev,
                                       meets_level=q["upper95"] <= lev + 1e-12, view_id=r["view_id"]))
    # level per swept value (same for all variants of one value)
    for r in recs_of(d, label, scen, col, VARIANTS[0]):
        lev = N_MC * level_of(r)
        if ylim[0] <= lev <= min(ylim[1], level_cut):
            ax.plot([r["pos"] - 0.45, r["pos"] + 0.45], [lev, lev], color=LEVEL_COLOUR, ls="--", lw=S.LW_REF,
                    zorder=2)
            drawn_level = True
    return drawn_level


def nonzero_rows(d, label, scen):
    rows = []
    rs = [r for r in d["boxes"] if r["label"] == label and r["scen"] == scen]
    for m in ALL_ROWS:
        if m in BOXM:
            if m == "q_unsafe" and all(not np.any(values_of(d, r, m)) for r in rs):
                continue
            rows.append(m)
        else:
            if any(r["counts"][m]["possible"] or r["counts"][m]["confirmed"] for r in rs):
                rows.append(m)
    return rows


LEVEL_HEADROOM = 1.18     # a level is drawn when it lies below this multiple of the row's largest bound


def row_limits(d, label, scen, metric, h_pt=None):
    """y range and ticks of a row. Box rows: the variants' 5th-95th percentile whiskers plus MARGIN of their span on
    each side. Count rows: 0 to the largest upper bound, and every level below LEVEL_HEADROOM times that bound (the
    levels drawn). Values close to lines or edges are hard to read, so at the row's height h_pt
    (pt; None in the first layout pass) the range also keeps whisker caps, level lines and markers (counts at 0
    included) CLEAR_PT from the frame. Returns (ylim, ticks, level_cut): levels up to level_cut are drawn."""
    rs = [r for r in d["boxes"] if r["label"] == label and r["scen"] == scen]
    bl = d["baselines"][f"{S.SETTINGS[label][0]:g}"][scen]
    refs = [] if metric == "q_cert" else [bl[m]["q_dep" if metric == "q_dep" else metric] for m in X.REF_BASE]
    refs = [v for v in refs if v is not None]
    if metric in CNTM:
        qs = [r["counts"][metric] for r in rs if r["counts"][metric]["applicable_n"]]
        top = max(q["applicable_n"] * (q["upper95"] or 0) for q in qs)
        top = max(top, 3.0)
        levels = [N_MC * level_of(r) for r in rs]
        level_cut = LEVEL_HEADROOM * top
        shown = [v for v in levels if v <= level_cut]
        cmin = min(q["possible"] for q in qs)
        cmax = max(q["possible"] for q in qs)
        if h_pt is None:
            hi = max([top * LEVEL_HEADROOM] + shown)
            lim = (-0.07 * hi, hi)
        else:
            tops = [(top + MARGIN * top, 0.0), (top, X.line_pad(LW_BOX)), (cmax, X.marker_pad(MS_MEAN))]
            tops += [(v, X.line_pad(S.LW_REF)) for v in shown]
            lim = X.limits_with_refs(tops, [(cmin, X.marker_pad(MS_MEAN))], refs, h_pt, S.LW_REF)
        ticks = [t for t in MaxNLocator(nbins=2, integer=True).tick_values(0, top)
                 if 0 <= t <= lim[0] + 0.85 * (lim[1] - lim[0])]
        return lim, ticks, level_cut
    st = [X.box_stats(values_of(d, r, metric)) for r in rs if len(values_of(d, r, metric))]
    lo, hi = min(s["p5"] for s in st), max(s["p95"] for s in st)
    if hi - lo < 1e-12:
        hi = lo + 1.0
    pad = MARGIN * (hi - lo)
    if h_pt is None:
        lim = (lo - pad, hi + pad)
    else:
        mlo, mhi = min(s["mean"] for s in st), max(s["mean"] for s in st)
        lim = X.limits_with_refs([(hi + pad, 0.0), (hi, X.line_pad(LW_BOX)), (mhi, X.marker_pad(MS_MEAN))],
                                 [(lo - pad, 0.0), (lo, X.line_pad(LW_BOX)), (mlo, X.marker_pad(MS_MEAN))], refs,
                                 h_pt, S.LW_REF)
    span = lim[1] - lim[0]
    ticks = [float(np.round(t, 6)) for t in MaxNLocator(nbins=4, steps=[1, 2, 2.5, 5, 10]).tick_values(*lim)
             if lim[0] <= t <= lim[1] - 0.1 * span]      # no tick at the top edge (it meets the row above)
    return lim, ticks, None


def legend_handles(a0, any_level, any_count):
    def box_handle(m):
        c, mk, _ = S.STYLE[m]
        return (Patch(facecolor=X.light(c), edgecolor=c, lw=LW_BOX),
                Line2D([], [], color="black", lw=LW_BOX, solid_capstyle="butt"),
                Line2D([], [], ls="", marker=mk, ms=MS_MEAN, mfc=c, mec="black", mew=0.6))
    handles = [box_handle(m) for m in VARIANTS]
    labels = [S.LABEL[m] for m in VARIANTS]
    for m in X.REF_BASE:
        h = Line2D([], [], color=S.STYLE[m][0], ls=X.REF_LS[m], lw=S.LW_REF)
        if m == "QT":
            h.set_path_effects(X.QT_EFFECT)
        handles.append(h)
        labels.append(S.LABEL[m])
    if any_level:
        handles.append(Line2D([], [], color=LEVEL_COLOUR, ls="--", lw=S.LW_REF))
        labels.append("level")
    handles.append(Patch(facecolor=BAND_COLOUR, edgecolor="none"))
    labels.append("setting")
    return handles, labels


def build(d, kind, label, scen, out: Path):
    """Two layout passes: the first measures each row's height in pt, the second sets the y ranges from it."""
    heights = _build(d, kind, label, scen, out, None)
    return _build(d, kind, label, scen, out, heights)


def _build(d, kind, label, scen, out: Path, row_pt):
    width_in, width_tex, txt = setup_text(kind)
    a0 = S.SETTINGS[label][0]
    rows = nonzero_rows(d, label, scen)
    blk = [R_CNT if m in CNTM else R_BOX for m in rows]
    hr = blk + [R_GAP] + blk
    nb = len(rows)
    # height: 150 pt for headers, labels and legend plus 46 pt per unit of row height in each band
    height_pt = min(HEIGHT_CAP, 150.0 + 2 * 46.0 * sum(blk))
    ncat = {c: max(r["pos"] for r in d["boxes"] if r["label"] == label and r["col"] == c) + 1 for c in COLS}
    fig = plt.figure(figsize=(width_in, height_pt / 72.27))
    # one gridspec per band (different column widths)
    outer = fig.add_gridspec(len(hr), 1, height_ratios=hr, hspace=0.0)
    axes = []
    for bi, cols in enumerate(BANDS):
        sub = outer[bi * (nb + 1):bi * (nb + 1) + nb, 0].subgridspec(nb, len(cols), height_ratios=blk,
                                                                   width_ratios=[ncat[c] for c in cols],
                                                                   hspace=HSPACE, wspace=WSPACE)
        a = np.empty((nb, len(cols)), dtype=object)
        for i in range(nb):
            for j in range(len(cols)):
                a[i, j] = fig.add_subplot(sub[i, j], sharex=a[0, j] if i else None, sharey=a[i, 0] if j else None)
        axes.append(a)
    name = f"sweepbox_{scen}_{label}"
    info = dict(figure=name, kind=kind, setting=label, tuple=S.SETTINGS[label][:4], scenario=scen,
                variants=[S.LABEL[m] for m in VARIANTS], rows=rows, bands=[list(b) for b in BANDS],
                override=dict(text_pt=txt, box_lw=LW_BOX, mean_marker_pt=MS_MEAN,
                              legend="one row above the panels, outside the axes"),
                whiskers="5th and 95th percentiles", boxes=[], counts=[], ranges={}, reference_lines=[],
                annotations=[], baselines=d["baselines"][f"{a0:g}"][scen])
    bl = d["baselines"][f"{a0:g}"][scen]
    any_level = False
    right = [dict(), dict()]
    for i, metric in enumerate(rows):
        lim, ticks, level_cut = row_limits(d, label, scen, metric, None if row_pt is None else row_pt[i])
        info["ranges"][metric] = dict(ylim=list(lim), ticks=ticks)
        for bi, cols in enumerate(BANDS):
            for j, col in enumerate(cols):
                ax = axes[bi][i, j]
                n = ncat[col]
                rs0 = recs_of(d, label, scen, col, VARIANTS[0])
                sp = [r["pos"] for r in rs0 if r["is_setting"]]
                for p in sp:
                    ax.axvspan(p - 0.5, p + 0.5, color=BAND_COLOUR, lw=0, zorder=0)
                if metric != "q_cert":
                    for m in X.REF_BASE:
                        v = bl[m]["q_dep" if metric == "q_dep" else metric]
                        if v is None:
                            continue
                        if lim[0] <= v <= lim[1]:
                            ln = ax.axhline(v, color=S.STYLE[m][0], ls=X.REF_LS[m], lw=S.LW_REF, zorder=2)
                            if m == "QT":
                                ln.set_path_effects(X.QT_EFFECT)
                            if j == 0:
                                info["reference_lines"].append(dict(row=metric, band=bi, baseline=m, value=v))
                                right[bi].setdefault(i, []).append((m, v, "line"))
                        elif j == 0:
                            right[bi].setdefault(i, []).append(
                                (m, v, "up" if v > lim[1] else "down"))
                if metric in CNTM:
                    any_level |= draw_counts(ax, d, label, scen, col, metric, lim, info, level_cut)
                else:
                    draw_boxes(ax, d, label, scen, col, metric, info)
                ax.set_xlim(-0.5, n - 0.5)
                labels_x = [r["x_label"] for r in rs0]
                ax.set_xticks(range(n), X.tick_labels(labels_x, sp[0] if sp else -1))
                ax.tick_params(axis="x", length=0, pad=3)
                ax.grid(False)
                ax.grid(True, which="major", axis="y")
                for p in range(1, n):      # thin separators between swept values
                    ax.axvline(p - 0.5, color="#DDDDDD", lw=0.6, zorder=0.5)
                ax.set_xlabel(XLAB[col])
                ylab = YLAB[metric]
                if kind == "report" and metric == "m_c":    # 11 pt sub/superscripts: two lines fit the row
                    ylab = "$m_{c,k}/$\n$c_{\\max}$"
                ax.set_ylabel(ylab, linespacing=1.0)
                if j > 0:
                    ax.yaxis.label.set_visible(False)
                    ax.tick_params(labelleft=False)
                if i != nb - 1:
                    ax.xaxis.label.set_visible(False)
                    ax.tick_params(labelbottom=False)
                if i == 0:
                    ax.set_title(header(label, col), fontsize=txt, pad=5, linespacing=1.15)
            axes[bi][i, 0].set_ylim(*lim)
            axes[bi][i, 0].set_yticks(ticks, [f"{t:g}" for t in ticks])
    for a in axes:
        fig.align_ylabels(a[:, 0])
    rlabels = [dict(), dict()]
    for bi in range(2):
        for i, items in right[bi].items():
            metric = rows[i]
            srt = sorted(items, key=lambda z: -z[1])
            rlabels[bi][i] = [(fig.text(0, 0, fmt_label(metric, m, v, where, lines=1), ha="left", va="center",
                                        fontsize=txt, linespacing=1.0), where, v) for m, v, where in srt]
            for (t, where, v), (m, _, _) in zip(rlabels[bi][i], srt):
                info["annotations"].append(dict(row=metric, band=bi, baseline=m, value=v, kind=where,
                                                text=t.get_text()))
    handles, labels = legend_handles(a0, any_level, any(m in CNTM for m in rows))
    ncol = len(labels)
    leg = fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=ncol,
                     frameon=False, fontsize=txt, handlelength=1.6, handleheight=0.9, handletextpad=0.35,
                     columnspacing=0.9, borderaxespad=0.0, borderpad=0.1,
                     handler_map={tuple: HandlerTuple(ndivide=1, pad=0.0)})
    layout(fig, axes, leg, rlabels, height_pt, outer, hr, nb)
    if row_pt is None:
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        heights = [axes[0][i, 0].get_window_extent(r).height * 72.0 / fig.dpi for i in range(nb)]
        plt.close(fig)
        return heights
    shifts = []
    for a in axes:
        shifts += X.separate_tick_labels(fig, a[-1, :])
    bad = S.text_overlaps(fig)
    if bad:
        raise ValueError(f"{name}: overlapping texts {bad[:6]}")
    pdf = fs.save(fig, out / name, png=False, strict=True)
    plt.close(fig)
    w, h = S.check_width(pdf, width_tex)
    if h > HEIGHT_CAP + 0.5:
        raise ValueError(f"height {h} pt above {HEIGHT_CAP}")
    info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=h, text_pt=txt,
                tick_label_shifts_pt=shifts,
                levels_exceeded=[c for c in info["counts"] if not c["meets_level"]])
    print(name, kind, f"{w} x {h} TeX pt", "rows", rows, "levels exceeded:", len(info["levels_exceeded"]),
          flush=True)
    return info


def fmt_label(metric, m, v, where, lines=2):
    val = {"m_c": f"{v:.2f}", "q_CBF": f"{v:.1f}", "q_unsafe": f"{v:.1f}", "T_mis": f"{v:.1f}",
           "q_dep": f"{int(round(v))}"}[metric]
    arrow = {"line": "", "up": r" $\blacktriangle$", "down": r" $\blacktriangledown$"}[where]
    return f"{S.LABEL[m]}{arrow}" + ("\n" if lines == 2 else " ") + val


def layout(fig, axes, leg, rlabels, height_pt, outer, hr, nb, gap_pt=4.0, band_gap_pt=7.0, iters=30):
    """Subplot box such that the y labels touch the left edge, the bottom labels the bottom edge, the baseline labels
    the right edge, and the column headers of band 1 sit gap_pt below the legend; the gap row between the bands is
    sized so that band 2's headers sit band_gap_pt below band 1's x labels."""
    fig.subplots_adjust(left=0.10, right=0.93, bottom=0.05, top=0.90)
    W, H = fig.get_size_inches()
    hr = list(hr)
    for _ in range(iters):
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        Hp = H * fig.dpi
        for bi in range(2):
            X.place_right_labels(fig, axes[bi], rlabels[bi])
        fig.canvas.draw()
        bb = fig.get_tightbbox(r)
        sp = fig.subplotpars
        dl, dr, db = -bb.x0 / W, (W - bb.x1) / W, -bb.y0 / H
        title_top = max(ax.title.get_window_extent(r).y1 for ax in axes[0][0]) / Hp
        leg_bottom = leg.get_window_extent(r).y0 / Hp
        dt = (leg_bottom - gap_pt / 72.0 / H) - title_top
        xl_bottom = min(ax.xaxis.get_tightbbox(r).y0 for ax in axes[0][-1])
        t2_top = max(ax.title.get_window_extent(r).y1 for ax in axes[1][0])
        need_px = band_gap_pt * fig.dpi / 72.0 - (xl_bottom - t2_top)
        unit_px = axes[0][0, 0].get_window_extent(r).height / hr[0]
        dg = need_px / unit_px
        if max(abs(dl), abs(dr), abs(db), abs(dt)) < 1e-4 and abs(need_px) < 0.3:
            break
        hr[nb] = max(0.05, hr[nb] + dg)
        outer.set_height_ratios(hr)
        fig.subplots_adjust(left=sp.left + dl, right=sp.right + dr, bottom=sp.bottom + db, top=sp.top + dt)
    else:
        raise ValueError("layout did not converge")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(S.KINDS), default="paper")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--reload", action="store_true")
    ap.add_argument("--load-only", action="store_true")
    ap.add_argument("--labels", nargs="+", default=None)
    ap.add_argument("--scens", nargs="+", default=list(SCENS))
    a = ap.parse_args()
    labels = a.labels or list(LABELS)
    d = load(labels, a.reload)
    if a.load_only:
        return
    out = a.out or HERE / ("out_paper" if a.kind == "paper" else "out_report")
    out.mkdir(parents=True, exist_ok=True)
    infos = []
    for label in labels:
        for scen in a.scens:
            infos.append(build(d, a.kind, label, scen, out))
    tag = "_".join(labels) + ("" if list(a.scens) == list(SCENS) else "_" + "_".join(a.scens))
    (out / f"sweeps_manifest_{a.kind}_box_ext_{tag}.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(B.CELLS): S.sha_file(B.CELLS),
                                                                 S.recorded_path(CACHE): S.sha_file(CACHE)},
        checks=d["checks"], figures=infos), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
