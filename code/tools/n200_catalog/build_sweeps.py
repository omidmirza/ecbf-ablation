#!/usr/bin/env python3
"""Sweep plots at a fixed setting.

sweepplot_<scen>_<label>.pdf (label Default, New1..New5 = T1..T5): columns = the swept parameter alpha (at the
setting's alpha0, delta, K), delta, K (categorical 10, 50, 100, 200, inf), b (budget normalizer), phi (bound factor;
L_mu = 0.0005 solid lines with markers, L_mu = 0.005 dashed with open markers; phi < 1 shaded: Assumption 2 may fail
with the reduced bound); rows = mean robustness margin m_c/c_max, CBF violations q+_CBF (%), safe-set violations
q+_unsafe (%), false certifications q_cert (count of N_MC = 200 with a confirmed or possible event), mission time
T_mis (s). Every other parameter stays at the setting. One line per variant (V-U, V-Delta, V-L, V-C) in
its colour and marker; the setting's own value: vertical dotted line; baselines c=0, c_max, ACI, QT at the setting's
alpha0 (matched_baselines_200, L_mu = 0.0005): horizontal reference lines (1.8 pt) where the metric applies
(not q_cert). Unfiltered control omitted. q rows: symlog, ticks 0, 0.1, 1, 10, 100.
sweepplot_T3_compact.pdf (--compact; not printed, the arXiv paper uses sweepbox_T3_compact.pdf of
build_sweeps_box.py): rows Corridor m_c/c_max, Corridor q+_CBF, Maze m_c/c_max, Maze q+_CBF at T3; columns alpha,
delta, K, b.

Data: sweep_cells.json (extract_sweep_cells.py) of N200_numbers/{final,new_anchors}/common/metrics.json: alpha,
delta, K from final main_grid_200; b from budget_normalizer_200 and phi x L_mu from
stress_drift_200 (final for the Default, new_anchors for T1-T5); baselines from final matched_baselines_200.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs

from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
CELLS = HERE / "sweep_cells.json"
KVALS = ["10", "50", "100", "200", "infinity"]
ROWS = ("m_c", "q_CBF", "q_unsafe", "q_cert", "T_mis")
YLAB = {"m_c": r"$\widehat m_c/c_{\max}$", "q_CBF": r"$\widehat q^+_{\mathrm{CBF}}$ (\%)",
        "q_unsafe": r"$\widehat q^+_{\mathrm{unsafe}}$ (\%)", "q_cert": r"$\widehat q_{\mathrm{cert}}$ (count)",
        "T_mis": r"$\widehat T_{\mathrm{mis}}$ (s)"}
COLS = ("alpha", "delta", "K", "b", "phi")
XLAB = {"alpha": r"$\alpha$", "delta": r"$\delta$", "K": r"$K$", "b": r"$b$", "phi": r"$\phi$"}


def _eq(have, want, key):
    if have is None:
        return False
    if key == "K" or isinstance(want, str):
        return str(have) == str(want)
    return abs(float(have) - float(want)) < 1e-12


class Data:
    def __init__(self):
        self.cells = json.loads(CELLS.read_text())["cells"]

    def _pick(self, fam, source, scen, **want):
        out = []
        for c in self.cells:
            if c["family"] != fam or (source and c["source"] != source):
                continue
            cf = c["configuration"]
            if cf.get("benchmark") != scen:
                continue
            if all(_eq(cf.get(k), v, k) for k, v in want.items()):
                out.append(c)
        return out

    def one(self, fam, source, scen, **want):
        r = self._pick(fam, source, scen, **want)
        if len(r) != 1:
            raise ValueError(f"{fam} {source} {scen} {want}: {len(r)} cells")
        if r[0]["n"] != 200:
            raise ValueError("n != 200")
        return r[0]

    def sweep(self, label, scen, col, method):
        """[(x, cell)] of the sweep `col` around setting `label` for `method` (and L_mu for phi: dict)."""
        a0, a, d, K, _ = S.SETTINGS[label]
        base = dict(alpha0=a0, alpha=a, delta=d, K=str(K), method=method, L_mu=0.0005)
        src = "final" if label == "Default" else "new_anchors"
        if col == "alpha":
            xs = [x for x in (0.05, 0.10, 0.15, 0.20, 0.25) if x < a0 - 1e-12]
            return [(x, self.one("main_grid_200", "final", scen, **dict(base, alpha=x))) for x in xs]
        if col == "delta":
            return [(x, self.one("main_grid_200", "final", scen, **dict(base, delta=x)))
                    for x in (0.025, 0.05, 0.10, 0.15, 0.20)]
        if col == "K":
            return [(i, self.one("main_grid_200", "final", scen, **dict(base, K=k))) for i, k in enumerate(KVALS)]
        if col == "b":
            return [(x, self.one("budget_normalizer_200", src, scen, **dict(base, budget_normalizer=x)))
                    for x in (2.11, 2.5, 3.0, 4.0)]
        if col == "phi":
            return {lm: [(x, self.one("stress_drift_200", src, scen, **dict(base, L_mu=lm, bound_scale=x,
                                                                         budget_normalizer=2.11)))
                         for x in (0.0, 0.25, 0.5, 0.75, 1.0)] for lm in (0.0005, 0.005)}
        raise ValueError(col)

    def baseline(self, scen, a0, method, L_mu=0.0005):
        r = [c for c in self._pick("matched_baselines_200", "final", scen, method=method, alpha0=a0, L_mu=L_mu)]
        if len(r) != 1:
            raise ValueError(f"baseline {scen} {a0} {method}: {len(r)}")
        return r[0]


def value(cell, row):
    if row == "q_cert":
        e = cell["events"].get("qcert")
        return None if not e or e.get("applicable_n") in (None, 0) else float(e["possible"])
    st = cell["stats"].get(row)
    if not st or not st.get("applicable") or st.get("mean") is None:
        return None
    v = float(st["mean"])
    return 100 * v if row in ("q_CBF", "q_unsafe") else v


def setting_x(label, col):
    a0, a, d, K, _ = S.SETTINGS[label]
    return {"alpha": a, "delta": d, "K": KVALS.index(str(K)), "b": 2.11, "phi": 1.0}[col]


def style_axis(ax, col, row, txt, label):
    if col == "alpha":
        xs = [x for x in (0.05, 0.10, 0.15, 0.20, 0.25) if x < S.SETTINGS[label][0] - 1e-12]
        ax.set_xlim(xs[0] - 0.025, xs[-1] + 0.025)
        tk = [0.05, 0.15, 0.25] if len(xs) == 5 else xs
        ax.set_xticks(tk, [f"{t:g}" for t in tk])
    elif col == "delta":
        ax.set_xscale("log")
        ax.set_xlim(0.025 / 1.3, 0.2 * 1.3)
        ax.set_xticks([0.025, 0.05, 0.1, 0.15, 0.2], ["0.025", "", "0.1", "", "0.2"])
        ax.xaxis.set_minor_locator(plt.NullLocator())
    elif col == "K":
        ax.set_xlim(-0.45, 4.45)
        ax.set_xticks(range(5), ["10", "50", "100", "", r"$\infty$"])   # 200 unlabelled: 100 and 200 collide
    elif col == "b":
        ax.set_xlim(2.11 - 0.2, 4.0 + 0.2)
        ax.set_xticks([2.11, 3, 4], ["2.11", "3", "4"])
    elif col == "phi":
        ax.set_xlim(-0.1, 1.1)
        ax.set_xticks([0, 0.5, 1], ["0", "0.5", "1"])
    if row in ("q_CBF", "q_unsafe"):
        ax.set_yscale("symlog", linthresh=0.1, linscale=0.6)
        ax.set_yticks([0, 0.1, 1, 10, 100], ["0", "0.1", "1", "10", "100"])
        ax.yaxis.set_minor_locator(plt.NullLocator())
    ax.grid(True, which="major", axis="both")


# V-U and V-Delta often coincide (same dashed style): V-U is drawn above V-Delta with its dashes shifted by one
# dash length, so that where the two curves coincide they alternate light and dark blue.
ZORDER = {"V-D": 3.0, "V-U": 3.2, "V-Cd": 3.4, "V-C": 3.6}


def line_style(m, dashed=False):
    if dashed:                              # phi column, L_mu = 0.005
        return {"V-U": (0, (3, 1.5)), "V-D": (2.25, (3, 1.5))}.get(m, (0, (3, 1.5)))
    return {"V-U": (0, (3.5, 3.5)), "V-D": (3.5, (3.5, 3.5))}.get(m, S.STYLE[m][2])


def draw_panel(ax, data, label, scen, col, row, txt, info):
    a0 = S.SETTINGS[label][0]
    ax.axvline(setting_x(label, col), color=S.REF_COLOUR, ls=":", lw=S.LW_REF, zorder=1, label="_nolegend_")
    if col == "phi":
        ax.axvspan(-0.1, 1.0, color="#E8E8E8", lw=0, zorder=0)
    if row != "q_cert":
        for m in S.BASELINES:
            v = value(data.baseline(scen, a0, m), row)
            if v is None:
                continue
            c, _, ls = S.STYLE[m]
            ax.axhline(v, color=c, ls=ls, lw=S.LW_REF, zorder=2, label="_nolegend_")
            info.setdefault("baselines", {}).setdefault(row, {})[m] = v
    for m in S.CERTIFIED:
        c, mk, ls = S.STYLE[m]
        sw = data.sweep(label, scen, col, m)
        series = sw.items() if isinstance(sw, dict) else [(0.0005, sw)]
        for lm, pts in series:
            xs = [x for x, cell in pts if value(cell, row) is not None]
            ys = [value(cell, row) for x, cell in pts if value(cell, row) is not None]
            dashed = lm == 0.005
            ax.plot(xs, ys, color=c, ls=line_style(m, dashed), lw=S.LW_DATA, marker=mk, ms=S.MS,
                    mfc="white" if dashed else c, mec=c if dashed else S.edge(c), mew=0.8 if dashed else S.MEW,
                    zorder=ZORDER[m], label="_nolegend_")
            info.setdefault("series", []).append(dict(method=m, col=col, row=row, L_mu=lm, x=xs, y=ys,
                                                      view_ids=[cell["view_id"] for _, cell in pts]))


def legend_handles(kind_phi=False):
    hv = [Line2D([], [], color=S.STYLE[m][0], ls=line_style(m), lw=S.LW_DATA, marker=S.STYLE[m][1], ms=S.MS,
                 mfc=S.STYLE[m][0], mec=S.edge(S.STYLE[m][0]), mew=S.MEW, label=S.LABEL[m]) for m in S.CERTIFIED]
    hb = [Line2D([], [], color=S.STYLE[m][0], ls=S.STYLE[m][2], lw=S.LW_REF, label=S.LABEL[m]) for m in S.BASELINES]
    return hv, hb


def set_limits(axes_row, row, data_max, compact=False):
    if row == "m_c":
        for ax in axes_row:
            ax.set_ylim(-0.06, 1.1)
            ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
    elif row in ("q_CBF", "q_unsafe"):
        for ax in axes_row:
            if compact:   # compact figure: every value of these rows is below 50 %
                if data_max >= 50:
                    raise ValueError("compact q row above 50 %")
                ax.set_yscale("symlog", linthresh=0.1, linscale=1.0)
                ax.set_yticks([0, 0.1, 1, 10], ["0", "0.1", "1", "10"])
                ax.yaxis.set_minor_locator(plt.NullLocator())
                ax.set_ylim(-0.02, 60)
            else:
                ax.set_ylim(-0.035, 160)
    elif row == "q_cert":
        top = max(4.0, data_max)
        top = 3.0 * top
        step = [s for s in (1, 2, 5, 10, 20, 25, 50, 100) if top / s <= 4][0]
        for ax in axes_row:
            ax.set_ylim(-0.06 * top, top)
            ax.set_yticks(np.arange(0, top + 1e-9, step))


def build(kind, out, labels, scens, compact):
    width_in, width_tex, txt = S.setup(kind)
    data = Data()
    manifest = []
    jobs = [(s, l) for l in labels for s in scens]
    if compact:
        jobs = [("compact", "New3")]
    for scen, label in jobs:
        if scen == "compact":
            rows = [("corridor", "m_c"), ("corridor", "q_CBF"), ("maze", "m_c"), ("maze", "q_CBF")]
            cols = ("alpha", "delta", "K", "b")
            height_pt = 258.0
            name = "sweepplot_T3_compact"
        else:
            rows = [(scen, r) for r in ROWS]
            cols = COLS
            height_pt = 425.0 if kind == "paper" else 465.0
            name = f"sweepplot_{scen}_{label}"
        fig, axes = plt.subplots(len(rows), len(cols), figsize=(width_in, height_pt / 72.27), sharex="col",
                                 sharey="row", squeeze=False)
        fig.subplots_adjust(wspace=0.12, hspace=0.10)
        info = dict(figure=name, setting=label, tuple=S.SETTINGS[label][:4], paired=S.SETTINGS[label][4],
                    rows=[f"{s}:{r}" for s, r in rows], cols=list(cols))
        for i, (sc, row) in enumerate(rows):
            dmax = 0.0
            for j, col in enumerate(cols):
                ax = axes[i, j]
                draw_panel(ax, data, label, sc, col, row, txt, info)
                style_axis(ax, col, row, txt, label)
                for ln in ax.get_lines():
                    y = np.asarray(ln.get_ydata(), dtype=float)
                    if len(y) and ln.get_transform() == ax.transData:
                        dmax = max(dmax, float(np.nanmax(y)))
                ax.set_xlabel(XLAB[col])
                ylab = YLAB[row] if scen != "compact" else (S.SCEN[sc] + "\n" + YLAB[row])
                ax.set_ylabel(ylab)
                if j > 0:
                    ax.yaxis.label.set_visible(False)
                if i < len(rows) - 1:
                    ax.xaxis.label.set_visible(False)
                ax.tick_params(labelsize=txt)
            if row in ("q_CBF", "q_unsafe"):   # baselines count too
                dmax = max([dmax] + list(info.get("baselines", {}).get(row, {}).values()))
            set_limits(axes[i], row, dmax, compact=scen == "compact")
        fig.align_ylabels(axes[:, 0])
        S.fill_box(fig, left=0.09, right=0.995, bottom=0.07, top=0.995)
        hv, hb = legend_handles()
        common = dict(handlelength=1.3, handletextpad=0.3, columnspacing=0.5, borderaxespad=0.2, borderpad=0.25,
                      labelspacing=0.2, framealpha=1.0)
        ref_lines = [ln for ax in axes.flat for ln in ax.get_lines() if ln.get_xydata().shape[0] == 2
                     and ln.get_linestyle() == ":" and ln.get_color() == S.REF_COLOUR]
        order = [axes[i, j] for i in range(len(rows)) for j in range(len(cols))]
        if scen != "compact":   # q_cert row first (no baseline lines), then the others
            qi = [r for _, r in rows].index("q_cert")
            order = [axes[qi, j] for j in range(len(cols))] + [a for a in order if a not in axes[qi]]
        ax_v, leg_v = S.place_legend(order, hv, txt, ncols=(2, 1), **common)
        if leg_v is None:
            ax_v, leg_v = S.place_legend(order, hv, txt, ncols=(2, 1), ignore=ref_lines, **common)
        if leg_v is None:
            raise ValueError(f"{name}: no free place for the variant legend")
        ax_v.add_artist(leg_v)
        ax_v.legend_ = None
        order_b = [a for a in order if a is not ax_v]
        # baselines legend: a panel without the q_cert restriction is fine; it must not cover any line
        ax_b, leg_b = S.place_legend(order_b, hb, txt, ncols=(2, 1), **common)
        if leg_b is None:
            ax_b, leg_b = S.place_legend(order_b, hb, txt, ncols=(2, 1), ignore=ref_lines, **common)
        if leg_b is None:
            raise ValueError(f"{name}: no free place for the baseline legend")
        bad = S.text_overlaps(fig)
        if bad:
            raise ValueError(f"{name}: overlapping texts {bad[:6]}")
        pdf = fs.save(fig, out / name, png=False, strict=True)
        plt.close(fig)
        w, h = S.check_width(pdf, width_tex)
        rowsn = [r for _, r in rows]
        info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=h, text_pt=txt,
                    legend_variants_panel=[rowsn[k // len(cols)] + "/" + cols[k % len(cols)]
                                           for k, a in enumerate(axes.flat) if a is ax_v],
                    legend_baselines_panel=[rowsn[k // len(cols)] + "/" + cols[k % len(cols)]
                                            for k, a in enumerate(axes.flat) if a is ax_b])
        cap = 430 if kind == "paper" else 470
        if compact and h > 260 + 0.5 or (not compact and h > cap + 0.5):
            raise ValueError(f"{name}: height {h} pt above the cap")
        manifest.append(info)
        print(name, kind, f"{w} x {h} TeX pt", "legends:", info["legend_variants_panel"], info["legend_baselines_panel"],
              flush=True)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(S.KINDS), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--labels", nargs="+", default=list(S.SETTINGS))
    ap.add_argument("--scens", nargs="+", default=list(S.SCEN))
    ap.add_argument("--compact", action="store_true")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = build(a.kind, a.out, a.labels, a.scens, a.compact)
    tag = "compact" if a.compact else "_".join(a.labels)
    (a.out / f"sweeps_manifest_{a.kind}_{tag}.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(CELLS): S.sha_file(CELLS)},
        figures=man), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
