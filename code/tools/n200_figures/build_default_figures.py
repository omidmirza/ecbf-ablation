#!/usr/bin/env python3
"""N200 figures of Appendix G (long version) and report Section 5 other than heatmaps and trajectory maps.

Figures are saved through the figure helper ieee_figure.py (strict save: vector PDF, zero margin, text >= body
size, line widths). Every figure is drawn at its target width: paper \\textwidth 505.89 pt (Appendix G is one column)
at the 10 pt body size, report \\textwidth 505.89 pt (geometry letterpaper, margin 0.75in) at the 11 pt body size.
V-L is drawn in #009E73; the T3 V-L theorem-check markers use it. All text is set at the body size (compact rule).

  default_margin_distributions.pdf          m_{c,k}/c_max per trajectory, default setting, Corridor | Maze
  distributions_default_<scen>_row2.pdf     q^+_{CBF,k} | q^+_{unsafe,k} per trajectory (%)
  distributions_default_<scen>_row3.pdf     p^max_k (alpha0 dotted) | h^min_k (m) per trajectory
  decisions_default_<scen>_K50.pdf          mean decision fractions q^+_nocert, q^+_nodep, q^+_use (four variants)
  theorem_checks_<scen>_New3.pdf            T3 (0.30, 0.20, 0.20, 50) V-L: events vs guarantee level at
                                            K in {50, inf} x L_mu in {0.0005, 0.005}, and risk-increase bound checks

Inputs (read only): default_assessments.json (load_default_assessments.py: per-trajectory saved assessments, hash-
checked against N200_INDEX.sqlite, means equal to the export) and N200_numbers/new_anchors/common/metrics.json
(detailed_200 at role 3 = T3; folder set by ECBF_N200_NUMBERS). Unfiltered control is omitted from plots (kept in tables).
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import sys
from pathlib import Path

if sys.platform == "darwin":   # MacTeX and Homebrew (latex, pdfinfo) when not on the PATH
    os.environ["PATH"] = "/Library/TeX/texbin:/opt/homebrew/bin:/usr/local/bin:" + os.environ.get("PATH", "")
sys.dont_write_bytecode = True
REPO = Path(__file__).resolve().parents[3]   # repository root (code/tools/n200_figures)
# The figure helper ieee_figure.py and the N200 summary exports are not part of this repository.
FIGURE_STYLE = Path(os.environ.get("ECBF_FIGURE_STYLE", REPO / "figure_style"))
NUM = Path(os.environ.get("ECBF_N200_NUMBERS", REPO / "N200_numbers"))
sys.path.insert(0, str(FIGURE_STYLE))


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


import ieee_figure as fs  # noqa: E402

import numpy as np  # noqa: E402
import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = Path(__file__).resolve().parent
DEFAULT_CACHE = HERE / "default_assessments.json"
NEW_ANCHORS = NUM / "new_anchors/common/metrics.json"
TEX_PT_PER_BP = 72.27 / 72.0
KINDS = {"paper": (505.89, 10.0), "report": (505.89, 11.0)}   # report geometry margin 0.75in

# the paper's method colours and markers
METHODS = ("Zero margin", "Conservative", "ACI", "QT", "V-U", "V-D", "V-Cd", "V-C")
LABEL = {"Zero margin": r"$c{=}0$", "Conservative": r"$c_{\max}$", "ACI": "ACI", "QT": "QT", "V-U": "V-U",
         "V-D": r"V-$\Delta$", "V-Cd": "V-L", "V-C": "V-C"}
STYLE = {"Zero margin": ("#7F7F7F", "v"), "Conservative": ("#B0B0B0", "s"), "ACI": ("#CC79A7", "D"),
         "QT": ("#F0E442", "P"), "V-U": ("#56B4E9", "o"), "V-D": ("#0072B2", "h"), "V-Cd": ("#009E73", "p"),
         "V-C": ("#E69F00", "^")}
CERTIFIED = ("V-U", "V-D", "V-Cd", "V-C")
VL = STYLE["V-Cd"][0]                       # T3 is a V-L setting: theorem-check markers in the V-L colour
SCEN = {"corridor": "Corridor", "maze": "Maze"}
MS, MEW, LW_BOX, LW_REF = 3.0, 0.5, 1.2, 1.8


def sha_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


class Builder:
    def __init__(self, kind: str, out: Path):
        self.kind = kind
        width_texpt, self.txt = KINDS[kind]
        self.width_in = width_texpt / 72.27              # figure width = template width in TeX pt
        self.content_bp = width_texpt / TEX_PT_PER_BP
        if not any(abs(self.width_in - w) < 1e-6 for w in fs.TEMPLATE_WIDTHS_IN):
            fs.TEMPLATE_WIDTHS_IN = tuple(fs.TEMPLATE_WIDTHS_IN) + (self.width_in,)
        if self.txt > 10:
            fs.MIN_TEXT_PT = 10.95
        self.out = out
        self.manifest = []

    # ------------------------------------------------------------------ helpers
    def _figure(self, height_pt, ncols, sharey=False, wspace=None, stack_height_pt=None):
        """Two panels side by side in the paper; stacked (shared x axis) in the report, where the 11 pt method
        and setting labels do not fit under half-width panels."""
        self.stacked = ncols == 2 and self.kind == "report"
        if self.stacked:
            fig, axes = fs.figure(width=self.width_in, height=(stack_height_pt or 1.75 * height_pt) / 72.0, nrows=2,
                                  ncols=1, sharex=True, layout="constrained")
        else:
            fig, axes = fs.figure(width=self.width_in, height=height_pt / 72.0, ncols=ncols, sharey=sharey,
                                  layout="constrained")
        t = self.txt
        plt.rcParams["text.latex.preamble"] += rf"\DeclareMathSizes{{{t:g}}}{{{t:g}}}{{{t:g}}}{{{t:g}}}"
        for ax in np.atleast_1d(axes):
            ax.tick_params(labelsize=t)
            ax.xaxis.label.set_size(t)
            ax.yaxis.label.set_size(t)
            ax.title.set_size(t)
        fig.get_layout_engine().set(w_pad=0.5 / 72, h_pad=0.5 / 72, hspace=0.04 if self.stacked else 0.02)
        if wspace is not None:
            fig.get_layout_engine().set(wspace=wspace)
        return fig, axes

    def _save(self, fig, name, info):
        # the content must span the template width: the figure is created 1/1.00375 wider than the content box
        # would need, and the tight bbox crops to the drawn content; check the result
        stem = self.out / name
        pdf = fs.save(fig, stem, png=False, strict=True)
        plt.close(fig)
        import shutil
        import subprocess
        size = subprocess.run([shutil.which("pdfinfo") or "pdfinfo", str(pdf)], capture_output=True, text=True).stdout
        w, h = [float(x) for x in size.split("Page size:")[1].split("pts")[0].split(" x ")]
        w_tex = w * TEX_PT_PER_BP
        target = KINDS[self.kind][0]
        if abs(w_tex - target) > 1.0:
            raise ValueError(f"{pdf.name}: width {w_tex:.2f} TeX pt, target {target}")
        info.update(pdf=recorded_path(pdf), sha256=sha_file(pdf), width_texpt=round(w_tex, 2), height_texpt=round(h * TEX_PT_PER_BP, 2),
                    text_pt=self.txt)
        self.manifest.append(info)
        print(pdf.name, f"{w_tex:.2f} x {h * TEX_PT_PER_BP:.2f} TeX pt", flush=True)

    def _fill_width(self, fig):
        """Freeze the constrained layout, then widen the outer axes so the drawn content spans the full width."""
        fig.canvas.draw()
        fig.set_layout_engine("none")
        W = fig.get_size_inches()[0]
        axes = sorted([a for a in fig.axes if a.get_visible()], key=lambda a: a.get_position().x0)
        for _ in range(3):
            fig.canvas.draw()
            bb = fig.get_tightbbox(fig.canvas.get_renderer())
            gl, gr = bb.x0 / W, (W - bb.x1) / W
            if abs(gl) < 1e-4 and abs(gr) < 1e-4:
                break
            if len(axes) == 2 and abs(axes[0].get_position().x0 - axes[1].get_position().x0) < 1e-6:   # stacked
                for a in axes:
                    q = a.get_position()
                    a.set_position([q.x0 - gl, q.y0, q.width + gl + gr, q.height])
                continue
            first, last = axes[0], axes[-1]
            p0 = first.get_position()
            first.set_position([p0.x0 - gl, p0.y0, p0.width + gl, p0.height])
            p1 = last.get_position()
            last.set_position([p1.x0, p1.y0, p1.width + gr, p1.height])
            if len(axes) == 2 and first is not last:   # keep the two panels equally wide
                a, b = first.get_position(), last.get_position()
                gap = b.x0 - a.x1
                w = (b.x1 - a.x0 - gap) / 2
                first.set_position([a.x0, a.y0, w, a.height])
                last.set_position([a.x0 + w + gap, b.y0, w, b.height])

    def _hide_top_xlabel(self, axes):
        if self.stacked:
            axes[0].xaxis.label.set_visible(False)

    def _strip_box(self, ax, data, methods, scale=1.0, na_y=None):
        """Jittered per-trajectory dots (left) and a box plot (right) per method; N-A where undefined."""
        rng = np.random.default_rng(0)          # jitter positions only
        for i, m in enumerate(methods):
            vals = np.array([v * scale for v in data[m] if v is not None], dtype=float)
            colour, marker = STYLE[m]
            if len(vals) == 0:
                ax.text(i, na_y, "N-A", ha="center", va="center", fontsize=self.txt)
                continue
            x = i - 0.17 + rng.uniform(-0.11, 0.11, len(vals))
            ax.plot(x, vals, linestyle="none", marker=marker, markersize=MS, markerfacecolor=colour,
                    markeredgecolor="black", markeredgewidth=MEW, alpha=0.75, label="_nolegend_")
            bp = ax.boxplot(vals, positions=[i + 0.22], widths=0.22, whis=1.5, showfliers=False, patch_artist=True,
                            manage_ticks=False)
            for part in ("whiskers", "caps", "medians"):
                for ln in bp[part]:
                    ln.set_linewidth(LW_BOX)
                    ln.set_color("black")
            for b in bp["boxes"]:
                b.set_facecolor("white")
                b.set_edgecolor("black")
                b.set_linewidth(LW_BOX)
        ax.set_xticks(range(len(methods)), [LABEL[m] for m in methods])
        ax.set_xlim(-0.6, len(methods) - 0.4)
        ax.tick_params(axis="x", length=0, pad=3)
        ax.grid(axis="x", visible=False)

    # ------------------------------------------------------------------ default figures
    def default_figures(self, cache: dict):
        cells = cache["cells"]

        def col(scen, key):
            return {m: [r[key] for r in cells[f"{scen}|{m}"]["values"]] for m in METHODS}

        # (1) normalized mean robustness margin, both scenarios
        fig, axes = self._figure(150, 2, sharey=True)
        for ax, scen in zip(axes, SCEN):
            self._strip_box(ax, col(scen, "m_c"), METHODS)
            ax.set_title(SCEN[scen], pad=3)
            ax.set_xlabel("Method")
            ax.set_ylim(-0.06, 1.08)
            ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
        axes[0].set_ylabel(r"$m_{c,k}/c_{\max}$")
        axes[1].set_ylabel(r"$m_{c,k}/c_{\max}$")
        if not self.stacked:
            axes[1].yaxis.label.set_visible(False)
        self._hide_top_xlabel(axes)
        self._fill_width(fig)
        self._save(fig, "default_margin_distributions",
                   dict(figure="default_margin_distributions", setting="(0.15,0.10,0.05,50)", metric="m_c",
                        methods=list(METHODS), source="default_assessments.json"))

        for scen in SCEN:
            # (2) CBF-violation and safe-set-violation fractions per trajectory (%)
            fig, axes = self._figure(150, 2)
            for ax, key, sym in ((axes[0], "q_CBF", r"$q^{+}_{\mathrm{CBF},k}$ (\%)"),
                                 (axes[1], "q_unsafe", r"$q^{+}_{\mathrm{unsafe},k}$ (\%)")):
                data = col(scen, key)
                hi = max(v * 100 for m in METHODS for v in data[m] if v is not None)
                self._strip_box(ax, data, METHODS, scale=100)
                ax.set_ylabel(sym)
                ax.set_xlabel("Method")
                ax.set_ylim(-0.05 * hi, 1.08 * hi)
            self._hide_top_xlabel(axes)
            self._fill_width(fig)
            self._save(fig, f"distributions_default_{scen}_row2",
                       dict(figure=f"distributions_default_{scen}_row2", metrics=["q_CBF", "q_unsafe"], unit="%"))

            # (3) estimated peak exceedance risk and minimum barrier value per trajectory
            fig, axes = self._figure(150, 2)
            data = col(scen, "p_max")
            self._strip_box(axes[0], data, METHODS, na_y=0.5)
            axes[0].axhline(0.15, color="black", linestyle=":", linewidth=LW_REF, label=r"$\alpha_0$", zorder=0.5)
            axes[0].set_ylabel(r"$\widehat p_k^{\max}$")
            axes[0].set_ylim(-0.04, 1.04)
            axes[0].set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
            axes[0].legend(loc="upper right", handlelength=1.6, fontsize=self.txt)
            na = [LABEL[m] for m in METHODS if all(v is None for v in data[m])]
            data = col(scen, "h_min")
            self._strip_box(axes[1], data, METHODS)
            axes[1].axhline(0.0, color="#555555", linestyle="-", linewidth=0.8, zorder=0.4)
            axes[1].set_ylabel(r"$h_k^{\min}$ (m)")
            for ax in axes:
                ax.set_xlabel("Method")
            self._hide_top_xlabel(axes)
            self._fill_width(fig)
            self._save(fig, f"distributions_default_{scen}_row3",
                       dict(figure=f"distributions_default_{scen}_row3", metrics=["p_max", "h_min"],
                            p_max_undefined_for=na, alpha0=0.15))

            # (4) decision fractions of the four certified variants (means over trajectories, %)
            ncol = 3 if self.kind == "paper" else 2     # the 11 pt legend needs two rows at the report width
            fig, ax = self._figure(112 if ncol == 3 else 130, 1)
            parts = (("q_nocert", "#B0B0B0", "", r"$\widehat q^{+}_{\mathrm{nocert}}$: no certified candidate"),
                     ("q_notrans", "#0072B2", "//", r"$\widehat q^{+}_{\mathrm{nodep}}$: certified, not deployable"),
                     ("q_use", "#E69F00", "..", r"$\widehat q^{+}_{\mathrm{use}}$: margin below $c_{\max}$"))
            means = {}
            for j, m in enumerate(CERTIFIED):
                left = 0.0
                for key, colour, hatch, _ in parts:
                    vals = [r[key] for r in cells[f"{scen}|{m}"]["values"]]
                    mu = 100 * float(np.mean(vals))
                    means[(m, key)] = mu
                    ax.barh(j, mu, left=left, height=0.62, color=colour, hatch=hatch, edgecolor="white", linewidth=0.6)
                    left += mu
                if abs(left - 100) > 1e-6:
                    raise ValueError(f"decision fractions of {scen} {m} sum to {left}")
            ax.set_yticks(range(4), [LABEL[m] for m in CERTIFIED])
            ax.set_ylim(-0.5, 5.6 if ncol == 3 else 6.9)   # empty band above the bars holds the legend
            ax.set_xlim(0, 100)
            ax.set_xticks([0, 20, 40, 60, 80, 100])
            ax.set_xlabel(r"Decisions (\%)")
            ax.set_ylabel("Method")
            ax.grid(axis="y", visible=False)
            handles = [Patch(facecolor=c, hatch=h, edgecolor="white", label=l) for _, c, h, l in parts]
            ax.legend(handles=handles, loc="upper center", ncol=ncol, fontsize=self.txt, handlelength=1.4,
                      columnspacing=0.8, handletextpad=0.3, borderaxespad=0.25)
            self._fill_width(fig)
            self._save(fig, f"decisions_default_{scen}_K50",
                       dict(figure=f"decisions_default_{scen}_K50", means_percent={f"{m}|{k}": round(v, 3)
                                                                                    for (m, k), v in means.items()}))

    # ------------------------------------------------------------------ theorem checks at T3 V-L
    def theorem_checks(self, anchors: dict):
        rows = {}
        for v in anchors["views"]:
            d = v["definition"]
            if d["family"] != "detailed_200" or "role3" not in (d.get("anchor_id") or ""):
                continue
            c = v["cells"][0]
            cf = c["configuration"]
            if cf["method"] != "V-Cd":
                continue
            key = (cf["benchmark"], str(cf["K"]), cf["L_mu"])
            if (cf["alpha0"], cf["alpha"], cf["delta"], cf["budget_normalizer"], cf["bound_scale"]) != (0.3, 0.2, 0.2, 2.11, 1.0):
                raise ValueError(f"unexpected T3 configuration {cf}")
            if c["n"] != 200 or key in rows:
                raise ValueError(f"T3 cell {key}")
            rows[key] = dict(view_id=v["view_id"], configuration_id=c["configuration_id"], cell=c)
        delta, b = 0.20, 2.11
        delta1 = delta / (2 * b * math.log(2) ** 2)
        groups = [("50", 0.0005), ("50", 0.005), ("infinity", 0.0005), ("infinity", 0.005)]
        glabel = {("50", 0.0005): r"$K{=}50$" "\n" r"$L_\mu{=}0.0005$", ("50", 0.005): r"$K{=}50$" "\n" r"$L_\mu{=}0.005$",
                  ("infinity", 0.0005): r"$K{=}\infty$" "\n" r"$L_\mu{=}0.0005$",
                  ("infinity", 0.005): r"$K{=}\infty$" "\n" r"$L_\mu{=}0.005$"}
        for scen in SCEN:
            fig, axes = self._figure(165, 2, stack_height_pt=265)
            ax = axes[0]
            info = dict(figure=f"theorem_checks_{scen}_New3", setting="T3 (0.30,0.20,0.20,50) V-L", b=b, phi=1.0,
                        delta=delta, delta1=delta1, groups=[])
            for gi, g in enumerate(groups):
                c = rows[(scen, *g)]["cell"]
                level = 100 * (delta1 if g[0] == "infinity" else delta)
                ax.plot([gi - 0.38, gi + 0.38], [level, level], color="#D55E00", linestyle="--", linewidth=LW_REF,
                        label="_nolegend_", zorder=1)
                entry = dict(K=g[0], L_mu=g[1], view_id=rows[(scen, *g)]["view_id"], level_percent=level)
                for off, fam, marker, face in ((-0.13, "qcert", "o", VL), (0.13, "qdep", "s", "white")):
                    e = c["trajectory_events"][fam]
                    n = e["applicable_n"]
                    count = e["possible"]                    # conservative count: confirmed + numerically possible
                    est, ub = 100 * count / n, 100 * e["upper95_possible"]
                    ax.errorbar([gi + off], [est], yerr=[[0], [ub - est]], fmt=marker, markersize=MS + 0.5,
                                markerfacecolor=face, markeredgecolor=VL, markeredgewidth=0.8, ecolor="black",
                                elinewidth=1.5, capsize=3, capthick=1.5, zorder=3)
                    entry[fam] = dict(confirmed=e["confirmed"], possible=count, n=n, upper95_percent=ub,
                                      meets_level=ub <= level + 1e-12)
                    if ub > level + 1e-12:
                        raise ValueError(f"{scen} {g} {fam}: bound {ub:.2f}% above level {level:.2f}%")
                info["groups"].append(entry)
            ax.set_xticks(range(4), [glabel[g] for g in groups])
            ax.set_xlim(-0.55, 3.55)
            ax.set_ylim(-0.6, 23.5)
            ax.set_ylabel(r"Trajectories (\%)")
            ax.set_xlabel("Test period and drift bound")
            ax.xaxis.label.set_visible(False)
            ax.tick_params(axis="x", length=0, pad=3)
            ax.grid(axis="x", visible=False)
            handles = [Line2D([], [], marker="o", linestyle="none", markersize=MS + 0.5, markerfacecolor=VL,
                              markeredgecolor=VL, label=r"false certification $\widehat q_{\mathrm{cert}}$"),
                       Line2D([], [], marker="s", linestyle="none", markersize=MS + 0.5, markerfacecolor="white",
                              markeredgecolor=VL, label=r"deployment $\widehat q_{\mathrm{dep}}$"),
                       Line2D([], [], color="#D55E00", linestyle="--", linewidth=LW_REF, label=r"level $\delta$ or $\delta_1$")]
            # empty band between the two levels (delta_1 ~ 9.9 %, delta = 20 %) holds the legend
            if self.stacked:   # full-width panel: the right half above delta_1 is empty
                ax.legend(handles=handles, loc="upper right", fontsize=self.txt, handlelength=1.5, borderaxespad=0.3)
            else:
                ax.legend(handles=handles, loc="center", fontsize=self.txt, handlelength=1.5,
                          bbox_to_anchor=(0.5, (15.0 + 0.6) / 24.1))

            # (b) risk-increase bound checks: excess counts at checked decisions
            ax = axes[1]
            for gi, g in enumerate(groups):
                c = rows[(scen, *g)]["cell"]
                ri = c["risk_increase_checks"]
                for off, lab, marker, face in ((-0.13, "direct", "D", VL), (0.13, "pairwise", "^", "white")):
                    r = ri[lab]
                    n = r["total_trajectories"]
                    k = r["any_excess"]
                    ub = 100 * (1 - 0.05 ** (1 / n)) if k == 0 else None
                    if k != 0:
                        raise ValueError(f"{scen} {g} {lab}: {k} trajectories with an excess; extend the plot")
                    ax.errorbar([gi + off], [100 * k / n], yerr=[[0], [ub]], fmt=marker, markersize=MS + 0.5,
                                markerfacecolor=face, markeredgecolor=VL, markeredgewidth=0.8, ecolor="black",
                                elinewidth=1.5, capsize=3, capthick=1.5, zorder=3)
                    info["groups"][gi][f"risk_increase_{lab}"] = dict(r, upper95_percent=ub)
                ax.text(gi, 2.05, f"{ri['direct']['assessed_decisions']:,}".replace(",", r"\,"), ha="center",
                        va="bottom", fontsize=self.txt)
            ax.text(-0.5, 2.75, r"checked decisions:", ha="left", va="bottom", fontsize=self.txt)
            ax.set_xticks(range(4), [glabel[g] for g in groups])
            ax.set_xlim(-0.55, 3.55)
            ax.set_ylim(-0.15, 4.1)
            ax.set_yticks([0, 1, 2])
            ax.set_ylabel(r"Trajectories (\%)")
            ax.set_xlabel("Test period and drift bound")
            ax.xaxis.label.set_visible(False)
            ax.tick_params(axis="x", length=0, pad=3)
            ax.grid(axis="x", visible=False)
            handles = [Line2D([], [], marker="D", linestyle="none", markersize=MS + 0.5, markerfacecolor=VL,
                              markeredgecolor=VL, label=r"$\widehat q_{\mathrm{A2}}$"),
                       Line2D([], [], marker="^", linestyle="none", markersize=MS + 0.5, markerfacecolor="white",
                              markeredgecolor=VL, label=r"$\widehat q_{\mathrm{pw}}$")]
            ax.legend(handles=handles, loc="upper right", ncol=2, fontsize=self.txt, handlelength=1.0,
                      borderaxespad=0.3, columnspacing=1.0)
            self._fill_width(fig)
            self._save(fig, f"theorem_checks_{scen}_New3", info)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(KINDS), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--only", choices=["default", "checks"], default=None)
    args = ap.parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)
    b = Builder(args.kind, args.out)
    if args.only in (None, "default"):
        b.default_figures(json.loads(DEFAULT_CACHE.read_text()))
    if args.only in (None, "checks"):
        b.theorem_checks(json.loads(NEW_ANCHORS.read_text()))
    manifest = dict(generator=recorded_path(Path(__file__)), generator_sha256=sha_file(Path(__file__).resolve()),
                    helper_sha256=sha_file(FIGURE_STYLE / "ieee_figure.py"), kind=args.kind,
                    inputs={recorded_path(DEFAULT_CACHE): sha_file(DEFAULT_CACHE), recorded_path(NEW_ANCHORS): sha_file(NEW_ANCHORS)},
                    figures=b.manifest)
    (args.out / f"figures_manifest_{args.kind}{'_' + args.only if args.only else ''}.json").write_text(
        json.dumps(manifest, indent=1) + "\n")


if __name__ == "__main__":
    main()
