#!/usr/bin/env python3
"""Per-trajectory distribution figures.

distributions_<scen>_<label>.pdf (label Default, New1..New5 = T1..T5): six panels in the layout of the paper's
distribution figure: mission time T_mis,k (s); mean robustness margin
m_c,k/c_max; CBF violations q+_CBF,k (%); safe-set violations q+_unsafe,k (%); estimated peak exceedance risk
p-hat_k^max with alpha0 dotted; minimum barrier value h_k^min (m). Per method (c=0, c_max, ACI, QT at the setting's
alpha0; V-U, V-Delta, V-L, V-C at the setting): jittered per-trajectory points (left) and a box plot (right; box =
quartiles, whiskers = 1.5 IQR, no fliers), N_MC = 200 trajectories each. Unfiltered control omitted (caption).
p-hat^max is undefined where no smaller-margin decision is executed (c_max): N-A.
Data: setting_assessments.json (load_setting_assessments.py: saved assessments, hash-checked against
N200_INDEX.sqlite, means equal to the N200 export). Saved through the figure helper (strict save), with the paper's
method colours.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import n200_style as S
from n200_style import plt, np, fs

from matplotlib.lines import Line2D  # noqa: E402
import matplotlib.colors as mcolors  # noqa: E402

HERE = Path(__file__).resolve().parent
CACHE = HERE / "setting_assessments.json"
METHODS = S.BASELINES + S.CERTIFIED
PANELS = (("T_mis", "Mission time", r"$T_{\mathrm{mis},k}$ (s)", 1.0),
          ("m_c", "Mean robustness margin", r"$m_{c,k}/c_{\max}$", 1.0),
          ("q_CBF", "CBF violations", r"$q^+_{\mathrm{CBF},k}$ (\%)", 100.0),
          ("q_unsafe", "Safe-set violations", r"$q^+_{\mathrm{unsafe},k}$ (\%)", 100.0),
          ("p_max", "Estimated peak exceedance risk", r"$\widehat p_k^{\max}$", 1.0),
          ("h_min", "Minimum barrier value", r"$h_k^{\min}$ (m)", 1.0))


def light(c, f=0.55):
    r, g, b = mcolors.to_rgb(c)
    return (1 - f * (1 - r), 1 - f * (1 - g), 1 - f * (1 - b))


def strip_box(ax, data, scale, txt, na_y):
    rng = np.random.default_rng(0)                  # jitter positions only
    stats = {}
    for i, m in enumerate(METHODS):
        vals = np.array([v * scale for v in data[m] if v is not None], dtype=float)
        colour, marker, _ = S.STYLE[m]
        if len(vals) == 0:
            ax.text(i, na_y, "N-A", ha="center", va="center", fontsize=txt)
            stats[m] = None
            continue
        x = i - 0.17 + rng.uniform(-0.12, 0.12, len(vals))
        ax.plot(x, vals, linestyle="none", marker=marker, markersize=S.MS, markerfacecolor=colour,
                markeredgecolor=S.edge(colour), markeredgewidth=S.MEW, alpha=0.8, label="_nolegend_", zorder=3)
        bp = ax.boxplot(vals, positions=[i + 0.23], widths=0.24, whis=1.5, showfliers=False, patch_artist=True,
                        manage_ticks=False, zorder=4)
        for part in ("whiskers", "caps", "medians"):
            for ln in bp[part]:
                ln.set_linewidth(S.LW_BOX)
                ln.set_color("black")
        for b in bp["boxes"]:
            b.set_facecolor(light(colour))
            b.set_edgecolor("black")
            b.set_linewidth(S.LW_BOX)
        q1, med, q3 = np.percentile(vals, [25, 50, 75])
        stats[m] = dict(n=int(len(vals)), mean=float(vals.mean()), median=float(med), q1=float(q1), q3=float(q3),
                        min=float(vals.min()), max=float(vals.max()))
    ax.set_xlim(-0.6, len(METHODS) - 0.4)
    ax.grid(axis="x", visible=False)
    return stats


def build(kind, out, labels, scens):
    width_in, width_tex, txt = S.setup(kind)
    cache = json.loads(CACHE.read_text())["cells"]
    manifest = []
    height_pt = 330.0 if kind == "paper" else 360.0
    for label in labels:
        a0 = S.SETTINGS[label][0]
        for scen in scens:
            name = f"distributions_{scen}_{label}"
            fig, axes = plt.subplots(3, 2, figsize=(width_in, height_pt / 72.27), sharex=True, squeeze=False)
            fig.subplots_adjust(wspace=0.25, hspace=0.30)
            info = dict(figure=name, setting=label, tuple=S.SETTINGS[label][:4], paired=S.SETTINGS[label][4],
                        alpha0=a0, methods=list(METHODS), n_mc=200, panels={})
            for k, (key, title, ylab, scale) in enumerate(PANELS):
                ax = axes.flat[k]
                data = {m: [r[key] for r in cache[f"{label}|{scen}|{m}"]["values"]] for m in METHODS}
                for m in METHODS:
                    if len(data[m]) != 200:
                        raise ValueError(f"{name} {m}: {len(data[m])} trajectories")
                info["panels"][key] = strip_box(ax, data, scale, txt, na_y=0.5)
                ax.set_title(title, loc="left", pad=2, fontsize=txt)
                ax.set_ylabel(ylab)
                ax.set_xlabel("Method")
                ax.set_xticks(range(len(METHODS)), [S.LABEL[m] for m in METHODS])
                ax.tick_params(axis="x", length=0, pad=3)
                if k < 4:
                    ax.xaxis.label.set_visible(False)
                if key in ("q_CBF", "q_unsafe"):
                    ax.set_yscale("symlog", linthresh=0.1, linscale=0.8)
                    ax.set_ylim(-0.03, 160)
                    ax.set_yticks([0, 0.1, 1, 10, 100], ["0", "0.1", "1", "10", "100"])
                    ax.yaxis.set_minor_locator(plt.NullLocator())
                if key in ("m_c", "p_max"):
                    ax.set_ylim(-0.06, 1.1)
                    ax.set_yticks([0, 0.5, 1], ["0", "0.5", "1"])
                if key == "p_max":
                    ax.axhline(a0, color=S.REF_COLOUR, ls=":", lw=S.LW_REF, zorder=2, label="_nolegend_")
                if key == "h_min":
                    ax.axhline(0.0, color="#555555", ls="-", lw=0.8, zorder=1, label="_nolegend_")
            S.fill_box(fig, left=0.08, right=0.995, bottom=0.08, top=0.97)
            h = [Line2D([], [], color=S.REF_COLOUR, ls=":", lw=S.LW_REF, label=r"$\alpha_0$")]
            ax_p = axes.flat[4]
            axl, leg = S.place_legend([ax_p], h, txt, ncols=(1,), handlelength=1.8, borderaxespad=0.25,
                                      framealpha=1.0)
            if leg is None:
                raise ValueError(f"{name}: no free place for the alpha0 legend")
            bad = S.text_overlaps(fig)
            if bad:
                raise ValueError(f"{name}: overlapping texts {bad[:6]}")
            pdf = fs.save(fig, out / name, png=False, strict=True)
            plt.close(fig)
            w, hh = S.check_width(pdf, width_tex)
            info.update(pdf=S.recorded_path(pdf), sha256=S.sha_file(pdf), width_texpt=w, height_texpt=hh, text_pt=txt,
                        legend_loc=leg._loc)
            manifest.append(info)
            print(name, kind, f"{w} x {hh} TeX pt", flush=True)
    return manifest


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--kind", choices=sorted(S.KINDS), required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--labels", nargs="+", default=list(S.SETTINGS))
    ap.add_argument("--scens", nargs="+", default=list(S.SCEN))
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    man = build(a.kind, a.out, a.labels, a.scens)
    (a.out / f"distributions_manifest_{a.kind}_{'_'.join(a.labels)}.json").write_text(json.dumps(dict(
        generator=S.recorded_path(Path(__file__)), generator_sha256=S.sha_file(Path(__file__).resolve()),
        style_sha256=S.sha_file(HERE / "n200_style.py"), inputs={S.recorded_path(CACHE): S.sha_file(CACHE)},
        figures=man), indent=1, default=str) + "\n")


if __name__ == "__main__":
    main()
