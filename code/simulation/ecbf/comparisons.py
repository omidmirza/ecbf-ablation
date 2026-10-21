"""Plots of fixed-setting distributions, event frequencies and decision shares from saved trajectory records."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from scipy.stats import beta

from .metrics import trajectory

WIDTH, HEIGHT = 245.7 / 72.27, 180 / 72.27
METHODS = ("c0", "cmax", "ACI", "QT", "V-U", "V-Delta", "V-L", "V-C")
CERTIFIED = ("V-U", "V-Delta", "V-L", "V-C")
LABELS = {"c0": r"$c=0$", "cmax": r"$c_{\max}$", "V-Delta": r"V-$\Delta$"}
COLOURS = {"c0": "#666666", "cmax": "#B0B0B0", "ACI": "#CC79A7", "QT": "#F0E442",
           "V-U": "#56B4E9", "V-Delta": "#0072B2", "V-L": "#000000", "V-C": "#E69F00"}
MARKERS = dict(zip(METHODS, ("v", "s", "D", "P", "o", "h", "p", "^")))
MODEL_FIELDS = ("benchmark", "T", "H", "c_max", "alpha0", "rho_buf", "L_mu", "beta_mode",
                "grid", "model_parameters")
CERTIFICATE_FIELDS = ("alpha", "delta", "kappa", "fraction", "policy", "budget_mode",
                      "budget_normalizer")
OMISSION = "Unfiltered is omitted to retain a readable comparison scale."
BOX_NOTE = ("Each point is one trajectory. Boxes show the median and interquartile range; "
            "whiskers extend to the most distant observations within 1.5 interquartile ranges.")
RISK_NOTE = ("Solid bars count confirmed events; the hatched extension includes unresolved events "
             r"and incomplete trajectories. Error bars give pointwise one-sided 95\% "
             "Clopper--Pearson upper bounds for the possible-event counts. "
             "Numerical risk enclosures include floating-point allowances and are not machine-certified.")
CAPTIONS = {
    "mission_time": r"Mission time $T_{\mathrm{mis},k}$, with the horizon assigned to unsuccessful trajectories.",
    "mean_margin": r"Mean applied margin $m_{c,k}/c_{\max}$.",
    "cbf_fraction": r"Fraction of decisions with a discrete CBF violation, $q^{+}_{\mathrm{CBF},k}$.",
    "unsafe_fraction": r"Fraction of decisions outside the safe set, $q^{+}_{\mathrm{unsafe},k}$.",
    "margin_use": r"Fraction of decisions using $c_t<c_{\max}$, $q^{+}_{\mathrm{use},k}$.",
    "peak_risk": (r"Maximum quadrature-estimated risk on decisions using $c_t<c_{\max}$, "
                  r"$\widehat p^{\max}_k$. Only trajectories with such decisions are included; "
                  r"an empty category means that the metric is undefined."),
    "first_use": (r"Time to the first decision using $c_t<c_{\max}$, $\tau^{\mathrm{use}}_k$. "
                  r"Only trajectories that use such a margin are included; an empty category is undefined."),
    "minimum_barrier": r"Minimum sampled barrier value along each trajectory.",
    "deployment_events": (r"Trajectory frequency of a decision with $p_t(c_t)>\alpha_0$, "
                          r"$\widehat q_{\mathrm{dep}}$. " + RISK_NOTE),
    "certification_events": (r"Trajectory frequency of a false certificate, $\widehat q_{\mathrm{cert}}$. "
                             + RISK_NOTE),
    "cbf_events": (r"Trajectory frequency of at least one discrete CBF violation, $\widehat q_{\mathrm{CBF}}$. "
                   r"Error bars give pointwise one-sided 95\% Clopper--Pearson upper bounds. "
                   r"This is an observed violation frequency, not the theorem's risk-event probability."),
    "contact_events": (r"Trajectory frequency of sampled geometric contact, $\widehat q_{\mathrm{contact}}$. "
                       r"Error bars give pointwise one-sided 95\% Clopper--Pearson upper bounds. "
                       r"This is an observed contact frequency, not the theorem's risk-event probability."),
    "decision_shares": (r"Decision shares averaged within each trajectory and then across trajectories: "
                        r"no certificate, a certificate that is not deployable, and use of $c_t<c_{\max}$."),
}


def _style():
    plt.rcParams.update({
        "text.usetex": True, "font.family": "serif", "font.serif": ["Computer Modern Roman"],
        "font.size": 10, "axes.labelsize": 10, "xtick.labelsize": 10, "ytick.labelsize": 10,
        "legend.fontsize": 10, "axes.linewidth": 1.2, "xtick.major.width": 1.2,
        "ytick.major.width": 1.2, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "axes.labelpad": 2, "legend.borderpad": 0.2, "legend.labelspacing": 0.15,
        "legend.handlelength": 1.2, "legend.handletextpad": 0.35,
        "text.latex.preamble": r"\usepackage{amsmath}\DeclareMathSizes{10}{10}{10}{10}",
        "pdf.fonttype": 42,
    })


def _load(data_dir, scenario):
    groups = {}
    for method in (*METHODS, "unfiltered"):
        folder = Path(data_dir) / scenario / method
        if not folder.exists():
            continue
        records = []
        for path in sorted(folder.glob("*.npz")):
            with np.load(path, allow_pickle=False) as arrays:
                rec = {key: arrays[key].item() if arrays[key].ndim == 0 else arrays[key].copy()
                       for key in arrays.files}
            if rec.get("benchmark") != scenario:
                raise ValueError(f"{path}: wrong scenario")
            if rec.get("public_method") != method:
                raise ValueError(f"{path}: wrong method")
            records.append(rec)
        if not records:
            raise ValueError(f"{folder}: no trajectory records")
        records.sort(key=lambda r: int(r["seed"]))
        seeds = [int(r["seed"]) for r in records]
        if len(seeds) != len(set(seeds)):
            raise ValueError(f"{folder}: duplicate seeds")
        fields = MODEL_FIELDS + (CERTIFICATE_FIELDS if method in CERTIFIED else ())
        for rec in records:
            for field in fields:
                if field not in rec or not np.array_equal(rec[field], records[0][field]):
                    raise ValueError(f"{folder}: inconsistent {field}")
        groups[method] = records
    if not groups:
        raise FileNotFoundError(f"No {scenario} records in {data_dir}; this function reads the saved trajectory "
                                "records (.npz) of the study, which are not part of this repository")
    first = next(iter(groups.values()))
    certified_first = next((groups[m][0] for m in CERTIFIED if m in groups), None)
    for method, records in groups.items():
        if [int(r["seed"]) for r in records] != [int(r["seed"]) for r in first]:
            raise ValueError(f"{scenario}: method seed sets differ")
        for field in MODEL_FIELDS:
            if not np.array_equal(records[0][field], first[0][field]):
                raise ValueError(f"{scenario}: methods have different {field}")
        if method in CERTIFIED:
            for field in CERTIFICATE_FIELDS:
                if not np.array_equal(records[0][field], certified_first[field]):
                    raise ValueError(f"{scenario}: variants have different {field}")
    return {m: groups[m] for m in METHODS if m in groups}


def _axes(methods, ylabel):
    fig = plt.figure(figsize=(WIDTH, HEIGHT))
    ax = fig.add_axes([0.23, 0.235, 0.755, 0.735])
    ax.set_xticks(range(len(methods)), [LABELS.get(m, m) for m in methods], rotation=38, ha="right")
    ax.set_xlim(-0.6, len(methods) - 0.4)
    ax.set_xlabel("Method")
    ax.set_ylabel(ylabel)
    ax.grid(axis="y", color="#D9D9D9", linewidth=1.2)
    ax.set_axisbelow(True)
    return fig, ax


def _save(fig, path):
    """Fit the complete vector content to one column, without scaling its text."""
    for _ in range(8):
        fig.canvas.draw()
        width = fig.get_tightbbox(fig.canvas.get_renderer()).width
        if abs(width - WIDTH) < 0.0001:
            break
        fig.set_size_inches(fig.get_figwidth() + WIDTH - width, fig.get_figheight(), forward=True)
    fig.savefig(path, bbox_inches="tight", pad_inches=0)
    plt.close(fig)


def _boxplot(groups, metrics, key, scale, ylabel, path):
    methods = list(groups)
    fig, ax = _axes(methods, ylabel)
    finite = []
    for i, method in enumerate(methods):
        values = np.asarray([r[key] for r in metrics[method]], float) * scale
        values = values[np.isfinite(values)]
        finite.extend(values)
        if not values.size:
            continue
        jitter = np.random.default_rng(1000 + i).uniform(-0.10, 0.10, len(values))
        ax.scatter(i - 0.14 + jitter, values, s=10, marker=MARKERS[method],
                   facecolor=COLOURS[method], edgecolor="white", linewidth=0.5, alpha=0.7, zorder=3)
        ax.boxplot([values], positions=[i + 0.16], widths=0.24, patch_artist=True, showfliers=False,
                   manage_ticks=False, boxprops=dict(facecolor=COLOURS[method], edgecolor="#182B38", linewidth=1.5),
                   medianprops=dict(color="#182B38", linewidth=1.8),
                   whiskerprops=dict(color="#182B38", linewidth=1.5),
                   capprops=dict(color="#182B38", linewidth=1.5))
    if not finite:
        ax.set_ylim(0, 1)
    else:
        ax.margins(y=0.10)
    if key == "p_max":
        alpha0 = float(next(iter(groups.values()))[0]["alpha0"])
        ax.axhline(alpha0, color="#444444", linestyle=":", linewidth=1.8, label=r"$\alpha_0$")
        lower, upper = ax.get_ylim()
        ax.set_ylim(lower, upper + 0.25 * (upper - lower))
        ax.legend(loc="upper right", frameon=False)
    _save(fig, path)


def _upper(count, total):
    return 1.0 if count == total else float(beta.ppf(0.95, count + 1, total - count))


def _eventplot(groups, event, ylabel, path):
    certified_only = event == "qcert"
    methods = [m for m in groups if not certified_only or m in CERTIFIED]
    if not methods:
        return False
    fig, ax = _axes(methods, ylabel)
    risk_event = event in ("qdep", "qcert")
    for i, method in enumerate(methods):
        records = groups[method]
        if risk_event:
            keys = [f"assessment_{event}_{part}" for part in ("confirmed", "possible")]
            if any(key not in rec for rec in records for key in keys):
                raise ValueError("Risk-event plots need the assessment saved in the trajectory records (.npz) of the "
                                 "study, which are not part of this repository")
            bounds = np.array([[rec[key] for key in keys] for rec in records], float)
            if not np.isfinite(bounds).all() or not np.isin(bounds, (0, 1)).all():
                raise ValueError(f"{method}: missing or invalid {event} classification")
            if np.any(bounds[:, 0] > bounds[:, 1]):
                raise ValueError(f"{method}: confirmed events exceed possible events")
            lower, upper = bounds.sum(axis=0).astype(int)
        else:
            flags = [np.any(rec["cbf"]) if event == "cbf" else np.any(np.asarray(rec["clearance"]) <= 0)
                     for rec in records]
            lower = upper = int(sum(flags))
        n = len(records)
        value, possible = 100 * lower / n, 100 * upper / n
        ax.bar(i, value, width=0.62, color=COLOURS[method], edgecolor="#182B38", linewidth=1.5, zorder=3)
        if upper > lower:
            ax.bar(i, possible - value, bottom=value, width=0.62, color=COLOURS[method],
                   edgecolor="#182B38", hatch="///", linewidth=1.5, alpha=0.45, zorder=3)
        bound = 100 * _upper(upper, n)
        ax.errorbar(i, possible, yerr=[[0], [max(0, bound - possible)]], fmt="none",
                    ecolor="#182B38", elinewidth=1.5, capsize=3, capthick=1.5, zorder=4)
    if risk_event:
        certified = next((groups[m][0] for m in methods if m in CERTIFIED), None)
        handles = [Patch(facecolor="#A0A0A0", edgecolor="#182B38", label="Confirmed"),
                   Patch(facecolor="#E0E0E0", edgecolor="#182B38", hatch="///", label="Possible")]
        if certified is not None and float(certified["fraction"]) >= 1 and float(certified["kappa"]) >= 1:
            delta = float(certified["delta"])
            ax.axhline(100 * delta, color="#444444", linestyle=":", linewidth=1.8, zorder=2)
            handles.append(Line2D([], [], color="#444444", linestyle=":", linewidth=1.8, label=r"$\delta$"))
        ax.legend(handles=handles, loc="upper left", ncol=len(handles), columnspacing=0.5, frameon=False)
    ax.set_ylim(-2, 133 if risk_event else 110)
    ax.set_yticks([0, 25, 50, 75, 100])
    _save(fig, path)
    return True


def _shares(groups, metrics, path):
    methods = [m for m in CERTIFIED if m in groups]
    if not methods:
        return False
    fig, ax = _axes(methods, r"Decisions (\%)")
    columns = ("q_nocert", "q_notrans", "q_use")
    values = []
    for method in methods:
        rows = np.asarray([[rec[key] for key in columns] for rec in metrics[method]], float)
        rows = rows[np.isfinite(rows).all(axis=1)]
        if not len(rows):
            raise ValueError(f"{method}: no nonempty trajectory for decision shares")
        if not np.allclose(rows.sum(axis=1), 1, atol=1e-12):
            raise ValueError(f"{method}: decision categories do not partition executed decisions")
        values.append(100 * rows.mean(axis=0))
    values = np.asarray(values)
    bottom = np.zeros(len(methods))
    styles = (("#B0B0B0", "///", "No certificate"), ("#E69F00", "", "Not deployable"),
              ("#56B4E9", "", "Margin use"))
    for j, (colour, hatch, label) in enumerate(styles):
        ax.bar(range(len(methods)), values[:, j], bottom=bottom, width=0.65, color=colour,
               edgecolor="white", linewidth=1.5, hatch=hatch, label=label, zorder=3)
        bottom += values[:, j]
    ax.set_ylim(0, 155)
    ax.set_yticks([0, 25, 50, 75, 100])
    ax.legend(loc="upper left", frameon=False)
    _save(fig, path)
    return True


def settings_caption(groups):
    """State the values saved in the generating records, not the current settings file."""
    certified = next((groups[m][0] for m in CERTIFIED if m in groups), None)
    record = certified if certified is not None else next(iter(groups.values()))[0]
    n = len(next(iter(groups.values())))
    levels = rf"$\alpha_0={record['alpha0']:g}$"
    if certified is not None:
        levels = rf"$(\alpha,\alpha_0,\delta)=({record['alpha']:g},{record['alpha0']:g},{record['delta']:g})$"
        policy = json.loads(record["policy"])
        if policy.get("period") is not None:
            levels += rf", $K={policy['period']}$"
        levels += rf", $b={record['budget_normalizer']:g}$, $\phi={record['fraction']:g}$"
    caption = (levels + rf", $L_\mu={record['L_mu']:g}$, $\rho_{{\mathrm{{buf}}}}={record['rho_buf']:g}$ m, "
               rf"$N_{{\mathrm{{MC}}}}={n}$ per method.")
    if certified is not None and (float(record["fraction"]) < 1 or float(record["kappa"]) < 1):
        caption += r" A bound factor $\phi<1$ is diagnostic; the guarantees do not apply."
    return caption


def make_comparisons(data_dir, scenario, out_dir):
    """Write vector PDFs and return their scientific captions; never modify raw runs."""
    groups = _load(data_dir, scenario)
    if not groups:
        raise ValueError("At least one filtered method is required for comparisons")
    _style()
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = {method: [trajectory(rec) for rec in records] for method, records in groups.items()}
    caption_settings = settings_caption(groups)
    captions = {}
    definitions = (
        ("mission_time", "T_mis", 1, r"$T_{\mathrm{mis},k}$ (s)"),
        ("mean_margin", "m_c", 1, r"$m_{c,k}/c_{\max}$ (1)"),
        ("cbf_fraction", "q_CBF", 100, r"$q^{+}_{\mathrm{CBF},k}$ (\%)"),
        ("unsafe_fraction", "q_unsafe", 100, r"$q^{+}_{\mathrm{unsafe},k}$ (\%)"),
        ("margin_use", "q_use", 100, r"$q^{+}_{\mathrm{use},k}$ (\%)"),
        ("peak_risk", "p_max", 1, r"$\widehat p^{\max}_k$ (1)"),
        ("first_use", "tau_use", 1, r"$\tau^{\mathrm{use}}_k$ (s)"),
        ("minimum_barrier", "h_min", 1000, r"$\min_t h(x_{k,t})$ (mm)"),
    )
    for name, key, scale, label in definitions:
        filename = f"{scenario}_{name}.pdf"
        _boxplot(groups, metrics, key, scale, label, out_dir / filename)
        captions[filename] = f"{scenario.capitalize()}: {CAPTIONS[name]} {BOX_NOTE} {caption_settings} {OMISSION}"
    for name, event, label in (
        ("deployment_events", "qdep", r"$\widehat q_{\mathrm{dep}}$ (\%)"),
        ("certification_events", "qcert", r"$\widehat q_{\mathrm{cert}}$ (\%)"),
        ("cbf_events", "cbf", r"$\widehat q_{\mathrm{CBF}}$ (\%)"),
        ("contact_events", "contact", r"$\widehat q_{\mathrm{contact}}$ (\%)"),
    ):
        filename = f"{scenario}_{name}.pdf"
        if _eventplot(groups, event, label, out_dir / filename):
            guide = ""
            if event in ("qdep", "qcert"):
                certified = next((groups[m][0] for m in CERTIFIED if m in groups), None)
                if certified is not None:
                    valid = float(certified["fraction"]) >= 1 and float(certified["kappa"]) >= 1
                    guide = (r" The dotted $\delta$ line is the theoretical bound of the variants, "
                             "conditional on the theorem assumptions; finite numerical checks do not prove those assumptions."
                             if valid else r" A bound factor $\phi<1$ is diagnostic; the guarantees do not apply.")
            captions[filename] = f"{scenario.capitalize()}: {CAPTIONS[name]}{guide} {caption_settings} {OMISSION}"
    filename = f"{scenario}_decision_shares.pdf"
    if _shares(groups, metrics, out_dir / filename):
        captions[filename] = f"{scenario.capitalize()}: {CAPTIONS['decision_shares']} {caption_settings}"
    return captions
