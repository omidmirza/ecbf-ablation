#!/usr/bin/env python3
"""Catalog fragments of the N200 tables.

Called only by ``build_tables.py --catalog`` after the edition fragments are written; the edition fragments and
tables_values.json.gz do not depend on this module. Writes

  appendices/generated_n200/catalog/{paper,report}/<name>.tex   (same layout engine and widths as the editions)
  appendices/generated_n200/catalog/CATALOG_TABLES.json         (manifest: type, parameters, data, restore line)
  tools/n200_tables/catalog_values.json.gz                      (every printed cell; check_tables.py re-derives it)

Fragments (both targets unless stated):
  sweep_{scen}_a{010..030}       full sweep: alpha -> delta -> K -> the four variants, all metrics, per alpha0
  compact_{scen}_a{010..030}     compact sweep: one row per (alpha, delta, K) with the best variant in its own
                                 column (rule below), m_c/c_max, q+_CBF, q+_unsafe, T_mis
  {delta,period,b}_all_{scen}    parameter tables, all four variants at the six settings (Default, S1-S5)
  phi_all_{scen}_{L0005,L005}    bound factor phi at one L_mu, all four variants, six settings
  lmu_all_{scen}                 L_mu in {0.0005, 0.005}, all four variants, six settings
  checks_all_{scen}_{L0005,L005} checks at K of the setting and K = infinity, all four variants
  detailed_events_all_{scen}     q_cert, q_dep at K and K = infinity x L_mu, all four variants
  topfive_events                 q_cert, q_dep of the default and S1-S5, all variants, Corridor and Maze
  default_{scen}, topfive_{scen}, topfive_full_{scen}   default and selected-setting tables (all variants)

Best variant of a compact cell (the bold rule of the key, applied to the four printed columns): the variant with
the most best means among m_c/c_max (lower), q+_CBF (lower), q+_unsafe (lower), T_mis (lower), ties
of a column all best; ties in the number of best means are broken by the smaller printed m_c/c_max, then
q+_unsafe, q+_CBF, T_mis, then the order V-U, V-Delta, V-L, V-C.
Restore line in an edition: ``\\ntable{../catalog/\\nassetkind/<name>}`` (n200_assets.tex prefixes
appendices/generated_n200/<kind>/), i.e. ``\\ntable{../catalog/paper/<name>}`` in the paper and
``\\ntable{../catalog/report/<name>}`` in the report. Layout tags start with Q, so they never meet an edition tag.
"""
from __future__ import annotations

import hashlib
import json
from decimal import Decimal

FIRST_CUT = 12   # value lines of a fragment's first page left for the catalog heading and the long caption
COMPACT_COLS = [["m_c"], ["q_CBF"], ["q_unsafe"], ["T_mis"]]
COMPACT_ORDER = ["m_c", "q_unsafe", "q_CBF", "T_mis"]   # tie-break order
COMPACT_KEY = (r"\textbf{Compact sweep.} One row per setting $(\alpha_0,\alpha,\delta,K)$ of the main grid. "
               r"Variant: the variant with the most best means in the four columns (best: the bold rule "
               r"of the key; ties of a column all count); $(k)$: the number $k$ of the four columns in which it has "
               r"the best mean; ties in $k$: smaller printed $\widehat m_c/c_{\max}$, then "
               r"$\widehat q_{\rm unsafe}^+$, $\widehat q_{\rm CBF}^+$, $\widehat T_{\rm mis}$, then the order "
               r"V-U, V-$\Delta$, V-L, V-C. $\pm$: one sample SD over the trajectories that enter the mean; "
               r"$^+$: decision fraction; $\downarrow$: lower is better.")


def compact_fragment(B, D, scen, alpha0, target, gid):
    lay = B.Layout(gid, [("a", r"$\alpha$"), ("d", r"$\delta$"), ("k", r"$K$")], COMPACT_COLS, 1,
                   ["1.6pt", "1.0pt", "0.6pt"], method_header=r"Variant $(k)$")
    alphas = [x / 100 for x in range(5, int(round(alpha0 * 100)), 5)]
    groups = []
    for a in alphas:
        for d in B.DELTAS:
            rows = []
            for K in B.KS:
                lv = (B.num(a), B.num(d), B.kfmt(K))
                disp = (rf"\rotatebox[origin=c]{{90}}{{$\alpha={lv[0]}$}}",
                        rf"\rotatebox[origin=c]{{90}}{{$\delta={lv[1]}$}}", f"${lv[2]}$")
                cfg = dict(alpha0=alpha0, alpha=a, delta=d, K=K, **B.BASE)
                found = {v: B.cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                         for v in B.VARIANTS}
                cand = [B.vrow({"": found[v]}, v, lv, disp, {}, (scen, a, d, K)) for v in B.VARIANTS]
                B.fill_rows(cand, [c[0] for c in COMPACT_COLS])
                score = []
                for i, r in enumerate(cand):
                    infos = {m: r["cells"][m][1] for m in COMPACT_ORDER}
                    wins = sum(1 for m in COMPACT_ORDER if infos[m].get("bold"))
                    tie = tuple(Decimal(infos[m]["mean_str"]) if infos[m].get("kind") == "meanSD" else Decimal("Infinity")
                                for m in COMPACT_ORDER)
                    score.append((-wins, tie, i))
                best = min(score)
                v = B.VARIANTS[best[2]]
                wins = -best[0]
                r = B.vrow({"": found[v]}, v, lv, disp,
                           dict(scenario=scen, table="compact", compact=True, wins=wins, **cfg), None)
                r["name"] = f"{B.VNAME[v]} ({wins})"
                rows.append(r)
            groups.append((0, len(rows), rows))
    allrows = [r for g in groups for r in g[2]]
    B.fill_rows(allrows, lay.specs)
    brows = B.baseline_rows(D, scen, alpha0, 0.0005)
    B.fill_rows(brows, lay.specs)
    frag = f"compact_{scen}_a{B.tag3(alpha0)}"
    budget = {"paper": (42, 48, 48), "report": (38, 50, 50)}[target]
    pages = B.paginate(groups, budget, True, B.boxes(target), same_sheet=False)
    return B.emit(lay, pages, frag, target, COMPACT_KEY,
                  lambda pg: B.settings_line(scen, alpha0=alpha0, b="2.11", phi="1", L_mu="0.0005"),
                  lambda k: (brows, B.settings_line(scen, tag="Baselines", alpha0=alpha0, L_mu="0.0005")))


def param_all_fragment(B, D, scen, kind, target, gid, L_mu=None):
    """All four variants at the six settings: kind delta / period (main grid),
    b, phi (one L_mu), lmu (anchor sweeps)."""
    cols, L = B.full_cols(target)
    pcol, vals = {"delta": (("p", r"$\delta$"), B.DELTAS), "period": (("p", r"$K$"), B.KS),
                  "b": (("p", r"$b$"), [2.11, 2.5, 3.0, 4.0]), "phi": (("p", r"$\phi$"), [0.0, 0.25, 0.5, 0.75, 1.0]),
                  "lmu": (("p", r"$L_\mu$"), [0.0005, 0.005])}[kind]
    lay = B.Layout(gid, [pcol], cols, L, ["0.6pt"], method_header="Method")   # block row per setting
    if kind == "delta":   # S1 and S2 share (alpha0, alpha, K): one block
        merged = [(0, [0]), (1, [1, 2]), (3, [3]), (4, [4]), (5, [5])]
        fields = ("alpha0", "alpha", "K")
    else:
        merged = [(i, [i]) for i in range(6)]
        fields = ("alpha0", "alpha", "delta") if kind == "period" else ("alpha0", "alpha", "delta", "K")
    groups = []
    for si, members in merged:
        name, t, _, _ = B.SIX[si]
        names = [B.SIX[m][0] for m in members]
        if kind == "delta" and len(members) > 1:   # S1 and S2 coincide once delta is dropped: one tuple
            assert len({B.tuple_plain(B.SIX[m][1], fields) for m in members}) == 1
            stext = (" and ".join(rf"\textbf{{{B.dn(B.SIX[m][0])}}}" for m in members)
                     + rf": {B.fields_sym(fields)} = {B.tuple_plain(t, fields)}")
        else:
            stext = B.block_text(name, t, fields=fields)   # the swept field is not printed
        whole = []
        for val in vals:
            rows = []
            for v in B.VARIANTS:
                if kind == "delta":
                    cfg = B.six_cfg(t, delta=val)
                    f = B.cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                    disp = f"${B.num(val)}$"
                elif kind == "period":
                    cfg = B.six_cfg(t, K=val)
                    f = B.cell_or_fail(D, "final", "main_grid_200", benchmark=scen, method=v, rho_buf=0.0, **cfg)
                    disp = f"${B.kfmt(val)}$"
                elif kind == "b":
                    f, cfg = B.anchor_cell(D, si, scen, v, "budget_normalizer_200", budget_normalizer=val)
                    disp = f"${B.b_value(val)}$"   # the default b prints as kappa
                elif kind == "phi":
                    f, cfg = B.anchor_cell(D, si, scen, v, "stress_drift_200", bound_scale=val, L_mu=L_mu)
                    disp = f"${val:g}$"
                else:
                    f, cfg = B.anchor_cell(D, si, scen, v, "stress_drift_200", L_mu=val)
                    disp = "$0.0005$" if val == 0.0005 else "$0.005$"
                r = B.vrow({"": f}, v, (str(val),), (disp,),
                           dict(scenario=scen, table=kind, setting="+".join(names), **cfg), (scen, si, str(val)))
                r["block_head"] = stext
                rows.append(r)
            if kind == "lmu" and target == "paper":
                whole += rows
            else:
                groups.append(((t["alpha0"],), len(rows) * L, rows))
        if whole:
            groups.append(((t["alpha0"],), len(whole) * L + 2, whole))   # + block row
    allrows = [r for g in groups for r in g[2]]
    B.fill_rows(allrows, lay.specs)
    frag = {"phi": f"phi_all_{scen}_{'L0005' if L_mu == 0.0005 else 'L005'}"}.get(kind, f"{kind}_all_{scen}")
    budget = {("lmu", "paper"): (32, 32, 32), ("lmu", "report"): (36, 36, 36)}.get((kind, target), B.B_MULTI[target])
    pages = B.paginate(groups, budget, True, B.boxes(target), same_sheet=False)
    lmd = None if L_mu is None else ("0.0005" if L_mu == 0.0005 else "0.005")
    bcache = {}

    def base_fn(k):
        a0 = k[0]
        if k not in bcache:
            if kind == "lmu":
                brows = []
                for lm, disp in ((0.0005, "0.0005"), (0.005, "0.005")):
                    brows += B.baseline_rows(D, scen, a0, lm, level=((disp,), (
                        rf"\rotatebox[origin=c]{{90}}{{$L_\mu={disp}$}}" if target == "paper" else f"${disp}$",)))
                bsets = B.settings_line(scen, tag="Baselines", alpha0=a0)
            else:
                brows = B.baseline_rows(D, scen, a0, L_mu if kind == "phi" else 0.0005)
                bsets = B.settings_line(scen, tag="Baselines", alpha0=a0, L_mu=lmd if kind == "phi" else "0.0005")
            B.fill_rows(brows, lay.specs)
            bcache[k] = (brows, bsets)
        return bcache[k]

    fixed = {"delta": dict(b="2.11", phi="1", L_mu="0.0005"), "period": dict(b="2.11", phi="1", L_mu="0.0005"),
             "b": dict(phi="1", L_mu="0.0005"), "phi": dict(b="2.11", L_mu=lmd), "lmu": dict(b="2.11", phi="1")}[kind]
    sets = lambda pg: B.settings_line(scen, **fixed)   # the settings are named in the block rows
    return B.emit(lay, pages, frag, target, None, sets, base_fn,
                  [pcol] if kind == "lmu" else (), box_alone=(kind == "lmu" and target == "report"))


def detailed_events_all_fragment(B, D, scen, target, gid):
    cols = [["qcert@K"], ["qdep@K"], ["qcert@inf"], ["qdep@inf"]]
    lay = B.Layout(gid, [("l", r"$L_\mu$")], cols, 1, ["0.6pt"], method_header="Method",
                   groups=[(B.two_lines(r"$K$ of the setting", r"(level $\delta$)"), 0, 2),
                           (B.two_lines(r"$K=\infty$", r"(level $\delta_1$)"), 2, 2)])
    groups = []
    for si in range(6):
        name, t, _, _ = B.SIX[si]
        bh = B.block_text(name, t, extra=rf"; $\delta_1={B.delta1_pct(t['delta'])}\%$", fields=B.NO_K)   # block row with delta_1, no K
        for lm in (0.0005, 0.005):
            rows = []
            for v in B.VARIANTS:
                fK, cfg = B.anchor_cell(D, si, scen, v, "detailed_200", L_mu=lm)
                fI, _ = B.anchor_cell(D, si, scen, v, "detailed_200", L_mu=lm, K="infinity")
                ld = "$0.0005$" if lm == 0.0005 else "$0.005$"
                rows.append(B.vrow({"K": fK, "inf": fI}, v, (str(lm),), (ld,),
                                   dict(scenario=scen, table="detailed_events", setting=name, **cfg), None))
                rows[-1]["block_head"] = bh
            groups.append((0, 4, rows))
    allrows = [r for g in groups for r in g[2]]
    B.fill_rows(allrows, lay.specs)
    budget = {"paper": (48, 48, 48), "report": (32, 40, 40)}[target]
    paired = [(0, 8, groups[i][2] + groups[i + 1][2]) for i in range(0, len(groups), 2)]
    pages = B.paginate(paired, budget, False, False, same_sheet=False)
    return B.emit(lay, pages, f"detailed_events_all_{scen}", target, None,
                  lambda pg: B.settings_line(scen, b="2.11", phi="1") + rf"\quad $\delta_1=\delta/(2\kappa\ln^2 2)$")


def checks_all_fragment(B, D, scen, L_mu, target, gid):
    lmd = "0.0005" if L_mu == 0.0005 else "0.005"
    frag = f"checks_all_{scen}_{'L0005' if L_mu == 0.0005 else 'L005'}"
    cols = [["qcert"], ["qdep"], ["rCBF"], ["runsafe"], ["qA2"], ["qpw"]]
    lay = B.Layout(gid, [("k", r"$K$")], cols, 1, ["0.6pt"], method_header="Method")   # block row per setting
    groups = []
    for si in range(6):
        name, t, _, _ = B.SIX[si]
        rows = []
        bh = B.block_text(name, t, fields=B.NO_K)   # rows at the setting's K and K = infinity, no K in the block row
        for K in (t["K"], "infinity"):
            for v in B.VARIANTS:
                f, cfg = B.anchor_cell(D, si, scen, v, "detailed_200", L_mu=L_mu, K=K)
                rows.append(B.vrow({"": f}, v, (str(K),), (f"${B.kfmt(K)}$",),
                                   dict(scenario=scen, table="checks", setting=name, **cfg), None))
                rows[-1]["block_head"] = bh
        groups.append(((t["alpha0"],), len(rows), rows))
    allrows = [r for g in groups for r in g[2]]
    B.fill_rows(allrows, lay.specs)
    sets = lambda pg: B.settings_line(scen, b="2.11", phi="1", L_mu=lmd)
    budget = {"paper": (40, 40, 40), "report": (24, 32, 40)}[target]
    pages = B.paginate(groups, budget, True, B.boxes(target), same_sheet=False)
    bc = {}

    def base_fn(k):
        if k not in bc:
            br = B.baseline_rows(D, scen, k[0], L_mu)
            B.fill_rows(br, lay.specs)
            B.apply_submission_ranges(br, scen, k[0])   # the submission's range cell (SUBMISSION_RANGES)
            bc[k] = (br, B.settings_line(scen, tag="Baselines", alpha0=k[0], L_mu=lmd))
        return bc[k]

    return B.emit(lay, pages, frag, target, None, sets, base_fn)


# ----------------------------------------------------------------------------------------------- manifest
SCN = {"corridor": "Corridor", "maze": "Maze"}


def describe(B, name):
    """(type, title, parameters) of a catalog fragment, for the manifest and the catalog documents."""
    p = name.split("_")
    if name.startswith(("sweep_", "compact_")):
        a0 = int(p[2][1:]) / 100
        kind = "sweep_full" if p[0] == "sweep" else "sweep_compact"
        what = "Full sweep" if p[0] == "sweep" else "Compact sweep (best variant)"
        return kind, f"{what}, {SCN[p[1]]}, $\\alpha_0={B.num(a0)}$", dict(scenario=p[1], alpha0=a0)
    if name.startswith("topfive_events"):
        return "events", "Default and selected settings: $\\widehat q_{\\rm cert}$, $\\widehat q_{\\rm dep}$, all variants", {}
    if name.startswith("detailed_events_all_"):
        return "checks", f"Detailed events, all variants, {SCN[p[3]]}", dict(scenario=p[3])
    if name.startswith("checks_all_"):
        lm = "0.0005" if p[3] == "L0005" else "0.005"
        return "checks", f"Checks, all variants, {SCN[p[2]]}, $L_\\mu={lm}$", dict(scenario=p[2], L_mu=float(lm))
    if name.startswith(("delta_all", "period_all", "b_all", "phi_all", "lmu_all")):
        sym = {"delta": "$\\delta$", "period": "$K$", "b": "$b$", "phi": "$\\phi$", "lmu": "$L_\\mu$"}[p[0]]
        extra = ""
        par = dict(scenario=p[2], parameter=p[0])
        if p[0] == "phi":
            lm = "0.0005" if p[3] == "L0005" else "0.005"
            extra = f", $L_\\mu={lm}$"
            par["L_mu"] = float(lm)
        return "parameter", f"Parameter {sym}, all variants, six settings, {SCN[p[2]]}{extra}", par
    if name.startswith("default_"):
        return "default_topfive", f"Default setting, all variants, all metrics, {SCN[p[1]]}", dict(scenario=p[1])
    if name.startswith("topfive_full_"):
        return "default_topfive", f"S1--S5, all variants, all metrics with $\\widehat q_{{\\rm dep}}$, $\\widehat q_{{\\rm cert}}$, {SCN[p[2]]}", dict(scenario=p[2])
    if name.startswith("topfive_"):
        return "default_topfive", f"Default and S1--S5, all variants, all metrics, {SCN[p[1]]}", dict(scenario=p[1])
    raise KeyError(name)


def build(B, D):
    """Returns {target: {name: tex}}; records cells in B.REC (the caller resets it first)."""
    frags = {t: {} for t in B.TARGETS}
    for target in B.TARGETS:
        T = target[0]
        F = frags[target]
        for scen in B.SCENS:
            s = scen[0]
            F[f"default_{scen}"] = B.default_fragment(D, scen, target, f"Qdf{s}{T}")
            F[f"topfive_{scen}"] = B.topfive_fragment(D, scen, target, f"Qtf{s}{T}")
            F[f"topfive_full_{scen}"] = B.topfive_full_fragment(D, scen, target, f"Qtg{s}{T}")
            for kind, tag in (("delta", "dl"), ("period", "pr"), ("b", "bn"), ("lmu", "lm")):
                F[f"{kind}_all_{scen}"] = param_all_fragment(B, D, scen, kind, target, f"Q{tag}{s}{T}")
            for lm, lt in ((0.0005, "L0005"), (0.005, "L005")):
                F[f"phi_all_{scen}_{lt}"] = param_all_fragment(B, D, scen, "phi", target, f"Qph{s}{lt[1:]}{T}", lm)
                F[f"checks_all_{scen}_{lt}"] = checks_all_fragment(B, D, scen, lm, target, f"Qck{s}{lt[1:]}{T}")
            F[f"detailed_events_all_{scen}"] = detailed_events_all_fragment(B, D, scen, target, f"Qde{s}{T}")
            for a0 in B.ALPHA0S:
                F[f"compact_{scen}_a{B.tag3(a0)}"] = compact_fragment(B, D, scen, a0, target, f"Qcp{s}{B.tag3(a0)}{T}")
                F[f"sweep_{scen}_a{B.tag3(a0)}"] = B.sweep_fragment(D, scen, a0, target, f"Qsw{s}{B.tag3(a0)}{T}")
        F["topfive_events"] = B.topfive_events_fragment(D, target, f"Qte{T}")
    return frags


def write(B, D, frags, entries):
    out = B.GEN / "catalog"
    widths = {"paper": "505.89pt \\textwidth, 10pt", "report": "505.89pt \\textwidth, 11pt Times"}
    for target in B.TARGETS:
        d = out / target
        d.mkdir(parents=True, exist_ok=True)
        for name, tex in frags[target].items():
            (d / f"{name}.tex").write_text(B.HEADER.format(target=target, width=widths[target]).replace(
                "do not edit by hand.", "catalog mode (--catalog); do not edit by hand.") + tex)
    # manifest: per fragment the data (files, families, view ids) and the restore line
    byf = {}
    for e in entries:
        k = (e["target"], e["fragment"])
        m = byf.setdefault(k, dict(views=set(), families=set(), srcs=set(), cells=0))
        m["cells"] += 1
        if e.get("view_id"):
            m["views"].add(e["view_id"])
            m["families"].add(e["family"])
            m["srcs"].add(e["src"])
    tables = []
    for target in B.TARGETS:
        for name in frags[target]:
            typ, title, par = describe(B, name)
            m = byf.get((target, name), dict(views=set(), families=set(), srcs=set(), cells=0))
            path = f"appendices/generated_n200/catalog/{target}/{name}.tex"
            tables.append(dict(name=name, target=target, type=typ, title=title, parameters=par, file=path,
                               sha256=hashlib.sha256((B.SRC / path).read_bytes()).hexdigest(),
                               pages=frags[target][name].count(r"\ntablebreak") + 1,
                               data=dict(files=sorted(B.rel_input(B.FILES[s]) for s in m["srcs"]),
                                         families=sorted(m["families"]), n_views=len(m["views"]),
                                         view_ids=sorted(m["views"])),
                               printed_cells=m["cells"],
                               restore=rf"\ntable{{../catalog/{target}/{name}}}"))
    man = dict(title="N200 catalog tables",
               generator="tools/n200_tables/build_tables.py --catalog (catalog_tables.py)",
               checked_by="tools/n200_tables/check_tables.py (catalog_values.json.gz)",
               restore_rule=r"\ntable{../catalog/\nassetkind/<name>} in any edition (n200_assets.tex)",
               compact_rule=COMPACT_KEY, inputs={k: B.rel_input(v) for k, v in B.FILES.items()}, sha256=D.sha,
               tables=tables)
    (out / "CATALOG_TABLES.json").write_text(json.dumps(man, indent=1))
    values = dict(generator="tools/n200_tables/build_tables.py --catalog", inputs={k: B.rel_input(v) for k, v in B.FILES.items()},
                  sha256=D.sha, fragments={t: sorted(frags[t]) for t in B.TARGETS}, cells=entries)
    B.write_gz(B.HERE / "catalog_values.json.gz", json.dumps(values, default=str, separators=(",", ":")))
    print(f"catalog: {len(entries)} printed cells recorded; fragments: paper {len(frags['paper'])}, "
          f"report {len(frags['report'])}; manifest {out / 'CATALOG_TABLES.json'}")
