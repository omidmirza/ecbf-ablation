# Table fragments: key and caption conventions

`build_tables.py` writes one LaTeX fragment per table to `appendices/generated_n200/{paper,report}/` of the paper
source tree; `check_tables.py` re-derives every printed value and checks the layout rules below.

## Table key

The symbols are defined once, in the boxed "Key to all tables" (`table_key`, report Section 2.1; `table_key_arxiv`,
a shorter key with event frequencies in percent). The Section VI tables (`sec6_short`, `sec6_combined`, `sec6_maze_arxiv`) carry a two-line key inside their
own box. No other fragment contains a key box; a caption adds at most one short note.

## Caption notes

| Fragment | Target | Caption note |
|---|---|---|
| `table_key`, `table_key_arxiv` | both / paper | none (boxed key only) |
| `sec6_short`, `sec6_maze_arxiv` | paper | none (two-line key inside); baselines at $\alpha_0=0.30$ |
| `sec6_combined` | paper | none (two-line key inside); Unfiltered, $c=0$, $c_{\max}$ once, ACI and QT at $\alpha_0=0.15$ and $0.30$ |
| `design_overview`, `arxiv_design_overview` | both / paper | none (comparison, varied values, fixed values) |
| `default_{scen}` | both | Default setting, four variants; baselines at $\alpha_0=0.15$. |
| `topfive_full_{scen}` | both | S1--S5, four variants; baselines at $\alpha_0=0.30$. |
| `delta_{scen}` | both | Paired variant; $\delta$ swept inside the tuple. |
| `period_{scen}` | both | Paired variant; $K$ swept inside the tuple. |
| `b_{scen}` | both | Paired variant; $b$ swept. |
| `phi_{scen}` | both | Paired variant; $\phi\times L_\mu$; baselines at both $L_\mu$; $\phi<1$: Assumption 2 may fail. |
| `detailed_events_{scen}` | both | Paired variant; $K=\infty$ columns: \checkmark against $\delta_1$. |
| `checks_{scen}_{L0005,L005}` | both | Paired variant at $K$ of the setting and $K=\infty$ (\checkmark against $\delta_1$). |
| `checks_{scen}_arxiv` | paper | Paired variant, S1--S5, both $L_\mu$, $K$ of the setting and $K=\infty$ (\checkmark against $\delta_1$); columns that are $0$ in every row are omitted. |
| `stress_{scen}_arxiv`, `stress_{scen}_cases` | paper / report | Paired variant, S1--S5, $\phi\times L_\mu$; $^\dagger$: bound scaled by $\phi<1$, Assumption 2 may fail. |
| `variation_{scen}` | both | $10^3\widehat V$, $100\widehat J$ over the first $N_{\rm var}=200$ decisions; last column: sizes of $\mathcal K_\upsilon$. |
| `n30_record_{scen}` | both | Settings selected on the first 30 trajectories; each row changes one of $b$, $\phi$, $L_\mu$ from the base. |

The paired variant of a setting is the variant reported for it: V-C at the default setting, S1, S2 and S5,
V-L at S3 and V-Delta at S4.

## Omitted columns

The arXiv checks and stress fragments omit every column that is 0 or N-A in every row. The omitted columns are listed
in the first line of the fragment after its header (`% OMITTED (0 or N-A in every row): ...`), so a caption can name
them.

## Layout

- Setting tuples are never split over lines. Parameter tables (delta, period, b, phi), `topfive_full` and
  `n30_record` print a full-width one-line block row per setting, for example
  "S3: V-L $(\alpha_0,\alpha,\delta,K) = (0.30, 0.20, 0.20, 50)$" (plus ", $L_\mu = \dots$" where $L_\mu$ varies);
  the rows give the swept value in their own column. The block row omits the swept parameter.
- Every table has a settings row, one header, block rows, centred names and values at the full line width, value
  columns of equal width when they fit, and settings and block rows padded 4 pt from the box rule.
- Baseline boxes are centred at the full line width. They use the grid of their table when no value column is
  dropped, and their own balanced grid otherwise.
- Rows of one line are stretched (`\arraystretch` between 0.93 and 1.4) to fill about 90% of the page; multi-page
  tables break with `\ntablebreak`.
- Appendix fragments are set at the body size; the Section VI tables keep the size of their `table*` float.

## Fragments

Both targets: `table_key`, `design_overview`, `default_{scen}`, `topfive_full_{scen}`, `delta_{scen}`,
`period_{scen}`, `b_{scen}`, `phi_{scen}`, `detailed_events_{scen}`, `checks_{scen}_{L0005,L005}`,
`variation_{scen}`, `n30_record_{scen}`, `star_maze`.

Paper only: `table_key_arxiv`, `arxiv_design_overview`, `sec6_short`, `sec6_combined`, `sec6_maze_arxiv`,
`main_maze_arxiv_v8` (arXiv Table I) and `main_corridor_arxiv_v8` (arXiv Table IV), `main_{maze,corridor}_extended_v8`,
`secondary_{scen}_arxiv`, `topfive_full_{scen}_arxiv_percent`, `topfive_full_{scen}_onepage`,
`variation_{scen}_arxiv`, `checks_{scen}_L0005_nodef`, `checks_{scen}_arxiv`, `stress_{scen}_arxiv`,
`stress_{scen}_arxiv_percent`, `stress_representative_arxiv_percent` (arXiv Table V) and `assumption_checks_arxiv` (a
file of `\providecommand` macros with the assumption-check totals and the main-grid event counts of arXiv App. F-D,
no table).

Report only: `stress_{scen}_cases` (the cells of `stress_{scen}_arxiv`).

Printed: in the arXiv paper `main_maze_arxiv_v8`, `main_corridor_arxiv_v8`, `stress_representative_arxiv_percent` and
`assumption_checks_arxiv`; in the report (report target) `table_key` (Section 2.1) and Tables 3-24: `design_overview`,
`default_{scen}`, `topfive_full_{scen}`, `star_maze`, `delta_{scen}`, `period_{scen}`, `b_{scen}`, `phi_{scen}`,
`stress_{scen}_cases`, `checks_{scen}_{L0005,L005}`, `variation_{scen}`. The paper source also includes the paper
fragments `table_key`, `table_key_arxiv`, `design_overview`, `default_{scen}`, `topfive_full_{scen}`,
`topfive_full_{scen}_onepage`, `star_maze`, `phi_{scen}`, `checks_{scen}_{L0005,L005}`, `variation_{scen}` and
`main_{maze,corridor}_extended_v8` (extended version). The other fragments are generated and checked but not printed.

`{scen}` is `maze` or `corridor`. S1--S5 are the selected settings (asset labels New1--New5 and T1--T5 in the
data files).
