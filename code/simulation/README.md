# Simulation code

This folder holds the simulation package `ecbf` that generated every trajectory of the ablation study, and a short
script `run_example.py` that runs one trajectory with the same settings. Relative to the package that produced the
data, this copy differs only in comments, docstrings, error messages, the caption notes of `comparisons.py`, the
unused `source` labels of the two scenes, the descriptive strings and two unused blocks (`robot`, `variants`) of the
Corridor scene file, three method names of the scene class (`obstacle_clearances`, `obstacle_kind`, `past_offset`),
and unused helper functions and one unused import that are left out; the computations are unchanged. Section 1 of
the ablation report (`pdf/ablation_report.pdf`) describes the implementation details that this code follows.

References in this folder: [1] is the arXiv paper, so [1, (57)] is its eq. (57) and [1, App. D-B] its Appendix D-B;
"report" is the ablation report, cited by its equation, section and table numbers.

## Modules

| File | Content |
| --- | --- |
| `ecbf/models.py` | Corridor and Maze dynamics (unicycle with speed and yaw-rate input in the Corridor, acceleration input in the Maze), the nominal inputs (goal tracking in the Corridor, pure pursuit on the reference path in the Maze), the CBF safety filter that returns the input $u_t(c)$ for every candidate margin and its feasibility (report Section 1.3.1), and the Jacobian bounds $K^z_x$, $K^z_u$ of the predicted nose point ([1, (63)], [1, (64)]). |
| `ecbf/geometry.py` | Scene geometry: obstacle clearances $\psi_i$ and $\psi = \min_i \psi_i$, signed distance to rectangles, the geometric check at the predicted nose point (including the rectangle condition, report eq. (21)), the set of effective radii, the exceedance fractions $P_R$ and $P_\infty$, the fraction bound $\bar P(c)$ and the sensitivity $L^{\mathrm{geom}}(c)$ ([1, (53)], [1, (54)], [1, (57)]). |
| `ecbf/certification/` | E-value tests and certification. `evidence.py`: mixture e-process over the gain set $\Lambda$ (report eq. (3)). `budgets.py`: level $\delta_r = \delta / [b (r+1) \ln^2(r+1)]$ with budget normalizer $b \ge 2.11$ (report eq. (4)). `policies.py`: test starts (one test for $K = \infty$, a new test every $K$ steps otherwise; the other policies are not used in the study). `monitor.py`: parallel tests, certificates and their monitoring periods $[\sigma_r, \tau_r(c))$. |
| `ecbf/core.py` | Certification monitor, deployment rule and supporting test index (`select_certified`, report eqs. (5) and (6)), the baselines ACI and QT ([1, (89)], [1, (90)]) and the method labels. |
| `ecbf/simulation.py` | Closed loop: disturbance distribution (report eq. (10)), the risk-increase terms of the variants (`Transfer.matrices`), and the record of one trajectory. |
| `ecbf/risk.py` | Estimated exceedance risk $\hat p_t(c)$ (report eq. (19)) by exact arc lengths in angle (report eq. (15)) and Gauss-Legendre quadrature in magnitude, with each magnitude interval split at the entry magnitudes of the obstacles (report eq. (14)). |
| `ecbf/risk_intervals.py` | Evaluation only: numerical risk bounds $p^-_t(c) \le p_t(c) \le p^+_t(c)$ (report Section 1.1.4), refined by bisection of the magnitude interval with the largest contribution to the width until a threshold is separated. No controller module uses it. |
| `ecbf/assessment.py` | Assessment of one trajectory: deployment events $p_t(c_t) > \alpha_0$ and false certificates, each as confirmed or possible (report Section 1.2.2), and the assumption checks (report Section 1.4). |
| `ecbf/metrics.py` | Per-trajectory metrics (report Section 1.2): mission time $T_{\mathrm{mis}}$, mean margin $c_t / c_{\max}$, goal, contact, infeasibility, first use of a margin below $c_{\max}$, clearance. |
| `ecbf/comparisons.py` | Plots of fixed-setting distributions and event frequencies from the saved trajectory records (`.npz`) of the study, which are not part of this repository (not needed to run trajectories). |
| `ecbf/*.json`, `ecbf/*.npy` | Corridor scene, Maze scene and Maze reference path. |

## Certificate bounds

`Transfer.matrices` in `simulation.py` returns, at decision $t$, the risk-increase term of the variant for every
candidate $c$ and every past step $s < t$, capped at 1. The certificate bound $d^V_t(c; \sigma, \tau)$ is its maximum
over $\sigma \le s < \tau$ (report eq. (8) and the table of variants above it; [1, (69)]-[1, (72)]). The terms use
the drift term $\Delta_{t,s} = \min\{1, L_\mu (t-s)\}$ (`eps`), the geometric correction $\eta_{t,s}(c)$ (`kap`) and
the uniform geometric correction $\bar\eta_{t,s} = \max_{c'} \eta_{t,s}(c')$, the fraction bound $\bar P(c)$ (`p0`),
the clearance term $G_{t,s}(c) = L^{\mathrm{geom}}(c) [g_t(c) - g_s]_+$ with the nose-point clearances
$g_t(c) = \psi(z_t(c))$ and $g_s = \psi(z_s)$, and the constants $L^{\mathrm{geom}}_x$, $L^{\mathrm{geom}}_u$
(report eq. (45)).

- V-U (`uniform`, label `V-U`): $L^{\mathrm{geom}}_x \|x_t - x_s\| + L^{\mathrm{geom}}_u \|u_t(c) - u_s\| +
  \Delta^{\mathrm{geom}}_t$ with $\Delta^{\mathrm{geom}}_t = \bar\eta_{t,s} + \Delta_{t,s}$; uniform check.
- V-Delta (`candresid`, label `V-D`): the same constants plus $\Delta^{\mathrm{geom}}_t(c) = \eta_{t,s}(c) +
  \Delta_{t,s}$; check for $c$.
- V-L (`commondelta`, label `V-Cd`): $\min\{\bar P(c), G_{t,s}(c)\} + \Delta^{\mathrm{geom}}_t$; uniform check.
- V-C (`candidate`, label `V-C`): $(1 - \Delta_{t,s}) \min\{\bar P(c), G_{t,s}(c) + \eta_{t,s}(c)\} + \Delta_{t,s}
  \bar P(c)$; check for $c$.

A failed required check sets the term to 1. For V-U and V-L a failed check of any candidate sets every term to 1,
and a past step with $g_s < g_{\min} = l$ sets its term to 1 for every variant. The bound factor $\phi$ scales
every term; the guarantees hold only for $\phi = 1$. `select_certified` in `core.py` applies the deployment rule
(report eq. (5)): it deploys the smallest candidate with a certificate whose bound $d^{(r)}_t(c)$, $\phi$ times the
largest term over the monitoring period $[\sigma_r, \tau_r(c))$, satisfies $\alpha + d^{(r)}_t(c) \le \alpha_0$,
and the margin $c_{\max}$ when no candidate qualifies. It also returns the supporting test index $r_t$ (report
eq. (6)).

## Requirements

Python 3.12 with NumPy and SciPy. Verified with Python 3.12.6, NumPy 1.26.4 and SciPy 1.17.0. `comparisons.py` also
needs Matplotlib (verified with 3.10.8) and a LaTeX installation, since it sets all figure text with `text.usetex`.
No GPU is used.

## Example

The script runs one trajectory from this folder and prints the mission time, the mean margin, the number of
certificates, the minimum clearance and the events.

```sh
python3 run_example.py --scenario corridor --variant V-C --alpha0 0.15 --alpha 0.10 \
    --delta 0.05 --K 50 --b 2.11 --phi 1.0 --L-mu 0.0005 --seed 970001
```

The arguments are the job fields of the study: scenario (`corridor`, `maze`), variant (`V-U`, `V-Delta`, `V-L`,
`V-C`, or a baseline `Unfiltered`, `Zero margin`, `Conservative`, `ACI`, `QT`), CBF-violation probability
$\alpha_0$, margin insufficiency $\alpha$, error confidence $\delta$, test period $K$ (an integer, or `infinity` for
one test), budget normalizer $b$, bound factor $\phi$, distribution-drift bound $L_\mu$ and seed. The baselines use
only $\alpha_0$, $L_\mu$ and the seed. As in the study, the script runs a variant with its package label (`V-U`,
`V-D`, `V-Cd`, `V-C`) followed by `+pK` for finite $K$ (budget mode `joint_log`, a new test every $K$ steps) or
`+once` for $K = \infty$ (budget mode `joint_log_once`, one test), and by `+f<phi>` when $\phi \ne 1$, for example
`V-C+p50`. With `--no-assessment` the risk events are not evaluated.

The command above reproduces the default setting for seed 970001 in
`data/per_trajectory/per_trajectory_values.csv` (mission time 35.6 s, 356 decisions, mean margin 0.7697, first
margin below $c_{\max}$ at 4.9 s, no event). The printed minimum clearance is the obstacle distance in metres, not
the column `h_min` of the per-trajectory data. A trajectory that does not reach the goal counts the full horizon in
$T_{\mathrm{mis}}$.

## Seeds

Each setting uses the common seeds 970001 to 970200, so $N_{\mathrm{MC}} = 200$ trajectories per setting. The seed
fixes both the initial state and the disturbance sequence, so every variant and baseline of a setting sees the same
disturbances.
