# Simulation code

This folder holds the simulation package `ecbf` that generated every trajectory of the ablation study, and a short
script `run_example.py` that runs one trajectory with the same settings. The package is the version used for the
study. For publication only comments, docstrings, two error messages and two descriptive strings of the Corridor scene
file were reworded; the computations are unchanged. Section 1 of the ablation report (`pdf/ablation_report.pdf`)
describes the implementation details that this code follows.

## Modules

| File | Content |
| --- | --- |
| `ecbf/models.py` | Corridor and Maze dynamics (unicycle, first-order and second-order input), the nominal controllers (goal tracking in the Corridor, pure pursuit on the reference path in the Maze), the CBF safety filter that returns the input $u_t(c)$ for every candidate margin and its feasibility, and the Jacobian bounds $K_x$, $K_u$ of the one-step prediction used by the uniform bounds. |
| `ecbf/geometry.py` | Scene geometry: obstacle clearance $\psi(z)$, signed distance to rectangles, the candidate conditions at the current point (including the rectangle condition), the set of offsets, and the candidate-dependent terms $p_0(c)$, $F$, $L(c)$ and the curvature term. |
| `ecbf/certification/` | E-value tests and certification. `evidence.py`: mixture e-process over the gain set $\Lambda$. `budgets.py`: allocated level $\delta_r = \delta / [b (r+1) \log^2(r+1)]$ with budget normalizer $b \ge 2.11$. `policies.py`: test starts (one test for $K = \infty$, a new test every $K$ steps otherwise). `monitor.py`: parallel tests, certificates and their monitoring periods $[\sigma_r(c), \tau_r(c))$. |
| `ecbf/core.py` | Certification monitor, selection of the smallest certified margin with $\alpha + d_t(c) \le \alpha_0$, where $d_t(c)$ is the maximum of the pairwise bound over the monitoring period, the baselines ACI and QT, and the method labels. |
| `ecbf/simulation.py` | Closed loop: disturbance distribution (two-component mixture with weight $(1+\sin(2 L_\mu t))/2$), the pairwise certificate bounds of the variants (`Transfer.matrices`), and the record of one trajectory. |
| `ecbf/risk.py` | Conditional risk $p_t(c)$ by exact arc lengths in angle and Gauss-Legendre quadrature in radius, with each radial interval split at the radii where an obstacle enters. |
| `ecbf/risk_intervals.py` | Evaluation only: lower and upper enclosures of $p_t(c)$ by radial bins, refined by bisection of the bin with the largest uncertainty until a threshold is separated. No controller module uses it. |
| `ecbf/assessment.py` | Assessment of one trajectory: deployment events $p_t(c_t) > \alpha_0$ and false certificates, each as confirmed or possible, and the checks of the assumptions. |
| `ecbf/metrics.py` | Per-trajectory metrics: mission time $T_{\mathrm{mis}}$, mean margin $c_t / c_{\max}$, goal, contact, infeasibility, first use of a margin below $c_{\max}$, clearance. |
| `ecbf/comparisons.py` | Plots of fixed-setting distributions and event frequencies (not needed to run trajectories). |
| `ecbf/*.json`, `ecbf/*.npy` | Corridor scene, Maze scene and Maze reference path. |

## Certificate bounds

`Transfer.matrices` in `simulation.py` returns the pairwise bound $D_{t,s}(c)$ for every candidate $c$ and every past
step $s < t$. Each variant uses the drift term $\varepsilon_{t,s} = \min\{1, L_\mu (t-s)\}$ for each $s$.

- V-U (`uniform`): uniform constants $L_x$, $L_u$ times the state and input distances, plus a common term
  $\Delta_{t,s}$.
- V-Delta (`candresid`, label `V-D`): uniform $L_x$, $L_u$, plus a candidate term $\Delta_{t,s}(c)$.
- V-L (`commondelta`, label `V-Cd`): candidate sensitivity $L(c)$ times the clearance gap, plus the common term
  $\Delta_{t,s}$.
- V-C (`candidate`): candidate sensitivity and candidate curvature term, mixed with $p_0(c)$ by the drift term.

A failed candidate condition sets the entry to 1. The bound factor $\phi$ scales the whole bound ($\phi = 1$ in the
valid setting). `select_certified` in `core.py` takes the maximum of $D_{t,s}(c)$ over the monitoring period of each
certificate and deploys the smallest candidate whose bound satisfies $\alpha + d_t(c) \le \alpha_0$. The margin
$c_{\max}$ is used when no candidate qualifies.

## Requirements

Python 3.12 with NumPy and SciPy. Verified with Python 3.12.6, NumPy 1.26.4 and SciPy 1.17.0. `comparisons.py` also
needs Matplotlib (verified with 3.10.8). No GPU is used.

## Example

The script runs one trajectory from this folder and prints the mission time, the mean margin, the number of
certificates, the minimum clearance and the events.

```sh
python3 run_example.py --scenario corridor --variant V-C --alpha0 0.15 --alpha 0.10 \
    --delta 0.05 --K 50 --b 2.11 --phi 1.0 --L-mu 0.0005 --seed 970001
```

The arguments are the job fields of the study: scenario (`corridor`, `maze`), variant (`V-U`, `V-Delta`, `V-L`,
`V-C`, or a baseline `Unfiltered`, `Zero margin`, `Conservative`, `ACI`, `QT`), $\alpha_0$, $\alpha$, $\delta$, $K$
(an integer, or `infinity` for one test), budget normalizer $b$, bound factor $\phi$, disturbance drift $L_\mu$ and
seed. The baselines use only $\alpha_0$, $L_\mu$ and the seed. With `--no-assessment` the risk events are not
evaluated. The command above reproduces the default setting for seed 970001 in `data/per_trajectory/per_trajectory_values.csv`
(mission time 35.6 s, 356 decisions, mean margin 0.7697, first margin below $c_{\max}$ at 4.9 s, no event). The printed
minimum clearance is the obstacle distance in metres, not the column `h_min` of the per-trajectory data. A trajectory that does not reach the goal counts the full horizon in
$T_{\mathrm{mis}}$.

## Seeds

Each setting uses the common seeds 970001 to 970200, so $N_{\mathrm{MC}} = 200$ trajectories per setting. The seed
fixes both the initial state and the disturbance sequence, so every variant and baseline of a setting sees the same
disturbances.
