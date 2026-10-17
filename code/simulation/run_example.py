"""Run one closed-loop trajectory and print its main metrics.

Example:
    python3 run_example.py --scenario corridor --variant V-C --alpha0 0.15 --alpha 0.10 \
        --delta 0.05 --K 50 --b 2.11 --phi 1.0 --L-mu 0.0005 --seed 970001
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ecbf.assessment import assess_record
from ecbf.core import make_method
from ecbf.simulation import simulate

# Report names of the certified variants and their labels in the package.
VARIANTS = {"V-U": "V-U", "V-Delta": "V-D", "V-L": "V-Cd", "V-C": "V-C"}
BASELINES = ("Unfiltered", "Zero margin", "Conservative", "ACI", "QT")
# Risk assessment used for every trajectory of the study.
ASSESSMENT = dict(mode="development", max_bins=128, risk_tol=1e-5)


def parse_K(text):
    return "infinity" if text in ("infinity", "inf") else int(text)


def make_job(args):
    """The job fields of the study (one row of the trajectory plan)."""
    method = VARIANTS.get(args.variant, args.variant)
    if method not in VARIANTS.values() and method not in BASELINES:
        raise ValueError("unknown variant: " + args.variant)
    job = dict(benchmark=args.scenario, method=method, alpha=args.alpha, alpha0=args.alpha0, delta=args.delta,
               K=parse_K(args.K), bound_scale=args.phi, budget_normalizer=args.b, rho_buf=0.0, L_mu=args.L_mu,
               seed=args.seed)
    if method in BASELINES:
        # The baselines do not use the certification settings.
        job.update(alpha=None, delta=None, K=None, bound_scale=None, budget_normalizer=None)
    return job


def run_job(job):
    """Simulate one trajectory as in the study and return the raw record."""
    common = dict(alpha0=job["alpha0"], model_over={"rho_buf": job["rho_buf"]}, beta_mode="sin",
                  L_mu=job["L_mu"], record_risk=True)
    if job["method"] in BASELINES:
        return simulate(job["benchmark"], make_method(job["method"]), job["seed"], **common)
    once = job["K"] == "infinity"
    # Label: variant, test schedule (one test, or a new test every K steps) and bound factor phi.
    label = job["method"] + ("+once" if once else "+p" + str(job["K"]))
    if job["bound_scale"] != 1.0:
        label += "+f" + str(job["bound_scale"])
    return simulate(job["benchmark"], make_method(label), job["seed"], alpha=job["alpha"], delta=job["delta"],
                    kappa=1.0, budget_mode="joint_log_once" if once else "joint_log",
                    budget_normalizer=job["budget_normalizer"], **common)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--scenario", choices=("corridor", "maze"), default="corridor")
    parser.add_argument("--variant", default="V-C", help="V-U, V-Delta, V-L, V-C or a baseline: " + ", ".join(BASELINES))
    parser.add_argument("--alpha0", type=float, default=0.15, help="risk level alpha_0")
    parser.add_argument("--alpha", type=float, default=0.10, help="test level alpha")
    parser.add_argument("--delta", type=float, default=0.05, help="error level delta")
    parser.add_argument("--K", default="50", help="test period K, or infinity for one test")
    parser.add_argument("--b", type=float, default=2.11, help="budget normalizer b (at least 2.11)")
    parser.add_argument("--phi", type=float, default=1.0, help="bound factor phi in [0, 1]")
    parser.add_argument("--L-mu", dest="L_mu", type=float, default=0.0005, help="disturbance drift L_mu")
    parser.add_argument("--seed", type=int, default=970001)
    parser.add_argument("--no-assessment", action="store_true", help="skip the evaluation of risk events")
    args = parser.parse_args()

    job = make_job(args)
    record = run_job(job)
    print("job:", json.dumps(job))
    print("termination:", record["termination"], "| decisions N =", int(record["N"]))
    if args.no_assessment:
        return
    out = assess_record(record, **ASSESSMENT)
    tau, mean_margin = out["tau_use"], out["m_c"]
    print(f"mission time T_mis = {out['T_mis']:.2f} s (goal reached: {out['goal']})")
    print("mean margin c_t/c_max =", "not defined" if mean_margin is None else f"{mean_margin:.4f}")
    print("first decision with c_t < c_max, tau_use =", "none" if tau is None or not math.isfinite(tau) else f"{tau:.2f} s")
    print(f"certificates issued = {out['n_certificates']}")
    print(f"minimum clearance = {out['clearance_min']:.4f} m | contact: {out['contact']} | infeasible: {out['infeas']}")
    if out["qdep_applicable"]:
        print(f"deployment event p_t(c_t) > alpha_0: confirmed {out['qdep_confirmed']}, possible {out['qdep_possible']}")
    if out["qcert_applicable"]:
        print(f"false certificate: confirmed {out['qcert_confirmed']}, possible {out['qcert_possible']}")


if __name__ == "__main__":
    main()
