"""Confirmed and possible events of one trajectory (report Section 1.2.2) and the assumption checks (report
Section 1.4); evaluation only, never used by the controller.

Confirmed means the numerical risk bounds, with their floating-point
allowances, establish the event. Possible includes every unresolved event.
These labels are not a claim of machine-verified interval arithmetic. A
failed/truncated trajectory retains prior confirmed events and contributes a
possible failure for its unobserved tail. Goal termination defines an absorbing
end of the mission with no further decisions; horizon termination is complete.
"""
from __future__ import annotations
import json
import math
import numpy as np
from .metrics import trajectory
from .models import make_model
from .risk_intervals import risk_interval
from .simulation import atom_level, beta_at


def _plain(value):
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (tuple, list, np.ndarray)):
        return [_plain(v) for v in value]
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (float, np.floating)):
        return float(value) if math.isfinite(float(value)) else None
    return value


def assess_record(rec, *, mode='development', model_over=None, max_bins=None,
                  risk_tol=1e-5):
    """Return scalar metrics and confirmed/possible trajectory-level events.

    mode='development' caps each radial calculation at 128 bins; 'refined' uses
    2048. Work-limit intervals remain possible failures, never false zeros.
    Model overrides must reproduce the generating simulator exactly. Only the
    two unicycle scenarios are supported. No input record is mutated.
    """
    if mode not in ('development', 'refined'):
        raise ValueError('mode must be development or refined')
    max_bins = (128 if mode == 'development' else 2048) if max_bins is None else int(max_bins)
    saved_over = rec.get('model_over', {})
    if isinstance(saved_over, str):
        saved_over = json.loads(saved_over)
    overrides = dict(saved_over or {})
    # rho_buf of the generating run (0 in the study, so h uses l alone); 0.05 is the model default
    overrides.setdefault('rho_buf', float(rec.get('rho_buf', .05)))
    if model_over is not None:
        for key, value in model_over.items():
            if key in overrides and overrides[key] != value:
                raise ValueError('assessment model overrides disagree with saved metadata')
            overrides[key] = value
    model = make_model(str(rec['benchmark']), **overrides)
    cmax, alpha, alpha0 = float(rec['c_max']), float(rec['alpha']), float(rec['alpha0'])
    if model.p.c_max != cmax:
        raise ValueError('assessment model support disagrees with stored c_max')
    N, H = int(rec['N']), int(rec['H'])
    x, u, margins = np.asarray(rec['x']), np.asarray(rec['u']), np.asarray(rec['c'])
    if len(x) != N+1 or len(u) != N or len(margins) != N:
        raise ValueError('incomplete decision record')
    grid = np.asarray(rec['grid'], float)
    if not np.array_equal(grid, model.grid):
        raise ValueError('assessment candidate grid disagrees with generating model')
    law = str(rec.get('beta_mode', 'sin'))
    if law not in ('sin', 'cos', 'const1') and atom_level(law) is None:
        raise ValueError('unrecognized disturbance law')
    lmu = float(rec.get('L_mu', 5e-4))
    atom = atom_level(law)
    z_cache, cache = {}, {}
    evaluations = interval_calls = 0
    numerical_errors = []
    def at(t, c, threshold):
        nonlocal evaluations, interval_calls
        key = (int(t), float(c))
        previous = cache.get(key)
        if previous is not None and (previous[1] < threshold or previous[0] > threshold
                                     or previous[2] == threshold):
            return previous[:2]
        try:
            if t not in z_cache:
                z_cache[t] = model.z(x[t], u[t])
            iv = risk_interval(model.scene, z_cache[t], c, beta_at(t, law, lmu), cmax,
                               threshold=threshold, max_bins=max_bins, tolerance=risk_tol,
                               atom_radius=None if atom is None else atom*cmax)
            interval_calls += 1
            evaluations += iv.evaluations
            lower, upper = iv.lower, iv.upper
        except (ValueError, ArithmeticError, FloatingPointError) as exc:
            lower, upper = 0.0, 1.0
            numerical_errors.append({'t': t, 'c': c, 'error': str(exc)})
        if previous is not None:
            lower, upper = max(lower, previous[0]), min(upper, previous[1])
            if lower > upper:
                lower, upper = 0.0, 1.0
                numerical_errors.append({'t': t, 'c': c, 'error': 'inconsistent enclosures'})
        cache[key] = (lower, upper, threshold)
        return lower, upper
    term = str(rec['termination'])
    missing_tail = term not in ('goal', 'horizon') or (term == 'horizon' and N < H)
    qdep_applicable = bool(np.all(np.isfinite(margins))) and not str(rec['method']).lower().startswith('unfiltered')
    certs = np.asarray(rec.get('certs', []), int).reshape(-1,4)
    qcert_applicable = str(rec['method']).startswith('V') or len(certs) > 0
    retained = []
    dep_confirmed = dep_possible = False
    dep_unresolved = dep_fail_n = 0
    # A quadrature value only sets assessment order; every negative decision is
    # excluded by an enclosure or the exact support shortcut.
    estimates = np.asarray(rec.get('p_dep', np.zeros(N)), float)
    order = np.argsort(np.nan_to_num(estimates, nan=1.0))[::-1]
    if qdep_applicable:
        for t in order:
            lo, hi = at(int(t), float(margins[t]), alpha0)
            confirmed, possible = lo > alpha0, hi > alpha0
            dep_confirmed |= confirmed
            dep_possible |= possible
            dep_fail_n += int(confirmed)
            dep_unresolved += int(possible and not confirmed)
            if possible:
                retained.append({'event': 'deployment', 't': int(t), 'c': float(margins[t]),
                                 'risk_lower': lo, 'risk_upper': hi, 'threshold': alpha0,
                                 'confirmed': confirmed, 'possible': possible})
    fc_confirmed = fc_possible = False
    fc_count = fc_possible_count = fc_unresolved = 0
    ptable = np.asarray(rec.get('p_applied', np.full((N,len(grid)), np.nan)), float)
    for j, r, s, e in certs:
        confirmed, possible, witness, checked = True, True, None, 0
        invalid = not (0 <= j < len(grid) and 0 <= s < e <= N)
        if invalid:
            confirmed = False
        else:
            # One verified upper<alpha is a sufficient non-false witness. It is
            # the earliest cheap route, but no quadrature value itself is evidence.
            times = np.arange(s,e)
            times = times[np.argsort(np.nan_to_num(ptable[s:e,j], nan=1.0))]
            for t in times:
                lo, hi = at(int(t), float(grid[j]), alpha)
                checked += 1
                if hi < alpha:
                    confirmed = possible = False
                    witness = int(t)
                    break
                if lo < alpha:
                    confirmed = False
        fc_confirmed |= confirmed
        fc_possible |= possible
        fc_count += int(confirmed)
        fc_possible_count += int(possible)
        fc_unresolved += int(possible and not confirmed)
        if possible:
            retained.append({'event': 'false_certificate', 'j': int(j), 'r': int(r),
                             'start': int(s), 'end': int(e), 'threshold': alpha,
                             'confirmed': confirmed, 'possible': possible,
                             'checked_risks': checked, 'invalid_window': invalid})
    feasible = np.asarray(rec.get('feasible_cmax', []), bool)
    loss = np.asarray(rec['Y'], float)
    a2 = np.asarray(rec.get('a2_excess', []), float)
    residual = np.asarray(rec.get('residual', []), float)
    # Baselines stop when their own requested margin is infeasible; V methods
    # and Conservative stop when the c_max fallback is infeasible. A generic
    # 'infeasible' label alone must not turn a baseline failure into evidence against Assumption 1.
    terminal_fallback_feasible = None
    terminal_filter_basis = 'not-an-infeasible-termination'
    own_margin_method = str(rec['method']) in ('ACI', 'QT', 'Zero margin')
    if term == 'infeasible':
        failed_margins = np.asarray(rec.get('failed_filter_margins', []), float).reshape(-1)
        failed_flags = np.asarray(rec.get('failed_filter_feasible', []), bool).reshape(-1)
        cmax_rows = np.flatnonzero(failed_margins == cmax)
        if len(failed_margins) == len(failed_flags) and cmax_rows.size:
            flags = failed_flags[cmax_rows]
            if np.all(flags == flags[0]):
                terminal_fallback_feasible = bool(flags[0])
                terminal_filter_basis = 'recorded-cmax-feasibility'
            else:
                terminal_filter_basis = 'inconsistent-recorded-cmax-flags'
        elif failed_margins.size or failed_flags.size:
            terminal_filter_basis = 'incomplete-recorded-terminal-flags'
        elif str(rec['method']).startswith('V') or str(rec['method']) == 'Conservative':
            # In this simulator these methods use own_ok=feas_cmax.
            terminal_fallback_feasible = False
            terminal_filter_basis = 'legacy-fallback-stopping-rule'
        else:
            terminal_filter_basis = 'legacy-terminal-fallback-unknown'
    executed_fallback_failure = bool(np.any(~feasible))
    terminal_fallback_failure = terminal_fallback_feasible is False
    assumption_status = {
        'global_fallback_feasibility_proved': False,
        'fallback_infeasibility_observed': executed_fallback_failure or terminal_fallback_failure,
        'executed_fallback_infeasibility_observed': executed_fallback_failure,
        'terminal_fallback_infeasibility_observed': terminal_fallback_failure,
        'terminal_fallback_feasible': terminal_fallback_feasible,
        'terminal_fallback_feasibility_unknown': term == 'infeasible' and terminal_fallback_feasible is None,
        'terminal_own_margin_infeasibility_observed': term == 'infeasible' and own_margin_method,
        'terminal_filter_assessment_basis': terminal_filter_basis,
        'support_excess_observed': bool(np.any(loss > cmax+1e-12)),
        'filter_residual_failure_observed': bool(np.any(residual < -1e-10)),
        'a2_quadrature_excess_observed': bool(np.any(a2[np.isfinite(a2)] > 1e-12)),
        'a2_quadrature_check_is_proof': False,
        'bound_multiplier_valid': float(rec.get('kappa', 1.0)) >= 1 and float(rec.get('fraction', 1.0)) >= 1,
        'risk_enclosure_machine_certified': False,
        'missing_tail': missing_tail,
        'termination': term,
        'goal_absorbing_no_further_decisions': term == 'goal',
        'numerical_error_count': len(numerical_errors),
    }
    out = trajectory(rec)
    out.update({
        'qdep_applicable': qdep_applicable, 'qcert_applicable': qcert_applicable,
        'qdep_confirmed': bool(dep_confirmed) if qdep_applicable else None,
        'qdep_possible': bool(dep_possible or missing_tail) if qdep_applicable else None,
        'qdep_prefix_possible': bool(dep_possible) if qdep_applicable else None,
        'qcert_confirmed': bool(fc_confirmed) if qcert_applicable else None,
        'qcert_possible': bool(fc_possible or missing_tail) if qcert_applicable else None,
        'qcert_prefix_possible': bool(fc_possible) if qcert_applicable else None,
        'n_dep_confirmed': dep_fail_n, 'n_false_confirmed': fc_count,
        'n_false_possible': fc_possible_count, 'n_certificates': len(certs),
        'risk_interval_unresolved': dep_unresolved+fc_unresolved,
        'deployment_intervals_unresolved': dep_unresolved,
        'certificate_events_unresolved': fc_unresolved,
        'risk_interval_calls': interval_calls, 'angular_evaluations': evaluations,
        'risk_assessment_mode': mode, 'risk_max_bins': max_bins,
        'numerical_standard': 'analytical radial enclosures with floating angular/antiderivative padding; not machine-certified',
        'assumption_status': assumption_status, 'retained_events': retained,
        'risk_numerical_errors': numerical_errors,
    })
    return _plain(out)
