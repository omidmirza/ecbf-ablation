"""Summable per-test error allocation, with test indices starting at one."""
from __future__ import annotations

import numpy as np

# smallest budget normalizer b: kappa = sum_{n>=2} 1/(n ln^2 n) ~ 2.1097, rounded up (report eq. (4))
JOINT_LOG_NORMALIZER = 2.11
BUDGET_MODES = ("legacy", "single", "joint_log", "joint_log_once")


def joint_log_delta(delta: float, index: int, budget_normalizer: float = JOINT_LOG_NORMALIZER) -> float:
    """Level delta_r = delta / [b (r+1) ln^2(r+1)] of test r >= 1, common to all candidates (report eq. (4));
    a larger budget normalizer b underspends delta."""
    if not np.isfinite(budget_normalizer) or budget_normalizer < JOINT_LOG_NORMALIZER:
        raise ValueError("budget_normalizer must be finite and at least 2.11")
    if not np.isfinite(delta) or not 0 < delta < 1:
        raise ValueError("delta must be finite and lie in (0, 1)")
    if isinstance(index, (bool, np.bool_)) or int(index) != index or index < 1:
        raise ValueError("test index must be a positive integer")
    n = int(index) + 1
    return float(delta / (budget_normalizer * n * np.log(n) ** 2))
