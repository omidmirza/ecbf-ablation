"""Summable per-test error allocation, with test indices starting at one."""
from __future__ import annotations

import numpy as np

JOINT_LOG_NORMALIZER = 2.11  # smallest supported budget normalizer b, above the fixed allocation normalizer
BUDGET_MODES = ("legacy", "single", "joint_log", "joint_log_once")


def joint_log_delta(delta: float, index: int, budget_normalizer: float = JOINT_LOG_NORMALIZER) -> float:
    """Common candidate budget for r >= 1; a larger budget normalizer conservatively underspends delta."""
    if not np.isfinite(budget_normalizer) or budget_normalizer < JOINT_LOG_NORMALIZER:
        raise ValueError("budget_normalizer must be finite and at least 2.11")
    if not np.isfinite(delta) or not 0 < delta < 1:
        raise ValueError("delta must be finite and lie in (0, 1)")
    if isinstance(index, (bool, np.bool_)) or int(index) != index or index < 1:
        raise ValueError("test index must be a positive integer")
    n = int(index) + 1
    return float(delta / (budget_normalizer * n * np.log(n) ** 2))
