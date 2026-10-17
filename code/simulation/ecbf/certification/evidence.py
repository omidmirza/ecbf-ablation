"""Mixture e-process updates."""
from __future__ import annotations

import numpy as np
from scipy.special import logsumexp


class GainSet:
    """The common gain set Lambda in (0, 1/(1-alpha)) and its log factors."""

    def __init__(self, gains, alpha: float):
        g = np.asarray(gains, dtype=float).ravel()
        if g.size == 0 or not 0 < alpha < 1:
            raise ValueError("need a nonempty gain set and 0 < alpha < 1")
        if np.any(g <= 0) or np.any(g >= 1 / (1 - alpha)):
            raise ValueError("every gain must lie in (0, 1/(1-alpha)) so that every factor is positive")
        self.gains = g
        self.alpha = float(alpha)
        self.log_up = np.log1p(g * alpha)            # log b when Z = 0
        self.log_down = np.log1p(-g * (1 - alpha))   # log b when Z = 1 (negative, finite)
        self.log_size = float(np.log(g.size))

    @property
    def size(self) -> int:
        return int(self.gains.size)

    def log_factor(self, z) -> np.ndarray:
        """log b^lambda for exceedance indicators z (shape (...,)); returns shape (..., |Lambda|)."""
        z = np.asarray(z, dtype=bool)
        return np.where(z[..., None], self.log_down, self.log_up)

    def log_correction_step(self, z_smallest: bool) -> np.ndarray:
        """log min{1, b^lambda(c^(1))}: zero without an exceedance of c^(1), log(1 - lambda (1 - alpha)) otherwise."""
        return self.log_down.copy() if bool(z_smallest) else np.zeros(self.size)

    def log_mixture(self, logw) -> np.ndarray:
        """log of |Lambda|^{-1} sum_lambda exp(logw) along the last axis."""
        return logsumexp(np.asarray(logw, dtype=float), axis=-1) - self.log_size

    def growth_rate(self, p: float) -> np.ndarray:
        """g_p(lambda) = (1-p) log(1 + lambda alpha) + p log(1 - lambda (1 - alpha)) (Proposition evidence-growth)."""
        return (1 - p) * self.log_up + p * self.log_down
