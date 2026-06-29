from __future__ import annotations

from dataclasses import dataclass
from itertools import product
from typing import Iterable, List

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class StateBound:
    state: str
    lower: float
    upper: float
    possible: bool
    required: bool
    excluded: bool
    fixed_positive: bool
    k_oxidised: int


def validate_marginals(marginals: Iterable[float], eps: float = 1e-12) -> np.ndarray:
    """
    Validate and return cysteine oxidation marginals as fractions in [0, 1].

    Values infinitesimally outside [0, 1] are snapped to the boundary to
    tolerate numerical noise. Values meaningfully outside [0, 1] raise an
    error rather than being silently clipped.
    """
    m = np.asarray(list(marginals), dtype=float)

    if m.ndim != 1:
        raise ValueError("Marginals must be a one-dimensional sequence.")

    if len(m) == 0:
        raise ValueError("At least one marginal is required.")

    if not np.all(np.isfinite(m)):
        raise ValueError("Marginals must be finite numeric values.")

    if np.any(m < -eps) or np.any(m > 1.0 + eps):
        raise ValueError("Marginals must lie in [0, 1].")

    m[m < eps] = 0.0
    m[m > 1.0 - eps] = 1.0

    return m


def enumerate_state_bounds(marginals: Iterable[float], eps: float = 1e-12) -> pd.DataFrame:
    """
    Enumerate observed-coordinate cysteine-redox substates and compute
    sharp theorem bounds for each state from first-order marginals.

    Scope
    -----
    These bounds are exact for the binary state space defined by the
    supplied marginals.

    Therefore:
      - if all protein cysteines are supplied, the bounds are over the full
        cysteine-redox oxiform space;
      - if only detected cysteines are supplied, the bounds are over the
        observed-coordinate projection and should not be interpreted as
        assignments of unmeasured cysteines.

    Theorem
    -------
    For a binary state s and site marginals m_j, define

        q_j(s) = m_j       if s_j = 1
        q_j(s) = 1 - m_j   if s_j = 0

    Then the sharp Frechet bounds for the probability assigned to state s are

        lower_s = max(0, sum_j q_j(s) - (R - 1))
        upper_s = min_j q_j(s)

    These bounds are sharp given only first-order marginals. CysNet uses them
    to determine whether each observed-coordinate substate is compatible,
    excluded, required or fixed.
    """
    m = validate_marginals(marginals, eps=eps)
    R = len(m)

    rows: List[StateBound] = []

    for bits in product([0, 1], repeat=R):
        bits_arr = np.asarray(bits, dtype=int)
        state = "".join(str(x) for x in bits)

        q = np.where(bits_arr == 1, m, 1.0 - m)

        lower = max(0.0, float(np.sum(q) - (R - 1)))
        upper = float(np.min(q))

        if abs(lower) < eps:
            lower = 0.0

        if abs(upper) < eps:
            upper = 0.0

        possible = upper > eps
        required = lower > eps
        excluded = not possible
        fixed_positive = required and abs(upper - lower) <= eps

        rows.append(
            StateBound(
                state=state,
                lower=lower,
                upper=upper,
                possible=possible,
                required=required,
                excluded=excluded,
                fixed_positive=fixed_positive,
                k_oxidised=state.count("1"),
            )
        )

    return pd.DataFrame([r.__dict__ for r in rows])


def classify_solution(bounds: pd.DataFrame, eps: float = 1e-12) -> str:
    """
    Classify whether the feasible observed-coordinate ensemble is exact or bounded.

    The feasible ensemble is exact if every state has identical lower and
    upper bounds. Otherwise, the data define a bounded non-singleton feasible
    set.
    """
    required_cols = {"state", "lower", "upper"}

    missing = required_cols.difference(bounds.columns)
    if missing:
        raise ValueError(f"Bounds table is missing required columns: {sorted(missing)}")

    widths = bounds["upper"].to_numpy(float) - bounds["lower"].to_numpy(float)

    if np.nanmax(np.abs(widths)) <= eps:
        return "exact_singleton"

    return "bounded"


def exact_distribution(bounds: pd.DataFrame, eps: float = 1e-12) -> pd.DataFrame:
    """
    Return the positive observed-coordinate substate distribution when exact.

    Raises
    ------
    ValueError
        If the feasible set is bounded rather than a singleton.
    """
    status = classify_solution(bounds, eps=eps)

    if status != "exact_singleton":
        raise ValueError("Solution is not exact; distribution remains bounded.")

    out = bounds[bounds["lower"] > eps].copy()
    out = out.rename(columns={"lower": "probability"})

    return out[["state", "probability", "k_oxidised"]].reset_index(drop=True)


def oxiform_union_bounds(marginals: Iterable[float], eps: float = 1e-12) -> tuple[float, float]:
    """
    Bounds for the fraction of molecules carrying at least one oxidised
    observed cysteine.

    This is the union of oxidised-coordinate events.

        lower = max_j m_j
        upper = min(1, sum_j m_j)

    These are the bounds used for oxiform-compatible copy-number scaling.
    """
    m = validate_marginals(marginals, eps=eps)

    lower = float(np.max(m))
    upper = float(min(1.0, np.sum(m)))

    return lower, upper


def fully_reduced_bounds(marginals: Iterable[float], eps: float = 1e-12) -> tuple[float, float]:
    """
    Bounds for the fraction of molecules reduced at all observed cysteines.

    This is the complement of the union of oxidised-coordinate events.

        lower = max(0, 1 - sum_j m_j)
        upper = 1 - max_j m_j
    """
    m = validate_marginals(marginals, eps=eps)

    lower = float(max(0.0, 1.0 - np.sum(m)))
    upper = float(1.0 - np.max(m))

    return lower, upper
