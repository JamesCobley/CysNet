import numpy as np
import pytest

from cysnet.theorem import (
    classify_solution,
    enumerate_state_bounds,
    exact_distribution,
    fully_reduced_bounds,
    oxiform_union_bounds,
)


def test_boundary_marginal_excludes_half_state_space():
    bounds = enumerate_state_bounds([0.5, 0.5, 0.0])

    excluded = set(bounds.loc[bounds["excluded"], "state"])
    expected = {"001", "011", "101", "111"}

    assert excluded == expected


def test_exact_00x_solution():
    bounds = enumerate_state_bounds([0.0, 0.0, 0.25])

    assert classify_solution(bounds) == "exact_singleton"

    dist = exact_distribution(bounds)
    got = dict(zip(dist["state"], dist["probability"]))

    assert np.isclose(got["000"], 0.75)
    assert np.isclose(got["001"], 0.25)
    assert len(got) == 2


def test_bounded_0xx_solution():
    bounds = enumerate_state_bounds([0.0, 0.25, 0.25])

    assert classify_solution(bounds) == "bounded"

    b = bounds.set_index("state")

    assert np.isclose(b.loc["000", "lower"], 0.50)
    assert np.isclose(b.loc["000", "upper"], 0.75)

    assert np.isclose(b.loc["001", "lower"], 0.00)
    assert np.isclose(b.loc["001", "upper"], 0.25)

    assert np.isclose(b.loc["010", "lower"], 0.00)
    assert np.isclose(b.loc["010", "upper"], 0.25)

    assert np.isclose(b.loc["011", "lower"], 0.00)
    assert np.isclose(b.loc["011", "upper"], 0.25)


def test_oxiform_union_bounds():
    lower, upper = oxiform_union_bounds([0.10, 0.20, 0.30])

    assert np.isclose(lower, 0.30)
    assert np.isclose(upper, 0.60)


def test_fully_reduced_bounds():
    lower, upper = fully_reduced_bounds([0.10, 0.20, 0.30])

    assert np.isclose(lower, 0.40)
    assert np.isclose(upper, 0.70)


def test_invalid_marginals_raise_error():
    with pytest.raises(ValueError):
        enumerate_state_bounds([-0.1, 0.5])

    with pytest.raises(ValueError):
        enumerate_state_bounds([0.2, 1.2])


def test_empty_marginals_raise_error():
    with pytest.raises(ValueError):
        enumerate_state_bounds([])
