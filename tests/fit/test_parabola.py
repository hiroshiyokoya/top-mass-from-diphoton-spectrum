"""TMDP.parabola_minimum: best fit from chi^2 on the template grid
(docs/REVIEW.md P2)."""
import os
import sys

import numpy as np
import pytest

ROOT_DIR = os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT_DIR, 'fit'))
import TMDP  # noqa: E402

GRID = np.arange(170.0, 176.5, 0.5)


def chi2_of(x, x0=173.27, sigma=0.8, c=95.0):
    return c + ((np.asarray(x) - x0) / sigma) ** 2


def test_exact_parabola():
    best, sigma, status = TMDP.parabola_minimum(GRID, chi2_of(GRID))
    assert status == 'ok'
    assert best == pytest.approx(173.27, abs=1e-9)
    assert sigma == pytest.approx(0.8, rel=1e-9)


def test_order_does_not_matter():
    perm = np.random.default_rng(0).permutation(len(GRID))
    best, sigma, status = TMDP.parabola_minimum(GRID[perm],
                                                chi2_of(GRID)[perm])
    assert status == 'ok' and best == pytest.approx(173.27, abs=1e-9)


def test_window_selects_points_near_minimum():
    # far from the minimum chi^2 deviates from a parabola; with window=4
    # only points within Delta chi^2 <= 4 (plus the nearest neighbours) count
    chi2 = chi2_of(GRID) + np.where(np.abs(GRID - 173.27) > 2.0,
                                    50.0 * (np.abs(GRID - 173.27) - 2.0), 0.0)
    best, sigma, status = TMDP.parabola_minimum(GRID, chi2)
    assert status == 'ok'
    assert best == pytest.approx(173.27, abs=1e-9)


def test_noisy_chi2():
    rng = np.random.default_rng(3)
    bests = [TMDP.parabola_minimum(GRID, chi2_of(GRID)
                                   + rng.normal(0, 0.3, len(GRID)))[0]
             for _ in range(200)]
    assert np.mean(bests) == pytest.approx(173.27, abs=0.05)


@pytest.mark.parametrize('x0', [169.0, 177.5])
def test_minimum_at_grid_edge(x0):
    best, sigma, status = TMDP.parabola_minimum(GRID, chi2_of(GRID, x0=x0))
    assert status == 'edge'
    assert best == (GRID[0] if x0 < GRID[0] else GRID[-1])
    assert np.isnan(sigma)
