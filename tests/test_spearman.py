import numpy as np
from scipy.stats import spearmanr
from time import time

from pymantra.network.enrichment.spearman import spearmans_correlation


def time_exec(func, *args, name="function", **kwargs):
    start = time()
    r = func(*args, **kwargs)
    end = time()
    print(f"{name} took {round(end - start, 3)} seconds")
    return r


def compare_results(r, p, cor, pval):
    assert np.abs(r[:10, 10:] - cor).sum() < 1e-10
    assert np.abs(p[:10, 10:] - pval).sum() < 1e-10


class TestSpearman:
    x = np.random.randn(5, 10)
    y = np.random.randn(5, 10)

    xnan = np.random.randn(5, 10)
    xnan[0, 0] = np.nan
    xnan[4, 5] = np.nan
    ynan = np.random.randn(5, 10)
    ynan[0, 0] = np.nan
    ynan[3, 9] = np.nan

    def test_2d(self):
        r, p = time_exec(spearmanr, self.x, self.y, name="scipy")
        cor, pval = time_exec(spearmans_correlation, self.x, self.y,
                              1, name="c++")
        compare_results(r, p, cor, pval)

    def test_2d_with_nan(self):
        r, p = time_exec(spearmanr, self.xnan, self.ynan,
                         nan_policy='omit', name="scipy")
        cor, pval = time_exec(spearmans_correlation, self.xnan, self.ynan,
                              1, name="c++")
        compare_results(r, p, cor, pval)

    def test_2d_omp(self):
        r, p = time_exec(spearmanr, self.x, self.y, name="scipy")
        cor, pval = time_exec(spearmans_correlation, self.x, self.y,
                              4, name="c++")
        compare_results(r, p, cor, pval)

    def test_2d_with_nan_omp(self):
        r, p = time_exec(spearmanr, self.xnan, self.ynan,
                         nan_policy='omit', name="scipy")
        cor, pval = time_exec(spearmans_correlation, self.xnan, self.ynan,
                              8, name="c++")
        compare_results(r, p, cor, pval)
