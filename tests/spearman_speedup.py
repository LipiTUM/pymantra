import timeit
from pymantra.network.enrichment import spearmans_correlation
from scipy.stats import spearmanr
import numpy as np


def scipy_fun():
    x = np.random.randn(10, 200)
    y = np.random.randn(10, 200)
    x[0, 10] = np.nan
    x[5, 25] = np.nan
    x[9, 104] = np.nan

    y[0, 10] = np.nan
    y[6, 25] = np.nan
    y[8, 104] = np.nan
    r, p = spearmanr(x, y, nan_policy='omit')


def cpp_fun():
    x = np.random.randn(10, 200)
    y = np.random.randn(10, 200)

    x[0, 10] = np.nan
    x[5, 25] = np.nan
    x[9, 104] = np.nan

    y[0, 10] = np.nan
    y[6, 25] = np.nan
    y[8, 104] = np.nan
    r, p = spearmans_correlation(x, y)


def omp_fun():
    x = np.random.randn(10, 200)
    y = np.random.randn(10, 200)

    x[0, 10] = np.nan
    x[5, 25] = np.nan
    x[9, 104] = np.nan

    y[0, 10] = np.nan
    y[6, 25] = np.nan
    y[8, 104] = np.nan
    r, p = spearmans_correlation(x, y, 4)


cpp_ = min(timeit.Timer(cpp_fun).repeat(repeat=10, number=1))
print(cpp_)
omp_ = min(timeit.Timer(omp_fun).repeat(repeat=10, number=1))
print(omp_)
scipy_ = min(timeit.Timer(scipy_fun).repeat(repeat=10, number=1))
print(scipy_, cpp_, omp_)
