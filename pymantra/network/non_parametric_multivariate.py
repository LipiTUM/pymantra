from typing import NamedTuple
import numpy as np
import pandas as pd
from scipy.stats import chi2


class MultivariateRankTest(NamedTuple):
    U2: float
    X: float
    df1: int
    df2: int
    pvalue: float


def multivariate_rank_test(X, Y, paired: bool = False) -> pd.DataFrame:
    r"""Multi-variate generalisation of the Wilcoxon signed-rank test

    A multivariate extension of the Wilcoxon signed-rank test as defined
    by Oja and Randles [Oja04].

    The test statistic :math:`U^2` is defined as

    .. math::
        U^2 = \frac{np}{4c^2_x} \lVert ave \{ S \left( A_x
            \left( X_i + X_j \right ) \right ) \} \rVert^2

    with :math:`U^2 \rightarrow \chi^2_p` (eq. 3 in [Oja04])

    The API of this function is purposely kept analogous to the
    py:func:`pingouin.multivariate_ttest` function in order to be able to use
    them interchangeably

    Parameters
    ----------
    X : np.ndarray

    Y : np.ndarray

    paired : bool, False

    Returns
    -------
    MultivariateRankTest

    References
    ----------
    .. [Oja04] Oja, H. and Randles, R.H., "Multivariate nonparametric tests"
    *Statistical Science*, 19(4), pp.598-605, 2004

    """
    raise NotImplementedError
