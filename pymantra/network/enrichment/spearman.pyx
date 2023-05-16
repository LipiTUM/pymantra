# distutils: language = c++

import numpy as np
from collections import namedtuple

from pymantra.network.enrichment.spearman cimport (
    spearmans_rank_1by1,
    spearmans_rank_2by2,
    HAS_PARALLEL
)

import warnings


def _check_type(x: np.ndarray, name: str):
    if x.dtype != 'float64':
        try:
            return x.astype(float)
        except ValueError:
            raise ValueError(
                f"Array {name} must have dtype 'float' or any dtype "
                f"that allows conversion to floats. Found {x.dtype}."
            )
    return x


def multi_threading():
    return HAS_PARALLEL()


SpearmansResults = namedtuple("SpearmansResults", ["correlation", "pvalue"])


def spearmans_correlation(
    x: np.ndarray, y: np.ndarray, n_threads: int = 1
):
    """Compute the spearman's correlation coefficient with a c++ backend.

    `nan` values are automatically ignored and `nan` is returned if less than
    three non-nan observations are available for a pair of features.

    If at least one 2D array is passed arrays will be returned, otherwise
    floats. If two 2D arrays of shape X x N and X x M are passed, the
    returned arrays will be of shape N x M.

    Parameters
    ----------
    x : np.ndarray
    y : np.ndarray
    n_threads : int

    Returns
    -------
    SpearmansResults
        NamedTuple of (1) correlations and (2) correlation-pvalues
    """
    xdim = x.ndim
    ydim = y.ndim
    if xdim > 2 or ydim > 2:
        raise ValueError(
            f"Only 1D and 2D arrays supported, but found {xdim}D and {ydim}D"
        )
    if n_threads > 1 and not HAS_PARALLEL():
        warnings.warn(
            "Parallel processing is only possible when compiled with "
            "OpenMP support", category=RuntimeWarning
        )

    # extension will crash if this is not checked
    x = _check_type(x, 'x')
    y = _check_type(y, 'y')

    if xdim == 2:
        if ydim == 1:
            cors = spearmans_rank_2by2(x.T.tolist(), [y.T.tolist()], n_threads)
        else:
            cors = spearmans_rank_2by2(x.T.tolist(), y.T.tolist(), n_threads)
        return SpearmansResults(np.array(cors[0]), np.array(cors[1]))

    elif ydim == 2:
        cors = spearmans_rank_2by2([x.T.tolist()], y.T.tolist(), n_threads)
        return SpearmansResults(np.array(cors[0]), np.array(cors[1]))

    return SpearmansResults(*spearmans_rank_1by1(x.T.tolist(), y.T.tolist()))
