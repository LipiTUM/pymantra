import warnings
import numpy as np
from typing import Tuple, Type, Union, List
# from skbio.diversity import alpha_diversity, beta_diversity
# from scipy.spatial.distance import pdist


def centred_log_ratio(data: np.ndarray) -> Tuple[Type[np.ndarray], np.ndarray]:
    r"""Centred log ratio transformation on a compositional matrix

    The clr transformation is defined as

    .. math::
        clr(x) & = &
            \left{
                log\left(\frac{x_1}{G(x)}\right),
                \dots,
                log\left(\frac{x_n}{G(x)}\right)
            \right} \\\\
            & = &
            \left{
                log(x_1) - log(G(x)), \dots, log(x_n) - log(G(x))
            \right}

    Note that prior to transformation all-zero features (components) are
    removed.

    Parameters
    ----------
    data : np.ndarray
        Compositional 2-dimensional matrix with features in rows and samples
        in columns

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        Tuple where the first element is the clr transformed version of data
        and the second element are the rows with zero counts that were removed

    examples
    --------
    >>> # Taken from the scikit-bio library
    >>> import numpy as np
    >>> x = np.array([[.1, .3, .4, .2], [.5, .25, .05, .1]]).T
    >>> centred_log_ratio(x)
    (array([[-0.79451346,  1.15129255],
           [ 0.30409883,  0.45814537],
           [ 0.5917809 , -1.15129255],
           [-0.10136628, -0.45814537]]), array([False, False, False, False]))
    """
    if data.ndim == 1:
        data = np.atleast_2d(data).T
    # filtering out rows with zero counts
    zero_rows = data.sum(axis=1) == 0
    if np.any(zero_rows):
        warnings.warn(
            "Rows with only zero values found. The following rows will be "
            f"removed: {np.where(zero_rows)[0]}"
        )
    non_zeros = data[~zero_rows, :]
    # checking for zero values to avoid undefined log
    if np.any(non_zeros == 0):
        raise ValueError(
            "Cells with zero counts found after removal of all zero features ("
            "components). Please modify your data to remove zero counts in "
            "non-constant features, e.g. by adding a pseudo count."
        )
    # ensuring equal row sums (i.e. all samples/compositions have the same
    # total)
    normed_mat = non_zeros / non_zeros.sum(axis=0)
    # log matrix
    log_mat = np.log(normed_mat)
    # geometric mean matrix
    geom_means = np.mean(log_mat, axis=0)
    # applying clr transformation
    clr_transformed = (log_mat - geom_means).squeeze()

    return clr_transformed, zero_rows


def rel_counts(data: np.ndarray) -> np.ndarray:
    rel = np.zeros(data.shape)
    for i in np.arange(data.shape[0]):
        rel[i, :] = data[i, :] / np.nansum(data[i, :])
    return rel


def filter_zero_fractions(
    data: np.ndarray,
    thresh: float = .5,
    return_data: bool = False
) -> Union[np.ndarray, List[int]]:
    to_retain = []
    for i in np.arange(data.shape[0]):
        if np.mean(data[i, :] == 0) < thresh:
            to_retain.append(i)
    if return_data:
        return data[to_retain, :]
    return to_retain
