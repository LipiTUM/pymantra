import numpy as np
from typing import Union, List, Tuple


def filter_zero_fractions(
    data: np.ndarray,
    thresh: float = .5,
    return_data: bool = False
) -> Union[np.ndarray, List[int]]:
    """Filter features by sparsity

    Parameters
    ----------
    data: np.ndarray
        Data to filter, features in rows
    thresh: float, 0.5
        Sparsity threshold to filter by. All features with at most this
        fraction of zero values will be retained
    return_data: bool, False
        Whether to return the filtered data or the filtering mask

    Returns
    -------
    Union[np.ndarray, List[int]]
        Filtered data or filtering mask
    """
    to_retain = []
    for i in np.arange(data.shape[0]):
        if np.mean(data[i, :] == 0) < thresh:
            to_retain.append(i)
    if return_data:
        return data[to_retain, :]
    return to_retain


def impute_metabolites(
    data: np.ndarray,
    frac: float = 2.,
    inplace: bool = False,
    return_mins: bool = False
) -> Union[None, np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Impute metabolites by a fraction of the feature-wise minimum

    Parameters
    ----------
    data: np.ndarray
        Data to impute
    frac: float, 2.0
        Fraction of the feature-wise minimum to impute. The minimum value will
        be multiplied by 1/`frac`
    inplace: bool, False
        Whether to impute inplace or return a new array
    return_mins: bool, False
        Whether to return values used for imputation

    Returns
    -------
    Union[None, np.ndarray, Tuple[np.ndarray, np.ndarray]]
        Imputed data if `inplace` is False, else None or Imputed data and
        row-wise imputation values if `return_mins` is True
    """
    if np.any(np.isnan(data)):
        raise ValueError(
            "No nan values are allowed in 'data'"
        )
    min_vals = np.zeros(data.shape[0])
    if inplace:
        for i in np.arange(data.shape[0]):
            zero_mask = data[i, :] == 0
            min_i = np.min(data[i, ~zero_mask]) / frac
            if np.any(zero_mask):
                data[i, zero_mask] = min_i
            min_vals[i] = min_i
        if return_mins:
            return min_vals
    else:
        data_ = data.copy()
        for i in np.arange(data.shape[0]):
            zero_mask = data[i, :] == 0
            min_i = np.min(data[i, ~zero_mask])/frac
            if np.any(zero_mask):
                data_[i, zero_mask] = min_i
            min_vals[i] = min_i
        if return_mins:
            return data_.astype(np.float64), min_vals
        return data_.astype(np.float64)


def apply_imputation(data: np.ndarray, min_vals: np.ndarray) -> np.ndarray:
    data_ = data.copy()
    for i in np.arange(data.shape[0]):
        data[i, data[i, :] == 0] = min_vals[i]
    return data_


def impute_relative(
    data: np.ndarray,
    inplace: bool = False,
    total_min: float = None,
    return_min: bool = False
) -> Union[None, np.ndarray, float, Tuple[np.ndarray, float]]:
    """Perform relative imputation

    Imputing relative data while maintaining constant sums over all samples
    Imputed value is defined as min(total_data)/#zero_in_sample

    Parameters
    ----------
    data: np.ndarray
        Data to impute
    inplace: bool, False
        Whether to impute inplace or return a new array
    total_min: float, optional
        If given, this value will be used to impute instead of computing the
        global minimum
    return_min: bool, False
        Whether to return the global minimum

    Returns
    -------
    Union[None, np.ndarray, float, Tuple[np.ndarray, float]]
        Imputed data if `inplace` is False, else None
    """
    zero_mask = data == 0
    zero_vals = np.sum(zero_mask, axis=0)
    if total_min is None:
        total_min = np.min(data[~zero_mask])
    if inplace:
        for i in np.arange(data.shape[1]):
            data[data[:, i] == 0, i] = total_min / zero_vals[i]
        if return_min:
            return total_min
    else:
        imputed = data.copy()
        for i in np.arange(data.shape[1]):
            imputed[imputed[:, i] == 0, i] = total_min/zero_vals[i]
        if return_min:
            return imputed.astype(np.float64), total_min
        return imputed.astype(np.float64)


def quotient_normalise(
    data: np.ndarray,
    return_dilutions: bool = False
) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
    """Perform quotient normalisation

    Parameters
    ----------
    data: np.ndarray
        Data to normalise
    return_dilutions: bool, False
        Whether to return sample-wise dilution values

    Returns
    -------
    Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]
        Normalised data or 2-tuple of normalised data (2D-array) and the
        sample-wise dilutions (1D-array)
    """
    # feature-wise median
    ref = np.nanmedian(data, axis=1)
    # dilution factors
    dilutions = np.empty(data.shape[1], dtype=np.float64)
    normed = np.empty(data.shape, dtype=np.float64)
    for i in np.arange(data.shape[1]):
        dilutions[i] = np.nanmedian(data[:, i]/ref)
        # applying dilutions
        normed[:, i] = data[:, i]/dilutions[i]
    if return_dilutions:
        return normed, dilutions
    return normed


def generalised_log(data: np.ndarray, epsilon: float = 1e-5) -> np.ndarray:
    r"""Apply the generalised log-transform

    The generalised log-transform is defined as

    .. math::
        log \frac{x + \sqrt[p]{x^p + \epsilon^p}}{2}

    and thus avoids log(0) situations.

    For this function :math:`log_2` is used.

    Parameters
    ----------
    data: np.ndarray
        Data to transform
    epsilon: float, 1e-5
        Constant to add

    Returns
    -------
    np.ndarray
        log-transformed data
    """
    return np.log2(
        (data + np.power(np.power(data, 10) + epsilon**10, 1/10)) / 2
    )


def scale_z_scores(data: np.ndarray, omit_nan: bool = True):
    """Apply z-score scaling

    Scale all features to zero mean and unit variance

    Parameters
    ----------
    data: np.ndarray
        Input data, features in rows, samples in columns
    omit_nan: bool, True
        Whether to omit or propagate nan values

    Returns
    -------
    np.ndarray
        Scaled data
    """
    if omit_nan:
        mean_fun = np.nanmean
        std_fun = np.nanstd
    else:
        mean_fun = np.mean
        std_fun = np.std
    centered = data.T - mean_fun(data, axis=1)
    return (centered / std_fun(centered, axis=0)).T
