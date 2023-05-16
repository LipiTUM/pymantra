import numpy as np
# scipy-backwards compatibility
try:
    from scipy.stats import median_abs_deviation as mad_fun
except ModuleNotFoundError:
    from scipy.stats import median_absolute_deviation as mad_fun


# NOTE: functions assume that features are in rows
def var_filter(data: np.ndarray, method: str) -> np.ndarray:
    if method == "nprsd":
        mad = mad_fun(data.T, nan_policy='omit')
        return mad / np.nanmedian(data, axis=1)
    if method == "rsd":
        return np.nanstd(data, axis=1) / np.nanmean(data, axis=1)
    # TODO: other options
    _implemented = {'nprsd', 'rsd'}
    raise ValueError(
        f"{method} is not implemented, please use one of "
        f"the following options: {', '.join(_implemented)}"
    )


def impute(data: np.ndarray) -> np.ndarray:
    imputed = data.copy()
    for i in np.arange(data.shape[0]):
        nans = np.nanmin(data[i, :])
        imputed[i, nans] = np.min(data[i, ~nans]) / 2
    return imputed


def generalised_log(data: np.ndarray, log_func=np.log2, c=1e-10) -> np.ndarray:
    return log_func(data + np.sqrt(data**2 + c))


def quotient_normalise(data: np.ndarray) -> np.ndarray:
    reference = np.nanmedian(data, axis=1)
    dilutions = np.nanmedian((data.T / reference).T, axis=0)
    return data / dilutions
