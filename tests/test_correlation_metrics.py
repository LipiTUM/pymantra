"""
Testing correlation metrics and covariance matrix estimators
for partial correlations
"""
import pathlib
import pandas as pd
import numpy as np
from pymantra.network.correlation_metrics import (
    pairwise_correlations, available_correlation_metrics
)


cov_estimators = [
    'schaefer_strimmer', 'rblw', 'ledoit_wolf',
    'graph_lasso', 'spearman'
]
# optimal shrinkage found by analytical estimator by corpcor::pcor.shrink
optimal_shrinkage = 0.372073
CONST = 1e-7

base = pathlib.Path(__file__).parent.absolute()
test_data = pd.read_csv(
    base / "test_data/correlation_test_data.csv", index_col=0
)
# NOTE: some covariance estimators cannot handle missing values
pcor_data = test_data.fillna(0)
pcor_output = pd.read_csv(
    base / "test_data/correlations_pcor.csv", index_col=0
)


# NOTE: does NOT contain output value tests, since only
# numpy/scipy/pandas functions are used
def test_correlation_functions():
    for corr_metric in available_correlation_metrics().difference('partial'):
        # numpy
        tmp_np = pairwise_correlations(test_data.to_numpy(), corr_metric)
        # pandas
        tmp_pd = pairwise_correlations(test_data, corr_metric)
        # NOTE: this does not check for correct values but only for the same
        # results with numpy and pandas
        assert np.any(abs(tmp_np - tmp_pd.to_numpy()) < CONST)
        assert tmp_np.shape == (test_data.shape[1], test_data.shape[1])


def test_partial_correlation_functions():
    # NOTE: some partial correlation estimators cannot handle zero-values
    for estimator in cov_estimators:
        if estimator == 'rblw':
            n_sample = pcor_data.shape[0]
            # numpy
            tmp_np = pairwise_correlations(
                np.cov(pcor_data.to_numpy().T), 'partial', n=n_sample,
                estimator=estimator
            )
            # pandas
            pd_cov = pd.DataFrame(np.cov(pcor_data.T), index=pcor_data.columns,
                                  columns=pcor_data.columns)
            tmp_pd = pairwise_correlations(
                pd_cov, 'partial', n=n_sample, estimator=estimator
            )
        else:
            # numpy
            tmp_np = pairwise_correlations(pcor_data.to_numpy(), 'partial',
                                           estimator=estimator)
            # pandas
            tmp_pd = pairwise_correlations(pcor_data, 'partial',
                                           estimator=estimator)
        # NOTE: this does not check for correct values but only for the same
        #       results with numpy and pandas
        assert np.any(abs(tmp_np - tmp_pd) < CONST)
        assert tmp_np.shape == (test_data.shape[1], test_data.shape[1])


def test_schaefer_strimmer_pcors():
    # nans = np.isnan(pcor_output.to_numpy())
    # numpy
    tmp_np = pairwise_correlations(pcor_data.to_numpy(), 'partial')
    assert np.all(abs(tmp_np - pcor_output.to_numpy()) < CONST)
    # pandas
    tmp_pd = pairwise_correlations(pcor_data, 'partial')
    assert np.all(abs(tmp_pd.to_numpy() - pcor_output.to_numpy()) < CONST)
