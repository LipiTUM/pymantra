"""
pytest tests for confounder correction
"""
import pathlib
import numpy as np
import pandas as pd
import statsmodels.api as sm

from pymantra.network.ld_estimation import confounder_correction


class TestCorrection:
    """
    The test class compares the output of the internal correction function with
    the outputs of the lmer function (R). Correction without random effects
    is not tested against any particular reference but only for shape, since
    it is very basic
    """
    total_data = sm.datasets.get_rdataset("dietox", "geepack").data
    total_data = total_data.loc[~total_data["Feed"].isna(), :]

    met_data = total_data[["Weight", "Feed"]]

    single_con_data = total_data[["Pig", "Time"]]
    double_con_data = total_data[["Pig", "Start", "Time"]]
    triple_con_data = total_data[["Pig", "Start", "Time", "Cu"]]

    base_path = pathlib.Path(__file__).parent
    single_lmer_reference = pd.read_csv(
        base_path / "test_data" / "single_correction.csv", index_col=0)
    double_lmer_reference = pd.read_csv(
        base_path / "test_data" / "double_correction.csv", index_col=0)
    triple_lmer_reference = pd.read_csv(
        base_path / "test_data" / "triple_correction.csv", index_col=0)

    def test_lmm_correction_single_random(self):
        corr_res = confounder_correction(
            self.met_data, self.single_con_data, random_effects=["Pig"])
        assert np.all(
            corr_res.values - self.single_lmer_reference.values < 1e-4)

    def test_lmm_correction_double_random(self):
        corr_res = confounder_correction(
            self.met_data, self.double_con_data,
            random_effects=["Pig", "Start"]
        )
        diff = corr_res.values[:, 0] - self.double_lmer_reference.values[:, 0]
        assert np.all(diff < 1e-4)

    def test_lmm_correction_triple_random(self):
        corr_res = confounder_correction(
            self.met_data, self.triple_con_data,
            random_effects=["Pig", "Start", "Cu"]
        )
        diff = corr_res.values[:, 0] - self.triple_lmer_reference.values[:, 0]
        assert np.all(diff < 1e-4)

    def test_lm_correction(self):
        corr_res = confounder_correction(
            self.met_data, self.single_con_data, random_effects=None)
        assert corr_res.shape == self.met_data.shape
