import pytest
from string import ascii_lowercase, ascii_uppercase
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from pymantra.datasets import example_graph
from pymantra.plotting import (
    plot_directed_graph, plot_undirected_graph, plot_reaction_association,
    residual_violinplot, plot_correlation_differences
)


N_REACTIONS = 10
N_MO = 5
N_SAMPLES = 4
N_CTRL = N_SAMPLES // 2

REACTIONS = [ascii_uppercase[i] for i in np.arange(N_REACTIONS)]
MO = [ascii_lowercase[i] for i in np.arange(N_MO)]
SAMPLES = [f"s{i}" for i in np.arange(N_SAMPLES)]


def _get_random_corrs(ni, nc, index, columns):
    return pd.DataFrame(np.random.randn(ni, nc), index=index, columns=columns)


def _get_random_pvals(size, index, low=None, columns=None):
    if low is None:
        pvals = np.random.uniform(size=size)
    else:
        pvals = np.random.uniform(low=low, size=size)

    if columns is None:
        return pd.Series(pvals, index)
    return pd.DataFrame(
        pvals.reshape(len(index), len(columns)), index=index, columns=columns)


class TestPlotting:
    np.random.seed(42)

    graph = example_graph()

    res = _get_random_corrs(N_SAMPLES, N_REACTIONS, SAMPLES, REACTIONS)
    mo = _get_random_corrs(N_MO, N_SAMPLES, MO, SAMPLES)
    corrs = _get_random_corrs(N_MO, N_REACTIONS, MO, REACTIONS)

    pvals = _get_random_pvals(N_REACTIONS, REACTIONS)
    unsig_pvals = _get_random_pvals(N_REACTIONS, REACTIONS, .1)

    groups = pd.Series(
        (["N"] * N_CTRL) + (["C"] * (N_SAMPLES - N_CTRL)), index=SAMPLES)

    group_corrs = {
        "N": _get_random_corrs(N_MO, N_REACTIONS, MO, REACTIONS),
        "C": _get_random_corrs(N_MO, N_REACTIONS, MO, REACTIONS)
    }
    group_pvals = {
        "N": _get_random_pvals(N_MO * N_REACTIONS, MO, columns=REACTIONS),
        "C": _get_random_pvals(N_MO * N_REACTIONS, MO, columns=REACTIONS)
    }

    def test_plot_directed_graph(self):
        plot_directed_graph(self.graph)
        plt.close()
        plot_directed_graph(self.graph, formula_as_reaction_label=True)
        plt.close()

    def test_plot_undirected_graph(self):
        plot_undirected_graph(self.graph)
        plt.close()
        plot_undirected_graph(self.graph, formula_as_reaction_label=True)
        plt.close()

    def test_plot_reaction_associations(self):
        plot_reaction_association(self.res, self.mo, self.corrs)
        plt.close()
        plot_reaction_association(self.res, self.mo, self.corrs, self.groups)
        plt.close()
        with pytest.raises(ValueError):
            plot_reaction_association(
                self.res.T, self.mo, self.corrs, self.groups)
            plt.close()
        with pytest.raises(ValueError):
            plot_reaction_association(
                self.res, self.mo, self.corrs.T, self.groups)
            plt.close()

        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=True,
            cluster=True, return_differences=True
        )
        plt.close()
        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=True,
            cluster=True, return_differences=False
        )
        plt.close()
        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=False,
            cluster=False, return_differences=True, reorder=False
        )
        plt.close()
        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=True,
            cluster=False, return_differences=False, reorder=False
        )
        plt.close()
        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=False,
            cluster=False, return_differences=True, reorder=True
        )
        plt.close()
        plot_correlation_differences(
            self.group_corrs, self.group_pvals, "N", "C", set_zero=True,
            cluster=False, return_differences=False, reorder=True
        )
        plt.close()

    def test_plot_residuals(self):
        residual_violinplot(self.res, self.groups, rotate_labels=True)
        plt.close()
        residual_violinplot(
            self.res, self.groups, pvalues=self.pvals, significance_only=True,
            drop_legend=True
        )
        plt.close()
        residual_violinplot(self.res, self.groups, pvalues=self.unsig_pvals)
        plt.close()
