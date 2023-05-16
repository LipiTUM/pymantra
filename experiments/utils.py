import warnings
import os
from argparse import ArgumentParser, Namespace
from warnings import warn
from typing import Dict, Optional, Callable, Union, Type, Tuple
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from scipy.stats import ranksums
from statsmodels.stats.multitest import multipletests
from skopt.space import Integer, Real, Categorical
from skopt.utils import use_named_args
from skopt import gp_minimize
from skopt.plots import plot_convergence
from sklearn.model_selection import (
    cross_val_score, RepeatedStratifiedKFold)
from sklearn.metrics import RocCurveDisplay, PrecisionRecallDisplay, auc
from sklearn import __version__ as sklearn_version

from pymantra import (
    MetaboliteLocalSearch, MultiOmicsLocalSearch, compute_reaction_estimates)
from pymantra.plotting import residual_violinplot


class CommandLineParser:
    """Just a helper class to avoid re-writing defaults for all experiments"""

    __slots__ = {"_parser", "_local_search_params", "parsed_arguments"}

    _parser: ArgumentParser
    _local_search_params: Dict[str, str]
    parsed_arguments: Optional[Namespace]

    def __init__(self, defaults: dict = None, *args, **kwargs):
        """
        Parameters
        ----------
        defaults : dict, optional
            Default arguments to overwrite pre-defined defaults
        kwargs
            Optional keyword arguments to pass to the argparse.ArgumentParser
            constructor
        """
        self._parser = ArgumentParser(*args, **kwargs)
        self._local_search_params = {
            "temperature": "temp",
            "l_min": "l_min",
            "l_max": "l_max",
            "max_iter": "max_iter",
            "min_reactions": "min_reactions",
            "p": "p"
        }
        self.parsed_arguments = None

        if defaults is None:
            defaults = {}
        # setting defaults
        # general settings for enrichment
        self.add_argument(
            "-n", "--n-threads", type=int,
            default=defaults.get("n_threads", 1),
            help="Number of threads to use when computing the enrichment"
        )
        self.add_argument(
            "-r", "--n-repeat", type=int,
            default=defaults.get("n_repeat", 5),
            help="Number of times to repeat the local search"
        )
        self.add_argument(
            "-c", "--combine-mode", type=str,
            default=defaults.get("combine_mode", "union"),
            choices=["union", "intersection"],
            help="Number of times to repeat the local search"
        )
        # local search-specific parameters
        self.add_argument(
            "-t", "--temperature", type=float,
            default=defaults.get("temperature", 10.),
            help="Initial simulated annealing temperature, exponentially"
                 "decaying every iteration. The higher `temp` the more likely "
                 "it is to a solution with a lower score at any iteration. "
                 "Intuitively more 'hops' will be performed at higher "
                 "temperature."
        )
        self.add_argument(
            "--l-min", type=int, default=defaults.get("l_min", 5),
            help="Minimal solution size"
        )
        self.add_argument(
            "--l-max", type=int, default=defaults.get("l_max", 10),
            help="Maximal solution size"
        )
        self.add_argument(
            "--max-iter", type=int, default=defaults.get("max_iter", 40),
            help="Maximum number of iterations before local search is "
                 "stopped, if the (sub)optimal results has not been found"
        )
        self.add_argument(
            "--min-reactions", type=int,
            default=defaults.get("min_reactions", 3),
            help="Minimum number of reactions to be contained in the solution"
        )
        self.add_argument(
            "-p", type=float, default=defaults.get("p", 2),
            help="Which Lp(Minkowski)-norm to use in the objective function."
        )
        self.add_argument(
            "--multipletest", default="bonferroni", type=str,
            help="Method for multipletest correction. Any option accepted by "
                 "statsmodels.stats.multitest.multipletests can be used."
        )

    def add_argument(self, *args, local_search_arg_name: str = None, **kwargs):
        """Add an argument to the parser

        Parameters
        ----------
        args
            Positional arguments to `parser.ArgumentParser.add_argument`
        local_search_arg_name: str, Optional
            Must be given if the argument should be returned when calling
            `get_local_search_kwargs`
        kwargs
            Keyword arguments to `parser.ArgumentParser.add_argument`

        Returns
        -------

        """
        arg = self._parser.add_argument(*args, **kwargs)
        if local_search_arg_name is not None:
            self._local_search_params[arg.dest] = local_search_arg_name

    def parse_args(self, *args, **kwargs):
        """Parse the command line arguments"""
        if self.parsed_arguments is None:
            self.parsed_arguments = self._parser.parse_args(*args, **kwargs)
        else:
            warn(
                "Arguments have already been parsed and stored in "
                "'parsed_arguments'"
            )

    def get_parsed_args(self):
        """Get all names of parsed arguments"""
        if self.parsed_arguments is None:
            warn(
                "Arguments have not been parsed yet. Please call `parse_args` "
                "first."
            )
        else:
            return self.parsed_arguments.__dict__

    def get_parsed_argnames(self):
        """Get all names of parsed arguments"""
        if self.parsed_arguments is None:
            warn(
                "Arguments have not been parsed yet. Please call `parse_args` "
                "first."
            )
        else:
            return set(self.parsed_arguments.__dict__.keys())

    def get_arg(self, arg_name: str, *args, **kwargs):
        """Get an argument by its name"""
        if self.parsed_arguments is None:
            warn(
                "Arguments have not been parsed yet. Will perform parsing now."
            )
            self.parsed_arguments = self._parser.parse_args(*args, **kwargs)
        return getattr(self.parsed_arguments, arg_name)

    def get_local_search_kwargs(self, *args, **kwargs):
        """Get all local search arguments as a dictionary"""
        if self.parsed_arguments is None:
            warn(
                "Arguments have not been parsed yet. Will perform parsing now."
            )
            self.parsed_arguments = self._parser.parse_args(*args, **kwargs)

        return {
            param_name: getattr(self.parsed_arguments, arg_name)
            for arg_name, param_name in self._local_search_params.items()
        }

    def set_defaults(self, **kwargs):
        """Set default values for parser"""
        self._parser.set_defaults(**kwargs)


def parse_to_graphml(graph: nx.Graph):
    """Format node and edge data for compatibility with graphml"""
    for node, attrs in graph.nodes(data=True):
        attrs_to_pop = [
            attr_key for attr_key, attr_val in attrs.items()
            if isinstance(attr_val, dict)
        ]
        for key in attrs_to_pop:
            graph.nodes[node].pop(key)

    for src, tgt, attrs in graph.edges(data=True):
        edge = (src, tgt)
        attrs_to_pop = [
            attr_key for attr_key, attr_val in attrs.items()
            if isinstance(attr_val, dict)
        ]
        for key in attrs_to_pop:
            graph.edges[edge].pop(key)


def load_data(folder: str, corrected_data: bool = False):
    """Loading data compute in the pre-processing scripts

    Parameters
    ----------
    folder: str
        Folder in which the pre-processed data is contained
    corrected_data: bool, False
        Whether to use confounder-corrected.

    Returns
    -------
    Tuple[pd.Series, pd.DataFrame, nx.Graph]
        sample groups, metabolite data, metabolite-reaction graph
    """
    sample_groups = pd.read_csv(
        os.path.join(folder, "sample_groups.csv"), index_col=0).iloc[:, 0]
    if corrected_data:
        corr_file = os.path.join(folder, "corrected_metabolome.csv")
        if os.path.exists(corr_file):
            metabolite_data = pd.read_csv(
                os.path.join(folder, "corrected_metabolome.csv"), index_col=0)
        else:
            warn("No corrected data available -- defaulting to uncorrected.")
            metabolite_data = pd.read_csv(
                os.path.join(folder, "metabolite_data.csv"), index_col=0)
    else:
        metabolite_data = pd.read_csv(
            os.path.join(folder, "metabolite_data.csv"), index_col=0)
    graph = nx.read_graphml(os.path.join(folder, "mantra_graph.graphml"))

    return sample_groups, metabolite_data, graph


# TODO: add reporting functions


def residual_pvalues(
    residuals: pd.DataFrame, groups: pd.Series,
    test_function: Callable = ranksums
) -> pd.Series:
    """Statistically testing the residual-distributions between groups

    Perform a statistical test (any binary test with the typical `scipy.stats`
    interface) to test whether the reaction-residuals between two groups
    are different.

    Parameters
    ----------
    residuals: pd.DataFrame
        Data frame containing reacion residual as returned by
        `pymantra.compute_reaction_estimates`
    groups: pd.Series
        Sample groups
    test_function: Callable, scipy.stats.ranksums
        Binary statistical test function, API must be the same as the
        `scipy.stats` function interfaces. Default is the wilcoxon rank-sum
        test.

    Returns
    -------
    pd.Series
        p-values in the same order as the columns of the input data frame
    """
    uni_groups = groups.unique()
    assert uni_groups.size == 2, \
        f"Exactly two groups must be provided but {uni_groups.size} found"

    control_mask = (groups == uni_groups[0]).to_numpy()
    pvals = np.zeros(residuals.shape[1])
    for i in range(residuals.shape[1]):
        _, pval = test_function(
            residuals.values[control_mask, i],
            residuals.values[~control_mask, i]
        )
        pvals[i] = pval

    return pd.Series(pvals, index=residuals.columns)


def residual_analysis(
    residuals: pd.DataFrame, sample_groups: pd.Series, base_file: str,
    correction: str, test_function: Callable = ranksums,
    plot: Tuple[plt.Figure, plt.axis] = None,
    violin_args: dict = None, return_qvals: bool = False, **kwargs
):
    """Analyse and plot reaction residuals

    Parameters
    ----------
    residuals: pd.DataFrame
        Residuals as computed in `pymantra.compute_reaction_estimates`
    sample_groups: pd.Series
        Sample group annotations
    base_file: str
        Base-file path to save data and plots. Should not contain any file
        endings.
    correction: str
        Multipletest-correction method. For available options see the
        documentation of `statsmodels.stats.multitest.multipletests`
    test_function: Callable, scipy.stats.ranksums
        Binary statistical test function, API must be the same as the
        `scipy.stats` function interfaces. Default is the wilcoxon rank-sum
        test.
    plot: Tuple[plt.Figure, plt.axis], optional
        Matplotlib figure and subplot axis to plot the violinplot onto
    violin_args: dict, optional
        Additional arguments to pass to :py:func:`seaborn.violinplot`.
        Must not contain 'x', 'y', 'hue', 'data' and 'ax'.
    return_qvals: bool, False
        If True, corrected pvalues are returned
    kwargs
        Optional keyword arguments to parse to `plt.subplots`
    """
    # computing p-values
    res_pvals = residual_pvalues(residuals, sample_groups, test_function)
    res_pvals.to_csv(f"{base_file}_pvalues.csv")

    # multipletest correction
    res_qvals = pd.Series(
        multipletests(res_pvals, method=correction)[1], index=res_pvals.index)
    res_qvals.to_csv(f"{base_file}_qvalues_{correction}.csv")

    # plotting distributions with q-values/significances
    save = True
    if plot is None:
        fig, ax = plt.subplots(
            figsize=kwargs.pop("figsize", (32, 12)), **kwargs)
    else:
        fig, ax = plot
        save = False
    residual_violinplot(
        residuals, sample_groups, res_qvals, ax=ax,
        **violin_args if violin_args is not None else {}
    )
    ax.set_xticks(
        ax.get_xticks(),
        [col.get_text().split(", ")[0] for col in ax.get_xticklabels()]
    )
    if save:
        fig.savefig(f"{base_file}_distributions.pdf")
        plt.close(fig)

    if return_qvals:
        return res_qvals


def correct_correlation_pvalues(
    pvals: np.ndarray, symmetric: bool = True, **kwargs
) -> np.ndarray:
    """Perform multipletest correction on correlation-matrix p-values

    Parameters
    ----------
    pvals: np.ndarray
        Symmetric matrix containing p-values for the respective pairwise
        correlation
    symmetric: bool, False
        Whether pvals is a symmetric matrix or whether row and column features
        are different
    kwargs
        Optional arguments passed to
        `statsmodels.stats.multitest.multipletests`

    Returns
    -------
    np.ndarray
        Symmetric matrix with corrected p-values. Same row/column order as the
        input matrix.
    """
    if symmetric:
        # correction is only happening on the lower triangular part of the
        # matrix
        triu_indices = np.triu_indices_from(pvals)
        flat_arr = pvals[triu_indices[0], triu_indices[1]]
    else:
        flat_arr = pvals.flatten()
    corr_arr = multipletests(flat_arr, **kwargs)[1]

    if symmetric:
        # storing the corrected p-values a new array
        corr_pvals = np.zeros(pvals.shape)
        corr_pvals[triu_indices] = corr_arr
        corr_pvals[np.tril_indices_from(corr_pvals)] = corr_pvals[
            np.triu_indices_from(corr_pvals)]
    else:
        corr_pvals = corr_arr.reshape(pvals.shape)

    return corr_pvals


def run_local_search(
    graph: nx.DiGraph,
    ls_class: Union[Type[MetaboliteLocalSearch], Type[MultiOmicsLocalSearch]],
    parser: CommandLineParser, file_base: str = None, plot: bool = True,
    **kwargs
) -> Union[MetaboliteLocalSearch, MultiOmicsLocalSearch]:
    """Run a local search and plot the results

    Parameters
    ----------
    graph: nx.Graph
        Metabolite-reaction or multi-omics reaction graph containing sample
        data (i.e. `add_reaction_estimates` and optionally
        `add_multiomics_associations` have been called)
    ls_class: Union[Type[MetaboliteLocalSearch], Type[MultiOmicsLocalSearch]]
        Local search class to use, `MetaboliteLocalSearch` if no multi-omics
        associations are included, else `MultiOmicsLocalSearch`
    parser: CommandLineParser
        Parser for command line arguments, importantly `parse_args` must be
        called prior to calling this function.
    file_base: str, optional
        Base file path to be extended inside the function. Should not contain
        any file endings (these will be added insided the function). Must be
        given when `plot` is True
    plot: bool, True
        Whether to plot and save score progression and subnetwork
    kwargs
        Optional keyword arguments for the __init__ functions of
        `MetaboliteLocalSearch` or `MultiOmicsLocalSearch`. Note that
        the arguments from `parser.get_local_search_kwargs` are parsed inside
        the function and should not be added here.

    Returns
    -------
    Union[MetaboliteLocalSearch, MultiOmicsLocalSearch]
        Local search object including the enrichment results
    """
    # initialize a new local search class
    lso = ls_class(
        graph, **kwargs, **parser.get_local_search_kwargs())
    # run n enrichment repeats with random seeds
    lso.run_repeated_local_search(
        parser.get_arg("n_repeat"), n_threads=parser.get_arg("n_threads"),
        min_comp_size=25
    )

    if plot:
        if file_base is None:
            raise ValueError(
                "'file_base' cannot be None when `plot` is set to True"
            )
        lso.plot_score_progression()
        plt.savefig(f"{file_base}_score_progressing.pdf")
        plt.close()

        if isinstance(lso, MetaboliteLocalSearch):
            lso.plot_subnetwork(graph)
        else:
            lso.plot_subnetwork(
                nx.get_node_attributes(graph, "node_type"), graph)
        plt.savefig(f"{file_base}_subnetwork.pdf")
        plt.close()

    return lso


def loss_name_from_sklearn_version():
    """Auxiliary function for scikit-learn version compatibility"""
    major, minor, _ = [int(x) for x in sklearn_version.split(".")]
    if major >= 1 and minor >= 2:
        # since version 1.2.0 => "log_loss"
        return "log_loss"
    # before 1.2.0 => "log"
    return "log"


def generate_search_space(*params_as_tuples):
    """Generate a scikit-optimize parameter search space"""
    search_space = []
    for param in params_as_tuples:
        if not isinstance(param[0], (float, int)):
            if len(param) == 2:
                search_space.append(Categorical(param[0], name=param[1]))
            else:
                search_space.append(
                    Categorical(param[0], name=param[1], **param[2:]))
        else:
            if isinstance(param[0], float):
                space = Real
            else:
                space = Integer
            if len(param) == 3:
                search_space.append(space(param[0], param[1], name=param[2]))
            else:
                search_space.append(
                    space(param[0], param[1], name=param[2], **param[3]))
    return search_space


def optimize_hyperparams(
    param_space: list, classifier, data: pd.DataFrame, groups: pd.Series,
    score: str, cv_args: dict = None, n_jobs: int = 1,
    summary: callable = np.nanmean, n_splits: int = 5, n_repeats: int = 5,
    **kwargs
):
    """Apply Bayesian Optimization with a defined search space and classifier
    """
    if cv_args is None:
        cv_args = {}
    else:
        n_splits = cv_args.pop("n_splits", n_splits)
        n_repeats = cv_args.pop("n_repeates", n_repeats)

    @use_named_args(param_space)
    def eval_space(**params):
        model = classifier(**params)
        cv = RepeatedStratifiedKFold(
            n_splits=n_splits, n_repeats=n_repeats, **cv_args)
        eval_res = cross_val_score(
            model, data, groups, cv=cv, n_jobs=n_jobs, scoring=score)
        return 1 - summary(eval_res)

    return gp_minimize(eval_space, param_space, **kwargs)


def plot_gp_optimization(results, file):
    """Plot the trajectory of the optimized objective (1 - AUC)"""
    figure, ax = plt.subplots(figsize=(16, 9))
    plot_convergence(results, ax=ax)
    figure.savefig(file)
    plt.close(figure)


def plot_evaluation(
    classifier, test_data, test_groups, ax, name=None, type_="ROC", **kwargs
):
    """Plot the evaluation of a classifier on an independent test"""
    if type_ == "ROC":
        roc = RocCurveDisplay.from_estimator(
            classifier, test_data, test_groups)
        roc.plot(ax=ax, name=name, **kwargs)
        ax.plot([0, 1], [0, 1], "--", c="tab:red")
        ax.set_title("ROC Curve")
        ax.set_xlabel("False Positive Rate")
        ax.set_ylabel("True Positive Rate")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
    else:
        pr = PrecisionRecallDisplay.from_estimator(
            classifier, test_data, test_groups)
        pr.plot(ax=ax, name=name, **kwargs)
        ax.plot([1, 0], [0, 1], "--", c="tab:red")
        ax.set_title("Precision Recall Curve")
        ax.set_xlabel("Recall")
        ax.set_ylabel("Precision")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)


def _plot_mean_curve(yvals, mean_x, aucs, curve_type, plot_ax):
    mean_y = np.mean(yvals, axis=0)
    if curve_type == "ROC":
        mean_y[-1] = 1.0
    else:
        mean_y[-1] = 0.0
    mean_auc = auc(mean_x, mean_y)
    std_auc = np.std(aucs)
    label = r"Mean %s (AUC = %0.2f $\pm$ %0.2f)" % \
            (curve_type, mean_auc, std_auc)
    plot_ax.plot(mean_x, mean_y, color="b", alpha=0.7, label=label)

    std_y = np.std(yvals, axis=0)
    y_upper = np.minimum(mean_y + std_y, 1)
    y_lower = np.maximum(mean_y - std_y, 0)
    plot_ax.fill_between(
        mean_x, y_lower, y_upper,
        color="tab:grey", alpha=0.4, label=r"$\pm$ 1 std. dev.",
    )
    if curve_type == "ROC":
        plot_ax.set(
            xlim=[0, 1], ylim=[0, 1],
            xlabel="False Positive Rate",
            ylabel="True Positive Rate",
            title="ROC curves"
        )
    else:
        plot_ax.set(
            xlim=[0, 1], ylim=[0, 1],
            xlabel="Recall",
            ylabel="Precision",
            title="Precision-Recall curves"
        )


def _combined_estimates(graph, data, groups, train, control_group, **kwargs):
    estimation_groups = pd.Series(
        [
            1 if i in train and groups[i] == control_group else 0
            for i in np.arange(data.shape[0])
        ],
        index=data.index
    )
    # computing reaction estimates on control group of training data
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return compute_reaction_estimates(
            graph, data, estimation_groups, control_group=1, **kwargs)
