"""
General options for improvement:
* different linear models (incl. regularized ones)
* different ways of correcting
"""
import sys
import warnings
import numpy as np
import pandas as pd
import networkx as nx
from tqdm import tqdm
from scipy.stats import f
from typing import Tuple, Dict, List, NamedTuple, Union, Callable
from sklearn.base import RegressorMixin
from sklearn.preprocessing import StandardScaler
# TODO: allow for ridge, lasso and/or elastic net?
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy.stats import ttest_ind, wilcoxon, kstest, norm
from multiprocessing import Process, Manager
with warnings.catch_warnings() as w:
    warnings.simplefilter("ignore")
    # outdated is a little annoying here
    from pingouin import multivariate_normality, multivariate_ttest
from statsmodels.regression.mixed_linear_model import MixedLM
from statsmodels.tools.sm_exceptions import ConvergenceWarning

from pymantra.network import reaction_graph_extraction
from pymantra.network.exceptions import R2Warning
from pymantra.network.non_parametric_multivariate import multivariate_rank_test


class LinearModel(NamedTuple):
    r"""
    Holding the results of a linear model between substrate
    and product intensities.

    substrates: List[str]
        Name of the substrates of the represented reaction
    products: List[str]
        Name of the products of the represented reaction
    model_fit : np.ndarray
        Normed residuals for all samples
    Model : RegressorMixin
        Underlying regression model
    nll : np.ndarray
        Negative log-likelihood aggregated over all products
        of the respective reaction
    r2 : np.ndarray
        :math:`R^2` per product
    f_statistics : np.ndarray
        f-statistics per product
    ftest_pvalue: np.ndarray
        f-test p-value per product
    active : bool
        Indicating whether the 'threshold' for an active
        reaction was met
    outlier_samples: list
        Outliers removed from the training data
    model_fit_type : str, "unset"
        The type of metric used for "model_fit". In most cases
        this would be "residuals", "normed_residuals" or
        "explained_variance"
    """
    substrates: list
    products: list

    model_fit: np.ndarray
    model: RegressorMixin
    nll: np.ndarray
    r2: np.ndarray
    f_statistics: np.ndarray
    ftest_pvalue: np.ndarray
    active: bool
    outlier_samples: list
    model_fit_type: str = "unset"


def _nll(X, Y, y_hat):
    # sum of squared distances
    ssd = np.sum((Y - y_hat) ** 2, axis=0)
    n_observ = X.shape[0] / 2
    # Negative log-likelihood
    # formula derived from the assumption of a scale of 1 and the log
    # pdf of a normal distribution
    return -n_observ * (np.log(2 * np.pi) + np.log(ssd / X.shape[0]) + 1)


def _r_square(X, Y, y_hat):
    # average residual sum of squares
    var_fit = np.sum((Y - y_hat) ** 2, axis=0) / X.shape[0]
    # average total sum of squares
    var_mean = np.sum((Y - np.mean(Y, axis=0)) ** 2, axis=0) / X.shape[0]
    return (var_mean - var_fit) / var_mean


def _test_explained_variances(var_ctr, var_case, alpha=.05):
    if var_ctr.ndim == 2 and var_ctr.shape[1] > 1:
        # multivariate case
        is_normal = multivariate_normality(var_ctr, alpha).normal and \
                        multivariate_normality(var_case, alpha).normal
        if is_normal:
            return multivariate_ttest(var_ctr, var_case)["p-val"]
        else:
            try:
                # TODO: implement this
                return multivariate_rank_test(var_ctr, var_case).pvalue
            except NotImplementedError:
                warnings.warn(
                    "Non-parametric multivariate test is currently not "
                    "implemented! Defaulting to parametric test, despite "
                    "non-normality."
                )
                return multivariate_ttest(var_ctr, var_case)["p-val"]
    # univariate case
    ctr_normal = kstest(var_ctr, norm.cdf).pvalue
    case_normal = kstest(var_case, norm.cdf).pvalue
    if ctr_normal and case_normal:
        return ttest_ind(var_ctr, var_case, nan_policy="omit")[1]
    else:
        return wilcoxon(var_ctr, var_case, nan_policy="omit"[1])


def _fstats(
    X, Y, y_hat, as_probability: bool = True,
    ssd_simple: np.ndarray = None, p_simple: int = None,
    norm_multiple: bool = False
):
    if ssd_simple is None:
        if norm_multiple:
            dist = np.linalg.norm(Y - np.mean(Y, axis=0), axis=1)
            ssd_simple = np.sum(dist ** 2, axis=0)
        else:
            ssd_simple = np.sum((Y - np.mean(Y, axis=0)) ** 2, axis=0)
        p_simple = 1
    elif p_simple is None:
        raise ValueError(
            "'p_simple' must be given, when 'ssd_simple' is given."
        )
    if norm_multiple:
        dist = np.linalg.norm(Y - y_hat, axis=1)
        ssd_fit = np.sum(dist ** 2, axis=0)
    else:
        ssd_fit = np.sum((Y - y_hat) ** 2, axis=0)
    df1 = X.shape[1] - p_simple
    df2 = X.shape[0] - X.shape[1]
    fstats = ((ssd_simple - ssd_fit) / df1) / (ssd_fit / df2)
    if as_probability:
        # probabilities from f-statistic
        return fstats, f.sf(fstats, df1, df2)
    return fstats


def _cooks_distance(y, yhat, X, p: [int, str] = 2):
    # for details on cooks' distance see
    # https://en.wikipedia.org/wiki/Cook%27s_distance
    leverage_mat = X * np.linalg.pinv(X).T
    leverage = leverage_mat.sum(axis=1)

    # NOTE: rank used instead of number of variables
    df = X.shape[0] - np.linalg.matrix_rank(X)
    res = _normed_residuals(y, yhat, p)
    # square root of MSE
    # NOTE: we use the p-norm here to get one distance per sample
    me = np.sqrt(res.dot(res) / df)
    # final distance calculation
    res_norm = res / (me * np.sqrt(1 - leverage))
    dist = (res_norm ** 2 / X.shape[1]) * (leverage / (1 - leverage))

    return dist, df


def _normed_residuals(y, y_hat, p=2, *args, **kwargs):
    """Compute the p-norm for each residual

    Parameters
    ----------
    y: np.ndarray
        True labels/target variables
    y_hat: np.ndarray
        Predicted labels/target variables
    p: Union[int, float, np.inf]
        p-norm to use
    args
        Just provided for compatibility with other function signatures
    kwargs
        Just provided for compatibility with other function signatures

    Returns
    -------
    np.ndarray
        Normed residuals per sample
    """
    if y.ndim > 1:
        return np.linalg.norm(y_hat - y, ord=p, axis=1) / y.shape[1]
    else:
        return y_hat - y


def _explained_variance(y, y_hat, *args, **kwargs):
    """Compute the explained variance per samples

    The explained variance is computed as the ratio of the residual sum of
    squares and the total sum of squares.

    Parameters
    ----------
    y: np.ndarray
        True labels/target variables
    y_hat: np.ndarray
        Predicted labels/target variables
    args
        Just provided for compatibility with other function signatures
    kwargs
        Just provided for compatibility with other function signatures

    Returns
    -------
    np.ndarray
        Explained variance per sample
    """
    # residual sum of squares
    res_ss = ((y - y_hat) ** 2).sum(axis=1)
    # total sum of squares
    total_ss = ((y - y.mean(axis=0)) ** 2).sum(axis=1)
    # explained variance as ratio between residual and total sum of squares
    expl_var = 1 - (res_ss / (total_ss + 1e-4))
    # preventing negative values since res / total is in (-inf, 1]
    expl_var[expl_var < 1e-4] = 1e-4
    return expl_var


def _compute_model(
    X: np.ndarray, Y: np.ndarray, residual_summary: Callable,
    model: RegressorMixin = None,
    substrates: list = None, products: list = None,
    p_tresh: float = 0.05, p: int = 2, simple_model: RegressorMixin = None,
    simple_data: np.ndarray = None, norm_multiple: bool = False,
    outlier_threshold: float = None, r2_threshold: float = .5,
    r2_warning: Callable = lambda: None, model_class=LinearRegression, **kwargs
) -> Union[LinearModel, None]:
    r"""

    The per-feature negative log-likelihood of the model is computed
    as the concentrated negative log-likelihood

    .. math::
        NLL(\theta) &= - \frac{N}{2} (log(2\pi) + log(\frac{SSR}{N} + 1)

    Parameters
    ----------
    X : np.ndarray
        Data of the explanatory variables. This reflects the substrate data,
        already corrected for all possible covariates.
    Y : np.ndarray
        Data of the target variables. This reflects the product data,
        already corrected for all possible covariates.
    residual_summary: Callable
        Function to compute a residual summary statistic
    model : LinearRegression, optional
        Linear model to use. Any class with a scikit-learn like API
        implementing `fit` and `predict` can be used here.
    substrates : list, optional
        Substrates to subset X by.
    products : list, optional
        Products to subset X by.
    p_tresh : float, default 0.05
        p-value threshold
    p : int, default 2
        lp-norm to use, only relevant when `residual_summary` is the p-norm
        function
    simple_model : LinearRegression, optional
        Simple baseline model to compare the actual model against for p-value
        computation
    simple_data : np.ndarray, optional
        Simple baseline data to compare the actual model against for p-value
        computation
    norm_multiple : bool, default False
        norm multivariable residuals
    outlier_threshold : float | bool, default None
        Cook's distance threshold for outlier detection. If `False`, no outlier
        removal will be performed
    r2_threshold : float, default 0.5
        Minimum :math:`R^2` the model needs to meet to not be
        discarded. If all models should be taken set to None.
    r2_warning: Callable, optional
        Function to raise a warning when the r2-threshold is not passed
    kwargs
        Optional keyword arguments for
        :py:class:`sklearn.linear_model.LinearRegression`

    Returns
    -------
    LinearModel | None
        Fitted linear model and negative log-likelihood as the mean NLL
        over all features. If model r2 restriction is set and not met
        None will be returned.
    """
    if model is None:
        model = model_class(**kwargs)
        model.fit(X, Y)
    # model prediction
    # e.g model.intercept_ + X @ model.coef_.T for a linear regression
    y_hat = model.predict(X)
    # cook's distance for outlier detection
    if isinstance(outlier_threshold, bool) and not outlier_threshold:
        has_outliers = False
    else:
        cook_dist, df = _cooks_distance(Y, y_hat, X, p)
        if outlier_threshold is None:
            # TODO: what's the better threshold here?
            # outlier_threshold = 4 / X.shape[0]
            outlier_threshold = f.sf(cook_dist, X.shape[1], df)
        outliers = cook_dist > outlier_threshold
        has_outliers = np.any(outliers)
    if has_outliers:
        xn = X[~outliers]
        yn = Y[~outliers]
        model.fit(xn, yn)
        y_hat = model.predict(xn)
    else:
        xn = X
        yn = Y
    # checking general model fit
    r2_val = r2_score(yn, y_hat)
    if r2_val < r2_threshold:
        r2_warning(r2_val)
        return
    # per feature nll => 'total' NLL is the mean over all features
    nll = _nll(xn, yn, y_hat)
    # TODO: we need to ensure that all values are above a certain
    #       threshold, since our requirement is that all products are
    #       predictable from the input
    # => compute all individual F-statistics + prob => test
    if simple_model:
        y_simple = simple_model.predict(simple_data[~outliers])
        ssd_simple = np.sum((yn - y_simple) ** 2, axis=0)
        if simple_model.coef_.ndim == 1:
            p_simple = simple_model.coef_.size
        else:
            p_simple = simple_model.coef_.shape[1]
        fstats, pvals = _fstats(xn, yn, y_hat, as_probability=True,
                                ssd_simple=ssd_simple,
                                p_simple=p_simple,
                                norm_multiple=norm_multiple)
        active = np.all(pvals < p_tresh)
    else:
        fstats, pvals = _fstats(xn, yn, y_hat, as_probability=True,
                                norm_multiple=norm_multiple)
        active = np.all(pvals < p_tresh)
    # TODO: additional test: hotelling
    # normed residuals normalized by the number of features
    res = np.zeros(X.shape[0])
    res[~outliers] = residual_summary(yn, y_hat, p=p)
    res[outliers] = np.nan

    return LinearModel(
        substrates, products, res, model, nll, _r_square(xn, yn, y_hat),
        fstats, pvals, active, list(np.where(outliers)[0])
    )


def _ld_models(
    reaction_nodes: Dict[str, Tuple[List[str], List[str]]],
    data: pd.DataFrame, residual_summary: Callable,
    case_data: pd.DataFrame = None, n_threads: int = 1,
    r2_warning: Callable = lambda x: lambda y: None,
    recompute_non_passing: bool = False, **kwargs
) -> Tuple[Dict[str, LinearModel], Dict[str, np.ndarray]]:
    if case_data is not None:
        def _model(
            sub_data, prod_data, reac, subs, prods, model_dict, res, inv=False
        ):
            if len(subs) == len(prods) and np.all(np.isin(subs, prods)):
                return
            lm = _compute_model(
                sub_data, prod_data,
                residual_summary,
                substrates=subs,
                products=prods,
                r2_warning=r2_warning(reac),
                **kwargs
            )
            if not lm:
                if not recompute_non_passing or inv:
                    return
                return _model(
                    case_data.loc[:, subs].values,
                    case_data.loc[:, prods].values,
                    reac, subs, prods, model_dict, res, True
                )
            use_data = data if inv else case_data
            model_dict[reac] = lm
            case_subs = use_data.loc[:, subs].values
            # NOTE: the RegressorMixin base class **always** implements a
            #       `predict` method
            case_pred = lm.model.predict(case_subs)
            normed_res = residual_summary(
                case_pred, use_data.loc[:, prods], p=kwargs.get("p", 2))
            res[reac] = normed_res
    else:
        def _model(sub_data, prod_data, reac, subs, prods, model_dict, _):
            if len(subs) == len(prods) and np.all(np.isin(subs, prods)):
                return
            model = _compute_model(
                sub_data, prod_data,
                residual_summary,
                substrates=subs,
                products=prods,
                r2_warning=r2_warning(reac),
                **kwargs
            )
            if model:
                model_dict[reac] = model
    if n_threads == 1:
        models = {}
        residuals = {}
        for reaction, (substrates, products) in reaction_nodes.items():
            _model(
                data.loc[:, substrates].values,
                data.loc[:, products].values,
                reaction, substrates,
                products, models, residuals
            )
    else:
        # FIXME: setting thread number correctly
        print("parallel option used")
        manager = Manager()
        residuals = manager.dict()
        models = manager.dict()
        # pool = manager.Pool()
        jobs = []
        for reaction, (substrates, products) in reaction_nodes.items():
            substrate_data = data.loc[:, substrates]
            product_data = data.loc[:, products]
            jobs.append(Process(
                target=_model,
                args=(substrate_data, product_data, reaction, substrates,
                      products, models, residuals)
            ))
        _ = []
        # FIXME: correct join operation
    return models, residuals


def _neighbour_models(
    graph: nx.DiGraph, initial_models: Dict[str, LinearModel],
    metabolome_data: pd.DataFrame, **kwargs
) -> Dict[str, LinearModel]:
    # extract reaction-reaction connection
    reaction_graph = reaction_graph_extraction(graph, include_attributes=False)
    # turn into adjacency list
    adjacency_list = {}
    for (src, tgt) in reaction_graph:
        adjacency_list.setdefault(src, set()).add(tgt)
        adjacency_list.setdefault(tgt, set()).add(src)
    neighbour_models = {}
    # compute neighbour-included models
    # TODO: parallelized option
    # TODO: what should we really 'correct' for here:
    #   * metabolites without direction to connection to products but to
    #     substrates
    #   * reaction-chain, i.e. mediating model A -> B -> C
    # TODO: outlier detection
    for reaction, model in initial_models.items():
        reaction_substrates = set(graph.predecessors(reaction))
        reaction_products = list(graph.successors(reaction))
        # NOTE: this cannot raise a KeyError unless there's an
        #       error in the package. In the case of incorrect
        #       node/edge type annotation reaction graph extraction
        #       will already throw an error
        neighbours = adjacency_list.get(reaction)
        # TODO: there have two options here
        #   1. correct data for neighbouring substrate
        #   2. 'joint' linear model and then check for coefficients of reaction
        # also note that we only use reaction substrates, but not the
        # products of adjacent reactions
        adj_substrates = set()
        if neighbours is not None:
            for neighbour in neighbours:
                if initial_models[neighbour].active:
                    adj_substrates.update(set(graph.predecessors(neighbour)))
            substrates = list(adj_substrates.union(reaction_substrates))
            substrate_data = metabolome_data.loc[:, substrates]
            product_data = metabolome_data.loc[:, reaction_products]
            # computing the linear model
            # TODO: pass the initial model (i.e. `model`) and compute
            #       f-statistics and p-value as the residual difference
            #       between the models

            sub_data = metabolome_data.loc[:, list(reaction_substrates)].values
            model = _compute_model(
                substrate_data, product_data,
                substrates=substrates,
                products=reaction_products,
                simple_model=model.model,
                simple_data=sub_data,
                **kwargs
            )
            # storing the results
            neighbour_models[reaction] = model
        else:
            neighbour_models[reaction] = LinearModel(
                [], [], np.array([]), LinearRegression(),
                np.array([]), np.array([]), np.array([]), np.ndarray([]),
                active=model.active, outlier_samples=[]
            )
    return neighbour_models


def confounder_correction(
    metabolome_data: pd.DataFrame, covariates: pd.DataFrame,
    random_effects: Union[str, List[str]] = None, **kwargs
) -> pd.DataFrame:
    """Correct for confounder using linear (mixed effect) models

    Perform a confounder correction by fitting a linear model or linear mixed
    effect model for each metabolite (as the target variable). The residual,
    representing the variance unexplained by the confounding factors, are
    returned as the confounder corrected data.

    For linear mixed effect models, only simple settings without interaction
    between variables are supported. If you have a more complex setup, we
    recommend correcting outside `pymantra` and using the corrected data
    for your analysis.

    Parameters
    ----------
    metabolome_data: pd.DataFrame
        Metabolomics data to correct with samples in rows and features in
        columns
    covariates: pd.DataFrame
        Confounder data with samples in rows and variables in columns
    random_effects: str | List[str]
        Random effect variables. If no random effects should be used (i.e.
        a 'classical' linear model instead of a linear mixed effect model) pass
        None.
    kwargs
        Optional keyword arguments to pass to :py:func:`MixedLM.from_formula`.

    Returns
    -------
    pd.DataFrame
        Confounder corrected data with the same shape, indices and columns as
        `metabolome_data`.
    """
    if random_effects is None:
        corr_model = LinearRegression()
        corr_model.fit(covariates, metabolome_data)
        pred_vals = corr_model.predict(metabolome_data.values)
        corr_residuals = pd.DataFrame(
            metabolome_data.values - pred_vals,
            index=metabolome_data.index, columns=metabolome_data.columns
        )
    else:
        # statsmodel cannot handle variables in formulas, which contain
        # whitespaces. Hence, we change feature names and reset them after
        # correction
        metabolites = metabolome_data.columns
        metabolome_data.columns = [
            f"s{i + 1}" for i in np.arange(metabolome_data.shape[1])]
        # same with covariate names
        covariates.columns = [x.replace(" ", "_") for x in covariates.columns]
        random_effects = [x.replace(" ", "_") for x in random_effects]

        lmm_data = pd.concat([metabolome_data, covariates], axis=1)

        fixed_effects = " + ".join(
            covariates.columns[~np.isin(covariates.columns, random_effects)])

        if isinstance(random_effects, str):
            random_effects = [random_effects]
        # this looks a bit hacky, but is the only way I manged to reproduce the
        # lmer behavior
        if len(random_effects) == 1:
            # single random effect variable
            def lmm(metabolite):
                model = MixedLM.from_formula(
                    f"{metabolite} ~ {fixed_effects}", lmm_data,
                    groups=lmm_data[random_effects[0]], **kwargs
                )
                return model.fit()
        else:
            # multiple independent random effects
            vcf = {
                rand_effect: f"0 + C({rand_effect})"
                for rand_effect in random_effects
            }
            lmm_data["Group"] = 1

            def lmm(metabolite):
                model = MixedLM.from_formula(
                    f"{metabolite} ~ {fixed_effects}", groups="Group",
                    vc_formula=vcf, data=lmm_data, **kwargs
                )
                return model.fit()

        corr_res = np.zeros(metabolome_data.shape)
        # TODO: parallelize and make more efficient
        with warnings.catch_warnings():
            # Parameter is often on the boundary
            warnings.simplefilter("ignore", ConvergenceWarning)
            desc = "Confounder Correction"
            for i, met in enumerate(tqdm(metabolome_data.columns, desc=desc)):
                try:
                    res = lmm(met)
                    corr_res[:, i] = res.resid
                # this is most likely a singular covariance matrix error
                except ValueError as err:
                    tb = sys.exc_info()[2]
                    raise ValueError(
                        f"Error correcting {met}: {err}"
                    ).with_traceback(tb)
        corr_residuals = pd.DataFrame(
            corr_res, index=metabolome_data.index, columns=metabolites)

    return corr_residuals


def _prep_data(
    graph, metabolome_data, covariates, scale, random_effects, **kwargs
):
    reaction_nodes = {}
    for node, ntype in nx.get_node_attributes(graph, 'node_type').items():
        if graph.nodes[node]['node_type'] == 'reaction':
            substrates = [
                sub for sub in graph.predecessors(node)
                if graph.nodes[sub]['node_type'] == 'metabolite'
            ]
            products = [
                prod for prod in graph.successors(node)
                if graph.nodes[prod]['node_type'] == 'metabolite'
            ]
            if substrates and products:
                reaction_nodes[node] = (substrates, products)

    # TODO: check if axis sizes match up
    if covariates is not None:
        corr_residuals = confounder_correction(
            metabolome_data, covariates, random_effects, **kwargs)
    else:
        corr_residuals = metabolome_data

    if scale:
        corr_residuals = pd.DataFrame(
            StandardScaler().fit_transform(corr_residuals.to_numpy()),
            index=corr_residuals.index, columns=corr_residuals.columns
        )

    return corr_residuals, reaction_nodes


def ld_estimation(
    graph: nx.DiGraph, metabolome_data: pd.DataFrame,
    covariates: pd.DataFrame = None, scale: bool = True,
    random_effects: Union[str, List[str]] = None, **kwargs
) -> Tuple[pd.DataFrame, Dict[str, LinearModel], Dict[str, LinearModel]]:
    warnings.warn(
        "'ld_estimation' is currently not supported! Please use "
        "'per_sample_ld_estimation'"
    )
    corr_residuals, reaction_nodes = _prep_data(
        graph, metabolome_data, covariates, scale, random_effects)

    # first step: compute linear models per reaction
    initial_models, _ = _ld_models(reaction_nodes, corr_residuals, **kwargs)
    # second step: including 1-hop neighbours, that are connected via 'active'
    #              reactions
    neighbour_models = _neighbour_models(graph, initial_models, corr_residuals,
                                         **kwargs)
    # extract the reaction values per sample
    initial_values = pd.DataFrame.from_dict(
        {reaction: pd.Series(model.model_fit)
         for reaction, model in initial_models.items()},
        orient="index"
    )
    neighbour_values = pd.DataFrame.from_dict(
        {reaction: pd.Series(model.model_fit)
         for reaction, model in neighbour_models.items()},
        orient="index"
    )
    res_diff = initial_values.abs() - neighbour_values.abs()
    return res_diff, neighbour_models, initial_models


def _scale_1D(x: np.ndarray, params: Tuple[float, float] = None):
    if params is None:
        params = (np.nanmean(x), np.nanstd(x))
        return (x - params[0]) / params[1], params
    return (x - params[0]) / params[1]


def per_sample_ld_estimation(
    graph: nx.DiGraph, metabolome_data: pd.DataFrame, groups: pd.Series,
    covariates: pd.DataFrame = None, compute_expl_var: bool = False,
    var_as_pval: bool = False, combined_models: bool = False,
    residual_summary: str = "expl_var", scale: bool = True,
    control_group: str = None, r2_threshold: float = .5,
    recompute_non_passing: bool = False, outlier_threshold: float = None,
    random_effects: Union[str, List[str]] = None, lmm_args: dict = None,
    verbose: bool = False, **kwargs
) -> Tuple[Dict[str, LinearModel],
           Dict[str, np.ndarray],
           Dict[str, pd.Series]]:
    r"""Compute linear reaction-models

    TODO

    Parameters
    ----------
    graph : nx.Graph
        Reaction-reaction graph
    metabolome_data : pd.DataFrame
        Metabolome data with samples in rows and metabolites in columns
    groups : pd.Series
        Sample group annotation
    covariates : pd.DataFrame, optional
        Confounder variables to correct for. Correction is done using a
        Linear Mixed Model. All variables (i.e. columns) not specified as
        random effect variables in `random_effects` are assumed to be fixed
        effects variables.
        Generally variables should be numerical (float or integer). If you
        have categorical data as strings you can use the `pandas.get_dummies`
        function to encode them as integers. Make sure to use `drop_first` to
        avoid introducing collinearity see
        (https://stackoverflow.com/questions/31498390/how-to-get-pandas-get-dummies-to-emit-n-1-variables-to-avoid-collinearity) # noqa: E501
        The correction currently only supports simple fixed and random effects
        inclusion. For more complex setups including factor interaction, it is
        recommended to do the correction beforehand and only pass the residuals
        to this function instead of the original metabolome data frame.
    compute_expl_var : bool, False
        Whether to return 1 - explained variance of the model or the residuals
    var_as_pval : bool, False
        Whether to return a p-value or a residual/explained variance value
    combined_models : bool, False
        Whether to compute the reference linear model based on both groups or
        only the 'control'
    residual_summary: str, "expl_var"
        Which method to use as residual summary statistic. Either "expl_var"
        for explained variance (RSS/TSS) or "norm" for p-norm
    scale : bool, True
        Whether to z-score scale metabolites
    control_group : str, optional
        Option to set which group should be viewed as the control.
        If None the first element in `groups` will be used.
    r2_threshold : float, .5
        Minimum :math:`R^2` value a method needs to achieve to be further
        considered
    recompute_non_passing: bool, False
        Whether to recompute models with case data that failed to pass the R2
        threshold for control data.
    outlier_threshold : float, optional
        Threshold to remove outliers by Cook's distance. If `None` a default
        on the basis of the survival function of a f-distribution is computed.
    random_effects: str | List[str], optional
        Random effects for confounder correction. If `covariates` is None this
        has no effect. Else, this specifies which column(s) of `covariates` to
        include as random effects, all other columns will be included as fixed
        effects. If this is None, all columns of `covariates` are assumed to be
        fixed effects.
    lmm_args: dict, optional
        Keyword arguments for :py:func:`MixedLM.from_formula`.
        Ignored unless `covariates` and `random_effects` are both not None.
    verbose: bool, False
        If True, warnings will be raised whenever a model does not pass the
        R2 filter
    kwargs
        Optional keyword arguments passed to model computation
        TODO

    Returns
    -------
    Tuple[Dict[str, LinearModel], Dict[str, np.ndarray], Dict[str, pd.Series]]
        Control models, case residuals and all scaled residuals per reaction.
        If `compute_expl_var` is set to True the last element will contain the
        explained variance instead of the residuals.
    """
    if groups.unique().size != 2:
        raise ValueError(
            "Exactly 2 groups must be given!"
        )
    if control_group is None:
        control_group = groups.unique()[0]
    group_mask = groups == control_group

    # data preparation includes
    # * extracting the substrates and products per reaction
    # * scaling
    # * (possibly) confounder correction
    corr_residuals, reaction_nodes = _prep_data(
        graph, metabolome_data, covariates, scale, random_effects,
        **(lmm_args or {})
    )

    if combined_models:
        warnings.warn(
            "Combined group computation currently not implemented")

    if residual_summary == "expl_var":
        res_sum = _explained_variance
    elif residual_summary == "norm":
        res_sum = _normed_residuals
    else:
        raise ValueError(
            "'residual_summary' must be 'expl_var' or 'norm', not "
            f"{residual_summary}"
        )

    # print-function lazy loading to avoid too many if-else repeats
    if verbose:
        def r2_warn_func(reac: str):
            def warn(r2: float):
                warnings.warn(
                    f"{reac} failed to pass R2 with {round(r2, 4)}",
                    category=R2Warning
                )
            return warn
    else:
        def r2_warn_func(reac: str):
            def warn(r2: float):
                return None
            return warn

    # first step: compute linear models per reaction for control samples
    control_models, case_residuals = _ld_models(
        reaction_nodes, corr_residuals.loc[group_mask, :], res_sum,
        case_data=corr_residuals.loc[~group_mask, :],
        r2_threshold=r2_threshold, r2_warning=r2_warn_func,
        recompute_non_passing=recompute_non_passing,
        outlier_threshold=outlier_threshold, **kwargs
    )

    if compute_expl_var:
        if var_as_pval:
            def assign_values(variances, var_arr, rea):
                variances[rea] = _test_explained_variances(
                    var_arr[group_mask], var_arr[~group_mask])
        else:
            def assign_values(variancess, var_arr, rea):
                variancess[rea] = var_arr
        # test case samples with multivariate test
        # => sample-wise
        explained_variances = {}
        for reaction, model in control_models.items():
            yhat = model.model.predict(
                corr_residuals.loc[groups.index, model.substrates])
            expl_var = _explained_variance(
                corr_residuals.loc[groups.index, model.products], yhat)
            assign_values(explained_variances, expl_var, reaction)
        return control_models, case_residuals, explained_variances

    # TODO: continue
    #   * z-scores for residuals
    #     (based on p-norm or is there a multidimensional one?)
    scaled_residuals = {}
    if var_as_pval:
        def transform_residuals(rea, mod):
            ypred = mod.model.predict(
                corr_residuals.loc[groups.index, model.substrates])
            res = res_sum(
                ypred,
                corr_residuals.loc[groups.index, model.products],
                kwargs.get("p", 2)
            )
            scaled_residuals[rea] = _test_explained_variances(
                res[group_mask], res[~group_mask])
    else:
        def transform_residuals(rea, mod):
            scaled_control, control_dist = _scale_1D(mod.model_fit)
            scaled_case = _scale_1D(case_residuals[rea], control_dist)
            scaled_residuals[rea] = pd.Series(
                np.hstack((scaled_control, scaled_case)),
                index=np.hstack((metabolome_data.index[group_mask],
                                 metabolome_data.index[~group_mask])),
            )
    # TODO: parallelized option
    for reaction, model in control_models.items():
        transform_residuals(reaction, model)
    return control_models, case_residuals, scaled_residuals
