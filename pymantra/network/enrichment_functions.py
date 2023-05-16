from warnings import warn
from typing import Set, Dict, List, Union
import pandas as pd
import networkx as nx

from pymantra.statics import Edge
from pymantra.network.ld_estimation import (
    per_sample_ld_estimation, confounder_correction)
from pymantra.network.reaction_associations import associate_multiomics_ld


def compute_reaction_estimates(
    graph: nx.DiGraph, metabolite_data: pd.DataFrame, sample_groups: pd.Series,
    covariates: pd.DataFrame = None,
    random_effects: Union[str, List[str]] = None, lmm_args: dict = None,
    residual_summary: str = "expl_var", return_all: bool = False,
    control_group: any = None, **kwargs
):
    """Generate reaction estimates

    Compute the linear-model estimates for a given graph and metabolomics data

    Parameters
    ----------
    graph : nx.DiGraph
        Metabolite-reaction graph. Metabolites need to be denoted as
        'metabolite' (via node attribute 'node_type') and reactions as
        'reaction'.
    metabolite_data : pd.DataFrame
        Metabolite data with samples in rows and metabolites in columns.
        Metabolite names need to match the metabolite node names in `graph` and
        indices need to match the indices of `sample_groups`.
    sample_groups : pd.Series
        Array indicating sample groups
    covariates : pd.DataFrame, optional
        Confounder variables to correct for. Correction is done using a
        Linear Mixed Model. All variables (i.e. columns) not specified as
        random effect variables in `random_effects` are assumed to be fixed
        effects variables.
        The correction currently only supports simple fixed and random effects
        inclusion. For more complex setups including factor interaction, it is
        recommended to do the correction beforehand and only pass the residuals
        to this function instead of the original metabolome data frame.
    random_effects: str | List[str], optional
        Random effects for confounder correction. If `covariates` is None this
        has no effect. Else, this specifies which column(s) of `covariates` to
        include as random effects, all other columns will be included as fixed
        effects. If this is None, all columns of `covariates` are assumed to be
        fixed effects.
    lmm_args: dict, optional
        Keyword arguments for :py:func:`MixedLM.from_formula`.
        Ignored unless `covariates` and `random_effects` are both not None.
    residual_summary: str, "expl_var"
        Which method to use as residual summary statistic. Either "expl_var"
        for explained variance (RSS/TSS) or "norm" for p-norm
    return_all : bool, False
        Whether to return all variables return by
        :py:func:`per_sample_ld_estimation` or only return the scaled residuals
    control_group : any, optional
        Name of the control group
    kwargs
        Keyword arguments. See :py:func:`per_sample_ld_estimation` for details.

    Returns
    -------
    Union[
        pd.DataFrame,
        Tuple[Dict[str, LinearModel], Dict[str, np.ndarray], pd.DataFrame]
    ]
        If `control_groups` is False a pd.DataFrame with samples as rows and
        reactions in columns. Else a 3-tuple as returned by
        :py:func:`per_sample_ld_estimation` only with `scaled_residuals` as
        a pd.DataFrame generated from the initially returned dictionary

    Examples
    --------
    >>> from pymantra.datasets import example_metabolome_enrichment_data
    >>> metabolite_data, sample_groups, graph = \
    ...     example_metabolome_enrichment_data()
    >>> compute_reaction_estimates(graph, metabolite_data, sample_groups)
    """
    if control_group is None:
        control_group = sample_groups[0]

    control_models, case_residuals, scaled_residuals = \
        per_sample_ld_estimation(
            graph, metabolite_data, sample_groups, control_group=control_group,
            covariates=covariates, random_effects=random_effects,
            lmm_args=lmm_args, residual_summary=residual_summary, **kwargs
        )

    res_df = pd.DataFrame.from_dict(scaled_residuals)
    res_df.fillna(value=0, inplace=True)

    if return_all:
        return control_models, case_residuals, res_df
    return res_df


def add_reaction_estimates(
    graph: nx.DiGraph, sample_groups: pd.Series = None,
    estimate_data: pd.DataFrame = None, metabolite_data: pd.DataFrame = None,
    control_group: any = None, return_estimates: bool = True, **kwargs
):
    """Add reaction estimates to a metabolite-reaction graph

    Add the linear model estimates to a given metabolite-reaction graph, either
    using pre-computed estimates or computing estimates via
    :py:func:`compute_reaction_estimates` and adding them directly.

    Parameters
    ----------
    graph : nx.DiGraph
        Metabolite-reaction graph. Metabolites need to be denoted as
        'metabolite' (via node attribute 'node_type') and reactions as
        'reaction'.
    sample_groups : pd.Series
        Array indicating sample group
    estimate_data : pd.DataFrame, optional
        Linear-model estimates per reaction and sample as generated by
        :py:func:`compute_reaction_estimates`. If None, estimates will be
        computed from `metabolite_data`, hence must be given.
    metabolite_data : pd.DataFrame, optional
        Metabolite data with samples in rows and metabolites in columns.
        If `estimate_data` is None this becomes a required parameter as it will
        be used to compute the reaction models.
        Metabolite names need to match the metabolite node names in `graph` and
        indices need to match the indices of `sample_groups`.
    control_group : any, optional
        Name of the control group
    return_estimates: bool, False
        Whether to return the linear model-base estimates computed in
        :py:func:`compute_reaction_estimates`
    kwargs
        Keyword arguments. See :py:func:`per_sample_ld_estimation` for details

    Examples
    --------
    >>> from pymantra.datasets import example_metabolome_enrichment_data
    >>> metabolite_data, sample_groups, graph = \
    ...     example_metabolome_enrichment_data()
    >>> residuals = \
    ...     compute_reaction_estimates(graph, metabolite_data, sample_groups)
    >>> add_reaction_estimates(graph, sample_groups, residuals)
    """
    if control_group is None:
        control_group = sample_groups.iloc[0]
    control_mask = sample_groups == control_group

    if estimate_data is None:
        if metabolite_data is None:
            raise ValueError(
                "'metabolite_data' must be given when 'estimate_data' is "
                "None"
            )
        estimate_data = compute_reaction_estimates(
            graph, metabolite_data, sample_groups, control_group=control_group,
            **kwargs
        )
    else:
        return_estimates = False

    control_data = estimate_data.loc[control_mask, :].to_dict('list')
    case_data = estimate_data.loc[~control_mask, :].to_dict('list')

    node_data = {
        reaction: {'0': control_data[reaction], '1': case_data[reaction]}
        for reaction in control_data.keys()
    }
    nx.set_node_attributes(graph, node_data, 'vec_group_data')

    if return_estimates:
        return estimate_data


def compute_multiomics_associations(
    residuals: pd.DataFrame, multi_omics: pd.DataFrame,
    sample_groups: pd.Series, covariates: pd.DataFrame = None,
    random_effects: Union[str, List[str]] = None,
    lmm_args: dict = None, **kwargs
):
    """Compute multi-omics associations with reaction estimates

    Compute the associations between multi-omics features and the residuals of
    reactions as estimated by the linear models.

    This is essentially a wrapper for :py:func:`associate_multiomics_ld` for
    interface consistency.

    Parameters
    ----------
    residuals : pd.DataFrame
        Linear model residuals matrix with samples in rows and reactions in
        columns
    multi_omics : pd.DataFrame
        Multi-omics measurements with samples in rows and multi-omics features
        in columns.
    sample_groups : pd.Series
        Array of sample groups
    covariates : pd.DataFrame, optional
        Confounder variables to correct for. Correction is done using a
        Linear Mixed Model. All variables (i.e. columns) not specified as
        random effect variables in `random_effects` are assumed to be fixed
        effects variables.
        The correction currently only supports simple fixed and random effects
        inclusion. For more complex setups including factor interaction, it is
        recommended to do the correction beforehand and only pass the residuals
        to this function instead of the original metabolome data frame.
    random_effects: str | List[str], optional
        Random effects for confounder correction. If `covariates` is None this
        has no effect. Else, this specifies which column(s) of `covariates` to
        include as random effects, all other columns will be included as fixed
        effects. If this is None, all columns of `covariates` are assumed to be
        fixed effects.
    lmm_args: dict, optional
        Keyword arguments for
        `statsmodels.regression.mixed_linear_model.MixedLM.from_formula`.
        Ignored unless `covariates` and `random_effects` are both not None.
    kwargs
        Keyword arguments passed to :py:func:`associate_multiomics_ld`

    Returns
    -------
    Tuple[Dict[str, pd.DataFrame], Dict[str, pd.DataFrame]]
        2-tuple where the first elements are correlations per group as a
        data frame of multi-omics x reaction and the second element are the
        correlation p-values per group in the same format as the correlations

    Examples
    --------
    >>> from pymantra.datasets import example_multiomics_enrichment_data
    >>> from pymantra import (
    ...     compute_reaction_estimates, compute_multiomics_associations)
    >>> metabolite_data, microbiome_data, sample_groups, graph = \
    ...     example_multiomics_enrichment_data()
    >>> residuals = \
    ...     compute_reaction_estimates(graph, metabolite_data, sample_groups)
    >>> compute_multiomics_associations(
    ...     residuals, microbiome_data, sample_groups)
    """
    if covariates is not None:
        corr_mo = confounder_correction(
            multi_omics, covariates, random_effects, **lmm_args
        )
        return associate_multiomics_ld(
            residuals, sample_groups, corr_mo, **kwargs)
    return associate_multiomics_ld(
        residuals, sample_groups, multi_omics, **kwargs)


def _add_multiomics_associations(
    graph: nx.DiGraph, associations: Dict[str, pd.DataFrame],
    sample_groups: pd.Series, edge_type: str, reaction_data: pd.DataFrame,
    multi_omics: pd.DataFrame, **kwargs
):
    if associations is None:
        associations = compute_multiomics_associations(
            reaction_data, sample_groups, multi_omics, **kwargs)
    corr_dict = {
        group: associations[group].to_dict()
        for group in sample_groups.unique()
    }

    edge_data = {
        edge: {
            group: corr_dict[group].get(edge[0], {}).get(edge[1], 0)
            for group in corr_dict.keys()
        }
        for edge in graph.edges if
        graph.edges[edge]['edge_type'] == edge_type
    }
    nx.set_edge_attributes(graph, edge_data, 'data')


def add_microbiome_associations(
    graph: nx.DiGraph, sample_groups: pd.Series,
    associations: Dict[str, pd.DataFrame] = None,
    residuals: pd.DataFrame = None, microbiome_data: pd.DataFrame = None,
    **kwargs
):
    """Add microbiome-reaction associations to a multi-omics graph

    Add the association estimates to a given multi-omics-reaction graph, either
    using pre-computed estimates or computing estimates via
    :py:func:`compute_multiomics_associations` and adding them directly.

    Parameters
    ----------
    graph : nx.DiGraph
        Metabolite-reaction graph containing additional reaction-organism
        connections. Usually when calling this function reaction estimates
        are already added to the graph.
    sample_groups : pd.Series
        Array indicating sample group
    associations : Dict[str, pd.DataFrame], optional
        Reaction-microbiome associations per group as generated by
        :py:func:`compute_multiomics_associations`
    residuals : pd.DataFrame
        Linear-model estimates per reaction and sample as generated by
        :py:func:`compute_reaction_estimates`. If `associations` is None
        this parameter is required.
    microbiome_data : pd.DataFrame
        Microbiome data with samples in rows and microbes in columns.
        If `associations` is None this becomes a required parameter.
        Microbe names need to match the organism node names in `graph` and
        indices need to match the indices of `sample_groups`.
    kwargs
        Keyword arguments to be passed to
        :py:func:`compute_multiomics_associations`

    Examples
    --------
    >>> from pymantra.datasets import example_multiomics_enrichment_data
    >>> from pymantra import (
    ...     compute_reaction_estimates, compute_multiomics_associations)
    >>> metabolite_data, microbiome_data, sample_groups, graph = \
    ...     example_multiomics_enrichment_data()
    >>> residuals = \
    ...     compute_reaction_estimates(graph, metabolite_data, sample_groups)
    >>> corrs, pvals = compute_multiomics_associations(
    ...     residuals, microbiome_data, sample_groups)
    >>> add_microbiome_associations(graph, sample_groups, corrs)
    """
    if associations is None and (microbiome_data is None or residuals is None):
        raise ValueError(
            "Either 'associations' or 'microbiome_data' and 'residuals' is "
            "required. See function documentation for more information."
        )
    _add_multiomics_associations(
        graph, associations, sample_groups, "REACTION_ORGANISM", residuals,
        microbiome_data, **kwargs
    )


def add_gene_associations(
    graph: nx.DiGraph, sample_groups: pd.Series,
    associations: Dict[str, pd.DataFrame] = None,
    residuals: pd.DataFrame = None, gene_data: pd.DataFrame = None, **kwargs
):
    """Add microbiome-reaction associations to a multi-omics graph

    Add the association estimates to a given multi-omics-reaction graph, either
    using pre-computed estimates or computing estimates via
    :py:func:`compute_multiomics_associations` and adding them directly.

    Parameters
    ----------
    graph : nx.DiGraph
        Metabolite-reaction graph containing additional reaction-organism
        connections. Usually when calling this function reaction estimates
        are already added to the graph.
    sample_groups : pd.Series
        Array indicating sample group
    associations : Dict[str, pd.DataFrame], optional
        Reaction-microbiome associations per group as generated by
        :py:func:`compute_multiomics_associations`
    residuals : pd.DataFrame
        Linear-model estimates per reaction and sample as generated by
        :py:func:`compute_reaction_estimates`. If `associations` is None
        this parameter is required.
    gene_data : pd.DataFrame
        Gene data with samples in rows and genes in columns.
        If `associations` is None this becomes a required parameter.
        Gene names need to match the gene node names in `graph` and
        indices need to match the indices of `sample_groups`.
    kwargs
        Keyword arguments to be passed to
        :py:func:`compute_multiomics_associations`
    """
    if associations is None and (residuals is None or gene_data is None):
        raise ValueError(
            "Either 'associations' or 'gene_data' is required. See function "
            "documentation for more information."
        )
    _add_multiomics_associations(
        graph, associations, sample_groups, "REACTION_GENE",
        residuals, gene_data, **kwargs
    )


def enrichment_pvalue(graph: nx.Graph, subnet: Set[Edge]) -> float:
    # to cpp?
    warn("pvalue computation currently not implemented")
    return -1.
