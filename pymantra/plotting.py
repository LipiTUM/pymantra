from typing import Union, Dict, List, Tuple
import warnings
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns

from pymantra.statics import (
    NODE_TYPE_LIST, EDGE_TYPE_LIST, DIRECT_EDGE_TYPE_LIST
)


NODE_COLOURS = dict(zip(
    NODE_TYPE_LIST, ["tab:orange", "tab:green", "tab:red", "tab:blue"]
))
NODE_SHAPES = dict(zip(
    NODE_TYPE_LIST, ["^", "s", "v", "o"]
))
EDGE_COLOURS = dict(zip(
    EDGE_TYPE_LIST,
    ["tab:grey", "tab:grey", "tab:grey", "tab:grey",
     "tab:grey", "tab:grey", "tab:grey", "tab:grey"]
))
EDGE_STYLES = dict(zip(
    EDGE_TYPE_LIST,
    ["solid", "solid", "dashed", "solid",
     "solid", "solid", "solid", "solid"]
))

DIRECT_EDGE_COLOURS = dict(zip(
    DIRECT_EDGE_TYPE_LIST,
    ["tab:grey", "tab:grey", "tab:grey", "tab:grey"]
))
DIRECT_EDGE_STYLES = dict(zip(
    DIRECT_EDGE_TYPE_LIST,
    ["solid", "solid", "solid", "dashed"]
))


def _legend(
    node_types: List[str], edge_types: List[str],
    node_colours: Dict[str, str] = None, node_shapes: Dict[str, str] = None,
    edge_colours: Dict[str, str] = None, edge_styles: Dict[str, str] = None,
    reaction_graph: bool = False, **kwargs
) -> List[Line2D]:
    if node_colours is None:
        node_colours = NODE_COLOURS
    if node_shapes is None:
        node_shapes = NODE_SHAPES
    if edge_colours is None:
        if reaction_graph:
            edge_colours = EDGE_COLOURS
        else:
            edge_colours = DIRECT_EDGE_COLOURS
    if edge_styles is None:
        if reaction_graph:
            edge_styles = EDGE_STYLES
        else:
            edge_styles = DIRECT_EDGE_STYLES
    handles = []
    for node_type in node_types:
        markersize = kwargs.pop("markersize", 12)
        handles.append(
            Line2D(
                [0], [0], color="w", marker=node_shapes[node_type],
                markerfacecolor=node_colours[node_type], label=node_type,
                markersize=markersize, **kwargs
            )
        )
    # placeholder between nodes and edges
    handles.append(Line2D([], [], color="none"))
    for edge_type in edge_types:
        lw = kwargs.pop("lw", 4)
        handles.append(
            Line2D(
                [0], [0], color=edge_colours.get(edge_type, 'tab:grey'),
                label=edge_type, linestyle=edge_styles.get(edge_type, 'solid'),
                lw=lw, **kwargs
            )
        )
    return handles


def _plot_graph(
    graph: Union[nx.Graph, nx.DiGraph],
    layout: callable = nx.kamada_kawai_layout,
    directed: bool = False, reaction_graph: bool = False,
    ax: plt.axis = None, reverse_directions: bool = False,
    formula_as_label: bool = False, show_labels: bool = True,
    edge_width: float = 1, **label_args
) -> plt.axis:
    if ax is None:
        ax = plt.subplots(figsize=(16, 9))
    positions = layout(graph)
    # nodes
    nodes_by_type = {}
    for node, node_type in nx.get_node_attributes(graph, "node_type").items():
        nodes_by_type.setdefault(node_type, set()).add(node)
    for node_type, nodes in nodes_by_type.items():
        nx.draw_networkx_nodes(
            graph, positions, nodelist=nodes,
            node_color=NODE_COLOURS[node_type],
            node_shape=NODE_SHAPES[node_type],
            edgecolors="tab:grey", ax=ax
        )
    # edges
    edges_by_type = {}
    for edge, edge_type in nx.get_edge_attributes(graph, "edge_type").items():
        edges_by_type.setdefault(edge_type, set()).add(edge)
    if reaction_graph:
        for edge_type, edges in edges_by_type.items():
            nx.draw_networkx_edges(
                graph, positions, edgelist=edges,
                edge_color=EDGE_COLOURS.get(edge_type, 'tab:grey'),
                style=DIRECT_EDGE_STYLES.get(edge_type, 'solid'),
                arrows=directed, ax=ax, width=edge_width
            )
    else:
        for edge_type, edges in edges_by_type.items():
            nx.draw_networkx_edges(
                graph, positions, edgelist=edges,
                edge_color=DIRECT_EDGE_COLOURS.get(edge_type, 'tab:grey'),
                style=DIRECT_EDGE_STYLES.get(edge_type, 'solid'),
                arrows=directed, ax=ax, width=edge_width
            )

    if show_labels:
        # node labels
        # => "Description" or "Formula" (reaction), "Name" (metabolite),
        #    "nodeLabel" (organism), "nodeLabel" (gene)
        label_key = {
            "reaction": "Formula" if formula_as_label else "Description",
            "metabolite": "Name",
            "organism": "nodeLabel",
            "gene": "nodeLabel"
        }
        labels = {
            node: node_data.get(
                label_key.get(node_data["node_type"], "nodeLabel"), node)
            for node, node_data in graph.nodes(data=True)
        }
        nx.draw_networkx_labels(
            graph, positions, labels=labels, ax=ax, **label_args)

    # legend
    # TODO: option to manually adapt colour scheme
    legend_handles = _legend(
        list(nodes_by_type.keys()), list(edges_by_type.keys()),
        reaction_graph=reaction_graph
    )
    ax.legend(handles=legend_handles)
    # axis formatting
    ax.grid('off')
    ax.axis('off')

    return ax


def plot_directed_graph(
    graph: nx.DiGraph,
    layout: callable = nx.kamada_kawai_layout,
    reaction_graph: bool = True,
    formula_as_reaction_label: bool = False,
    ax: plt.axis = None,
    **label_args
) -> plt.axis:
    """
    Plotting a directed network obtained from the database

    Parameters
    ----------
    graph: nx.DiGraph
        networkx DiGraph object, that has to contain node and edge type for
        each element labelled as <node/edge>_type
    layout: callable, default nx.kamada_kawai_layout
        Function to compute node positions for a :obj:`nx.DiGraph`, that has to
        be compatible with networkx draw_* functions . Default is
        :meth:`nx.kamada_kawai_layout`. If you want to use any of the networkx
        functions with particular parameter settings please use lambda
        functions to set the parameters.
    reaction_graph: bool, default True
        Whether the graph enter contains reaction nodes or not. This should be
        equivalent to whether the input graph was generated with
        :meth:`NetworkGenerator.get_reaction_subgraph` or
        :meth:`NetworkGenerator.get_subgraph`.
    formula_as_reaction_label: bool, False
        If True reaction formulas will be used as node labels, otherwise the
        reaction description as provided in the source database.
    ax: plt.axis, Optional
        Axis to plot the network onto.
    label_args:
        Keyword arguments for `networkx.draw_networkx_labels`

    Examples
    --------
    >>> from pymantra import datasets
    >>> mantra_graph = datasets.example_graph()
    >>> plot_directed_graph(mantra_graph)

    Returns
    -------
    plt.axis
        Axis containing the network plot including legend
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 9))
    return _plot_graph(
        graph, layout, directed=True, reaction_graph=reaction_graph, ax=ax,
        formula_as_label=formula_as_reaction_label, **label_args
    )


def plot_undirected_graph(
    graph: nx.Graph,
    layout: callable = nx.kamada_kawai_layout,
    reaction_graph: bool = False,
    formula_as_reaction_label: bool = False,
    ax: plt.axis = None,
    **label_args
) -> plt.axis:
    """
    Plotting an undirected network obtained from the database

    Parameters
    ----------
    graph: nx.Graph
        networkx Graph object, that has to contain node and edge type for each
        element labelled as node-/edge-type
    layout: callable, default nx.kamada_kawai_layout
        Function to compute node positions for a :obj:`nx.DiGraph`, that has to
        be compatible with networkx draw_* functions . Default is
        :meth:`nx.kamada_kawai_layout`. If you want to use any of the networkx
        functions with particular parameter settings pleas use lambda functions
        to set the parameters.
    reaction_graph: bool, default False
        Whether the graph enter contains reaction nodes or not. This should be
        equivalent to whether the input graph was generated with
        :meth:`NetworkGenerator.get_reaction_subgraph` or
        :meth:`NetworkGenerator.get_subgraph`.
        Generally, It is recommended to use :py:meth:`~plot_directed_graph`
        when plotting a reaction graph.
    formula_as_reaction_label: bool, False
        If True reaction formulas will be used as node labels, otherwise the
        reaction description as provided in the source database.
    ax: plt.axis, Optional
        Axis to plot the network onto.
    label_args:
        Keyword arguments for `networkx.draw_networkx_labels`

    Examples
    --------
    >>> from pymantra import datasets
    >>> mantra_graph = datasets.example_graph()
    >>> plot_undirected_graph(mantra_graph)

    Returns
    -------
    plt.axis
        Axis containing the network plot including legend
    """
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 9))
    return _plot_graph(
        graph, layout, directed=False, reaction_graph=reaction_graph, ax=ax,
        formula_as_label=formula_as_reaction_label, **label_args
    )


def _remove_zero_features(df: pd.DataFrame):
    zero_idxs = df.index[(df == 0).all(axis=1)]
    zero_cols = df.columns[(df == 0).all(axis=0)]
    return df.drop(index=zero_idxs, columns=zero_cols)


def _prep_correlation_dfs(
    corrs: Dict[str, pd.DataFrame], ref_group: str, tgt_group: str,
    pvals: pd.Series, set_zero: bool, thresh: float, strip_column_names: bool
):
    ref_df = corrs[ref_group].copy()
    tgt_df = corrs[tgt_group].copy()
    if set_zero:
        ref_df[pvals[ref_group] > thresh] = 0
        tgt_df[pvals[tgt_group] > thresh] = 0

    if strip_column_names:
        col_map = {col: col.split(", ")[0] for col in ref_df.columns}
        ref_df.rename(columns=col_map, inplace=True)
        tgt_df.rename(columns=col_map, inplace=True)

    return ref_df, tgt_df


def _plot_clustmap(data, cluster, reorder, return_df, ax, **kwargs):
    if cluster:
        if return_df:
            return data, sns.clustermap(data, **kwargs)
        return sns.clustermap(data, **kwargs)
    else:
        if reorder:
            clust = sns.clustermap(data)
            rows = clust.dendrogram_row.reordered_ind
            cols = clust.dendrogram_col.reordered_ind
            if return_df:
                return data, \
                    sns.heatmap(data.iloc[rows, cols], ax=ax, **kwargs)
            return sns.heatmap(data.iloc[rows, cols], ax=ax, **kwargs)
        if return_df:
            return data, sns.heatmap(data, ax=ax, **kwargs)
        return sns.heatmap(data, ax=ax, **kwargs)


def plot_correlation_averages(
    corrs: Dict[str, pd.DataFrame], pvals: Dict[str, pd.DataFrame],
    ref_group: str, tgt_group: str, set_zero: bool = True, thresh: float = .05,
    cluster: bool = False, ax: plt.axis = None, reorder: bool = False,
    strip_column_names: bool = False, return_averages: bool = False,
    remove_all_zeros: bool = False, **kwargs
) -> Union[Union[sns.matrix.ClusterGrid, plt.axis],
           Tuple[pd.DataFrame, Union[sns.matrix.ClusterGrid, plt.axis]]]:
    """Plot the average multi-omics correlation over multiple sample groups

    Parameters
    ----------
    corrs: Dict[str, pd.DataFrame]
        Group-wise correlations as returned by
        `pymantra.compute_multiomics_associations`
    pvals: Dict[str, pd.DataFrame]
        Group-wise correlation p-values as returned by
        `pymantra.compute_multiomics_associations`
    ref_group: str
        Name of the reference group.
    tgt_group: str
        Name of the target group. The correlation values of this group will be
        subtracted from the correlations of `ref_group`
    set_zero: bool, False
        Whether to set correlations with a p-value > some threshold (see
        `thresh`) to zero
    thresh: float, .05
        p-value cutoff above which all correlations will be set to zero. Only
        relevant if `set_zero` is True
    cluster: bool, False
        Whether to cluster features
    ax: plt.axis, optional
        Axis to plot onto. Only used if `cluster` if False
    reorder: bool, False
        Whether to reorder columns and rows by clustering. This essentially
        means plotting a clustermap but leaving away the clustering trees and
        only reordering. Only relevant if `cluster` is False.
    strip_column_names: bool, False
        Whether to strip the column names (i.e. reaction labels) to only retain
        the first reaction annotation
    return_averages: bool, False
        Whether to return the correlation difference data frame
    remove_all_zeros: bool, False
        Whether to remove features with no significant associations. Only
        used if `pvals` is given.
    kwargs
        Keyword arguments to parse to `seaborn.heatmap` (`cluster` is False)
        or `seaborn.clustermap`

    Returns
    -------
    Union[Union[sns.matrix.ClusterGrid, plt.axis],
          Tuple[pd.DataFrame, Union[sns.matrix.ClusterGrid, plt.axis]]]
        A seaborn `ClusterGrid` object, if `cluster` is True, a matplotlib axis
        otherwise

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
    ...     residuals, microbiome_data, sample_groups, comparison=("0", "1"))
    >>> diff_associations, clust_map = plot_correlation_averages(
    ...     corrs, pvals, "0", "1", cluster=True)
    """
    ref_df, tgt_df = _prep_correlation_dfs(
        corrs, ref_group, tgt_group, pvals, set_zero, thresh,
        strip_column_names
    )
    mean_data = pd.DataFrame(
        np.mean([ref_df, tgt_df], axis=0), columns=ref_df.columns,
        index=ref_df.index
    )
    if remove_all_zeros:
        mean_data = _remove_zero_features(mean_data)
    return _plot_clustmap(
        mean_data, cluster, reorder, return_averages, ax,
        cmap=kwargs.pop("cmap", "vlag"), vmin=kwargs.pop("vmin", -1),
        vmax=kwargs.pop("vmax", 1),
        xticklabels=kwargs.pop("xticklabels", True),
        yticklabels=kwargs.pop("yticklabels", True), **kwargs
    )


def plot_correlation_differences(
    corrs: Dict[str, pd.DataFrame], pvals: Dict[str, pd.DataFrame],
    ref_group: str, tgt_group: str, set_zero: bool = True, thresh: float = .05,
    cluster: bool = False, ax: plt.axis = None, reorder: bool = False,
    strip_column_names: bool = False, return_differences: bool = False,
    remove_all_zeros: bool = False, **kwargs
) -> Union[Union[sns.matrix.ClusterGrid, plt.axis],
           Tuple[pd.DataFrame, Union[sns.matrix.ClusterGrid, plt.axis]]]:
    """Plot the multi-omics correlation differences between two sample groups

    Parameters
    ----------
    corrs: Dict[str, pd.DataFrame]
        Group-wise correlations as returned by
        `pymantra.compute_multiomics_associations`
    pvals: Dict[str, pd.DataFrame]
        Group-wise correlation p-values as returned by
        `pymantra.compute_multiomics_associations`
    ref_group: str
        Name of the reference group.
    tgt_group: str
        Name of the target group. The correlation values of this group will be
        subtracted from the correlations of `ref_group`
    set_zero: bool, False
        Whether to set correlations with a p-value > some threshold (see
        `thresh`) to zero
    thresh: float, .05
        p-value cutoff above which all correlations will be set to zero. Only
        relevant if `set_zero` is True
    cluster: bool, False
        Whether to cluster features
    ax: plt.axis, optional
        Axis to plot onto. Only used if `cluster` if False
    reorder: bool, False
        Whether to reorder columns and rows by clustering. Only relevant if
        `cluster` is False,
    strip_column_names: bool, False
        Whether to strip the column names (i.e. reaction labels) to only retain
        the first reaction annotation
    return_differences: bool, False
        Whether to return the correlation difference data frame
    remove_all_zeros: bool, False
        Whether to remove features with no significant associations. Only
        used if `pvals` is given.
    kwargs
        Keyword arguments to parse to `seaborn.heatmap` (`cluster` is False)
        or `seaborn.clustermap`

    Returns
    -------
    Union[Union[sns.matrix.ClusterGrid, plt.axis],
          Tuple[pd.DataFrame, Union[sns.matrix.ClusterGrid, plt.axis]]]
        A seaborn `ClusterGrid` object, if `cluster` is True, a matplotlib axis
        otherwise

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
    ...     residuals, microbiome_data, sample_groups, comparison=("0", "1"))
    >>> diff_associations, clust_map = plot_correlation_differences(
    ...     corrs, pvals, "0", "1", cluster=True, return_differences=True)
    """
    ref_df, tgt_df = _prep_correlation_dfs(
        corrs, ref_group, tgt_group, pvals, set_zero, thresh,
        strip_column_names
    )
    diff_data = ref_df - tgt_df
    if remove_all_zeros:
        diff_data = _remove_zero_features(diff_data)
    return _plot_clustmap(
        diff_data, cluster, reorder, return_differences, ax,
        cmap=kwargs.pop("cmap", "vlag"), vmin=kwargs.pop("vmin", -2),
        vmax=kwargs.pop("vmax", 2),
        xticklabels=kwargs.pop("xticklabels", True),
        yticklabels=kwargs.pop("yticklabels", True), **kwargs
    )


def plot_reaction_association(
    residuals: pd.DataFrame, omics_data: pd.DataFrame,
    corrs: pd.DataFrame, groups: pd.Series = None, top_n: int = 3,
    axes: List[plt.axis] = None, pal: Dict[str, str] = None
):
    """Plot the highest correlating associations between reactions and
    non-metabolomic omics-data

    Parameters
    ----------
    residuals: pd.DataFrame
        Residual values as computed by `pymantra.compute_reaction_estimates`
    omics_data: pd.DataFrame
        Multi-omics data associated with `residuals`. Rows must be features
        (e.g. microbial species or transcript abundances) and columns samples.
    corrs: pd.DataFrame
        Correlations between omics-data (rows) and reaction residuals (columns)
    groups: pd.Series
        Sample groups, index must be the same as the indices in `residuals_`
    top_n: int, 3
        Number of associations to plot. If `axes` is given, len(axes) will be
        used instead.
    axes: List[plt.axis], optional
        (Flat) List of axes to plot onto. If `None`, `top_n` axes will be
        created inside the function. If you want to have a multi-row layout,
        we recommend using `plt.subplots` with the `ncols`/`nrows` parameters
        and then passing the axes with `axes.flatten()`.
    pal: Dict[str, str], optional
        Colour palette to use. Only use if `groups` is given. The keys must be
        the group names that appear in `groups` and the values must be strings
        defining colors in a matplotlib-compatible format.

    Returns
    -------
    List[plt.axis]
        List of axes containing the association plots

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
    ...     residuals, microbiome_data, sample_groups, comparison=("0", "1"))
    >>> diff_associations, clust_map = plot_correlation_differences(
    ...     corrs, pvals, "0", "1", cluster=True, return_differences=True)
    >>> plot_reaction_association(
    ...     residuals, microbiome_data, corr_associations, sample_groups)
    """
    # shape checking
    if residuals.shape[0] != omics_data.shape[1]:
        raise ValueError(
            "'residuals' and 'microbiome_data' need to have matching samples")

    n, m = residuals.shape[1], omics_data.shape[0]
    if n != corrs.shape[1] or m != corrs.shape[0]:
        raise ValueError(
            "'corrs' must have multi-omics data in rows and residuals in"
            f"columns. Expected {m} rows, found {corrs.shape[0]}.\n"
            f"Expected {n} columns, found {corrs.shape[1]}."
        )

    # (N * M) x 2 array giving the indices of correlations sorted by absolute
    # correlation values/differences
    top_corrs = np.dstack(
        np.unravel_index(
            np.argsort(-corrs.abs().to_numpy().ravel()), corrs.shape)
    ).squeeze()

    if axes is None:
        _, axes = plt.subplots(ncols=top_n, figsize=(16, 9))
        if top_n == 1:
            axes = [axes]

    if groups is None:
        # plotting without group coloring
        for i, ax in enumerate(axes):
            yidx, xidx = top_corrs[i, :]
            ax.scatter(residuals.iloc[:, xidx], omics_data.iloc[yidx, :])
            ax.set_xlabel(corrs.columns[xidx])
            ax.set_ylabel(corrs.index[yidx])
    else:
        # plotting with group coloring
        masks = {group: groups == group for group in groups.unique()}
        if pal is None:
            pal = dict(zip(
                masks.keys(), plt.rcParams["axes.prop_cycle"].by_key()['color']
            ))
        for i, ax in enumerate(axes):
            for group, mask in masks.items():
                yidx, xidx = top_corrs[i, :]
                ax.scatter(
                    residuals.loc[mask, :].iloc[:, xidx],
                    omics_data.loc[:, mask].iloc[yidx, :],
                    c=pal.get(group, "tab:grey"), label=group
                )
                ax.set_xlabel(residuals.columns[xidx])
                ax.set_ylabel(omics_data.index[yidx])
        # legend is only called shown on the last plot
        plt.legend()

    return axes


def _residual_to_seaborn_long(
    residuals_: pd.DataFrame, groups: pd.Series, id_vars: List[str] = None,
):
    """Turn a wide residual data frame into a long format"""
    residuals = residuals_.copy()
    residuals["Group"] = groups
    if id_vars is None:
        id_vars = ["Group"]
    return residuals.melt(
        id_vars=id_vars, value_vars=residuals.columns,
        value_name="Residual", var_name="Reaction"
    )


def _move_legend(
    ax: plt.axis, legend_position: Tuple[str, Tuple[float, float]],
    drop_legend: bool = False
):
    """Move the legend of a seaborn plot"""
    if drop_legend:
        ax.legend_.remove()
    elif legend_position is None:
        sns.move_legend(ax, "lower right", bbox_to_anchor=(1.05, 0))
    else:
        sns.move_legend(
            ax, legend_position[0], bbox_to_anchor=legend_position[1])


def _rotate_xticks(
    ax: plt.axis, labels: List[str], rotate_labels: bool
):
    """Format the xticklabels of a reaction plot and optionally rotate"""
    if rotate_labels:
        ax.set_xticks(ax.get_xticks(), labels, rotation=90, ha="center")
    else:
        ax.set_xticks(ax.get_xticks(), labels)


def residual_violinplot(
    residuals_: pd.DataFrame, groups: pd.Series, pvalues: pd.Series = None,
    ax: plt.axis = None, fontsize: int = 12, significance_only: bool = True,
    legend_position: Tuple[str, Tuple[float, float]] = None,
    drop_legend: bool = False, rotate_labels: bool = False,
    plot_significant_features: bool = False, thresh: float = .05, **kwargs
) -> plt.axis:
    """Plot the residual values per group as a violinplot

    Parameters
    ----------
    residuals_: pd.DataFrame
        Residual values as computed by `pymantra.compute_reaction_estimates`
    groups: pd.Series
        Sample groups, index must be the same as the indices in `residuals_`
    pvalues: pd.Series, optional
        p-values to be added to the plot
    ax: plt.axis, optional
        Axis to plot onto. If `None`, a new figure will be initialised inside
        the function
    fontsize: int, 7
        Font size for p-value annotations
    significance_only: bool, True
        If True and `pvalues` is not None, the p-values will be shown as
        asterisks or 'ns' indicating different levels of significance. If False
        rounded p-values are displayed.
        If `pvalues` is None, this has no effect.
    legend_position: Tuple[str, Tuple[float, float]], optional
        Where to position the legend. First element is a string compatible with
        the typical matplotlib legend location names, the second is a 2-tuple
        giving the x- and y-axis offset (via `bbox_to_anchor`)
    drop_legend: bool, False
        If set to True legend drawing will be suppressed
    rotate_labels: bool, False
        Whether to rotate the x-axis ticklabels by 90 degrees or not
    plot_significant_features: bool, False
        Whether to show only significant features. Significance is defined via
        the `thresh` parameter. Can only be used if `pvalues` is given.
    thresh: float, .05
        p-value threshold used to subset if `plot_significant_features` is True
    kwargs
        Optional keyword arguments to `seaborn.violinplot`

    Returns
    -------
    plt.axis
        Matplotlib axis used for plotting

    Examples
    --------
    >>> from pymantra.datasets import example_multiomics_enrichment_data
    >>> from pymantra import (
    ...     compute_reaction_estimates, compute_multiomics_associations)
    >>> metabolite_data, microbiome_data, sample_groups, graph = \
    ...     example_multiomics_enrichment_data()
    >>> residuals = \
    ...     compute_reaction_estimates(graph, metabolite_data, sample_groups)
    >>> residual_violinplot(residuals, sample_groups)
    """
    # refactoring residuals data into long-format with group annotation
    sig_features = residuals_.columns
    if plot_significant_features:
        if pvalues is None:
            warnings.warn(
                "'plot_significant' can only be used when 'pvalues' is set")
        else:
            sig_features = pvalues.index[pvalues < thresh]
            if sig_features.empty:
                warnings.warn(
                    "No significant features detected, with threshold "
                    f"{thresh}. Defaulting to top-10 lowest p-values."
                )
                sig_features = pvalues.index[np.argsort(pvalues)[:10]]
    res_long = _residual_to_seaborn_long(
        residuals_.loc[:, sig_features], groups)

    # actual plotting
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 9))
    sns.violinplot(
        data=res_long, x="Reaction", y="Residual", hue="Group", ax=ax,
        inner=kwargs.pop("inner", "quart"), split=kwargs.pop("split", True),
        saturation=kwargs.pop("saturation", 1), **kwargs
    )

    # removing legend or moving outside plot to see all significance
    # annotations
    _move_legend(ax, legend_position, drop_legend)

    # get the maximum value on the y-axis to locate the p-value annotation
    _, ymax = ax.get_ylim()
    ymax -= ymax / 10

    # p-value annotation
    if pvalues is not None:
        if significance_only:
            def get_annotation(reac):
                pval = pvalues[reac]
                if pval < 0.05:
                    if pval < 0.01:
                        if pval < 0.005:
                            if pval < 0.001:
                                return "****"
                            else:
                                return "***"
                        else:
                            return "**"
                    else:
                        return "*"
                else:
                    return "ns"
        else:
            def get_annotation(reac):
                return str(round(pvalues[reac], ndigits=3))

        # annotate the p-values on top
        for ticklabel in ax.get_xticklabels():
            reaction = ticklabel.get_text()
            xpos = ticklabel.get_position()[0]
            ax.text(
                xpos, ymax, get_annotation(reaction), ha="center", va="bottom",
                fontsize=fontsize
            )

    # xtick rotation while cutting merged reactions into one node label for
    # readability
    _rotate_xticks(ax, sig_features, rotate_labels)

    return ax
