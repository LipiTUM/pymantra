from string import ascii_lowercase
from typing import Tuple, List, Iterable
import numpy as np
import pandas as pd
# from umap import UMAP
import networkx as nx
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import seaborn as sns
from sklearn.decomposition import PCA
from scipy.spatial.distance import pdist
from skbio.stats.ordination import OrdinationResults, pcoa

from pymantra.plotting import (
    _residual_to_seaborn_long, NODE_SHAPES, NODE_COLOURS)


def _plot_embedding(
    data, groups, pcs, rel_var_ex, ax, file, embed_name, show=False,
    show_legend=False, cmap=None, **kwargs
):
    if ax is None:
        fig, ax = plt.subplots(figsize=(16, 9))
    if groups is None:
        ax.scatter(
            data[:, pcs[0]],
            data[:, pcs[1]],
        )
    else:
        unique_groups = np.unique(groups)

        if cmap is None:
            if unique_groups.size < 11:
                cmap_name = "tab10"
            else:
                cmap_name = "tab20"
            cols = [
                plt.get_cmap(cmap_name)(i)
                for i in np.arange(unique_groups.size)
            ]
            cmap = dict(zip(unique_groups, cols))

        for i, group in enumerate(unique_groups):
            mask = groups == group
            ax.scatter(
                data[mask, pcs[0] - 1], data[mask, pcs[1] - 1], c=cmap[group],
                label=group
            )
        if show_legend:
            ax.legend()
    if rel_var_ex is None:
        ax.set_xlabel(f"{embed_name}{pcs[0]}")
        ax.set_ylabel(f"{embed_name}{pcs[1]}")
    else:
        ax.set_xlabel(
            f"{embed_name}{pcs[0]} "
            f"({round(rel_var_ex[pcs[0] - 1], ndigits=4) * 100}%)"
        )
        ax.set_ylabel(
            f"{embed_name}{pcs[1]} "
            f"({round(rel_var_ex[pcs[1] - 1], ndigits=4) * 100}%)"
        )
    if file is None:
        if show:
            plt.show()
    else:
        plt.savefig(file, **kwargs)


def mds(
    data: np.ndarray,
    dist: str = "braycurtis",
    return_plot: bool = False,
    plot_as_svg: bool = True,
    **kwargs
) -> OrdinationResults:
    """Compute a MDS

    Parameters
    ----------
    data: np.ndarray
        Data with samples in columns
    dist: str, "braycurtis"
        Distance to base the MDS on
    return_plot: bool, False
        Whether to return the plot
    plot_as_svg: bool, False
        Whether to return the plot as svg (or png). Only relevant if
        `return_plot` is True
    kwargs
        Optional keyword arguments to parse to `skbio.stats.ordination.pcoa`
    """
    # generate distance matrix
    dist_mat = pdist(data.T, dist)
    # compute pcoa from distances
    pcoa_res = pcoa(dist_mat, **kwargs)
    if return_plot:
        if plot_as_svg:
            return pcoa_res.svg
        return pcoa_res.png
    return pcoa_res


def plot_mds(
    data: OrdinationResults,
    pcs: Tuple[int, int] = (1, 2),
    sample_groups: np.ndarray = None,
    ax: plt.axis = None,
    file: str = None,
    show_legend: bool = True,
    **kwargs
):
    """Plot a MDS

    Parameters
    ----------
    data: OrdinationResults
        Result of an MDS computation as computed by `mds`
    pcs: Tuple[int, int], (1, 2)
        Which latent dimensions to plot
    sample_groups: np.ndarray, optional
        Sample group annotations
    ax: plt.axis, optional
        Axis to plot onto. If not given a new figure will be initialised.
    file: str, optional
        File path to save the plot to. If none, `plt.show` will be called
    show_legend: bool, True
        Whether to draw the color legend. Only used if `sample_groups` is not
        None
    kwargs
        Optional keyword arguments for `plt.savefig`. Only relevant if `file`
        is given
    """
    _plot_embedding(
        data.samples.values, sample_groups, pcs, data.proportion_explained,
        ax, file, "PCo", show_legend=show_legend, **kwargs
    )


def plot_pca(
    data: np.ndarray,
    pcs: Tuple[int, int] = (1, 2),
    sample_groups: np.ndarray = None,
    ax: plt.axis = None,
    file: str = None,
    show_legend: bool = False,
    **kwargs
):
    """Plot a PCA

    Parameters
    ----------
    data: np.ndarray
        Sample data to reduce via PCA. Samples must be in columns.
    pcs: Tuple[int, int], (1, 2)
        Which latent dimensions to plot
    sample_groups: np.ndarray, optional
        Sample group annotations
    ax: plt.axis, optional
        Axis to plot onto. If not given a new figure will be initialised.
    file: str, optional
        File path to save the plot to. If none, `plt.show` will be called
    show_legend: bool, True
        Whether to draw the color legend. Only used if `sample_groups` is not
        None
    kwargs
        Optional keyword arguments for `plt.savefig`. Only relevant if `file`
        is given
    """
    prcomp = PCA()
    embedding = prcomp.fit_transform(data.T)
    _plot_embedding(
       embedding, sample_groups, pcs, prcomp.explained_variance_ratio_,
       ax, file, "PC", show_legend=show_legend, **kwargs
    )


def reaction_jitterplot(
    residuals_: pd.DataFrame, groups: pd.Series, subtype: pd.Series,
    figsize: Tuple[int, int] = (16, 9), **kwargs
) -> sns.FacetGrid:
    """Plot the residual values per reaction as a jitterplot with subtype
    annotation

    Parameters
    ----------
    residuals_: pd.DataFrame
        Residual values as computed by `pymantra.compute_reaction_estimates`
    groups: pd.Series
        Sample groups, index must be the same as the indices in `residuals_`
    subtype: pd.Series
        Subtype annotation for each sample
    figsize: Tuple[int, int], (16, 9)
        Figure size in inches
    kwargs
        Optional keyword arguments to `seaborn.violinplot`
    """
    # refactoring residuals data into long-format with group annotation
    res_df = residuals_.copy()
    res_df["Subtype"] = subtype
    res_long = _residual_to_seaborn_long(
        res_df, groups, ["Group", "Subtype"])

    # actual plotting
    catp = sns.catplot(
        data=res_long, x="Group", y="Residual", hue="Subtype", col="Reaction",
        **kwargs
    )
    catp.set_titles("{col_name}")
    catp.figure.set_size_inches(*figsize)

    return catp


def residual_pca(
    metabolite_data: pd.DataFrame, residuals: pd.DataFrame, groups: pd.Series,
    file: str = None, plot: Tuple[plt.Figure, List[plt.axis]] = None,
    pca_kwargs: dict = None, **kwargs
):
    """Plot the comparison between original metabolite and the residual data

    Parameters
    ----------
    metabolite_data: pd.DataFrame
        Metabolite data with samples in columns
    residuals: pd.DataFrame
        Residuals as computed in `pymantra.compute_reaction_estimates`
    groups: pd.Series
        Sample group annotations
    file: str
        File path to save the plot in
    plot: Tuple[plt.Figure, List[plt.axis]], optional
        Matplotlib figure and subplot axes to plot onto
    pca_kwargs: dict, optional
        PCA plotting keyword arguments
    kwargs
        Optional keyword arguments to parse to `plt.subplots`
    """
    if plot is None:
        fig, axes = plt.subplots(ncols=2, **kwargs)
    else:
        fig, axes = plot
    if pca_kwargs is None:
        pca_kwargs = {}
    plot_pca(
        metabolite_data, sample_groups=groups[metabolite_data.columns],
        ax=axes[0], show_legend=False, **pca_kwargs
    )
    plot_pca(
        residuals.T, sample_groups=groups[residuals.index], ax=axes[1],
        **pca_kwargs
    )
    plt.legend()
    plt.tight_layout()
    if file is not None:
        fig.savefig(file)
        plt.close(fig)


def paper_theme():
    # TODO: matplotlib theme, especially regarding default fontsizes,
    #       linewidths and padding
    pass


def metabolome_layout():
    """Generate the layout for the metabolomics base figures"""
    fig = plt.figure(figsize=(24, 18))
    axes = np.array(
        [
            [
                # a
                plt.subplot2grid((8, 2), (0, 0), colspan=2, rowspan=3),
                # d
                plt.subplot2grid((8, 2), (5, 0), colspan=2, rowspan=3)
            ],
            [
                # b
                plt.subplot2grid((8, 2), (3, 0), rowspan=2),
                # c
                plt.subplot2grid((8, 2), (3, 1), rowspan=2),
            ]
        ]
    )
    return fig, axes


def annotate_subplots(*axes, fontsize=30):
    """Helper function to annotate subplots with lower case letters"""
    for ax, label in zip(axes, ascii_lowercase):
        ax.text(
            0., .95, label, transform=ax.transAxes, fontsize=fontsize,
            verticalalignment="top", horizontalalignment="left",
            fontweight="bold",
            bbox=dict(facecolor="none", edgecolor="none", pad=0)
        )


def _draw_labels(graph, positions, ax, formula_as_label):
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
    nx.draw_networkx_labels(graph, positions, labels=labels, ax=ax)


def _edge_style_legend(edge_styles: dict):
    return [
        plt.Line2D(
            [0, 1], [0, 1], linestyle=style, color="k", label=edge_type)
        for edge_type, style in edge_styles.items()
    ]


def _shape_legend(node_types: Iterable, uniform_color: bool = True):
    if uniform_color:
        return [
            plt.Line2D(
                [0], [0], color="w", marker=NODE_SHAPES[node_type],
                markerfacecolor="tab:grey", label=node_type,
                markersize=21
            )
            for node_type in node_types
        ]
    return [
        plt.Line2D(
            [0], [0], color="w", marker=NODE_SHAPES[node_type],
            markerfacecolor=NODE_COLOURS[node_type], label=node_type,
            markersize=21
        )
        for node_type in node_types
    ]


def _add_colorbar(cmap: str, vmin: float, vmax: float, ax):
    scaled_map = plt.cm.ScalarMappable(
        cmap=cmap, norm=plt.Normalize(vmin=vmin, vmax=vmax))
    scaled_map._A = []

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="2.5%", pad=0.1)
    plt.colorbar(scaled_map, cax=cax)


def plot_graph_by_degree(
    graph: nx.Graph, ax: plt.axis = None, cmap: str = None,
    edge_width: float = 1, show_labels: bool = True,
    formula_as_label: bool = True, **node_args
):
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 9))
    positions = nx.kamada_kawai_layout(graph)

    # workaround for networkx's super laborious degree interface
    degree_list = [graph.degree[n] for n in graph.nodes]
    vmin = np.min(degree_list)
    vmax = np.max(degree_list)

    # nodes by type with degree colouring
    nodes_by_type = {}
    for node, node_type in nx.get_node_attributes(graph, "node_type").items():
        nodes_by_type.setdefault(node_type, set()).add(node)
    for node_type, nodes in nodes_by_type.items():
        nx.draw_networkx_nodes(
            graph, positions, nodelist=nodes,
            node_color=[graph.degree(node) for node in nodes],
            node_shape=NODE_SHAPES[node_type], edgecolors="tab:grey",
            ax=ax, cmap=cmap, vmin=vmin, vmax=vmax, **node_args
        )

    # edges
    edges_by_type = {}
    for edge, edge_type in nx.get_edge_attributes(graph, "edge_type").items():
        edges_by_type.setdefault(edge_type, set()).add(edge)
    for edge_type, edges in edges_by_type.items():
        nx.draw_networkx_edges(
            graph, positions, edgelist=edges,
            edge_color='tab:grey',
            style='solid',
            arrows=isinstance(graph, nx.DiGraph),
            ax=ax, width=edge_width
        )

    # node labels
    if show_labels:
        _draw_labels(graph, positions, ax, formula_as_label)

    # legend and colorbar
    _add_colorbar(cmap, vmin, vmax, ax)
    ax.legend(handles=_shape_legend(nodes_by_type.keys()))
    ax.grid("off")
    ax.axis("off")

    return ax


def plot_corr_db_graph(
    g: nx.Graph, pos: dict = None, cmap: str = "vlag", ax: plt.axis = None,
    edge_kwargs: dict = None, indicate_egde_source: bool = False
):
    if ax is None:
        _, ax = plt.subplots(figsize=(16, 9))
    if pos is None:
        pos = nx.kamada_kawai_layout(g)

    # nodes
    reaction_nodes = []
    microbe_nodes = []
    for node, node_type in nx.get_node_attributes(g, "node_type").items():
        if node_type == "reaction":
            reaction_nodes.append(node)
        else:
            microbe_nodes.append(node)
    nx.draw_networkx_nodes(
        g, pos, nodelist=reaction_nodes, node_color=NODE_COLOURS["reaction"],
        node_shape=NODE_SHAPES["reaction"], ax=ax
    )
    nx.draw_networkx_nodes(
        g, pos, nodelist=microbe_nodes, node_color=NODE_COLOURS["organism"],
        node_shape=NODE_SHAPES["organism"], ax=ax
    )

    # edges
    if indicate_egde_source:
        db_edges = []
        db_weights = []
        non_data_edges = []
        corr_edges = []
        corr_weights = []
        for edge, edge_type in nx.get_edge_attributes(g, "source").items():
            weight = g.edges[edge].get("corr")
            if edge_type == "db":
                if weight is None:
                    non_data_edges.append(edge)
                else:
                    db_edges.append(edge)
                    db_weights.append(weight)
            else:
                corr_edges.append(edge)
                corr_weights.append(weight)

        vmin = np.nanmin(db_weights + corr_weights)
        vmax = np.nanmax(db_weights + corr_weights)

        nx.draw_networkx_edges(
            g, pos, edgelist=db_edges, edge_color=db_weights,
            edge_cmap=plt.get_cmap(cmap), edge_vmax=vmax, edge_vmin=vmin,
            ax=ax, **edge_kwargs if edge_kwargs is not None else {}
        )
        nx.draw_networkx_edges(
            g, pos, edgelist=corr_edges, edge_color=corr_weights,
            edge_cmap=plt.get_cmap(cmap), edge_vmax=vmax, edge_vmin=vmin,
            style="--", ax=ax, **edge_kwargs if edge_kwargs is not None else {}
        )
        nx.draw_networkx_edges(
            g, pos, edgelist=non_data_edges, edge_color="tab:grey", ax=ax,
            **edge_kwargs if edge_kwargs is not None else {}
        )
    else:
        weights = [g.edges[edge].get("corr") for edge in g.edges]
        vmin = np.nanmin(weights)
        vmax = np.nanmax(weights)
        nx.draw_networkx_edges(
            g, pos, edge_color=weights, edge_cmap=plt.get_cmap(cmap),
            edge_vmax=vmax, edge_vmin=vmin, ax=ax,
            **edge_kwargs if edge_kwargs is not None else {}
        )

    # node labels
    _draw_labels(g, pos, ax, True)

    # legend
    # colorbar
    _add_colorbar(cmap, vmin, vmax, ax)
    # node shapes and types, edge styles
    legend_elements = _shape_legend(["reaction", "organism"], False)
    legend_elements += _edge_style_legend({
        "Database": "-", "Correlation-only": "--"})

    ax.legend(handles=legend_elements)
    # removing grid and axis
    ax.grid("off")
    ax.axis("off")
