import numpy as np
import pandas as pd
import networkx as nx


from pymantra.network.exceptions import NodeTypeError


def _hypergraph_incidence_matrix(
    graph: nx.DiGraph
) -> pd.DataFrame:
    r"""
    Turn the metabolic network incidence matrix into a hypergraph
    incidence matrix where the edges (i.e. the matrix columns)
    represent the reactions. Consequently, a 1 represent a product
    relation and a -1 a substrate relation.

    The hypergraph incidence matrix :math:`I_{HG}` is calculated as:

    .. math::
        I_{HG} \coloneqq I_M \cdot \left|I_R^T\right|

    where :math:`I_M` is the subset of the incidence matrix containing
    all metabolite nodes and :math:`I_R` is the subset of the incidence
    matrix containing all reaction nodes.


    Parameters
    ----------
    graph: nx.DiGraph
        Metabolic network as a directed graph (:obj:`nx.DiGraph`)

    Returns
    -------
    pd.DataFrame
        Hypergraph incidence matrix of shape (m, r) where m is the number
        of metabolite nodes in `graph` and r is the number of reaction nodes
        in `graph`. The matrix does not have column sums of 0 unlike a 'real'
        incidence matrix

        *Note*: if `graph` is an undirected graph all values will be positive.
    """
    # extracting relevant information from the graph
    node_arr = np.array(graph.nodes)
    node_types = np.array(
        list(nx.get_node_attributes(graph, 'node_type').values())
    )
    if node_types.size == 0:
        raise NodeTypeError(
            "No node type annotations were found in 'graph'. Please make sure "
            "node types are assigned as attributes ('node_type') for each "
            "node."
        )
    # TODO: make dynamic to change when `database.statics` changes
    metabolite_nodes = node_types == 'metabolite'
    reaction_nodes = node_types == 'reaction'
    # computing hypgergraph incidence matrix
    i_g = nx.incidence_matrix(graph, oriented=isinstance(graph, nx.DiGraph))
    i_hg = i_g[metabolite_nodes, :] @ np.abs(i_g[reaction_nodes, :].T)
    return pd.DataFrame(i_hg.toarray(), index=node_arr[metabolite_nodes],
                        columns=node_arr[reaction_nodes])


def io_estimation(
    graph: nx.DiGraph, metabolome_data: pd.DataFrame,
    groups: pd.Series, control_group: any = None
):
    unique_groups = groups.unique()
    if unique_groups.size != 2:
        raise ValueError(
            "Exactly 2 different groups are allowed in 'groups'"
        )
    if control_group is None:
        control_group = unique_groups[0]
    control_mask = groups == control_group

    incidence = _hypergraph_incidence_matrix(graph)
    incidence = incidence.loc[metabolome_data.columns, :]

    substrate_incidence = incidence.copy()
    substrate_incidence.loc[substrate_incidence == -1] = 0

    reaction_differences = metabolome_data.values @ incidence.values
    substrate_sums = metabolome_data.values  @ substrate_incidence.values

    initial = pd.DataFrame(
        reaction_differences / substrate_sums,
        index=metabolome_data.index, columns=incidence.columns
    )

    initial_changes = \
        initial.loc[control_mask, :] - initial.loc[~control_mask, :]
    print(initial_changes)

    # TODO
    balanced = np.array([])

    return initial, balanced
