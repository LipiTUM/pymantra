import os
import pathlib
from typing import Union, Tuple
import networkx as nx


def reduce_reaction_nodes(graph: nx.DiGraph):
    """Clean-up a metabolite-reaction graph

    This functions merges reaction which have the same substrates and products
    (e.g. because they are coming from different organisms) and removes
    reaction nodes which have only one participant, i.e. transport reactions.

    All modifications are done **inplace**.

    This utility function is only recommended if you did not use the pymantra
    database module to generate your network, as the already applies this
    function internally.

    Parameters
    ----------
    graph : nx.DiGraph
        Directed metabolite-reaction graph, optionally with multi-omics
        reaction connections
    """
    # gathering reaction node information
    reaction_nodes = {}
    to_remove = set()
    for node, node_type in nx.get_node_attributes(graph, "node_type").items():
        if node_type == "reaction":
            subs = tuple(sorted(graph.predecessors(node)))
            prod = tuple(sorted(graph.successors(node)))
            if subs == prod or not subs or not prod:
                to_remove.add(node)
            else:
                # reaction nodes are stored by key, which will be identical
                # if reaction nodes have the exact same set of substrates
                # and products
                reaction_nodes.setdefault((subs, prod), set()).add(node)
    # remove nodes with no substrates or products or where substrates and
    # products are the same
    graph.remove_nodes_from(to_remove)
    # merge reaction nodes with the same substrates and products
    nx.set_node_attributes(graph, values=False, name="merged")
    for nodes in reaction_nodes.values():
        relabel = {}
        if len(nodes) > 1:
            node_list = list(nodes)
            relabel[node_list[0]] = ", ".join(node_list)
            for i in range(1, len(nodes)):
                nx.contracted_nodes(
                    graph, node_list[0], node_list[i], copy=False)
        nx.relabel_nodes(graph, relabel, copy=False)
        nx.set_node_attributes(
            graph, values={k: True for k in relabel.keys()}, name="merged")
    return graph


def read_env(filepath: Union[str, pathlib.Path], add_to_env: bool = True):
    """Process the content of an .env file

    Variables will either be added to the current environment or returned as a
    dictionary

    Parameters
    ----------
    filepath: Union[str, pathlib.Path]
        Path to the .env file
    add_to_env: bool, False
        If True the variables will be added to the global environment, else
        they will be returned as a dictionary

    Returns
    -------
    Union[None, Dict[str, str]]
        None if `add_to_env` is False, otherwise a dictionary containing the
        variable name/value pairs
    """
    env_vars = {}
    with open(filepath, "r") as file:
        for line in file:
            try:
                key, val = line.strip().split("=")
                env_vars[key] = val
            except ValueError:
                raise ValueError(
                    "Invalid .env file. Please check against the .env.template"
                )
    if add_to_env:
        for key, var in env_vars.items():
            os.environ[key] = var
    else:
        return env_vars


def get_auth_from_env(filepath: Union[str, pathlib.Path]) -> Tuple[str, str]:
    """Get the neo4j username and password from a .env file

    Variables should be named "NEO4J_USER" and "NEO4J_PASSWORD". For a .env
    template use the .env.template in the root of the pymantra project folder

    Parameters
    ----------
    filepath: Union[str, pathlib.Path]
        Path to the .env file

    Returns
    -------
    Tuple[str, str]
    """
    env_data = read_env(filepath, False)
    return env_data["NEO4J_USER"], env_data["NEO4J_PASSWORD"]
