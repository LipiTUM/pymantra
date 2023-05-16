import pathlib
from typing import Tuple
from string import ascii_lowercase
import pandas as pd
import networkx as nx


_DATA_PATH = pathlib.Path(__file__).parent.absolute() / "data"


def example_graph() -> nx.DiGraph:
    """Get an example mantra graph

    Returns
    -------
    nx.DiGraph
        Example graph containing metabolites, reactions, organisms and genes
    """
    return nx.read_graphml(_DATA_PATH / "example_graph.graphml")


def _data_graph(with_organism=False) -> nx.DiGraph:
    metabolite_nodes = [x for x in ascii_lowercase[:6]]
    reaction_nodes = ["A", "B", "C"]
    substrate_edges = [
        ("a", "A"), ("c", "B"), ("e", "C"), ("d", "C")]
    product_edges = [
        ("A", "b"), ("A", "c"), ("B", "d"), ("B", "e"), ("C", "f")]

    graph = nx.DiGraph()
    for node in metabolite_nodes:
        graph.add_node(node, node_type="metabolite")
    for node in reaction_nodes:
        graph.add_node(node, node_type="reaction")

    for src, tgt in substrate_edges:
        graph.add_edge(src, tgt, edge_type="SUBSTRATE")
    for src, tgt in product_edges:
        graph.add_edge(src, tgt, edge_type="PRODUCT")

    if with_organism:
        organism_nodes = ["O1", "O2", "O3"]
        organism_edges = [
            ("A", "O1"), ("A", "O2"), ("B", "O2"), ("B", "O3"), ("C", "O3")]

        for node in organism_nodes:
            graph.add_node(node, node_type="organism")
        for src, tgt in organism_edges:
            graph.add_edge(src, tgt, edge_type="REACTION_ORGANISM")

    return graph


def example_metabolome_enrichment_data(
) -> Tuple[pd.DataFrame, pd.Series, nx.DiGraph]:
    """Get example data containing metabolome data, sample groups and a graph

    Returns
    -------
    Tuple[pd.DataFrame, pd.Series, nx.DiGraph]
        3-tuple of metabolome data, sample groups and a corresponding
        'metabolic' graph
    """
    met_data = pd.read_csv(_DATA_PATH / "metabolite_data.csv", index_col=0)
    sample_names = [f"S{i}" for i in range(7)]
    sample_groups = pd.Series(
        ["0", "0", "0", "0", "1", "1", "1"], index=sample_names)
    return met_data, sample_groups, _data_graph()


def example_multiomics_enrichment_data(
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.Series, nx.DiGraph]:
    """Get example data containing metabolome and microbiome data, sample
    groups and a graph

    Returns
    -------
    Tuple[pd.DataFrame, pd.Series, nx.DiGraph]
        4-tuple of metabolome data, microbiome data, sample groups and a
        corresponding 'metabolic' graph
    """
    met_data = pd.read_csv(_DATA_PATH / "metabolite_data.csv", index_col=0)
    mic_data = pd.read_csv(_DATA_PATH / "microbiome_data.csv", index_col=0)
    sample_names = [f"S{i}" for i in range(7)]
    sample_groups = pd.Series(
        ["0", "0", "0", "0", "1", "1", "1"], index=sample_names)
    return met_data, mic_data, sample_groups, _data_graph(True)
