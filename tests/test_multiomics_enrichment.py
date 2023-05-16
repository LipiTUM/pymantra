import pathlib
from string import ascii_lowercase
import pandas as pd
import networkx as nx
import pytest

from pymantra import (
    compute_reaction_estimates, add_reaction_estimates,
    compute_multiomics_associations, add_microbiome_associations,
    add_gene_associations
)
from pymantra.network import MultiOmicsLocalSearch


def test_microbiome_enrichment():
    # generate a random bipartite metabolic network
    metabolite_nodes = [x for x in ascii_lowercase[:6]]
    reaction_nodes = ["A", "B", "C"]
    organism_nodes = ["O1", "O2", "O3"]
    substrate_edges = [
        ("a", "A"), ("c", "B"), ("e", "C"), ("d", "C")]
    product_edges = [
        ("A", "b"), ("A", "c"), ("B", "d"), ("B", "e"), ("C", "f")]
    organism_edges = [
        ("A", "O1"), ("A", "O2"), ("B", "O2"), ("B", "O3"), ("C", "O3")]

    graph = nx.DiGraph()
    for node in metabolite_nodes:
        graph.add_node(node, node_type="metabolite")
    for node in reaction_nodes:
        graph.add_node(node, node_type="reaction")
    for node in organism_nodes:
        graph.add_node(node, node_type="organism")

    for src, tgt in substrate_edges:
        graph.add_edge(src, tgt, edge_type="SUBSTRATE")
    for src, tgt in product_edges:
        graph.add_edge(src, tgt, edge_type="PRODUCT")
    for src, tgt in organism_edges:
        graph.add_edge(src, tgt, edge_type="REACTION_ORGANISM")

    # generate random data
    sample_names = [f"S{i}" for i in range(7)]
    sample_groups = pd.Series(
        ["0", "0", "0", "0", "1", "1", "1"], index=sample_names)

    base_path = pathlib.Path(__file__).parent.absolute()
    metabolite_data = pd.read_csv(
        base_path / "test_data/metabolite_data.csv", index_col=0)
    organism_data = pd.read_csv(
        base_path / "test_data/microbiome_data.csv", index_col=0)

    # compute models per reaction
    residuals = compute_reaction_estimates(
        graph, metabolite_data, sample_groups)
    add_reaction_estimates(graph, sample_groups, residuals)

    corrs, pvals = compute_multiomics_associations(
        residuals, organism_data, sample_groups, comparison=("0", "1"))
    add_microbiome_associations(graph, sample_groups, corrs)

    mmls = MultiOmicsLocalSearch(
        graph, "organism", 10., 1e-4, 5, 10, 10, 2)

    mmls.run_local_search(n_threads=1, min_comp_size=4)
    mmls.run_local_search(n_threads=4, min_comp_size=4)


def test_gene_enrichment():
    # generate a random bipartite metabolic network
    metabolite_nodes = [x for x in ascii_lowercase[:6]]
    reaction_nodes = ["A", "B", "C"]
    gene_nodes = ["G1", "G2", "G3"]
    substrate_edges = [
        ("a", "A"), ("c", "B"), ("e", "C"), ("d", "C")]
    product_edges = [
        ("A", "b"), ("A", "c"), ("B", "d"), ("B", "e"), ("C", "f")]
    gene_edges = [
        ("A", "G1"), ("A", "G2"), ("B", "G2"), ("B", "G3"), ("C", "G3")]

    graph = nx.DiGraph()
    for node in metabolite_nodes:
        graph.add_node(node, node_type="metabolite")
    for node in reaction_nodes:
        graph.add_node(node, node_type="reaction")
    for node in gene_nodes:
        graph.add_node(node, node_type="gene")

    for src, tgt in substrate_edges:
        graph.add_edge(src, tgt, edge_type="SUBSTRATE")
    for src, tgt in product_edges:
        graph.add_edge(src, tgt, edge_type="PRODUCT")
    for src, tgt in gene_edges:
        graph.add_edge(src, tgt, edge_type="REACTION_GENE")

    # generate random data
    sample_names = [f"S{i}" for i in range(7)]
    sample_groups = pd.Series(
        ["0", "0", "0", "0", "1", "1", "1"], index=sample_names)

    base_path = pathlib.Path(__file__).parent.absolute()
    metabolite_data = pd.read_csv(
        base_path / "test_data/metabolite_data.csv", index_col=0)
    gene_data = pd.read_csv(
        base_path / "test_data/microbiome_data.csv", index_col=0)
    gene_data.columns = [f"G{i + 1}" for i in range(gene_data.shape[1])]

    # compute models per reaction
    residuals = compute_reaction_estimates(
        graph, metabolite_data, sample_groups)
    add_reaction_estimates(graph, sample_groups, residuals)

    corrs, pvals = compute_multiomics_associations(
        residuals, gene_data, sample_groups)

    with pytest.raises(ValueError):
        add_gene_associations(graph, sample_groups)
    add_gene_associations(graph, sample_groups, corrs)

    mmls = MultiOmicsLocalSearch(
        graph, "gene", 10., 1e-4, 5, 10, 10, 2)

    mmls.run_repeated_local_search(5, n_threads=1, min_comp_size=4)
    mmls.run_repeated_local_search(5, n_threads=4, min_comp_size=4)

    mmls.plot_subnetwork(graph)
