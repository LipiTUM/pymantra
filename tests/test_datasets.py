from pymantra.datasets import (
    example_graph, example_metabolome_enrichment_data,
    example_multiomics_enrichment_data
)


def test_example_graph():
    _ = example_graph()


def test_metabolome_example():
    _ = example_metabolome_enrichment_data()


def test_multiomics_example():
    _ = example_multiomics_enrichment_data()
