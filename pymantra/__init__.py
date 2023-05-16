import warnings
import importlib.metadata
from pymantra.network import (
    # classes
    LocalSearch,
    MetaboliteLocalSearch,
    MultiOmicsLocalSearch,
    EnrichmentResults,
    # functions
    reaction_graph_extraction,
    per_sample_ld_estimation,
    compute_reaction_estimates,
    compute_multiomics_associations,
    add_reaction_estimates,
    add_microbiome_associations,
    add_gene_associations,
    confounder_correction,
    # exceptions
    NodeTypeError
)
from pymantra.network.enrichment.spearman import multi_threading


if not multi_threading():
    warnings.warn(
        "mantra was compiled without OpenMP support. Multi-threaded local "
        "search will not be available.", category=RuntimeWarning
    )


__all__ = [
    # classes
    "LocalSearch",
    "MetaboliteLocalSearch",
    "MultiOmicsLocalSearch",
    "EnrichmentResults",
    # functions
    "reaction_graph_extraction",
    "per_sample_ld_estimation",
    "compute_reaction_estimates",
    "compute_multiomics_associations",
    "add_reaction_estimates",
    "add_microbiome_associations",
    "add_gene_associations",
    "confounder_correction",
    # exceptions
    "NodeTypeError"
]


__version__ = importlib.metadata.version("pymantra")


if __name__ == "__main__":
    from argparse import ArgumentParser
    import pytest

    parser = ArgumentParser()
    parser.add_argument(
        "--test", action='store_true', default=False, type=bool,
        help="Run all test files for the module"
    )

    args = parser.parse_args()
    if args.test:
        test_return = pytest.main(["-x", "tests"])
        if test_return == pytest.ExitCode.TESTS_FAILED:
            print(
                "\033([1;31m Failed tests report \033[0m"
            )
            exit(1)
        else:
            print(
                "\033([1;31m Tests passed \033[0m"
            )
            exit(0)

    # TODO: add options to install neo4j and sql databases via
    #       python -m pymantra ...
