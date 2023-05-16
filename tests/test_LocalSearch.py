"""
Test suite for local search classes
"""
import pytest
import os
import pathlib
import numpy as np
import networkx as nx
from string import ascii_lowercase, ascii_uppercase

from pymantra.network import (
    MetaboliteLocalSearch, MultiOmicsLocalSearch, reaction_graph_extraction)
from pymantra.network import EnrichmentResults, RepeatedEnrichmentResults
from pymantra.network.enrichment.LocalSearch import EnrichmentBaseResults


class TestLocalSearch:
    np.random.seed(42)

    reaction_values = [
        {"x": -0.42, "y": .15},
        {"x": .76, "y": -3.},
        {"x": 3.34, "y": -3.7},
        {"x": .78, "y": -.36},
        {"x": -.32, "y": .02},
        {"x": .16, "y": .37},
        {"x": -.3, "y": .04}
    ]
    reaction_graph_edges = {
        (0, 1), (0, 2), (0, 5), (1, 2),
        (2, 3), (2, 5), (3, 4), (5, 6)
    }

    ext_reaction_values = [
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": 0, "y": 0},
        {"x": -0.42, "y": .15},
        {"x": .76, "y": -3.},
        {"x": 3.34, "y": -3.7},
        {"x": .78, "y": -.36}
    ]
    ext_reaction_graph_edges = {
        (0, 14), (0, 15), (1, 13), (1, 14), (2, 12),
        (2, 14), (3, 15), (4, 15), (5, 12), (5, 15),
        (6, 14), (6, 15), (7, 13), (7, 14), (8, 13),
        (9, 13), (10, 12), (11, 12)
    }

    # bipartite metabolite-reaction graph
    ext_reaction_graph = nx.Graph()
    for i, values in enumerate(ext_reaction_values):
        if i < 12:
            ext_reaction_graph.add_node(ascii_lowercase[i], data=values,
                                        node_type="metabolite")
        else:
            ext_reaction_graph.add_node(ascii_lowercase[i], data=values,
                                        node_type="reaction")
    for edge in ext_reaction_graph_edges:
        ext_reaction_graph.add_edge(ascii_lowercase[edge[0]],
                                    ascii_lowercase[edge[1]],
                                    edge_type="reaction_reaction",
                                    data={'x': 0, 'y': 0})

    # multi-omics graphs
    microbiome_edges = [
        (20, 12, {"x": -.42, "y": .15}),
        (20, 13, {"x": .76, "y": -1.}),
        (21, 12, {"x": -.32, "y": .45}),
        (21, 14, {"x": -.42, "y": .55}),
        (22, 13, {"x": .78, "y": -.36}),
        (23, 13, {"x": .34, "y": .37}),
        (23, 14, {"x": .24, "y": .47}),
        (23, 15, {"x": .30, "y": .17}),
    ]
    mo_graph = ext_reaction_graph.copy()
    for edge in microbiome_edges:
        if edge[0] not in mo_graph.nodes:
            mo_graph.add_node(edge[0], node_type="organism")
        mo_graph.add_edge(edge[0], edge[1], **edge[2])

    gene_edges = [
        (12, {"x": -.42, "y": .15}),
        (12, {"x": .76, "y": -1.}),
        (12, {"x": .34, "y": .37}),
        (12, {"x": .78, "y": -.36})
    ]

    def test_reaction_graph(self):
        computed_edges = reaction_graph_extraction(self.ext_reaction_graph)
        expected_solution = {('m', 'o'), ('m', 'p'), ('n', 'o'), ('o', 'p')}
        assert computed_edges == expected_solution, \
            "Incorrect edge extraction"

    def test_reaction_only(self):
        reaction_graph = nx.Graph()
        for i, values in enumerate(self.reaction_values):
            values = {
                group: [x + .001, x - .01, x + 0.05]
                for group, x in values.items()
            }
            reaction_graph.add_node(ascii_lowercase[i], vec_group_data=values,
                                    node_type="reaction")

        for edge in self.reaction_graph_edges:
            reaction_graph.add_edge(
                ascii_lowercase[edge[0]], ascii_lowercase[edge[1]],
                edge_type="reaction_reaction", data={'x': 0, 'y': 0}
            )

        expected_solutions = [{"b", "c", "d"}, {"a", "b", "c", "d"}]

        mlso = MetaboliteLocalSearch(
            network=reaction_graph, temp=20, delta_min=1e-10,
            l_min=3, l_max=5, max_iter=2, is_reaction_graph=True,
            p=1, objective_function="ld_reactions"
        )

        # single-thread
        mlso.run_local_search(n_threads=1, min_comp_size=5)
        # print(mlso.solution)
        assert mlso.solution[0].subgraph in expected_solutions

        # multi-thread
        mlso.run_local_search(n_threads=2, min_comp_size=5)
        assert mlso.solution[0].subgraph in expected_solutions

        mlso.score_final_solution(["x", "y"])

    def test_from_metabolite_reaction_graph(self):
        expected_solutions = [{"n", "o", "p"}, {"m", "o", "p"}]
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )
        # single-thread
        mlso.run_local_search(n_threads=1, min_comp_size=4)
        assert mlso.solution[0].subgraph in expected_solutions

    def test_from_metabolite_reaction_graph_parallel(self):
        expected_solutions = [{"n", "o", "p"}, {"m", "o", "p"}]
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )
        # multi-thread
        mlso.run_local_search(n_threads=2, min_comp_size=4)
        assert mlso.solution[0].subgraph in expected_solutions

    def test_repeated_search(self):
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )
        # single-thread
        mlso.run_repeated_local_search(3, n_threads=1, min_comp_size=4)
        # plotting
        with pytest.raises(ValueError):
            mlso.plot_subnetwork()
        mlso.plot_subnetwork(
            node_types=nx.get_node_attributes(
                self.ext_reaction_graph, "node_type")
        )

        mlso.plot_subnetwork(self.ext_reaction_graph.to_directed())

    def test_metabolomics_results(self):
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )
        with pytest.raises(ValueError) as ve:
            _ = mlso.objective_values
            print(ve)
        # single-thread
        mlso.run_local_search(n_threads=1, min_comp_size=4)
        # plotting
        mlso.plot_subnetwork(self.ext_reaction_graph.to_directed())
        mlso.plot_score_progression()
        # scores
        assert mlso.score_progression is not None
        # objective values
        assert mlso.objective_values is not None
        # properties
        assert isinstance(mlso.converged, bool)
        assert isinstance(mlso.l_min, int)
        assert isinstance(mlso.l_max, int)
        assert isinstance(mlso.min_reactions, int)
        assert isinstance(mlso.temp, float)
        assert isinstance(mlso.max_iter, int)
        assert isinstance(mlso.p, float)

    def test_setter_getter(self):
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )

        with pytest.raises(ValueError) as negative_error:
            mlso.set_l_max(-4)
            print(f"Expected exception: {negative_error}")
        with pytest.raises(ValueError) as non_numeric_error:
            mlso.set_l_max("x")
            print(f"Expected exception: {non_numeric_error}")

        mlso.set_temp(25.)
        assert mlso.temp == 25.

        mlso.set_l_min(2)
        assert mlso.l_min == 2

        mlso.set_l_max(10)
        assert mlso.l_max == 10

        mlso.set_min_reactions(4)
        assert mlso.min_reactions == 4

        mlso.set_max_iter(70)
        assert mlso.max_iter == 70

        mlso.set_delta_min(.4)
        assert mlso.delta_min == .4

        mlso.set_p(2)
        assert mlso.p == 2

        mlso.set_seed("a", 5)

    @staticmethod
    def _random_graph(node_labels=ascii_lowercase):
        graph = nx.erdos_renyi_graph(10, .4, directed=True)
        nx.set_node_attributes(
            graph, values={i: "reaction" for i in range(10)}, name="node_type")
        return nx.relabel_nodes(graph, {i: node_labels[i] for i in range(10)})

    def test_multicomponent_plotting(self):
        mlso = MetaboliteLocalSearch(
            self.ext_reaction_graph, temp=20, delta_min=1e-10, l_min=3,
            l_max=5, max_iter=10, is_reaction_graph=False,
            p=1, objective_function="metabolic_reactions"
        )
        # NOTE: this not how it should be used in real-world scenarios
        #       just a simplification to not have to build a more complex
        #       example network
        mlso._solution = [
           RepeatedEnrichmentResults(
               {"a", "b", "c"}, [.3, .4], [[.2, .25, .3], [.27, .35, .4]])
        ]
        mlso._score_progression = [
            np.abs(np.random.randn(5)), np.abs(np.random.randn(5))]

        mlso.plot_subnetwork(
            nx.compose(
                self._random_graph(), self._random_graph(ascii_uppercase)))
        mlso.plot_score_progression()

    def test_multiomics_result_plotting(self):
        mlso = MultiOmicsLocalSearch(
            self.ext_reaction_graph, "organism", temp=20, delta_min=1e-10,
            l_min=3, l_max=5, max_iter=10, min_reactions=2,
            is_reaction_graph=False, p=1
        )
        # single-thread
        mlso.run_local_search(n_threads=1, min_comp_size=4)
        mlso.plot_subnetwork(self.ext_reaction_graph.to_directed())

    def test_microbiome(self):
        # TODO
        pass

    def test_transcriptome(self):
        # TODO
        pass

    def test_min_reactions(self):
        # TODO
        pass


class TestReporting:
    data_path = pathlib.Path(__file__).parent.absolute() / "test_data"

    single_results = EnrichmentResults({"a", "b", "c"}, 1.0, True)
    single_result_file = data_path / "single_results.json"

    repeated_results = RepeatedEnrichmentResults(
        {"a", "b", "c"}, [1.0, 2.0, 3.0], [[1, 2, 3], [4, 5, 6]])
    repeated_result_file = data_path / "repeated_results.json"

    @staticmethod
    def _write_read_check(file: pathlib.Path, results: EnrichmentBaseResults):
        # write results data
        class_type = type(results)
        results.to_json(file)
        # read data back in
        read_results = class_type.from_json(file)
        # remove file since we don't need to store it
        os.remove(file)
        # check equality of original and written/read objects
        assert read_results == results, \
            f"'{class_type.__name__}' changed after writing and reading " \
            "back in"

    def test_single_enrichment_result(self):
        self._write_read_check(self.single_result_file, self.single_results)

    def test_repeated_enrichment_result(self):
        self._write_read_check(
            self.repeated_result_file, self.repeated_results)
