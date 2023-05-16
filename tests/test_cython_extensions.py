from string import ascii_lowercase
import networkx as nx

from pymantra.network.enrichment.LSO.lso import LocalSearchOptimization


def _generate_example_graph():
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

    return ext_reaction_graph


class TestLSO:
    graph = _generate_example_graph()

    temp = 20
    delta_min = 1e-5
    l_min = 4
    l_max = 5
    max_iter = 20
    min_reactions = 4
    p = 1

    lso = LocalSearchOptimization(
        graph, temp, delta_min, l_min, l_max, max_iter, "metabolic_reactions",
        p, min_reactions
    )

    def test_getters(self):
        assert self.lso.temp() == self.temp
        assert self.lso.delta_min() == self.delta_min
        assert self.lso.l_min() == self.l_min
        assert self.lso.l_max() == self.l_max
        assert self.lso.max_iter() == self.max_iter
        assert self.lso.min_reactions() == self.min_reactions
        assert self.lso.p() == self.p

    def test_setters(self):
        self.lso.set_temp(30.)
        assert self.lso.temp() == 30.

        self.lso.set_delta_min(1e-7)
        assert self.lso.delta_min() == 1e-7

        self.lso.set_l_min(3)
        assert self.lso.l_min() == 3

        self.lso.set_l_max(7)
        assert self.lso.l_max() == 7

        self.lso.set_max_iter(25)
        assert self.lso.max_iter() == 25

        self.lso.set_min_reactions(3)
        assert self.lso.min_reactions() == 3

        self.lso.set_p(2)
        assert self.lso.p() == 2
