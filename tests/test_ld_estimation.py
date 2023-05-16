"""
pytest tests for estimation reaction values
based on the 'local dependence' method
"""
import json
import pathlib
import numpy as np
import pandas as pd
import networkx as nx
from sklearn.linear_model import LinearRegression
from pymantra.network.ld_estimation import (
    _nll, _explained_variance, _compute_model, per_sample_ld_estimation)


def sim_data(n=100, m=10):
    x = np.array([np.random.normal(size=n) for i in range(m)]).T
    y = (x + np.random.normal(size=m)) * np.abs(np.random.normal())
    return x, y


X, Y = sim_data()


def test_nll():
    lr = LinearRegression()
    lr.fit(X, Y)
    nll = _nll(X, Y, lr.predict(X))
    print(nll)


def test_model_computation():
    _compute_model(X, Y, _explained_variance)


class TestLDEstimation:
    base_path = pathlib.Path(__file__).parent.absolute()
    test_data = base_path / 'test_data'

    node_file = open(test_data / 'test_metabolite_nodes.json', 'r')
    nodes = json.load(node_file)

    edge_file = open(test_data / 'test_metabolite_edges.json', 'r')
    edges = json.load(edge_file)

    graph = nx.DiGraph()
    for node, data in nodes.items():
        graph.add_node(node, **data)
    for edge, data in edges.items():
        src, tgt = edge.split("__")
        graph.add_edge(src, tgt, **data)

    metabolites = pd.read_csv(test_data / 'test_metabolite_data.csv')
    sample_groups = pd.Series(
        np.random.randint(0, 2, size=metabolites.shape[0]),
        dtype=str
    )
    sample_groups.index = metabolites.index

    def test_ld_estimation(self):
        # NOTE: currently disabled, since we do not use this function
        # linear models from simulated data and network
        # res = ld_estimation(self.graph, self.metabolites)
        # print(res)
        pass

    def test_per_sample_ld_estimation(self):
        for expl_var in [False, True]:
            for as_pval in [False, True]:
                _, _, scaled_res = per_sample_ld_estimation(
                    self.graph, self.metabolites, self.sample_groups,
                    control_group="0", compute_expl_var=expl_var,
                    var_as_pval=as_pval
                )
        print(scaled_res)

    def test_parallel_per_sample_ld_estimation(self):
        for expl_var in [False, True]:
            for as_pval in [False, True]:
                models, case_res, scaled_res = per_sample_ld_estimation(
                    self.graph, self.metabolites, self.sample_groups,
                    control_group="0", n_threads=2, compute_expl_var=expl_var,
                    var_as_pval=as_pval
                )
