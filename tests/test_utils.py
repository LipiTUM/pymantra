import os
import pathlib

from pymantra.datasets import example_graph
from pymantra.database import reduce_reaction_nodes, read_env


def test_reduce_reaction_nodes():
    graph = example_graph()
    reduce_reaction_nodes(graph)


def test_env_reading():
    read_env(pathlib.Path(__file__).parent.parent.absolute() / ".env.template")
    assert os.getenv("NEO4J_USER") == "neo4j"
    assert os.getenv("NEO4J_PASSWORD") == "123"
