import networkx as nx
import matplotlib.pyplot as plt

from pymantra.statics import EDGE_BY_NODE_TYPE
from pymantra.database import reduce_reaction_nodes
from pymantra.plotting import plot_directed_graph

metabolite_nodes = {f"m{i}" for i in range(6)}
reaction_nodes = {f"r{i}" for i in range(3)}
organism_nodes = {f"o{i}" for i in range(4)}

edges = {
    # metabolite - reaction/ reaction - metabolite edges
    ("m0", "r0"), ("r0", "m1"), ("r0", "m2"),
    ("m1", "r1"), ("r1", "m3"), ("r1", "m4"),
    ("m2", "r2"), ("r2", "m5"),
    # organism - reaction edges
    ("o0", "r0"), ("o1", "r0"),
    ("o0", "r1"), ("o2", "r1"), ("o3", "r1"),
    ("o1", "r2"), ("o2", "r2"), ("o3", "r2")
}

graph = nx.DiGraph()

for node in metabolite_nodes:
    graph.add_node(node, node_type="metabolite")
for node in reaction_nodes:
    graph.add_node(node, node_type="reaction")
for node in organism_nodes:
    graph.add_node(node, node_type="organism")

for src, tgt in edges:
    src_type = graph.nodes[src]["node_type"]
    tgt_type = graph.nodes[tgt]["node_type"]
    graph.add_edge(
        src, tgt, edge_type=EDGE_BY_NODE_TYPE[(src_type, tgt_type)])

graph = reduce_reaction_nodes(graph)

plot_directed_graph(graph)
plt.show()
