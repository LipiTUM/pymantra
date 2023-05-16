import pickle
import pathlib
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pymantra.database import (
    NetworkGenerator, APINetworkGenerator, get_auth_from_env)
from pymantra.plotting import plot_directed_graph, plot_undirected_graph

from utils import parse_to_graphml


AUTH = get_auth_from_env("../.env")

try:
    generator = NetworkGenerator("bolt://127.0.0.1:7687", AUTH)
except ConnectionError:
    print(
        "Connection to local database failed. Defaulting to online.",
        "If you want to use a local database, adapt the URI and/or .env file"
    )
    generator = APINetworkGenerator()

base_path = pathlib.Path(__file__).parent.absolute() / "Metabolome"
data = pd.read_csv(base_path / "metabolite_data.csv", index_col=0)

metabolites = set(data.index)
edges = generator.get_reaction_subgraph(
    set(), set(), metabolites, reaction_organism=("Abbreviation_KEGG", "hsa"))
graph = generator.as_networkx(
    edges=edges, reaction_subgraph=True)

pickle.dump(graph, open(base_path / "mantra_graph.pickle", "wb"))

# dictionaries are not handled in graphml => removing
parse_to_graphml(graph)
nx.write_graphml(graph, base_path / "mantra_graph.graphml")

# plotting
fig, ax = plt.subplots(figsize=(16, 9))
plot_directed_graph(graph, ax=ax, font_size=5)
plt.tight_layout()
fig.savefig(base_path / "mantra_graph_directed.pdf")
plt.close(fig)

fig, ax = plt.subplots(figsize=(16, 9))
plot_undirected_graph(graph, ax=ax, font_size=5)
plt.tight_layout()
fig.savefig(base_path / "mantra_graph.pdf")
plt.close(fig)
