import os
import pickle
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pymantra.database import (
    NetworkGenerator, APINetworkGenerator, get_auth_from_env,
    reduce_reaction_nodes
)
from pymantra.namemapping import metaboanalyst_name_mapping
from pymantra.plotting import plot_directed_graph, plot_undirected_graph

from utils import parse_to_graphml
from mapping_utils import get_metabolite_cluster, generate_map, apply_mapping
from preprocessing_utils import generalised_log

AUTH = get_auth_from_env("../.env")

os.chdir("MetabolomeMicrobiome")

# ============ #
# Name Mapping #
# ============ #
metabolomics_data = pd.read_csv("mapped_scaled_metabolites.csv", index_col=0)

metabolite_meta = pd.read_excel(
    "41564_2018_306_MOESM3_ESM.xlsx", sheet_name="dataset_s1", index_col=0,
    skiprows=[0]
)
metabolite_meta['Exact Match to Standard (* = isomer family)'] = [
    str(x).replace("*", "")
    for x in metabolite_meta['Exact Match to Standard (* = isomer family)']
]

# we use metaboanalyst to first get the HMDB and KEGG IDs for all metabolites
name_map_file = "name_map.csv"
if not os.path.exists(name_map_file) or not os.path.isfile(name_map_file):
    name_map = metaboanalyst_name_mapping(
        np.unique(
            metabolite_meta['Exact Match to Standard (* = isomer family)'])
    )
    name_map.to_csv(name_map_file)
else:
    name_map = pd.read_csv(name_map_file, index_col=0)

# mapping from HMDB to internal IDs
base_mapping = {}
for i in np.arange(name_map.shape[0]):
    for id_ in str(name_map['HMDB'][i]).split(";"):
        for franz_id in get_metabolite_cluster(name_map['Query'][i],
                                               metabolite_meta):
            if isinstance(name_map['HMDB'][i], str):
                base_mapping[id_] = franz_id

hmdb_ids = [
    x for id_ in name_map.loc[~name_map['HMDB'].isna(), 'HMDB']
    for x in id_.split(";")
]
mantra_mapping = generate_map("id_mapping.pickle", hmdb_ids)

# formatting output and renaming data
mapping_tuples = set()
for met, mantra_ids in mantra_mapping.items():
    met_id = base_mapping.get(met)
    if met_id is None:
        continue
    if mantra_ids:
        if isinstance(mantra_ids, str):
            mapping_tuples.add((met_id, mantra_ids))
        elif len(mantra_ids) == 1:
            mapping_tuples.add((met_id, mantra_ids[0]))
        elif len(mantra_ids) > 1:
            for sub_id in mantra_ids:
                mapping_tuples.add((met_id, sub_id))

n_nan = sum([len(x) == 0 for x in mantra_mapping.values()])
print(
    f"{n_nan} metabolites could not be mapped to mantra IDs\n"
    f"Continuing network generation with {metabolomics_data.shape[0] - n_nan} "
    f"metabolites."
)

scaled_metabolome = apply_mapping(metabolomics_data, mapping_tuples)
scaled_metabolome.to_csv("metabolite_data.csv")

quotient_metabolome = pd.read_csv("quotient_normalised.csv", index_col=0)
mapped_quotient = apply_mapping(quotient_metabolome, mapping_tuples)
generalised_log(mapped_quotient).to_csv("log_norm_metabolite_data.csv")

validation_metabolome = pd.read_csv(
    "metabolome_scaled_validation.csv", index_col=0)
apply_mapping(validation_metabolome, mapping_tuples).to_csv(
    "validation_metabolite_data.csv")


# ================== #
# Network Generation #
# ================== #
metabolites = set(scaled_metabolome.index)

try:
    generator = NetworkGenerator("bolt://127.0.0.1:7687", AUTH)
except ConnectionError:
    print(
        "Connection to local database failed. Defaulting to online.",
        "If you want to use a local database, adapt the URI and/or .env file"
    )
    generator = APINetworkGenerator()

edge_file = "edges.pickle"
metabolomics_network_file = "metabolomics_network"
if os.path.exists(f"{metabolomics_network_file}.pickle"):
    metabolomics_graph = pickle.load(
        open(f"{metabolomics_network_file}.pickle", "rb"))
else:
    if os.path.exists(edge_file):
        edges = pickle.load(open(edge_file, "rb"))
    else:
        edges = generator.get_reaction_subgraph(set(), set(), metabolites)
        pickle.dump(edges, open(edge_file, "wb"))
    metabolomics_graph = generator.as_networkx(
        edges=edges, reaction_subgraph=True, reduce=False)

    to_delete = [
        node for node in metabolomics_graph.nodes
        if metabolomics_graph.nodes[node]["node_type"] == "metabolite" and
        node not in scaled_metabolome.index
    ]
    metabolomics_graph.remove_nodes_from(to_delete)
    to_delete = [
        node for node in metabolomics_graph.nodes
        if metabolomics_graph.degree(node) == 0
    ]
    if to_delete:
        metabolomics_graph.remove_nodes_from(to_delete)

    metabolomics_graph = reduce_reaction_nodes(metabolomics_graph)

    pickle.dump(
        metabolomics_graph, open(f"{metabolomics_network_file}.pickle", "wb"))

    parse_to_graphml(metabolomics_graph)
    nx.write_graphml(
        metabolomics_graph, f"{metabolomics_network_file}.graphml")

matched = [
    node in metabolites for node, node_type in
    nx.get_node_attributes(metabolomics_graph, 'node_type').items()
    if node_type == "metabolite"
]
print(f"{sum(matched)} metabolites out of {len(matched)} nodes measured")

# plotting
fig, ax = plt.subplots(figsize=(16, 9))
plot_directed_graph(metabolomics_graph, ax=ax, font_size=5, show_labels=False)
plt.tight_layout()
fig.savefig("mantra_graph_directed.pdf")
plt.close(fig)

fig, ax = plt.subplots(figsize=(16, 9))
plot_undirected_graph(
    metabolomics_graph, ax=ax, font_size=5, show_labels=False)
plt.tight_layout()
fig.savefig("mantra_graph.pdf")
plt.close(fig)
