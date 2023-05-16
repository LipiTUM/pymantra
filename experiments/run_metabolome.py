"""
Run the enrichment with metabolome data only
"""
import os
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns

from pymantra import compute_reaction_estimates, add_reaction_estimates
from pymantra.plotting import (
    plot_directed_graph, plot_undirected_graph, NODE_COLOURS)

from utils import CommandLineParser, load_data, residual_analysis
from visualisation_utils import (
    residual_pca, annotate_subplots, metabolome_layout, plot_graph_by_degree)

colours = {
    "Control": plt.get_cmap("Set1")(0), "Tumor": plt.get_cmap("Set1")(1)}

# setting up argument parser
parser = CommandLineParser()

parser.parse_args()
correction = parser.get_arg("multipletest")

# loading data
base_folder = "Metabolome"
sample_groups, metabolite_data, graph = load_data(base_folder)

# supplementary figure: degree distributions
degrees_by_type = {"metabolite": [], "reaction": []}
n_types = {"metabolite": 0, "reaction": 0}
for node, degree in graph.degree:
    node_type = graph.nodes[node]["node_type"]
    degrees_by_type[node_type].append(degree)
    n_types[node_type] += 1

print("Found:")
for node_type, n_type in n_types.items():
    print(f"{n_type} {node_type}")

fig, (ax1, ax2) = plt.subplots(
    figsize=(16, 18), nrows=2, gridspec_kw={'height_ratios': [2 / 3, 1 / 3]})

plot_graph_by_degree(
    graph, ax1, show_labels=False, cmap="Spectral", node_size=100)

node_types = []
degrees = []
for node_type, nt_degs in degrees_by_type.items():
    node_types += len(nt_degs) * [node_type]
    degrees += nt_degs

sns.kdeplot(
    data=pd.DataFrame({"Degree": degrees, "NodeType": node_types}),
    x="Degree", hue="NodeType", ax=ax2, fill=True,
    palette={nt: NODE_COLOURS[nt] for nt in n_types.keys()}
)
ax2.set_xlabel("Degree")
ax2.set_ylabel("Frequency")

plt.tight_layout()

fig.savefig("Supplement/NodeDegrees.pdf")

top_degree_table = []
for node in graph.nodes:
    if graph.nodes[node]["node_type"] == "metabolite" and \
            graph.degree(node) > 10:
        top_degree_table.append(
            {
                "Label": graph.nodes[node]["nodeLabel"],
                "Name": graph.nodes[node]["Name"],
                "Degree": graph.degree(node),
            }
        )
pd.DataFrame(top_degree_table).to_latex(
    "Supplement/NodeDegrees.tex", index=False)

# figure layout
fig, axes = metabolome_layout()

plot_directed_graph(
    graph, ax=axes[0, 0], font_size=5, show_labels=False, edge_width=.5)

# compute and add reaction estimates
res_file = os.path.join(base_folder, "residuals.csv")
if os.path.exists(res_file):
    residuals = pd.read_csv(res_file, index_col=0)
else:
    residuals = compute_reaction_estimates(
        graph, metabolite_data.T, sample_groups, control_group="Control",
        verbose=True, recompute_non_passing=True
    )
    residuals.to_csv(res_file)

add_reaction_estimates(graph, sample_groups, residuals)

# plotting PCA for original metabolomics and residual data
residual_pca(
    metabolite_data, residuals, sample_groups,
    os.path.join(base_folder, "residual_pca.pdf"),
    plot=(fig, axes[1, :]), pca_kwargs={"cmap": colours}
)

# plotting p-values for residuals
res_qvals = residual_analysis(
    residuals, sample_groups, os.path.join(base_folder, "residual"),
    correction, plot=(fig, axes[0, 1]),
    violin_args={"palette": colours, "drop_legend": True,
                 "plot_significant_features": True},
    return_qvals=True
)

# res_dict = {
#     "Name": [],
#     "Formula": [],
#     "Reaction Node": [],
# }
# for idx in res_qvals.index:
#     reactions = idx.split(", ")
#     node_name = reactions[0]
#     for reaction in idx.split(", "):
#         # TODO: get these
#         res_dict["Name"].append()
#         res_dict["Formula"].append()
#         res_dict["Reaction Node"].append(node_name)
# # supplementary table mapping reaction names
# pd.DataFrame(res_dict).to_latex("Supplement/residual_table.tex", index=False)
# supplementary table describing reactions
node_df = pd.DataFrame.from_dict(
    {
        node: data for node, data in graph.nodes(data=True)
        if data["node_type"] == "reaction"
    },
    orient="index"
).drop(columns=["vec_group_data"])
node_df.to_csv("Supplement/reaction_information.csv")

# TODO: include *all* reactions
cols_oi = ["Formula", "Description", "ReconMap3ID", "ReactomeID", "KeggID"]
node_df.loc[res_qvals.index, cols_oi].to_csv("Supplement/residual_table.csv")

# adding subplot annotation (a-d)
annotate_subplots(axes[0, 0], axes[1, 0], axes[1, 1], axes[0, 1])

fig.savefig("Figures/MetabolomeFigure.pdf")
plt.close(fig)

# supplementary figures
fig, ax = plt.subplots(figsize=(16, 9))

sig_nodes = set()
for node in res_qvals.index:
    sig_nodes = sig_nodes.union(set(nx.ego_graph(graph, node).nodes))
subgraph = sorted(nx.connected_components(
    graph.subgraph(sig_nodes).to_undirected()), key=len, reverse=True)[0]
plot_undirected_graph(graph.subgraph(subgraph), ax=ax)

plt.tight_layout()
fig.savefig(os.path.join("Supplement", "SignificantNetwork.pdf"))

fig, ax = plt.subplots(figsize=(16, 9))
plot_directed_graph(graph, ax=ax, show_labels=False)
fig.savefig("Metabolome/unlabeled_graph.pdf")
plt.close(fig)
