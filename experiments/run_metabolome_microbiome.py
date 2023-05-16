"""
Run the enrichment with paired metabolome and microbiome data
"""
import os
import pickle

import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pymantra import (
    compute_reaction_estimates, compute_multiomics_associations,
    add_reaction_estimates, add_microbiome_associations, MetaboliteLocalSearch
)
from pymantra.plotting import plot_correlation_differences

from utils import (
    CommandLineParser, load_data, residual_analysis, run_local_search)
from visualisation_utils import (
    residual_pca, metabolome_layout, annotate_subplots)


colours = {"Control": "#1b9e77ff", "CD": "#cf45b3ff"}


# setting up argument parser
parser = CommandLineParser()

parser.add_argument(
    "--use-corrected", action="store_true",
    help="Path to a file containing confounder information. All variables in "
         "this file will be used to correct for if a file path is passed. "
         "By default all variables are assumed to be fixed effects. If any "
         "variables are supposed to be random effects, use --random-effects "
         "to specify them"
)
parser.set_defaults(use_corrected=False)
parser.parse_args()

correction = parser.get_arg("multipletest")

# figure layout
# metabolomics figure
met_fig, met_axes = metabolome_layout()

# loading data
# data is already processed (in pre_process_metabolome_microbiome.py)
base_folder = "MetabolomeMicrobiome"
sample_groups, metabolite_data, graph = load_data(
    base_folder, corrected_data=parser.get_arg("use_corrected"))
if parser.get_arg("use_corrected"):
    microbiome_data = pd.read_csv(
        os.path.join(base_folder, "corrected_microbiome.csv"), index_col=0)
else:
    microbiome_data = pd.read_csv(
        os.path.join(base_folder, "microbiome_data.csv"), index_col=0)

# only compare CD to control
samples_oi = sample_groups.index[sample_groups != "UC"]

sample_groups = sample_groups[samples_oi]
metabolite_data = metabolite_data.loc[:, samples_oi]
microbiome_data = microbiome_data.loc[:, samples_oi]

# loading pre-computed graph (generate_metabolome_microbiome_network.py)
graph = pickle.load(
    open(os.path.join(base_folder, "metabolomics_network.pickle"), "rb"))


# compute and add reaction estimates
res_file = os.path.join(base_folder, "residual.csv")
if os.path.exists(res_file):
    residuals = pd.read_csv(res_file, index_col=0)
else:
    models, _, residuals = compute_reaction_estimates(
        graph, metabolite_data.T, sample_groups, control_group="Control",
        recompute_non_passing=True, return_all=True
    )
    residuals.to_csv(res_file)

# adding estimates to graph
add_reaction_estimates(graph, sample_groups, residuals)

# counting node types
reaction_nodes = 0
for node, node_type in nx.get_node_attributes(graph, "node_type").items():
    if node_type == "reaction":
        reaction_nodes += 1
print(f"Found: {reaction_nodes} reaction nodes")

# plotting PCA for original metabolomics and residual data
residual_pca(
    metabolite_data, residuals, sample_groups,
    os.path.join(base_folder, "residual_pca.pdf"),
    plot=(met_fig, met_axes[1, :]), pca_kwargs={"cmap": colours}
)

# plotting p-values for residuals
violin_args = {
    "palette": colours,
    "drop_legend": True,
    "plot_significant_features": True,
    "thresh": .05
}
residual_analysis(
    residuals, sample_groups, os.path.join(base_folder, "residual"),
    correction, plot=(met_fig, met_axes[0, 1]), violin_args=violin_args
)
res_qvals = pd.read_csv(
    f"{base_folder}/residual_qvalues_{correction}.csv", index_col=0)

# compute and add microbiome-reaction associations
corr_associations, pvals = compute_multiomics_associations(
    residuals, microbiome_data.T, sample_groups, comparison=("Control", "CD"))
add_microbiome_associations(graph, sample_groups, corr_associations)

pickle.dump(
    corr_associations,
    open(f"{base_folder}/corr_associations.pickle", "wb")
)
pickle.dump(pvals, open(os.path.join(base_folder, "corr_pvals.pickle"), "wb"))

# Optional: plotting the associations between microbiome data and residuals
# plot_reaction_association(
#     residuals, microbiome_data, corr_diffs, sample_groups, axes=mo_axes[1:])

# correlation difference heatmaps plot separate from paper figure
clust_map = plot_correlation_differences(
    corr_associations, pvals, "CD", "Control", cluster=True, cmap="viridis")
clust_map.savefig(os.path.join(base_folder, "correlation_clustermap.pdf"))
plt.close(clust_map.figure)

# generate local search object and run enrichment
np.random.seed(42)
lso = run_local_search(
    graph, MetaboliteLocalSearch, parser,
    os.path.join(base_folder, "local_search"), omics="organism",
    delta_min=1e-4, plot=False
)
lso.plot_subnetwork(graph, ax=met_axes[0, 0])
lso.solution[0].to_json(os.path.join(base_folder, "enrichment_results.json"))

# finishing metabolite figure
annotate_subplots(
    met_axes[0, 0], met_axes[1, 0], met_axes[1, 1], met_axes[0, 1])
plt.tight_layout()
# => Figures not MetabolomeMicrobiome!
met_fig.savefig("Figures/FranzosaMetabolomeFigure.pdf")
plt.close(met_fig)
