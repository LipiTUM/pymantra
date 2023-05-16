import networkx as nx
import matplotlib.pyplot as plt

from pymantra.network import (
    compute_reaction_estimates, add_reaction_estimates,
    compute_multiomics_associations, add_microbiome_associations,
    MultiOmicsLocalSearch
)
from pymantra.datasets import example_multiomics_enrichment_data


# loading example data
metabolite_data, microbiome_data, sample_groups, graph = \
    example_multiomics_enrichment_data()


# compute and add reaction estimates
residuals = compute_reaction_estimates(graph, metabolite_data, sample_groups)
add_reaction_estimates(graph, sample_groups, residuals)

# compute and add microbiome-reaction associations
corr_associations, pvals = compute_multiomics_associations(
    residuals, microbiome_data, sample_groups, comparison=("0", "1"))
add_microbiome_associations(graph, sample_groups, corr_associations)

# generate local search object and run
# NOTE: these are randomly chose parameters
mo_lso = MultiOmicsLocalSearch(graph, "organism", 10., 1e-4, 4, 10, 10, 2)
mo_lso.run_local_search(n_threads=1, min_comp_size=4)

# report the results
print(f"Local Search solution: {mo_lso.solution}")
# saving results
mo_lso.solution[0].to_json("example_multiomics_solution.json")

mo_lso.plot_score_progression()
plt.show()

mo_lso.plot_subnetwork(graph)
plt.show()

# option omitting metabolite nodes
mo_lso.plot_subnetwork(node_types=nx.get_node_attributes(graph, "node_type"))
plt.show()
