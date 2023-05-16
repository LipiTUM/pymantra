import networkx as nx
import matplotlib.pyplot as plt

from pymantra.network import (
    compute_reaction_estimates, add_reaction_estimates,
    MetaboliteLocalSearch
)
from pymantra.datasets import example_metabolome_enrichment_data


# loading example data
metabolite_data, sample_groups, graph = example_metabolome_enrichment_data()

# compute and add reaction estimates #
residuals = compute_reaction_estimates(graph, metabolite_data, sample_groups)
add_reaction_estimates(graph, sample_groups, residuals)

# generate local search object and run
# NOTE: these are randomly chose parameters
m_lso = MetaboliteLocalSearch(graph, 10., 1e-4, 2, 10, 10, 2)
m_lso.run_local_search(n_threads=1, min_comp_size=2)

# report the results
print(f"Local Search solution: {m_lso.solution}")
# saving results
m_lso.solution[0].to_json("example_metabolome_solution.json")

m_lso.plot_score_progression()
plt.show()

m_lso.plot_subnetwork(graph)
plt.show()

# option omitting metabolite nodes
m_lso.plot_subnetwork(node_types=nx.get_node_attributes(graph, "node_type"))
plt.show()
