import matplotlib.pyplot as plt

from pymantra.network import (
    compute_reaction_estimates, compute_multiomics_associations,)
from pymantra.datasets import example_multiomics_enrichment_data
from pymantra.plotting import (
    plot_directed_graph, residual_violinplot, plot_correlation_differences,
    plot_reaction_association
)


# loading example data
metabolite_data, microbiome_data, sample_groups, graph = \
    example_multiomics_enrichment_data()

# we set a custom color palette for the plots
colors = {
    "Control": plt.get_cmap("Set2")(0), "CD": plt.get_cmap("Set2")(2)}


# plotting the full example graph
fig, ax = plt.subplots(figsize=(16, 9))
plot_directed_graph(graph, ax=ax)
plt.show()
plt.close(fig)

# compute reaction estimates
residuals = compute_reaction_estimates(graph, metabolite_data, sample_groups)
# plot the residuals per group as violins
residual_violinplot(
    residuals, sample_groups, palette=colors, plot_significant_features=False)
plt.show()

# compute microbiome-reaction associations
corr_associations, pvals = compute_multiomics_associations(
    residuals, microbiome_data, sample_groups, comparison=("0", "1"))
# plot the difference in associations between groups in a heatmap
# `return_differences` also returns the correlation differences instead of
# just clust_map
# `reorder` reorders rows and columns according to hierarchical clustering but
# without showing the cluster tree (to show use `cluster=True`)
diff_associations, clust_map = plot_correlation_differences(
    corr_associations, pvals, "0", "1", reorder=True, return_differences=True,
    set_zero=False
)
plt.show()
# plot the reaction/microbe pairs with the highest difference in correlations
# between groups as scatter plots colored by group
plot_reaction_association(
    residuals, microbiome_data, corr_associations, sample_groups, pal=colors)
plt.show()
