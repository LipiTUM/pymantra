import pandas as pd

from pymantra.network import compute_reaction_estimates, add_reaction_estimates
from pymantra.datasets import example_metabolome_enrichment_data


# loading example data and graph
metabolite_data, sample_groups, graph = example_metabolome_enrichment_data()
control_mask = sample_groups == '0'

# define some confounding variables
confounders = pd.DataFrame(
    {
        "site": [],
        "gender": [],
        "age": []
    }
)

# compute and add reaction estimates with confounding factors
residuals = compute_reaction_estimates(
    graph, metabolite_data, sample_groups, covariates=confounders,
    random_effects=[""]
)
add_reaction_estimates(graph, sample_groups, residuals)
