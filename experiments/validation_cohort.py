import pathlib
import json
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier

from utils import (
    generate_search_space, optimize_hyperparams, plot_gp_optimization,
    plot_evaluation, _combined_estimates
)


def load_groups(file):
    """Load sample groups without UC and return them as a pd.Series"""
    groups = pd.read_csv(file, index_col=0)["Diagnosis"]
    return groups[groups != "UC"]


data_path = pathlib.Path(__file__).parent / "MetabolomeMicrobiome"
fig_path = pathlib.Path(__file__).parent / "Figures"
suppl_path = pathlib.Path(__file__).parent / "Supplement"

enrichment_subgraph = json.load(
    open(data_path / "enrichment_results.json", "r"))
graph = pickle.load(open(data_path / "metabolomics_network.pickle", "rb"))

met_data = pd.read_csv(data_path / "metabolite_data.csv", index_col=0)
sample_groups = load_groups(data_path / "sample_groups.csv")
int_groups = pd.Series(
    [int(group == "CD") for group in sample_groups], index=sample_groups.index
)
met_data = met_data.loc[:, sample_groups.index]

validation_met_data = pd.read_csv(
    data_path / "validation_metabolite_data.csv", index_col=0)
validation_groups = load_groups(data_path / "validation_groups.csv")
validation_met_data = validation_met_data.loc[:, validation_groups.index]


residuals = _combined_estimates(
    graph, pd.concat([met_data, validation_met_data], axis=1).T,
    pd.concat([sample_groups, validation_groups]), range(sample_groups.size),
    "Control", recompute_non_passing=True
)
discovery_data = residuals.loc[met_data.columns, :]
validation_data = residuals.loc[validation_groups.index, :]


# optimize hyperparameters on discovery cohort
rf_optim_file = data_path / "rf_optim_params.pickle"
if not rf_optim_file.exists():
    rf_search_space = generate_search_space(
        (20, 150, "n_estimators"),
        (10, discovery_data.shape[0] / 2, "max_depth"),
        (["gini", "entropy"], "criterion"),
        ([True, False], "oob_score")
    )
    rf_opt_params = optimize_hyperparams(
        rf_search_space, RandomForestClassifier, discovery_data, int_groups,
        "roc_auc", n_jobs=4
    )
    print(f"Optimized AUC value {1 - rf_opt_params.fun}")
    plot_gp_optimization(
        rf_opt_params, data_path / "RFOptimizationResults.pdf")

    pickle.dump(rf_opt_params.x, open(rf_optim_file, "wb"))


# total model
# train classifier with optimal parameters from CV on entire discovery cohort
# NOTE: we decided to *NOT* use the optimized version => remove
rf_opt_settings = pickle.load(
    open(data_path / "rf_optim_params.pickle", "rb"))
total_model = RandomForestClassifier(random_state=123)
total_model.fit(discovery_data, sample_groups)

met_model = RandomForestClassifier(random_state=123)
met_model.fit(met_data.T, sample_groups)

# roc => main figure
fig, ax = plt.subplots(figsize=(12, 12))
plot_evaluation(
    total_model, validation_data, validation_groups, ax,
    "Validation Cohort", type_="PR", linewidth=3
)
plt.tight_layout()
fig.savefig(fig_path / "validation_model.pdf")
plt.close(fig)

# roc + pr curves
fig, axes = plt.subplots(figsize=(24, 12), ncols=2)

plot_evaluation(
    total_model, validation_data, validation_groups, axes[0],
    "RandomForest with mantra"
)
plot_evaluation(
    met_model, validation_met_data.T, validation_groups, axes[0],
    "RandomForest with metabolite data"
)

plot_evaluation(
    total_model, validation_data, validation_groups, axes[1],
    "RandomForest with mantra", type_="PR"
)
plot_evaluation(
    met_model, validation_met_data.T, validation_groups, axes[1],
    "RandomForest with metabolite data", type_="PR"
)

plt.tight_layout()
fig.savefig(suppl_path / "validation_model.pdf")
plt.close(fig)
