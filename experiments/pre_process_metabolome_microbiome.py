import os
import json
import pandas as pd
import numpy as np
from skbio.stats.composition import clr
import matplotlib.pyplot as plt
from sklearn.preprocessing import StandardScaler

from pymantra.network import confounder_correction

from preprocessing_utils import (
    filter_zero_fractions,
    impute_metabolites,
    apply_imputation,
    impute_relative,
    quotient_normalise,
    generalised_log
)
from visualisation_utils import plot_pca


os.chdir("MetabolomeMicrobiome")


meta_columns = [
    "SRA_metagenome_name", "Age", "Diagnosis", "Fecal.Calprotectin",
    "antibiotic", "immunosuppressant", "mesalamine", "steroids"
]

# ==================
# microbe processing
# ==================
otu_data = pd.read_excel(
    "41564_2018_306_MOESM6_ESM.xlsx", header=0, index_col=0, skiprows=[0])
# separating discovery and validation cohort
validation_cols = [
    col for col in otu_data.columns if col.startswith("Validation")]
validation = otu_data.loc[
    ~np.isin(otu_data.index, meta_columns), validation_cols]
otu_data.loc["Diagnosis", validation_cols].to_csv("validation_groups.csv")

otu_data = otu_data.loc[:, ~np.isin(otu_data.columns, validation.columns)]

# metadata
sample_meta = otu_data.loc[meta_columns, :].T
sample_meta.to_csv("sample_meta.csv")
sample_meta["Diagnosis"].to_csv("sample_groups.csv")

# OTU data
otu_data = otu_data.loc[~np.isin(otu_data.index, meta_columns), :]
otu_data.to_csv("OTU_data.csv")
# filtering out too sparse microbes
otu_idxs = filter_zero_fractions(otu_data.values)
otu_imputed, otu_min = impute_relative(
    otu_data.values[otu_idxs, :], return_min=True)
# clr transformation
otu_clr = clr(otu_imputed.T).T
# storing processed data
pd.DataFrame(
    otu_clr, index=otu_data.index[otu_idxs],
    columns=otu_data.columns
).to_csv("microbiome_data.csv")

# processing validation cohort according to discovery cohort parameters
pd.DataFrame(
    clr(impute_relative(
        validation.values[otu_idxs, :], total_min=otu_min).T).T,
    columns=validation.columns, index=validation.index[otu_idxs]
).to_csv("microbiome_validation.csv")


# =====================
# metabolome processing
# =====================
metabolome_data = pd.read_excel(
    "41564_2018_306_MOESM4_ESM.xlsx", index_col=0,
    skiprows=[0, 2, 3, 4, 5, 6, 7, 8]
)
# metabolome_data = pd.read_csv("metabolome_data.csv", index_col=0)

validation_met = metabolome_data.loc[:, validation_cols]
prism_cols = metabolome_data.columns[
    ~np.isin(metabolome_data.columns, validation_cols)]
metabolome_data = metabolome_data.loc[:, prism_cols]

# metabo_mds_raw = mds(metabolome_data)
# plot_mds(
#     metabo_mds_raw, sample_groups=sample_meta["Diagnosis"],
#     file="raw_metabolome_mds.pdf"
# )
# plt.close()

# TODO: remove
# # filtering metabolites with too many zero elements
# metabo_idxs = filter_zero_fractions(metabolome_data.values)
metabo_idxs = np.arange(metabolome_data.shape[0])
# imputing zeros with half minimum
metabolome_imp, imp_vals = impute_metabolites(
    metabolome_data.values, return_mins=True)
# quotient normalisation
metabo_quotient = pd.DataFrame(
    quotient_normalise(metabolome_imp[metabo_idxs, :]),
    index=metabolome_data.index[metabo_idxs], columns=metabolome_data.columns
)
metabo_quotient.to_csv("quotient_normalised.csv")

# =========================
# metabolite classification
# =========================
# log-transform and scaling
log_data = pd.DataFrame(
    generalised_log(metabo_quotient.values),
    index=metabo_quotient.index, columns=metabo_quotient.columns
)
scaler = StandardScaler()
scaled_data = pd.DataFrame(
    scaler.fit_transform(log_data.T).T, index=log_data.index,
    columns=log_data.columns
)

fig, ax = plt.subplots(figsize=(16, 9))
plot_pca(
    scaled_data.values, sample_groups=sample_meta["Diagnosis"], ax=ax)
ax.set_title("Quotient normalised Metabolome")
plt.savefig("processing_comparison.pdf")
plt.close(fig)

# saving
metabo_quotient.to_csv("mapped_metabolites_quotient.csv")
scaled_data.to_csv("mapped_scaled_metabolites.csv")

# processing validation cohort with discovery cohort parameters
quotient_validation = pd.DataFrame(
    quotient_normalise(apply_imputation(
        validation_met.values[metabo_idxs, :], imp_vals)),
    columns=validation_met.columns, index=validation_met.index[metabo_idxs]
)
log_validation = pd.DataFrame(
    generalised_log(quotient_validation.values),
    index=quotient_validation.index, columns=quotient_validation.columns
)
scaled_validation = pd.DataFrame(
    scaler.transform(log_validation.T).T, index=log_validation.index,
    columns=log_validation.columns
)
scaled_validation.to_csv("metabolome_scaled_validation.csv")

# =====================
# confounder correction
# =====================
fixed_effects = [
    "Age", "antibiotic", "immunosuppressant", "mesalamine", "steroids"]
sample_meta[fixed_effects] = sample_meta[fixed_effects].astype("category")

metabo_corrected = confounder_correction(
    scaled_data, sample_meta[fixed_effects], random_effects=None)
metabo_corrected.to_csv("corrected_metabolome.csv")

microbe_corrected = confounder_correction(
    otu_clr, sample_meta[fixed_effects], random_effects=None)
metabo_corrected.to_csv("corrected_microbiome.csv")

# ======================
# species classification
# ======================
species = {
    otu: {"Genus": otu.split("_")[0], "Species": otu.split("_")[1]}
    for otu in otu_data.index
}
json.dump(
    species, open("species_map.json", "w"), indent=4)
