import os
import json
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt

from pymantra.namemapping import NameMapper

from preprocessing_utils import scale_z_scores
from visualisation_utils import plot_pca


os.chdir("Metabolome")

# loading metadata and extracting metabolomics-defined subtypes identified
# in the original paper
metadata = pd.read_excel(
    "FUSCCTNBC.xlsx", sheet_name="Table S5", skiprows=[0], index_col=0)
metadata.index = [idx + "_T" for idx in metadata.index]

subtypes = metadata["Metabolomic_subtype"]
subtypes.to_csv("subtypes.csv")

# loading metabolite abundance -- these do not require additional
# processing apart from scaling and mapping the IDs
data = pd.read_excel(
    "FUSCCTNBC.xlsx", sheet_name="Table S1", skiprows=[0], index_col=1
).drop(columns=["Peak"])

sample_groups = pd.Series(
    ["Control" if col.endswith("N") else "Tumor" for col in data.columns],
    index=data.columns
)
sample_groups.to_csv("sample_groups.csv")

# extracting database IDs
metabolite_meta = pd.read_excel(
    "FUSCCTNBC.xlsx", sheet_name="Table S3", skiprows=[0], index_col=2)

db_map = metabolite_meta["HMDB ID"].to_dict()
for met in metabolite_meta.index[metabolite_meta["HMDB ID"].isna()]:
    kegg = metabolite_meta.loc[met, "KEGG ID"]
    if isinstance(kegg, str):
        db_map[met] = kegg
    elif not np.isnan(metabolite_meta.loc[met, "ChEBI ID"]):
        db_map[met] = str(int(metabolite_meta.loc[met, "ChEBI ID"]))


# mapping IDs to mantra IDs
mapping_file = "mantra_id_map.json"
if not os.path.exists(mapping_file):
    mapper = NameMapper()
    mantra_map = {}
    for metabolite, id_ in tqdm(db_map.items()):
        mantra_ids = []
        if isinstance(id_, str):
            if id_.startswith("HMDB"):
                mantra_ids = mapper.map_id(id_, "hmdb", "internal")
            elif id_.startswith("C"):
                mantra_ids = mapper.map_id(id_, "kegg", "internal")
        elif not np.isnan(id_):
            mantra_ids = mapper.map_id(id_, "chebi", "internal")
        mantra_map[metabolite] = mantra_ids

    json.dump(mantra_map, open(mapping_file, "w"))
else:
    mantra_map = json.load(open(mapping_file, "r"))

# drop all metabolites that could not be matched to mantra
n_mapped = 0
to_drop = set()
for met in data.index:
    if not mantra_map[met]:
        to_drop.add(met)
    else:
        n_mapped += len(mantra_map[met])

data.drop(index=list(to_drop), inplace=True)

# rename the data frame to mantra IDs and save
# NOTE: all metabolites have exactly one match, so we can just take [0] safely
data.index = [mantra_map[metabolite][0] for metabolite in data.index]
# imputing missing data with 1/2 the row-wise minimum
mins = np.nanmin(data, axis=1) / 2
for i in np.arange(data.shape[0]):
    data.values[i, np.isnan(data.values[i, :])] = mins[i]
data.to_csv("metabolite_data.csv")

# plot the PCA
fig, ax = plt.subplots(figsize=(16, 9))
plot_pca(
    scale_z_scores(data.values), sample_groups=sample_groups, ax=ax,
    file="data_pca.pdf"
)
plt.close(fig)
