import os
import pickle
import numpy as np
import pandas as pd
from tqdm import tqdm

from pymantra.namemapping import NameMapper


def get_metabolite_cluster(name, meta):
    matches = meta['Exact Match to Standard (* = isomer family)'] == name
    return meta.index[matches].tolist()


def generate_map(file: str, ids):
    if not os.path.exists(file):
        name_mapper = NameMapper()
        mapping = {
            str(idx): name_mapper.map_id(idx, "hmdb", "internal")
            if str(idx).startswith("HMDB") else None
            for idx in ids
        }
        pickle.dump(mapping, open(file, "wb"))
        return mapping
    return pickle.load(open(file, "rb"))


def apply_mapping(data: pd.DataFrame, mapping_tups: set):
    index = []
    mapped_data = np.zeros((len(mapping_tups), data.shape[1]))
    for i, (idx, int_idx) in enumerate(tqdm(mapping_tups)):
        rep_data = data.loc[idx, :]
        if rep_data.values.ndim > 1:
            mapped_data[i, :] = rep_data.sum()
        else:
            mapped_data[i, :] = rep_data
        index.append(int_idx)
    mapped_df = pd.DataFrame(
        mapped_data, index=index,
        columns=data.columns
    )
    return mapped_df.groupby(mapped_df.index).sum()
