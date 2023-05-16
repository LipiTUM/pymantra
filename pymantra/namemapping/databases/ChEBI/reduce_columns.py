import pandas as pd


accession_cols = [1, 2, 4]
all_chebi = pd.read_csv("database_accession.tsv", sep="\t")
all_chebi.iloc[:, accession_cols].to_csv(
    "accession_chebi.tsv", sep="\t", index=False, header=False)
