import pathlib
import pandas as pd

from pymantra.namemapping import ChEBIQuery


class TestChEBI:
    db = ChEBIQuery()

    base_path = pathlib.Path(__file__).parent.absolute() / "test_data"
    ids = pd.read_csv(base_path / "test_ids.tsv", sep="\t")

    mappings = {
        inchi: chebi for inchi, chebi in zip(ids["inchi"], ids["chebi"])}

    def test_inchi2chebi(self):
        for key, mapping in self.mappings.items():
            self.db.inchi_to_chebi(key)

    def test_chebi2inchi(self):
        for key, mapping in self.mappings.items():
            self.db.inchi_to_chebi(mapping)
