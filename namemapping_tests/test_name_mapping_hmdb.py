import pathlib
import pandas as pd

from pymantra.namemapping import HMDBQuery


class TestHMDB:
    db = HMDBQuery()

    base_path = pathlib.Path(__file__).parent.absolute() / "test_data"
    ids = pd.read_csv(base_path / "test_ids.tsv", sep="\t")

    mappings = {chebi: hmdb for chebi, hmdb in zip(ids["chebi"], ids["hmdb"])}

    def test_get_column(self):
        for key, mapping in self.mappings.items():
            self.db.get_column("chebi", "hmdb", key)

    def test_get_multiple_columns(self):
        for key, mapping in self.mappings.items():
            self.db.get_multiple_columns("chebi", ["hmdb", "pubchem"], key)

    def test_taxonomy_from_foreign_id(self):
        for key, mapping in self.mappings.items():
            self.db.taxonomy_from_foreign_id("chebi", "sub_class", key)

    def test_multi_taxonomy_from_foreign_id(self):
        for key, mapping in self.mappings.items():
            self.db.multi_taxonomy_from_foreign_id(
                "chebi", ["class", "sub_class"], key)
