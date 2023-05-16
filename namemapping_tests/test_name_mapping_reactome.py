import pathlib
import pandas as pd

from pymantra.namemapping import ReactomeQuery


class TestReactome:
    db = ReactomeQuery()

    base_path = pathlib.Path(__file__).parent.absolute() / "test_data"
    ids = pd.read_csv(base_path / "test_ids.tsv", sep="\t")

    mappings = {
        str(int(chebi)): reactome for chebi, reactome in
        zip(ids["chebi"], ids["reactome"])
    }

    def test_chebi(self):
        for key, mapping in self.mappings.items():
            self.db.chebi_to_reactome(key)

    def test_reactome2chebi(self):
        for key, mapping in self.mappings.items():
            self.db.chebi_to_reactome(mapping)

    def test_reactome2ncbi(self):
        for key in self.mappings.values():
            self.db.reactome_to_ncbi(key)

    def test_ncbi2reactome(self):
        # TODO
        pass
