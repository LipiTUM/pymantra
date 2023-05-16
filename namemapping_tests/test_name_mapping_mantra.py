import pathlib
import pandas as pd

from pymantra.namemapping import MantraDBQuery


class TestMantra:
    db = MantraDBQuery()

    base_path = pathlib.Path(__file__).parent.absolute() / "test_data"
    ids = pd.read_csv(base_path / "test_ids.tsv", sep="\t")

    chebi_mappings = {
        mantra: str(int(chebi)) for mantra, chebi
        in zip(ids["abbreviation"], ids["chebi"])
    }
    hmdb_mappings = {
        mantra: hmdb for mantra, hmdb
        in zip(ids["abbreviation"], ids["hmdb"])
    }
    kegg_mappings = {
        mantra: kegg for mantra, kegg
        in zip(ids["abbreviation"], ids["kegg"])
    }
    reactome_mappings = {
        mantra: reactome for mantra, reactome
        in zip(ids["abbreviation"], ids["reactome"])
    }
    vmh_mappings = {
        mantra: vmh for mantra, vmh
        in zip(ids["abbreviation"], ids["vmh"])
    }

    def test_chebi2internal(self):
        for chebi, internal in self.chebi_mappings.items():
            self.db.chebi_to_internal(chebi)

    def test_internal2chebi(self):
        for chebi, internal in self.chebi_mappings.items():
            self.db.internal_to_chebi(internal)

    def test_hmdb2internal(self):
        for hmdb, internal in self.hmdb_mappings.items():
            self.db.hmdb_to_internal(hmdb)

    def test_internal2hmdb(self):
        for hmdb, internal in self.hmdb_mappings.items():
            self.db.internal_to_hmdb(internal)

    def test_kegg2internal(self):
        for kegg, internal in self.kegg_mappings.items():
            self.db.kegg_to_internal(kegg)

    def test_internal2kegg(self):
        for kegg, internal in self.kegg_mappings.items():
            self.db.internal_to_kegg(internal)

    def test_reactome2internal(self):
        for reactome, internal in self.reactome_mappings.items():
            self.db.reactome_to_internal(reactome)

    def test_internal2reactome(self):
        for reactome, internal in self.reactome_mappings.items():
            self.db.internal_to_reactome(internal)

    def test_vmh2internal(self):
        for vmh, internal in self.vmh_mappings.items():
            self.db.vmh_to_internal(vmh)

    def test_internal2vmh(self):
        for vmh, internal in self.vmh_mappings.items():
            self.db.internal_to_vmh(internal)
