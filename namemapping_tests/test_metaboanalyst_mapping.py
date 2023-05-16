import pathlib
import pandas as pd

from pymantra.namemapping import metaboanalyst_name_mapping


class TestMetaboanlystMapping:

    base_path = pathlib.Path(__file__).parent.absolute() / "test_data"
    ma_mapping = pd.read_csv(base_path / "metaboanalyst_ids.csv", index_col=0)

    def test_from_name(self):
        mapping = metaboanalyst_name_mapping(self.ma_mapping["Query"])
        assert all(
            mapping["HMDB"] == self.ma_mapping.loc[mapping.index, "HMDB"])

    def test_from_hmdb(self):
        metaboanalyst_name_mapping(
            self.ma_mapping["HMDB"], source_type="hmdb")
        # NOTE: there seems to be a bug in the metaboanalyst API
        # assert all(
        #     mapping.loc[self.ma_mapping["HMDB"], "KEGG"] == \
        #         self.ma_mapping["KEGG"]
        # )

    def test_from_kegg(self):
        metaboanalyst_name_mapping(
            self.ma_mapping["KEGG"], source_type="kegg")
        # NOTE: there seems to be a bug in the metaboanalyst API
        # assert all(
        #     mapping.loc[self.ma_mapping["KEGG"], "HMDB"] == \
        #         self.ma_mapping.loc[mapping.index, "HMDB"]
        # )
