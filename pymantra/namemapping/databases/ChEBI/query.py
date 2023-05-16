import pathlib
import warnings
from sqlite3 import ProgrammingError
from typing import List
from pymantra.namemapping.databases.sqlite_base import SQLiteBase, unique_list


class ChEBIQuery(SQLiteBase):
    """Query class to map from and to ChEBI"""
    def __init__(self, *args, **kwargs):
        db_file = pathlib.Path(__file__).parent.absolute() / "chebi.db"
        super(ChEBIQuery, self).__init__(db_file, *args, **kwargs)
        self.id_columns = {'kegg_id', 'chebi_id', 'inchi', 'foreign_id'}

    def active_connection(self) -> bool:
        try:
            self.execute_query("select kegg_id from chebi limit 1")
            return True
        except ProgrammingError as err:
            warnings.warn(f"Reactome database connection unavailable.\n{err}")
            return False

    @unique_list
    def _accession_query(
        self, src_id: str, src_col: str, tgt_col: str
    ) -> List[str]:
        if src_col not in self.id_columns:
            raise ValueError(
                f"'{src_col}' is not a supported ID type for the ChEBI query "
                "class, please use on of the following: "
                f"{', '.join(self.id_columns)}"
            )
        if tgt_col not in self.id_columns:
            raise ValueError(
                f"'{tgt_col}' is not a supported ID type for the ChEBI query "
                "class, please use on of the following: "
                f"{', '.join(self.id_columns)}"
            )
        query = \
            f"select {tgt_col} from accessions where {src_col}=? " \
            f"and source_type='KEGG COMPOUND'"
        return self._process_single_results(
            self.execute_query(query, (src_id,)))

    def kegg_to_chebi(self, kegg_id: str):
        return self._accession_query(kegg_id, "foreign_id", "chebi_id")

    def chebi_to_kegg(self, chebi_id: str):
        return self._accession_query(chebi_id.replace("CHEBI:", ""),
                                     "chebi_id", "foreign_id")

    @unique_list
    def inchi_to_chebi(self, inchi: str) -> List[str]:
        """Map from Inchi to ChEBI"""
        return self._process_single_results(
            self._query_from_to_("inchi", inchi, "inchi", "chebi_id"))

    @unique_list
    def chebi_to_inchi(self, chebi_id: str) -> List[str]:
        """Map from ChEBI to Inchi """
        return self._process_single_results(
            self._query_from_to_(
                "inchi", chebi_id.replace("CHEBI:", ""), "chebi_id", "inchi"))
