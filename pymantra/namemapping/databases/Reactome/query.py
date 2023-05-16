import pathlib
import warnings
from sqlite3 import ProgrammingError
from pymantra.namemapping.databases.sqlite_base import SQLiteBase, unique_list


class ReactomeQuery(SQLiteBase):
    """Queries to map Reactome IDs to ChEBI and NCBI"""
    def __init__(self, *args, **kwargs):
        db_file = pathlib.Path(__file__).parent.absolute() / "reactome.db"
        super(ReactomeQuery, self).__init__(db_file, *args, **kwargs)
        self.id_columns = {'reactome_id', 'chebi_id', 'ncbi_id'}

    def active_connection(self) -> bool:
        try:
            self.execute_query("select reactome_id FROM reactome LIMIT 1")
            return True
        except ProgrammingError as err:
            warnings.warn(f"Reactome database connection unavailable.\n{err}")
            return False

    @unique_list
    def reactome_to_chebi(self, reactome_id: str):
        """Map a Reactome ID to ChEBI"""
        return self._process_single_results(
            self._query_from_to_(
                "chebi", reactome_id, "reactome_id", "chebi_id")
        )

    @unique_list
    def chebi_to_reactome(self, chebi_id: str):
        """Map a ChEBI ID to Reactome"""
        return self._process_single_results(
            self._query_from_to_("chebi", chebi_id.replace("CHEBI:", ""),
                                 "chebi_id", "reactome_id")
        )

    @unique_list
    def reactome_to_ncbi(self, reactome_id: str):
        """Map a Reactome ID to NCBI"""
        return self._process_single_results(
            self._query_from_to_("ncbi", reactome_id, "reactome_id", "ncbi_id")
        )

    @unique_list
    def ncbi_to_reactome(self, ncbi_id: str):
        """Map a NCBI ID to Reactome"""
        return self._process_single_results(
            self._query_from_to_("ncbi", ncbi_id, "ncbi_id", "reactome_id")
        )
