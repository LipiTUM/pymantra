import pathlib
import warnings
from sqlite3 import ProgrammingError
from typing import List, Union
from pymantra.namemapping.databases.sqlite_base import SQLiteBase, unique_list


class MantraDBQuery(SQLiteBase):
    """Query class to map mantra IDs from and to other databases"""
    def __init__(self, *args, **kwargs):
        db_file = pathlib.Path(__file__).parent.absolute() / "mantra.db"
        super(MantraDBQuery, self).__init__(db_file, *args, **kwargs)
        self.id_columns = {
            'internal_name', 'kegg', 'reactome', 'chebi', 'hmdb', 'vmh'}

    def active_connection(self) -> bool:
        try:
            self.execute_query("select kegg from mmnet limit 1")
            return True
        except ProgrammingError as err:
            warnings.warn(f"mantra database connection unavailable.\n{err}")
            return False

    def _query(self, id_, src_column, tgt_column) -> List[str]:
        if src_column not in self.id_columns:
            raise ValueError(
                f"'{src_column}' is not a supported ID type for the mantra DB "
                "query class, please use on of the following: "
                f"{', '.join(self.id_columns)}"
            )
        if tgt_column not in self.id_columns:
            raise ValueError(
                f"'{tgt_column}' is not a supported ID type for the mantra DB "
                "query class, please use on of the following: "
                f"{', '.join(self.id_columns)}"
            )
        query = f"select {tgt_column} from metabolites where {src_column}=?"
        return self._process_single_results(self.execute_query(query, (id_,)))

    @unique_list
    def _id_to_internal_(
        self, id_: str, column: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        if not single_value:
            return self._query(id_, column, "internal_name")
        res = self._query(id_, column, "internal_name")
        if res:
            return res[0]
        return None

    def kegg_to_internal(
        self, kegg_id: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        """Map a KEGG ID to a mantra ID"""
        return self._id_to_internal_(kegg_id, 'kegg', single_value)

    def reactome_to_internal(
        self, reactome_id: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        """Map a Reactome ID to a mantra ID"""
        return self._id_to_internal_(reactome_id, 'reactome', single_value)

    def chebi_to_internal(
        self, chebi_id: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        """Map a ChEBI ID to a mantra ID"""
        return self._id_to_internal_(chebi_id, 'chebi', single_value)

    def hmdb_to_internal(
        self, hmdb_id: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        """Map a HMDB ID to a mantra ID"""
        return self._id_to_internal_(hmdb_id, 'hmdb', single_value)

    def vmh_to_internal(
        self, vmh_id: str, single_value: bool = False
    ) -> Union[List[str], str, None]:
        """Map a VMH ID to a mantra ID"""
        return self._id_to_internal_(vmh_id, 'vmh', single_value)

    @unique_list
    def internal_to_kegg(self, internal_name: str) -> List[str]:
        """Map a mantra ID to KEGG"""
        return self._query(internal_name, "internal_name", "kegg")

    @unique_list
    def internal_to_reactome(self, internal_name: str) -> List[str]:
        """Map a mantra ID to Reactome"""
        return self._query(internal_name, "internal_name", "reactome")

    @unique_list
    def internal_to_chebi(self, internal_name: str) -> List[str]:
        """Map a mantra ID to ChEBI"""
        return self._query(internal_name, "internal_name", "kegg")

    @unique_list
    def internal_to_hmdb(self, internal_name: str) -> List[str]:
        """Map a mantra ID to HMDB"""
        return self._query(internal_name, "internal_name", "hmdb")

    @unique_list
    def internal_to_vmh(self, internal_name: str) -> List[str]:
        """Map a mantra ID to VMH"""
        return self._query(internal_name, "internal_name", "vmh")
