import pathlib
import warnings
from sqlite3 import ProgrammingError
from typing import List, Union, Set
from pymantra.namemapping.databases.sqlite_base import (
    SQLiteBase, unique_list, UnknownMappingError
)
from .xml_processing import ID_TAGS, TAXONOMY_TAGS


ID_COLUMNS = set(ID_TAGS.values())
ID_COLUMNS.add("hmdb")
TAXONOMY_COLUMNS = set(TAXONOMY_TAGS.values())
TAXONOMY_COLUMNS.add("hmdb")


def _single_query_table(from_: str, to_: str) -> str:
    """
    Finding the correct table to query from for a single
    column query
    TODO
    Parameters
    ----------
    from_
    to_

    Returns
    -------

    """
    if from_ in ID_COLUMNS and to_ in ID_COLUMNS:
        table = "ids"
    elif from_ in TAXONOMY_COLUMNS and to_ in TAXONOMY_COLUMNS:
        table = "taxonomy"
    else:
        raise UnknownMappingError(f"No mapping possible from {from_} to {to_}")
    return table


def _multi_query_table(from_: str, to_: List[str]) -> str:
    """
    Finding the correct table to query from for a multi-column
    query
    TODO
    Parameters
    ----------
    from_
    to_

    Returns
    -------

    """
    if all(t_ in ID_COLUMNS for t_ in to_) and from_ in ID_COLUMNS:
        table = "ids"
    elif all(t_ in TAXONOMY_COLUMNS for t_ in to_) and from_ in \
            TAXONOMY_COLUMNS:
        table = "taxonomy"
    else:
        raise UnknownMappingError(f"No mapping possible from {from_} to {to_}")
    return table


class HMDBQuery(SQLiteBase):
    """Query class to map IDs between databases"""
    def __init__(self, *args, **kwargs):
        db_file = pathlib.Path(__file__).parent.absolute() / "hmdb.db"
        super(HMDBQuery, self).__init__(db_file, *args, **kwargs)
        self.id_columns = ID_COLUMNS

    def active_connection(self) -> bool:
        try:
            self.execute_query("select name FROM hmdb LIMIT 1")
            return True
        except ProgrammingError as err:
            warnings.warn(f"HMDB database connection unavailable.\n{err}")
            return False

    @unique_list
    def get_column(self, from_: str, to_: str, match: str) -> List[str]:
        """Map from one database to another or to a HMDB hierarchy"""
        return self._process_single_results(
            self._query_from_to_(
                _single_query_table(from_, to_), match, from_, to_)
        )

    @unique_list
    def get_multiple_columns(
        self, from_: str, to_: List[str], match: str
    ) -> List[tuple]:
        """Map from one database to multiple other databases or HMDB
        hierarchies
        """
        return self._query_from_to_(
            _multi_query_table(from_, to_), match, from_, ','.join(to_))

    def taxonomy_from_foreign_id(
        self, id_colum: str, taxonomy_column: str,
        foreign_id: str
    ) -> List[str]:
        """Get a HMDB hierarchy entry from a non-HMDB ID"""
        hmdb_id = self.get_column(id_colum, 'hmdb', foreign_id)
        if len(hmdb_id) == 1:
            return self.get_column('hmdb', taxonomy_column, hmdb_id[0])
        else:
            ids = set()
            for id_ in hmdb_id:
                ids = ids.union(
                    set(self.get_column('hmdb', taxonomy_column, id_)))
            return list(ids)

    def multi_taxonomy_from_foreign_id(
        self, id_colum: str, taxonomy_columns: List[str],
        foreign_id: str
    ) -> Union[List[tuple], Set[tuple]]:
        """Get a HMDB entries for multiple hierarchies from a non-HMDB ID"""
        hmdb_id = self.get_column(id_colum, 'hmdb', foreign_id)
        if len(hmdb_id) == 1:
            return self.get_multiple_columns(
                'hmdb', taxonomy_columns, hmdb_id[0])
        else:
            return {
                self.get_multiple_columns('hmdb', taxonomy_columns, id_)
                for id_ in hmdb_id
            }
