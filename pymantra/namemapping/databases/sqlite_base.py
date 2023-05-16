from abc import ABC, abstractmethod
from typing import Set, Union, List
import sqlite3
import pathlib
import traceback


def unique_list(func):
    def unique(*args, **kwargs):
        qres = func(*args, **kwargs)
        if isinstance(qres, list):
            return list(set(qres))
        else:
            return qres
    return unique


class SQLiteBase(ABC):
    """Base class handling sqlite database connections"""
    db: sqlite3.Connection
    id_columns: Set[str]

    def __init__(self, file: Union[str, pathlib.Path], *args, **kwargs):
        try:
            if isinstance(file, pathlib.Path):
                self.db = sqlite3.connect(
                    f"{file.as_uri()}?mode=ro", uri=True, *args, **kwargs)
            else:
                self.db = sqlite3.connect(
                    f"file:{file}?mode=ro", uri=True, *args, **kwargs)
        except sqlite3.OperationalError:
            raise ConnectionError(
                f"Database path '{file}' is either invalid or does not "
                "contain an sqlite database"
            )

    def __del__(self):
        self.close()

    # for usage with with(subclass()) as ...
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        self.close()
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_value, tb)

    @property
    def cursor(self) -> sqlite3.Cursor:
        return self.db.cursor()

    @abstractmethod
    def active_connection(self) -> bool:
        """
        Verifying the connection to the database

        Returns
        -------
        bool
            True if connection is established, else False
        """
        raise NotImplementedError

    def close(self):
        """
        Closing the database connection. Should always be called
        once the object instance is not needed anymore
        """
        # will be False if database login is incorrect at init
        if hasattr(self, "db"):
            self.db.close()

    def execute_query(
        self, query: str, parameters: Union[list, tuple] = None
    ) -> List[tuple]:
        if parameters:
            cursor = self.db.execute(query, parameters)
        else:
            cursor = self.db.execute(query)
        return cursor.fetchall()

    def _query_from_to_(
        self, table: str, src_id: str, src_col: str, tgt_col: str
    ) -> List[tuple]:
        """
        Auxiliary function to query with a column based on the value
        of a second column

        Parameters
        ----------
        table: str
        src_id: str
        src_col: str
        tgt_col: str

        Returns
        -------
        List[tuple]

        """
        # NOTE: checks for names needs to be done outside!
        query = f"select {tgt_col} from {table} where {src_col}=?"
        return self.execute_query(query, (src_id,))

    @staticmethod
    def _process_single_results(results: List[tuple]) -> List[str]:
        return [
            result[0] for result in results
        ]


class UnknownMappingError(ValueError):
    def __init__(self, message: str):
        self.message = message
