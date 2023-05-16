import warnings
import traceback
from typing import Dict, Tuple, List, Set, Union
from collections import namedtuple
import numpy as np
import pandas as pd
from .databases.HMDB.query import HMDBQuery
from .databases.Reactome.query import ReactomeQuery
from .databases.ChEBI.query import ChEBIQuery
from .databases.mantra_db.query import MantraDBQuery
from .databases.sqlite_base import UnknownMappingError


MultiIdMapping = namedtuple("MultiIdMatch", "from_ to_")
MetaboliteIdentification = namedtuple(
    "MetaboliteIdentification", "kegg reactome"
)


class NameMapper:
    """Metabolite ID mapping class

    Mapping between HMDB, ChEBI, NCBI, KEGG, Reactome, InChI, SMILES, VMH and
    the internal database IDs. The sources of mapping are coming from the
    resources themselves.

    Parameters
    ----------
    hmdb : HMDBQuery
        Interface to query from or to HMDB IDs
    chebi : ChEBIQuery
        Interface to query from or to ChEBI IDs
    reactome : ReactomeQuery
        Interface to query from or to Reactome IDs
    mantra_db : MantraDBQuery
        Interface to query from or to mantra-internal IDs
    query_functions : Dict[Tuple[str, str], callable]
        Dictionary mapping (source ID type, target ID type) to the correct
        mapping function
    """
    hmdb: HMDBQuery
    chebi: ChEBIQuery
    reactome: ReactomeQuery
    mantra_db: MantraDBQuery
    query_functions: Dict[Tuple[str, str], callable]

    def __init__(self, **sqlite3_args):
        """Construct NameMaper instance

        Parameters
        ----------
        sqlite3_args:
            Optional keywora arguments to be passed to
            :py:func:`sqlite3.connect` for all database connections
        """
        # setting up database connections
        self.hmdb = HMDBQuery(**sqlite3_args)
        self.chebi = ChEBIQuery(**sqlite3_args)
        self.reactome = ReactomeQuery(**sqlite3_args)
        self.mantra_db = MantraDBQuery(**sqlite3_args)
        # used to match query to function and database
        self._query_functions = {
            ("reactome", "chebi"): self.reactome.reactome_to_chebi,
            ("chebi", "reactome"): self.reactome.chebi_to_reactome,
            ("reactome", "ncbi"): self.reactome.reactome_to_ncbi,
            ("ncbi", "reactome"): self.reactome.ncbi_to_reactome,
            ("chebi", "kegg"): self.chebi.chebi_to_kegg,
            ("kegg", "chebi"): self.chebi.kegg_to_chebi,
            ("chebi", "inchi"): self.chebi.chebi_to_inchi,
            ("inchi", "chebi"): self.chebi.inchi_to_chebi,
            ("reactome", "kegg"): self.reactome_to_kegg,
            ("kegg", "reactome"): self.kegg_to_reactome,
            ("reactome", "inchi"): self.reactome_to_inchi,
            ("inchi", "reactome"): self.inchi_to_reactome,
            # hmdb queries
            ("hmdb", "inchi"): lambda x: self.hmdb.get_column(
                "hmdb", "inchi", x),
            ("inchi", "hmdb"): lambda x: self.hmdb.get_column(
                "inchi", "hmdb", x),
            ("hmdb", "chebi"): lambda x: self.hmdb.get_column(
                "hmdb", "chebi", x),
            ("chebi", "hmdb"): lambda x: self.hmdb.get_column(
                "chebi", "hmdb", x),
            ("hmdb", "pubchem"): lambda x: self.hmdb.get_column(
                "hmdb", "pubchem", x),
            ("pubchem", "hmdb"): lambda x: self.hmdb.get_column(
                "pubchem", "hmdb", x),
            ("hmdb", "smiles"): lambda x: self.hmdb.get_column(
                "hmdb", "smiles", x),
            ("smiles", "hmdb"): lambda x: self.hmdb.get_column(
                "smiles", "hmdb", x),
            ("hmdb", "kegg"): lambda x: self.hmdb.get_column(
                "hmdb", "kegg", x),
            ("kegg", "hmdb"): lambda x: self.hmdb.get_column(
                "kegg", "hmdb", x),
            ("hmdb", "vmh"): lambda x: self.hmdb.get_column(
                "hmdbid", "vmh", x),
            ("vmh", "hmdb"): lambda x: self.hmdb.get_column(
                "vmh", "hmdb", x),
            # mapping from and to internal ids
            ("inchi", "internal"): self.inchi_to_internal,
            ("kegg", "internal"): self.mantra_db.kegg_to_internal,
            ("reactome", "internal"): self.mantra_db.reactome_to_internal,
            ("chebi", "internal"): self.mantra_db.chebi_to_internal,
            ("hmdb", "internal"): self.mantra_db.hmdb_to_internal,
            ("vmh", "internal"): self.mantra_db.vmh_to_internal,
            ("internal", "kegg"): self.mantra_db.internal_to_kegg,
            ("internal", "reactome"): self.mantra_db.internal_to_reactome,
            ("internal", "chebi"): self.mantra_db.internal_to_chebi,
            ("internal", "hmdb"): self.mantra_db.internal_to_hmdb,
            ("internal", "vmh"): self.mantra_db.internal_to_vmh,
        }

    def close(self):
        """
        Ensuring all database connections are closed
        """
        if hasattr(self, 'hmdb'):
            self.hmdb.close()
        if hasattr(self, 'chebi'):
            self.chebi.close()
        if hasattr(self, 'reactome'):
            self.reactome.close()
        if hasattr(self, 'mantra_db'):
            self.mantra_db.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        if exc_type is not None:
            traceback.print_exception(exc_type, exc_val, exc_tb)

    def __del__(self):
        pass
        # self.close()

    @property
    def conversion_options(self):
        """List all conversion options

        Returns
        -------
        List[str]
            List of conversion options
        """
        return list(self._query_functions.keys())

    def print_conversion_options(self):
        """Print all conversion options"""
        for option in self._query_functions.keys():
            print(f"{option[0]} to {option[1]}")

    @staticmethod
    def _intermediate_mapping_(id_: str, inter: callable, fin: callable):
        int_id = inter(id_)
        ids = []
        for id_ in int_id:
            mapped_ids = fin(id_)
            for iid in mapped_ids:
                ids.append(iid)
        return ids

    def reactome_to_kegg(self, reactome_id: str) -> List[str]:
        """Map Reactome ID to KEGG ID(s)

        Parameters
        ----------
        reactome_id: str
            Reactome ID

        Returns
        -------
        List[str]
            List of mapped KEGG IDs
        """
        return self._intermediate_mapping_(
            reactome_id, self.reactome.reactome_to_chebi,
            self.chebi.chebi_to_kegg
        )

    def kegg_to_reactome(self, kegg_id: str) -> List[str]:
        """Map KEGG ID to Reactome ID

        Parameters
        ----------
        kegg_id: str
            KEGG ID

        Returns
        -------
        List[str]
            List of mapped Reactome IDs
        """
        return self._intermediate_mapping_(
            kegg_id, self.chebi.kegg_to_chebi,
            self.reactome.chebi_to_reactome
        )

    def reactome_to_inchi(self, reactome_id: str) -> List[str]:
        """Map Reactome ID to InCHI key

        Parameters
        ----------
        reactome_id: str
            Reactome ID

        Returns
        -------
        List[str]
            List of mapped InCHI keys
        """
        return self._intermediate_mapping_(
            reactome_id, self.reactome.reactome_to_chebi,
            self.chebi.chebi_to_inchi
        )

    def inchi_to_reactome(self, inchi_id: str) -> List[str]:
        """Map InCHI key to Reactome

        Parameters
        ----------
        inchi_id: str
            InCHI ID

        Returns
        -------
        List[str]
            Mapped Reactome ID
        """
        return self._intermediate_mapping_(
            inchi_id, self.chebi.inchi_to_chebi,
            self.reactome.chebi_to_reactome
        )

    def inchi_to_internal(self, inchi_id: str) -> List[str]:
        """Map InCHI key to mantra-internal ID

        Parameters
        ----------
        inchi_id : str
            InCHI key

        Returns
        -------
        List[str]
            Mapped internal ID
        """
        return [
            intern for reactome_id in self.inchi_to_reactome(inchi_id)
            for intern in self.mantra_db.reactome_to_internal(reactome_id)
        ]

    def map_id(
        self, id_: str, id_type: str,
        map_to: Union[str, List[str]],
        **kwargs
    ) -> Union[List[str], List[tuple], Set[tuple]]:
        """Wrapper for name mapping functions. Takes an ID from a supported
        database and converts it to corresponding identifiers from other
        databases.

        Parameters
        ----------
        id_: str
            ID to map
        id_type: str
            database/ID type from which `id_` is originating
        map_to: Union[str, List[str]]
            String or list of strings specifying to which databases `id_`
            should be mapped to

        Returns
        -------
        Union[List[str], List[tuple], Set[tuple]]
            List of strings is `map_to` is a string, where each element
            represents a match with the target database.
            List or set of tuples, if `map_to` is a list. Each element
            represents a mapping, where the first element is the mapped ID and
            the second element is the database from which this ID is coming

        Examples
        --------
        >>> from pymantra.namemapping import NameMapper
        >>> name_map = NameMapper()
        >>>
        >>> hmdb_ids = ["HMDB0003255", "HMDB0001051", "HMDB0006404"]
        >>> hmdb_mapping = [name_map.map_id(id_, "hmdb", "internal")]
        >>>
        >>> kegg_ids = ["C00317", "C02154", "C05274"]
        >>> for id_ in kegg_ids:
        >>>     print(name_map.map_id(id_, "kegg", "internal"))
        """
        id_type = id_type.lower()
        map_to = map_to.lower()
        if isinstance(map_to, str):
            query_fun = self._query_functions.get((id_type, map_to))
            if query_fun:
                return query_fun(id_, **kwargs)
            try:
                return self.hmdb.taxonomy_from_foreign_id(id_type, map_to, id_)
            except UnknownMappingError:
                warnings.warn(f"No mapping found from {id_type} to {map_to}")
                return []
        try:
            return self.hmdb.multi_taxonomy_from_foreign_id(
                id_type, map_to, id_
            )
        except UnknownMappingError:
            warnings.warn(f"No mapping found from {id_type} to {map_to}")
            return []

    def map_to_many(
        self, ids: Dict[str, List[str]], map_to: Union[str, List[str]],
        remove_duplicates: bool = True
    ) -> Union[Dict[str, List[List[str]]],
               Dict[str, List[str]],
               Dict[str, List[Tuple[str, str]]]]:
        """Mapping multiple entries of multiple databases onto multiple
        other databases.

        Parameters
        ----------
        ids: Dict[str, List[str]]
            All ids to query (values) by id type (keys)
        map_to: Union[str, List[str]]
            ID type to map
        remove_duplicates: bool, default True
            If True only the first match will be returned if multiple matches
            are found

        Returns
        -------
        Union[Dict[str, List[List[str]]],
              Dict[str, List[str]],
              Dict[str, List[Tuple[str, str]]]]
        """
        mapping = {}
        if isinstance(map_to, str):
            map_to = map_to.lower()
            for src_type, ids in ids.items():
                src_type = src_type.lower()
                qfun = self._query_functions.get((src_type, map_to))
                if qfun is None:
                    warnings.warn(
                        f"Mapping from {src_type} to {map_to} is not "
                        f"implemented. "
                        f"For a full list of available options call "
                        f"`NameMapper.print_conversion_options`"
                    )
                mapping[src_type] = []
                for id_ in ids:
                    mapped = qfun(id_)
                    if mapped and remove_duplicates:
                        mapping[src_type].append(mapped[0])
                    else:
                        mapping[src_type].append(mapped)
        else:
            for src_type, ids in ids.items():
                src_type = src_type.lower()
                mappable_tgts = [
                    tgt_type for tgt_type in map_to
                    if self._query_functions.get((src_type, tgt_type.lower()))
                ]
                if not mappable_tgts:
                    warnings.warn(
                        f"Mapping from {src_type} to any of the id types "
                        f"specified in `map_to` is not implemented. "
                        f"For a full list of available options call "
                        f"`NameMapper.print_conversion_options`"
                    )
                mapping[src_type] = []
                match = None
                for id_ in ids:
                    for tgt_type in mappable_tgts:
                        mapped = self._query_functions.get(
                            (src_type, tgt_type.lower()))(id_)
                        if mapped:
                            if remove_duplicates:
                                match = MultiIdMapping(id_, mapped[0])
                            else:
                                match = [MultiIdMapping(id_, mapped_)
                                         for mapped_ in mapped]
                        break
                    if match:
                        mapping[src_type].append(match)
                    else:
                        if remove_duplicates:
                            mapping[src_type].append(
                                MultiIdMapping(None, None)
                            )
                        else:
                            mapping[src_type].append(
                                [MultiIdMapping(None, None)]
                            )
                    match = None
        return mapping

    def _to_kegg_reactome(
        self, options: Union[dict, np.ndarray], dbs: pd.Index = None
    ) -> MetaboliteIdentification:
        kegg = None
        reactome = None
        if isinstance(options, dict):
            iterator = options.items()
        else:
            iterator = zip(dbs, options)
        for db, id_ in iterator:
            tmp_kegg = self.map_id(id_, db, "kegg")
            if tmp_kegg:
                kegg = tmp_kegg[0]
            tmp_reactome = self.map_id(id_, db, "reactome")
            if tmp_reactome:
                reactome = tmp_reactome[0]
            if kegg and reactome:
                break
        return MetaboliteIdentification(kegg, reactome)

    def map_data(
        self, data: Union[pd.DataFrame, List[Dict[str, str]]],
        remove_na: bool = False
    ) -> Union[pd.DataFrame, List[MetaboliteIdentification]]:
        """Maps multiple metabolites with one or multiple database identifiers
        to KEGG and Reactome IDs

        Parameters
        ----------
        data: Union[pd.DataFrame, List[Dict[str, str]]
            database identifiers to map.
            If a pandas DataFrame, each row represents a metabolite
            and each column a database (i.e. each cell is a database
            identifier) If a list, each list entry represents a metabolite
            and each key - value pair a database - ID pair.

        remove_na: bool, default False
            If True metabolites for which no match was found are removed.
            This option is only considered when `data` is a DataFrame.

        Returns
        -------
        Union[pd.DataFrame, List[MetaboliteIdentification]]
            If `data` is a DataFrame, so is the output. Columns in this case
            are 'kegg' and 'reactome' and index is the same as data.index.
            Else a list of 2-tuples, where each tuple represents the database
            ids for KEGG (0, kegg) and Reactome (1, reactome) found for the
            respective input.
        """
        # TODO: reporting for source id/db of matches?
        if isinstance(data, pd.DataFrame):
            # NOTE: this could also be done using iterrows using the same
            #       scheme as the list case. Avoided for performance reasons.
            mapped = [
                self._to_kegg_reactome(
                    data.values[i, :], data.columns)._asdict()
                for i in range(data.shape[0])
            ]
            mapped = pd.DataFrame(mapped, index=data.index)
            if remove_na:
                kegg_mask = (mapped['kegg'] != 'None')
                reactom_mask = (mapped['reactome'] != 'None')
                mapped = mapped.loc[kegg_mask | reactom_mask, :]
        else:
            mapped = [self._to_kegg_reactome(row) for row in data]
        return mapped
