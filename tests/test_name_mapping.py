"""
pytest file to test whether ID mapping between different database IDs works
correctly and whether mapping to internal names is working as expected
"""
import pathlib
from pymantra.namemapping import NameMapper, MultiIdMapping
from typing import Dict, List
import pandas as pd

name_mapper = NameMapper()

test_file = pathlib.Path(__file__).parent.absolute() / "test_data/test_ids.tsv"
test_ids = pd.read_csv(test_file, sep="\t")
test_ids["chebi"] = test_ids["chebi"].astype(int).astype(str)


def ids_matched(control: pd.Series, test: list) -> List[bool]:
    if len(test) != control.size:
        return [False]
    return [
        control[i] in test[i] for i in range(control.size)
    ]


def many_to_many_matches(
    control: pd.DataFrame, test: Dict[str, List[MultiIdMapping]]
) -> Dict[str, List[bool]]:
    matches = {}
    for db, mappings in test.items():
        matches[db] = []
        if isinstance(mappings[0], list):
            # many_to_many with duplicates/multiple matches
            for i, maps in enumerate(mappings):
                tests = []
                for mapping in maps:
                    if mapping.from_ is None or mapping.to_ is None:
                        if control[db][i] is None:
                            # this avoids false assertions
                            tests.append(True)
                        else:
                            tests.append(False)
                    else:
                        row = control[db] == mapping.from_
                        tests.append(
                            (control.loc[row, :] == mapping.to_).sum(
                                axis=1).values > 0)
                matches[db].append(any(tests))
        else:
            # many_to_many without duplicates/multiple matches
            for i, mapping in enumerate(mappings):
                if mapping.from_ is None or mapping.to_ is None:
                    if control[db][i] is None:
                        # this avoids false assertions
                        matches[db].append(True)
                    else:
                        matches[db].append(False)
                else:
                    row = control[db] == mapping.from_
                    matches[db].append(
                        (control.loc[row, :] == mapping.to_).sum() > 0)
    return matches


# base function tests
def test_base_functions():
    with NameMapper() as nm:
        pass

    nm = NameMapper()
    _ = nm.conversion_options
    nm.print_conversion_options()


# =========================== #
# Tests: database => internal #
# These functions are mostly  #
# auxiliary ones.             #
# =========================== #
# NOTE: these are automatically also tests for `map_id`
def test_kegg_to_internal():
    mapping = [name_mapper.map_id(id_, "kegg", "internal")
               for id_ in test_ids['kegg']]
    assert all(ids_matched(test_ids['abbreviation'], mapping))


def test_reactome_to_internal():
    mapping = [name_mapper.map_id(id_, "reactome", "internal")
               for id_ in test_ids['reactome']]
    assert all(ids_matched(test_ids['abbreviation'], mapping))


def test_chebi_to_internal():
    _ = [
        name_mapper.map_id(id_, "chebi", "internal")
        for id_ in test_ids['chebi']
    ]


def test_hmdb_to_internal():
    mapping = [name_mapper.map_id(id_, "hmdb", "internal")
               for id_ in test_ids['hmdb']]
    assert all(ids_matched(test_ids['abbreviation'], mapping))


def test_vmh_to_internal():
    mapping = [name_mapper.map_id(id_, "vmh", "internal")
               for id_ in test_ids['vmh']]
    assert all(ids_matched(test_ids['abbreviation'], mapping))


# =========================== #
# Tests: many to many mapping #
# =========================== #
def test_map_to_many_with_duplicates():
    ids = {
        "hmdb": ["HMDB0003255", "HMDB0001051", "HMDB0006404"],
        "kegg": ["C00317", "C02154", "C05274"]
    }
    _ = name_mapper.map_to_many(ids, "internal", False)
    _ = name_mapper.map_to_many(ids, ["internal", "chebi"], False)

    _ = name_mapper.map_to_many(ids, "internal", True)
    _ = name_mapper.map_to_many(ids, ["internal", "chebi"], True)


# ====================================== #
# Tests: any\internal => kegg + reactome #
# This is also the main function used in #
# the pipeline                           #
# ====================================== #
def test_map_data():
    test_columns = ['chebi', 'hmdb', 'vmh', 'inchi']

    df = test_ids.copy(deep=True)
    df.index = df['abbreviation']
    mapped = name_mapper.map_data(df.loc[:, test_columns])
    assert all(mapped == df.loc[:, ['kegg', 'reactome']])

    name_mapper.map_data(df.loc[:, test_columns], remove_na=True)
