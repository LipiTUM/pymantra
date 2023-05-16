try:
    from lxml import etree
    lxml_import = True
except ImportError:
    lxml_import = False
from typing import Tuple, Dict
from tqdm import tqdm


SUB_TAGS = {"secondary_accessions", "synonyms"}


ID_TAGS = {
    "secondary_accessions": "secondary_accession",
    "name": "metabolite_name",
    "synonyms": "synonyms",
    "inchikey": "inchi_key",
    "inchi": "inchi",
    "smiles": "smiles",
    "cas_registry_number": "cas_registry",
    "iupac_name": "iupac",
    "kegg_id": "kegg",
    "vmh_id": "vmh",
    "pubchem_compound_id": "pubchem",
    "foodb_id": "foodb",
    "chebi_id": "chebi"
}


TAXONOMY_TAGS = {
    "direct_parent": "direct_parent",
    "kingdom": "kingdom",
    "super_class": "super_class",
    "class": "class",
    "sub_class": "sub_class",
    "molecular_framework": "molecular_framework"
}


RELEVANT_TAGS = set(ID_TAGS.keys()).union(set(TAXONOMY_TAGS.keys()))
RELEVANT_TAGS.add("accession")
RELEVANT_TAGS.add("synonym")
RELEVANT_TAGS.add("metabolite")


def parse_tag(tag: str):
    return tag.replace("{http://www.hmdb.ca}", "")


if lxml_import:
    def parse_xml_file(
        file: str
    ) -> Tuple[Dict[str, Dict[str, str]], Dict[str, Dict[str, str]]]:
        iterator = etree.iterparse(file, events=('start', 'end'))
        _, root = next(iterator)
        ids = {}
        taxonomy = {}
        current_accession = None
        current_tag = None
        accession_data = {}
        for event, element in tqdm(iterator, desc="Parsing", total=147798450):
            tag = parse_tag(element.tag)
            if tag not in RELEVANT_TAGS:
                continue
            if event == "start":
                if current_tag not in SUB_TAGS:
                    if tag in SUB_TAGS:
                        current_tag = tag
                    elif tag == "accession":
                        current_accession = element.text
                    elif tag in TAXONOMY_TAGS.keys() or tag in ID_TAGS.keys():
                        accession_data[tag] = element.text
                elif accession_data.get(current_tag):
                    accession_data[current_tag] += f",{element.text}"
                else:
                    accession_data[current_tag] = element.text
            else:
                if tag == "metabolite":
                    if current_accession is not None:
                        # this means we reached the end of the metabolite entry
                        ids[current_accession] = {
                            sql_key: accession_data.get(hmdb_key, "")
                            for hmdb_key, sql_key in ID_TAGS.items()
                        }
                        taxonomy[current_accession] = {
                            sql_key: accession_data.get(hmdb_key, "")
                            for hmdb_key, sql_key in TAXONOMY_TAGS.items()
                        }
                        current_accession = None
                        current_tag = None
                        accession_data = {}
                elif tag in SUB_TAGS:
                    current_tag = None
            element.clear()
        return ids, taxonomy
else:
    def parse_xml_file(file: str):
        raise ImportError(
            "'lxml' package could not be import but is required parse HMDB "
            "data. Please install the package to generate the HMDB mapping "
            "table."
        )


def _write_table_(
    file: str, tags: Dict[str, str], id_data: Dict[str, Dict[str, str]],
    delim: str = "\t"
):
    with open(file, 'w') as outfile:
        for accession, accession_data in tqdm(id_data.items(),
                                              desc="Writing..."):
            outfile.write(accession)
            for column in tags.values():
                outfile.write(delim)
                cell = accession_data.get(column, 'None')
                if cell:
                    if len(str(cell)) > 6000:
                        cell = str(cell)[0:6000]
                    outfile.write(cell)
                else:
                    outfile.write('None')
            outfile.write("\n")


def write_id_table(
    file: str, id_data: Dict[str, Dict[str, str]], delim: str = "\t"
):
    _write_table_(file, ID_TAGS, id_data, delim)


def write_taxonomy_table(
    file: str, id_data: Dict[str, Dict[str, str]], delim: str = "\t"
):
    _write_table_(file, TAXONOMY_TAGS, id_data, delim)


if __name__ == "__main__":
    # NOTE: the .xml file can be downloaded from https://hmdb.ca/downloads
    id_data, tax_data = parse_xml_file("hmdb_metabolites.xml")
    write_id_table("id_table.tsv", id_data)
    write_taxonomy_table("taxonomy_table.tsv", tax_data)
