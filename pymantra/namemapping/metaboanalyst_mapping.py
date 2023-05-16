import json
import requests
import pandas as pd


MA_URL = "http://api.xialab.ca/mapcompounds"
HEADERS = {
    "Content-Type": "application/json",
    "cache-control": "no-cache"
}


def metaboanalyst_name_mapping(
    metabolites, source_type: str = "name"
) -> pd.DataFrame:
    """Convert a metabolite name or database ID to other IDs

    Convert between different metabolite naming conventions and database ID
    types to obtain a unified or mappable format.

    Parameters
    ----------
    metabolites: Iterable[str]
        An iterable of metabolite strings to query
    source_type: str, "name"
        What type of identifier `metabolites` specifies. This must be
        compatible with the Metaboanalyst API. Currently available options are
        (adapted from https://www.metaboanalyst.ca/docs/APIs.xhtml):
        "name" - Compound name (e.g., 1,3-Diaminopropane)
        "hmdb" - Human Metabolome Database (e.g., HMDB0000002)
        "pubchem" - PubChem Substance and Compound databases(e.g., 428)
        "chebi" - Chemical Entities of Biological Interest(e.g., 15725)
        "metlin" - Metabolite and Chemical Entity Database (e.g., 5081)
        "kegg" - KEGG COMPOUND Database (e.g., C00986)


    Returns
    -------
    pd.DataFrame
        DataFrame where the indices are the strings given in `metabolites` and
        the columns are different identifier types
    """
    dstr = "{\n\t\"queryList\": \"" \
           f"{';'.join(metabolites)}\",\n\t\"inputType\": " \
           f"\"{source_type}\"\n}}"
    ma_query = requests.request(
        "POST", MA_URL, data=dstr.encode("utf-8"), headers=HEADERS)
    if not ma_query.ok:
        raise ValueError(
            f"Unexpected error in metaboanalyst query: {ma_query.reason}")
    ma_mapping = pd.DataFrame.from_dict(json.loads(ma_query.text))
    ma_mapping.index = metabolites
    return ma_mapping
