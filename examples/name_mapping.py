import json
import pathlib

from pymantra.namemapping import metaboanalyst_name_mapping, NameMapper


# just a list of metabolite names, you can find the file here:
# TODO: add github file link
metabolites = json.load(
    open(pathlib.Path(__file__).parent.absolute() / "metabolites.json", "r"))
print(metabolites)

# getting database IDs from metabolite names
# can be skipped, if at least one ID type per metabolite is known
# this might take a few second to run
name_map = metaboanalyst_name_mapping(metabolites)

# we use HMDB IDs to map to mantra IDs
mapper = NameMapper()
mantra_ids = {
    hmdb_id: mapper.map_id(hmdb_id, "hmdb", "internal")
    for hmdb_id in name_map["HMDB"]
}
print(mantra_ids)
