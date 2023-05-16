from typing import NamedTuple


class Edge(NamedTuple):
    """`Edge` data type as a 2-tuple representing source and target node

    Attributes
    ----------
    source: str
        Source node name
    target: str
        Target node name
    """
    source: str
    target: str


# Convention:
#
# * edge types should have the same spelling as node types, but
#   in capital letters.
# * order of node types in EDGE_TYPE_LIST indicates edge direction except
#   for metabolite - reaction connections which can go both ways
# * node types should not contain any underscores
#
# To change the name of a relation in the database the only object that
# should be edited is 'EDGE_TYPE_LIST'
#
# The order EDGE_TYPE_LIST and NODE_TYPE_LIST should always stay the same.
#
EDGE_TYPES = {"GENE_ORGANISM", "REACTION_GENE", "REACTION_ORGANISM",
              "SUBSTRATE", "PRODUCT"}
EDGE_TYPE_LIST = ["GENE_ORGANISM", "REACTION_GENE", "REACTION_ORGANISM",
                  "SUBSTRATE", "PRODUCT"]

NODE_TYPES = {"reaction", "organism", "gene", "metabolite"}
NODE_TYPE_LIST = ["reaction", "organism", "gene", "metabolite"]

NODE_ATTRIBUTES = dict(zip(
    NODE_TYPE_LIST, [
        # TODO: adapt 'ID' attributes to match edge table requirements
        # reaction
        {
            "nodeLabel": "abbreviation",
            "Description": "description",
            "Formula": "formula",
            "KeggID": "kegg_id",
            "RheaID": "rhea",
            "ReconMap3ID": "reconMap3",
            "ReactomeID": "reactome_id",
            "isHuman": "isHuman",
            "isMicrobe": "isMicrobe"
        },
        # organism
        {
            "nodeLabel": "species",
            "Abbreviation_KEGG": "abbreviations_kegg",
            "Family": "family",
            "Genus": "genus",
            "Species": "species",
            "KeggID": "ids_kegg",
            "VMH": "organisms_vmh"
        },
        # gene
        {
            "nodeLabel": "gene_number"
        },
        # metabolite
        {
            "nodeLabel": "abbreviation",
            "Name": "fullName",
            "IUPAC": "iupac",
            "KeggID": "KEGG",
            "PubChemID": "pubchem",
            "ChEBI_ID": "CHEBI",
            "ReactomeID": "reactome_id",
            "HMDBID": "hmdb_id",
            "FooDBID": "foodb_id",
            "BiGGID": "bigg_id",
            "NCBIID": "ncbi",
            "PDBID": "pdb_id",
            "UniprotID": "uniprot",
            "EnsemblID": "ensembl",
            "miRBaseID": "mirbase"
        }
    ]
))

EDGE_ATTRIBUTES = dict(zip(
    EDGE_TYPE_LIST, [
        # organism_gene
        {
            'enrichment': 'gene_number',
            'tgt': 'species'
        },
        # gene_reaction
        {
            'enrichment': 'reaction_abbreviation',
            'tgt': 'gene_number'
        },
        # reaction_organism
        {
            'enrichment': 'reaction_abbreviation',
            'tgt': 'species'
        },
        # substrate
        {},
        # product
        {}
    ]
))

# NOTE: do NOT change! This ensures that we don't need to adapt
#       Neo4jGenerator and network Generator, even if global
#       type names are adapted
EDGE_TYPE_NAMES = dict(zip(["organism_gene", "gene_reaction",
                            "organism_reaction",
                            "substrate",
                            "product"],
                           EDGE_TYPE_LIST))
NODE_TYPE_NAMES = dict(zip(["reaction", "organism", "gene", "metabolite"],
                           NODE_TYPE_LIST))


EDGE_BY_NODE_TYPE = dict(zip(
    [
        # gene => organism
        (NODE_TYPE_NAMES['gene'], NODE_TYPE_NAMES['organism']),
        # reaction => gene
        (NODE_TYPE_NAMES['reaction'], NODE_TYPE_NAMES['gene']),
        # reaction => organism
        (NODE_TYPE_NAMES['reaction'], NODE_TYPE_NAMES['organism']),
        # substrate
        (NODE_TYPE_NAMES['metabolite'], NODE_TYPE_NAMES['reaction']),
        # product
        (NODE_TYPE_NAMES['reaction'], NODE_TYPE_NAMES['metabolite']),
    ],
    EDGE_TYPE_LIST
))


NODE_TYPES_BY_EDGE = dict(zip(
    EDGE_TYPE_LIST,
    [
        # gene => organism
        (NODE_TYPE_NAMES['gene'], NODE_TYPE_NAMES['organism']),
        # reaction => gene
        (NODE_TYPE_NAMES['reaction'], NODE_TYPE_NAMES['gene']),
        # reaction => organism
        (NODE_TYPE_NAMES['organism'], NODE_TYPE_NAMES['reaction']),
        # substrate
        (NODE_TYPE_NAMES['metabolite'], NODE_TYPE_NAMES['reaction']),
        # product
        (NODE_TYPE_NAMES['reaction'], NODE_TYPE_NAMES['metabolite']),
    ]
))


NODE_FILES = dict(zip(
    NODE_TYPE_LIST,
    ['reactions', 'organisms', 'human_catalysts', 'metabolites']
))
EDGE_FILES = dict(zip(
    EDGE_TYPE_LIST,
    ['catalyst_organism', 'reaction_catalyst', 'reaction_organism',
     'substrate_relations', 'product_relations']
))


# ================================================= #
# Edge types for subgraphs EXCLUDING reaction nodes #
# ================================================= #
DIRECT_EDGE_TYPES = {
    'metabolite_metabolite', 'gene_metabolite', 'organism_gene',
    'organism_metabolite'
}

DIRECT_EDGE_TYPE_LIST = [
    'metabolite_metabolite', 'gene_metabolite', 'organism_gene',
    'organism_metabolite'
]

DIRECT_EDGE_TYPE_NAMES = dict(zip(
    ['metabolite_metabolite', 'gene_metabolite',
     'organism_gene', 'organism_metabolite'],
    DIRECT_EDGE_TYPE_LIST
))
DIRECT_NODE_TYPES_BY_EDGE = dict(zip(
    DIRECT_EDGE_TYPE_LIST,
    [
        (NODE_TYPE_NAMES['metabolite'], NODE_TYPE_NAMES['metabolite']),
        (NODE_TYPE_NAMES['gene'], NODE_TYPE_NAMES['metabolite']),
        (NODE_TYPE_NAMES['organism'], NODE_TYPE_NAMES['gene']),
        (NODE_TYPE_NAMES['organism'], NODE_TYPE_NAMES['metabolite']),
    ]
))


# characters that raise an exception in neo4j if part of node label
forbidden_node_characters = ["'", "\"", "{", "}", "[", "]"]
