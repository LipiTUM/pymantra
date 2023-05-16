from pymantra.database import NetworkGenerator, reduce_reaction_nodes

species = {
    "Acidipropionibacterium virtanenii", "Acidothermus cellulolyticus",
    "Azoarcus sp.", "Candidatus Pantoea", "Intestinibaculum porci",
    "Microbacterium wangchenii", "Permianibacter aggregans",
    "Plasmodium falciparum 3D7", "Pseudolysobacter antarcticus",
    "Pseudomonas viridiflava"
}
metabolites = {
    "6pgc", "6pgl", "coF420", "coF420h", "g6p", "h", "h2o", "nad", "nadp",
    "nadph"
}

generator = NetworkGenerator("bolt://127.0.0.1:7687", ("<user>", "<password>"))

edges = generator.get_reaction_subgraph(
    species, set(), metabolites, reaction_organism=("Abbreviation_KEGG", "hsa")
)

network = generator.as_networkx(
    edges=edges, reaction_subgraph=True)
network = reduce_reaction_nodes(network)
