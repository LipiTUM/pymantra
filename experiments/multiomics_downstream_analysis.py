import pathlib
import json
import pickle
import requests
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from pymantra.database import NetworkGenerator, get_auth_from_env
from pymantra.network import compute_multiomics_associations
from pymantra.plotting import (
    plot_correlation_differences, plot_correlation_averages, _plot_clustmap,
    _remove_zero_features
)

from visualisation_utils import annotate_subplots, plot_corr_db_graph


base_path = pathlib.Path(__file__).parent.parent
data_path = base_path / "experiments" / "MetabolomeMicrobiome"
figure_path = base_path / "experiments" / "Figures"
suppl_path = base_path / "experiments" / "Supplement"


def _reaction_organism_query(reactions, organisms):
    netgen = NetworkGenerator(
        "bolt://127.0.0.1:7687", get_auth_from_env(base_path / ".env"))
    query = \
        "MATCH (r:reaction)-[e:REACTION_ORGANISM]-(o:organism) " \
        "WHERE r.nodeLabel in $reactions and o.nodeLabel in $organisms " \
        "RETURN r, o"
    # TODO: split reactions into their individual parts => merge back together
    #       after query
    unlisted_nodes = set()
    node_mapping = {}
    for node in reactions:
        orig_nodes = node.split(",")
        unlisted_nodes = unlisted_nodes.union(orig_nodes)
        for on in orig_nodes:
            node_mapping[on] = node
    edges = netgen.run(
        query, reactions=list(unlisted_nodes), organisms=list(organisms)
    ).values()
    return edges, node_mapping


def requery_organisms(g: nx.Graph, organisms):
    reactions = [
        node for node, node_type in nx.get_node_attributes(g, "node_type")
        if node_type == "reaction"
    ]
    orgs = {org: " ".join(org.split(" ")[:2]) for org in organisms}
    return _reaction_organism_query(reactions, list(orgs.values()))


def get_metabolic_connections(reactions, organisms):
    edges, node_mapping = _reaction_organism_query(reactions, organisms)

    g = nx.Graph()
    g.add_nodes_from(reactions, node_type="reaction")
    g.add_nodes_from(organisms, node_type="organism")

    db_organisms = {}
    for (reaction, organism) in edges:
        db_organisms[organism.get("nodeLabel")] = organism.get("VMH")
        reac_label = reaction.get("nodeLabel")
        g.add_edge(
            node_mapping.get(reac_label, reac_label),
            organism.get("nodeLabel"),
            source="db"
        )

    remove_organisms = set()
    remove_reactions = set()
    for node, degree in g.degree:
        if degree == 0:
            if g.nodes[node]["node_type"] == "organism":
                match = False
                for org, names in db_organisms.items():
                    if node in names:
                        for _, tgt in g.edges(node):
                            g.add_edge(node, tgt)
                        match = True
                        break
                if not match:
                    short_node = " ".join(node.split(" ")[:2])
                    if short_node in g.nodes:
                        for _, tgt in g.edges(node):
                            g.add_edge(node, tgt)
                    else:
                        remove_organisms.add(node)
            else:
                remove_reactions.add(node)
    g.remove_nodes_from(remove_reactions)
    g.remove_nodes_from(remove_organisms)

    return g


def add_edge_data(df: pd.DataFrame, g: nx.Graph):
    g_ = g.copy()
    for org in df.index:
        org_sub = org.replace("_", " ")
        for reaction in df.columns:
            val = df.loc[org, reaction]
            if (org_sub, reaction) in g.edges:
                if val == 0:
                    g_.remove_edge(org_sub, reaction)
                else:
                    g_.edges[(org_sub, reaction)]["corr"] = val
                    g_.edges[(org_sub, reaction)]["weight"] = abs(val)
            elif val != 0:
                if org_sub not in g_.nodes:
                    g_.add_node(org_sub, node_type="organism")
                if reaction not in g_.nodes:
                    g_.add_node(reaction, node_type="reaction")
                g_.add_edge(
                    org_sub, reaction, corr=val, weight=abs(val),
                    source="corr"
                )

    g_.remove_nodes_from([
        node for node in g_.nodes
        if node not in df.columns and node.replace(" ", "_") not in df.index
    ])
    g_.remove_nodes_from(list(nx.isolates(g_)))

    return g_


def format_yticklabels(axis, remove_specification: bool = False):
    if remove_specification:
        return axis.set_yticklabels([
            " ".join(tick.get_text().split("_")[:2])
            for tick in axis.get_yticklabels()
        ])
    return axis.set_yticklabels([
        tick.get_text().replace("_", " ") for tick in axis.get_yticklabels()])


def query_species_enzymes(species_name: str, enzymes: list):
    # get species ID on uniprot
    species_id_query = requests.get(
        "https://rest.uniprot.org/taxonomy/search",
        params={"query": species_name}
    )
    if species_id_query.status_code != 200:
        raise ValueError(f"ID query failed for {species_name}")

    matches = [
        match["taxonId"]
        for match in json.loads(species_id_query.text)["results"]
        if match["scientificName"] == species_name
    ]
    if len(matches) == 0:
        matches = [
            match["taxonId"]
            for match in json.loads(species_id_query.text)["results"]
            if species_name in match["scientificName"]
        ]
        if len(matches) == 0:
            print(f"No results found for {species_name}")
            return [0] * len(enzymes)
        else:
            print(f"Results for {species_name} only under requery")
    if len(matches) > 1:
        print(
            f"Multiple results found for {species_name}: "
            f"{', '.join(matches)}.\nConcatenating results for all IDs."
        )

    # query enzymes using organism ID and ec codes
    match_counts = []
    for ec_code in enzymes:
        enzyme_matches = set()
        for species_id in matches:
            species_enzyme_query = requests.get(
                "https://rest.uniprot.org/uniprotkb/search",
                params={
                    "query": f"(taxonomy_id:{species_id}) AND (ec:{ec_code})"}
            )
            if species_enzyme_query.status_code != 200:
                raise ValueError(f"Enzyme query failed for {ec_code}")
            ec_matches = {
                entry["primaryAccession"]
                for entry in json.loads(species_enzyme_query.text)["results"]
            }
            enzyme_matches = enzyme_matches.union(ec_matches)
        match_counts.append(len(enzyme_matches))
    return match_counts


# main
mo_fig, mo_axes = plt.subplots(figsize=(48, 36), ncols=2, nrows=2)

# plotting all reaction - microbe correlations
corr_associations = pickle.load(
    open(data_path / "corr_associations.pickle", "rb"))
pvals = pickle.load(open(data_path / "corr_pvals.pickle", "rb"))

n_sig_cd = (pvals["CD"] < .05).sum().sum()
n_sig_ctrl = (pvals["Control"] < .05).sum().sum()
n_sig_shared = ((pvals['Control'] < .05) & (pvals['CD'] < .05)).sum().sum()
print(
    f"{n_sig_ctrl} significant correlations in Control\n"
    f"{n_sig_cd} significant correlations in CD\n"
    f"{n_sig_shared} shared significant correlations"
)

corr_diffs, _ = plot_correlation_differences(
    corr_associations, pvals, "CD", "Control", ax=mo_axes[0, 0], reorder=True,
    strip_column_names=True, return_differences=True, cmap="vlag",
    remove_all_zeros=True
)
mo_axes[0, 0].set_xticklabels([
    tick.get_text().split(", ")[0] for tick in mo_axes[0, 0].get_xticklabels()
])
format_yticklabels(mo_axes[0, 0])
mo_axes[0, 0].set_title("Correlation Differences")

corr_means, _ = plot_correlation_averages(
    corr_associations, pvals, "CD", "Control", ax=mo_axes[0, 1], reorder=True,
    strip_column_names=True, return_averages=True, cmap="vlag",
    remove_all_zeros=True
)
mo_axes[0, 1].set_xticklabels([
    tick.get_text().split(", ")[0] for tick in mo_axes[0, 1].get_xticklabels()
])
format_yticklabels(mo_axes[0, 1])
mo_axes[0, 1].set_title("Correlation Averages")


graph_file = data_path / "multiomics_graph.pickle"
if graph_file.exists():
    graph = pickle.load(open(graph_file, "rb"))
else:
    graph = get_metabolic_connections(
        list(corr_associations["CD"].columns),
        [idx.replace("_", " ") for idx in corr_associations["CD"].index]
    )
    pickle.dump(graph, open(graph_file, "wb"))

full_mean_graph = add_edge_data(corr_means, graph)
full_diff_graph = add_edge_data(corr_diffs, graph)

plot_corr_db_graph(
    full_diff_graph, ax=mo_axes[1, 0], edge_kwargs={"width": 2.0})
plot_corr_db_graph(
    full_mean_graph, ax=mo_axes[1, 1], edge_kwargs={"width": 2.0})

# finishing multi-omics figure
annotate_subplots(*mo_axes.flatten())
plt.tight_layout()
mo_fig.savefig(suppl_path / "TotalCorrelationPlots.pdf")
plt.close(mo_fig)


# plotting only top-correlations with reactions in enrichment solution
subgraph_reactions = json.load(
    open(data_path / "enrichment_results.json", "r"))["subgraph"]
res = pd.read_csv(data_path / "residual.csv", index_col=0)
microbiome_data = pd.read_csv(data_path / "microbiome_data.csv", index_col=0)
sample_groups = pd.read_csv(
    data_path / "sample_groups.csv", index_col=0)["Diagnosis"]

# TODO: plot average and different correlations for best subgraph
#   => also correlation network (edge style => whether edge exists in database)
best_corrs, best_pvals = compute_multiomics_associations(
    res.loc[:, subgraph_reactions], microbiome_data.T, sample_groups,
    comparison=("Control", "CD")
)

best_fig, axes = plt.subplots(figsize=(32, 14), ncols=2, nrows=1)

# heatmaps
best_corr_diffs, diff_clustmap = plot_correlation_differences(
    best_corrs, best_pvals, "CD", "Control", reorder=True,
    strip_column_names=True, return_differences=True, cmap="vlag",
    ax=axes[0], remove_all_zeros=True
)
format_yticklabels(axes[0])
# axes[0, 0].set_title("Correlation Difference")

# best_corr_means, mean_clustmap = plot_correlation_averages(
#     best_corrs, best_pvals, "CD", "Control", reorder=True,
#     strip_column_names=True, return_averages=True, cmap="vlag",
#     ax=axes[0, 1], remove_all_zeros=True
# )
# format_yticklabels(axes[0, 1])
# axes[0, 1].set_title("Correlation Average")
#
# # graphs
# mean_graph = add_edge_data(best_corr_means, graph)
diff_graph = add_edge_data(best_corr_diffs, graph)

plot_corr_db_graph(diff_graph, ax=axes[1], edge_kwargs={"width": 2.0})
# plot_corr_db_graph(mean_graph, ax=axes[1, 1], edge_kwargs={"width": 2.0})

plt.tight_layout()
best_fig.savefig(figure_path / "FranzosaBestCorrelation.pdf")
plt.close(best_fig)

# multiple test corrected
significant_species = {}
fig, axes = plt.subplots(figsize=(32, 28), ncols=2, nrows=2)
for group, ax in zip(["Control", "CD"], axes):
    corrs = best_corrs[group].copy()
    corrs[best_pvals[group] > .05] = 0
    corrs = _remove_zero_features(corrs)
    print(
        f"{(corrs != 0).sum().sum()} significant assocations with "
        f"{corrs.shape[0]} microbial species in {group}"
    )
    significant_species[group] = [idx.replace("_", " ") for idx in corrs.index]

    # cluster map
    _plot_clustmap(
        corrs, False, True, False, ax[0], cmap="vlag", vmin=-1, vmax=1)
    format_yticklabels(ax[0])
    ax[0].set_title(group)
    # graph
    group_graph = add_edge_data(corrs, graph)
    plot_corr_db_graph(group_graph, ax=ax[1], edge_kwargs={"width": 2.0})

annotate_subplots(*axes.flatten())

plt.tight_layout()
fig.savefig(figure_path / "MultiOmicsCorrelations.pdf")
plt.close(fig)


# uniprot queries
ec_codes = ["1.4.1.14", "2.6.1.57", "2.6.1.88", "2.6.1.42"]
group_enzyme_matches = {
    group: {
        species: query_species_enzymes(species, ec_codes)
        for species in group_species
    }
    for group, group_species in significant_species.items()
}
enzyme_match_dfs = []
for group, group_matches in group_enzyme_matches.items():
    df = pd.DataFrame(group_matches).T
    df.columns = ec_codes
    df["SignificanceGroup"] = group
    enzyme_match_dfs.append(df)

pickle.dump(enzyme_match_dfs, open(data_path / "enzyme_matches.pickle", "wb"))
pd.concat(enzyme_match_dfs).to_csv(suppl_path / "enzyme_matches.csv")
