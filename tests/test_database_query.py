import os
import pathlib
import warnings
import networkx as nx
from pymantra.statics import Edge
from pymantra.database import read_env
from pymantra.plotting import (
    plot_undirected_graph,
    plot_directed_graph
)
from pymantra.database.query import NetworkGenerator
from pymantra.statics import (
    NODE_TYPE_NAMES,
    EDGE_TYPE_NAMES,
    DIRECT_EDGE_TYPE_NAMES
)


read_env(pathlib.Path(__file__).parent.parent / ".env")
AUTH = (os.getenv("NEO4J_USER"), os.getenv("NEO4J_PASSWORD"))
PORT = "bolt://127.0.0.1:7687"


def _check_db_availability(func):
    """Helper decorator to avoid crashing when database is unavailable"""
    def _has_connection(*args):
        if args[0].__test__:
            return func(*args)
        return None

    return _has_connection


# ==================================================== #
# pytest test suites to check whether subgraph queries #
# work as expected                                     #
# ==================================================== #
class TestNeo4j:
    try:
        gen = NetworkGenerator(PORT, AUTH)
        __test__ = True
    except ConnectionError:
        warnings.warn(
            "neo4j service not running or invalid credentials, skipping "
            "database tests"
        )
        # avoiding tests to crash and exit immediately
        __test__ = False

    sub_organisms = {'hsa'}
    sub_genes = {'G6PDH', 'PGLS'}
    sub_metabolites = {'G6P', '6PGL', '6PG'}
    expected_reactions = {'G6P=>6PGL', '6PGL=>6PG'}
    expected_reaction_edges = {
        EDGE_TYPE_NAMES['organism_gene']: {Edge('hsa', 'G6PDH'),
                                           Edge('hsa', 'PGLS')},
        EDGE_TYPE_NAMES['gene_reaction']: {Edge('G6PDH', 'G6P=>6PGL'),
                                           Edge('PGLS', '6PGL=>6PG')},
        EDGE_TYPE_NAMES['organism_reaction']: {Edge('hsa', 'G6P=>6PGL'),
                                               Edge('hsa', '6PGL=>6PG')},
        EDGE_TYPE_NAMES['substrate']: {Edge('G6P', 'G6P=>6PGL'),
                                       Edge('6PGL', '6PGL=>6PG')},
        EDGE_TYPE_NAMES['product']: {Edge('G6P=>6PGL', '6PGL'),
                                     Edge('6PGL=>6PG', '6PG')}
    }
    expected_edges = {
        DIRECT_EDGE_TYPE_NAMES['organism_metabolite']: {
            Edge('hsa', 'G6P'), Edge('hsa', '6PGL'), Edge('hsa', '6PG')
        },
        DIRECT_EDGE_TYPE_NAMES['gene_metabolite']: {
            Edge('PGLS', '6PG'), Edge('PGLS', '6PGL'),
            Edge('G6PDH', '6PGL'), Edge('G6PDH', 'G6P')
        },
        DIRECT_EDGE_TYPE_NAMES['organism_gene']: {
            Edge('hsa', 'PGLS'), Edge('hsa', 'G6PDH')
        },
        DIRECT_EDGE_TYPE_NAMES['metabolite_metabolite']: {
            Edge('G6P', '6PGL'), Edge('6PGL', '6PG')
        }
    }

    @_check_db_availability
    def test_base_functions(self):
        self.gen.valid_port()
        self.gen.active_connection()

        with NetworkGenerator(PORT, AUTH) as generator:
            generator.active_connection()
        NetworkGenerator(PORT, AUTH).close()

    # @_check_db_availability
    # def test_db_stats(self):
    #     _ = self.gen.nodes
    #     _ = self.gen.nodes
    #     _ = self.gen.edges
    #     _ = self.gen.edges
    #     _ = self.gen.neo4j_nodes

    @_check_db_availability
    def test_reaction_organism_subset(self):
        subgraph = self.gen.get_reaction_subgraph(
            set(), set(), self.sub_metabolites,
            reaction_organism=("Abbreviation_KEGG", "hsa")
        )
        for edge_type, edges_ in subgraph.items():
            matches = [
                edge in self.expected_reaction_edges[edge_type]
                for edge in edges_
            ]
            assert all(matches)

    @_check_db_availability
    def test_reaction_subgraph_wgenes(self):
        """
        Test whether the subgraph including reaction nodes is
        turning out as expected if genes are included
        """
        subgraph = self.gen.get_reaction_subgraph(
            self.sub_organisms, self.sub_genes, self.sub_metabolites)
        for edge_type, edges_ in subgraph.items():
            matches = [edge in self.expected_reaction_edges[edge_type]
                       for edge in edges_]
            if not all(matches):
                print(
                    edge_type, edges_, self.expected_reaction_edges[edge_type],
                    sep="\n"
                )
            assert all(matches)

    @_check_db_availability
    def test_reaction_subgraph_wogenes(self):
        """
        Test whether the subgraph including reaction nodes is
        turning out as expected if genes are not included
        """
        subgraph = self.gen.get_reaction_subgraph(
            self.sub_organisms, set(), self.sub_metabolites)
        for edge_type, edges_ in subgraph.items():
            matches = [edge in self.expected_reaction_edges[edge_type]
                       for edge in edges_]
            if not all(matches):
                print(
                    edge_type, edges_, self.expected_reaction_edges[edge_type],
                    sep="\n"
                )
            assert all(matches)

    @_check_db_availability
    def test_subgraph_wgenes(self):
        """
        Test whether the subgraph omitting reaction nodes is
        turning out as expected if genes are included
        """
        subgraph = self.gen.get_subgraph(
            self.sub_organisms, self.sub_genes, self.sub_metabolites)
        for edge_type, edges_ in subgraph.items():
            matches = [
                edge in self.expected_edges[edge_type] for edge in edges_]
            if not all(matches):
                print(
                    edge_type, "\n", edges_, "\n",
                    self.expected_edges[edge_type]
                )
            assert all(matches)

    @_check_db_availability
    def test_subgraph_wogenes(self):
        """
        Test whether the subgraph omitting reaction nodes is
        turning out as expected if genes are not included
        """
        subgraph = self.gen.get_subgraph(
            self.sub_organisms, set(), self.sub_metabolites)
        for edge_type, edges_ in subgraph.items():
            matches = [
                edge in self.expected_edges[edge_type] for edge in edges_]
            if not all(matches):
                print(
                    edge_type, "\n", edges_, "\n",
                    self.expected_edges[edge_type]
                )
            assert all(matches)

    @_check_db_availability
    def test_to_networkx(self):
        """
        Converting the database graph to a networkx object.
        Checks are whether nodes and edges are the same and whether
        they are assigned the correct types.
        """
        nx_graph = self.gen.as_networkx(
            {
                NODE_TYPE_NAMES['metabolite']: self.sub_metabolites,
                NODE_TYPE_NAMES['organism']: self.sub_organisms,
                NODE_TYPE_NAMES['gene']: self.sub_genes,
            },
            self.expected_edges
        )
        for node, node_type in \
                nx.get_node_attributes(nx_graph, 'node_type').items():
            if node_type == NODE_TYPE_NAMES['metabolite']:
                assert node in self.sub_metabolites
            elif node_type == NODE_TYPE_NAMES['organism']:
                assert node in self.sub_organisms
            elif node_type == NODE_TYPE_NAMES['gene']:
                assert node in self.sub_genes
            else:
                raise ValueError(f"Unknown node type {node_type} for {node}")
        for edge, edge_type in \
                nx.get_edge_attributes(nx_graph, 'edge_type').items():
            assert edge in self.expected_edges[edge_type]

    @_check_db_availability
    def test_to_networkx_reaction(self):
        """
        Converting the database graph to a networkx object.
        Checks are whether nodes and edges are the same and whether
        they are assigned the correct types.
        """
        nx_graph = self.gen.as_networkx(
            {
                NODE_TYPE_NAMES['metabolite']: self.sub_metabolites,
                NODE_TYPE_NAMES['organism']: self.sub_organisms,
                NODE_TYPE_NAMES['gene']: self.sub_genes,
                NODE_TYPE_NAMES['reaction']: set()
            },
            self.expected_reaction_edges, reaction_subgraph=True
        )
        for node, node_type in \
                nx.get_node_attributes(nx_graph, 'node_type').items():
            if node_type == NODE_TYPE_NAMES['metabolite']:
                assert node in self.sub_metabolites
            elif node_type == NODE_TYPE_NAMES['organism']:
                assert node in self.sub_organisms
            elif node_type == NODE_TYPE_NAMES['gene']:
                assert node in self.sub_genes
            else:
                raise ValueError(f"Unknown node type {node_type} for {node}")
        for edge, edge_type in \
                nx.get_edge_attributes(nx_graph, 'edge_type').items():
            assert edge in self.expected_reaction_edges[edge_type]

    def test_as_networkx_options(self):
        nodes = {
            NODE_TYPE_NAMES['metabolite']: self.sub_metabolites,
            NODE_TYPE_NAMES['organism']: self.sub_organisms,
            NODE_TYPE_NAMES['gene']: self.sub_genes,
            NODE_TYPE_NAMES['reaction']: set()
        }
        self.gen.as_networkx(nodes=nodes, edges=None, include_attributes=False)
        self.gen.as_networkx(
            nodes=nodes, edges=None, include_attributes=False,
            reaction_subgraph=True
        )
        self.gen.as_networkx(
            nodes=None, edges=self.expected_reaction_edges,
            include_attributes=False, reaction_subgraph=True
        )
        self.gen.as_networkx(
            nodes=None, edges=self.expected_reaction_edges,
            include_attributes=True, reaction_subgraph=False
        )

    @_check_db_availability
    def test_undirected_plotting(self):
        """
        Testing plot functions (networkx objects only)
        """
        nx_graph = self.gen.as_networkx(
            edges=self.gen.get_subgraph(
                self.sub_organisms, set(), self.sub_metabolites)
        )
        plot_undirected_graph(nx_graph)

    @_check_db_availability
    def test_directed_plotting(self):
        """
        Testing plot functions (networkx objects only)
        """
        nx_graph = self.gen.as_networkx(
            edges=self.gen.get_reaction_subgraph(
                self.sub_organisms, self.sub_genes, self.sub_metabolites),
            reaction_subgraph=True
        )
        plot_directed_graph(nx_graph, reaction_graph=True)

    @_check_db_availability
    def test_getters(self):
        # get all IDs of a specific node type
        _ = self.gen.metabolite_ids
        _ = self.gen.organism_ids
        _ = self.gen.gene_ids
        _ = self.gen.reaction_ids
        # get node type by database index
        self.gen.get_node_by_id(0)
        self.gen.get_node_by_id(0, True)
        # get stats
        _ = self.gen.n_nodes
        _ = self.gen.n_edges
        # getting database attributes of a node
        self.gen.get_node_attributes("M03072")
        # getting neighbours of a particular node
        self.gen.get_node_neighbours("M03072")
        self.gen.get_node_neighbours("M03072", "metabolite")
        self.gen.get_node_neighbours("M03072", "metabolite", as_strings=True)
        # getting edges of a particular node
        self.gen.get_node_edges("M03072")
        self.gen.get_node_edges("M03072", "metabolite")
