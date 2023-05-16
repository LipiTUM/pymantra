import warnings
from neo4j import Result
from neo4j.graph import Node, Relationship
from typing import (
    Dict, Set, List, Tuple, Union, Optional)
import networkx as nx

from pymantra.database.exceptions import IncorrectNodeType, IncorrectEdgeType
from pymantra.database.base import Neo4jBaseConnector
from pymantra.database.utils import reduce_reaction_nodes
from pymantra.statics import (
    Edge,
    NODE_TYPES, EDGE_TYPES,
    NODE_TYPE_NAMES, EDGE_TYPE_NAMES,
    NODE_TYPES_BY_EDGE, DIRECT_EDGE_TYPE_NAMES, DIRECT_NODE_TYPES_BY_EDGE
)


def _check_edge_types(edge_types):
    for edge_type in edge_types:
        if edge_type not in EDGE_TYPES:
            raise IncorrectNodeType(
                f"Edge type '{edge_type}' is not supported. Please use only "
                f"the following: {', '.join(EDGE_TYPES)}"
            )


def _check_node_types(node_types):
    for node_type in node_types:
        if node_type not in NODE_TYPES:
            raise IncorrectNodeType(
                f"Node type '{node_type}' is not supported. Please use only "
                f"the following: {', '.join(NODE_TYPES)}"
            )


class NetworkGenerator(Neo4jBaseConnector):
    """
    Querying a neo4j database containing the reference network generated with
    the Neo4jGenerator class.

    Most query functions depend on the requirements test with the `Verifier`
    class. To ensure that all functions work as expected, only databases test
    for their correctness should be used.
    """
    def __init__(
        self, uri: str, auth: Union[Tuple[str, str], None] = None, **kwargs
    ):
        """
        Parameters
        ----------
        uri : str
            database uri

        auth : Tuple[str, str], Optional
            database credentials int the form of (user, password).
            If database is not secured pass None.

        Raises
        ------
        ConnectionError
            If database is not reachable with the given parameters
        """
        super(NetworkGenerator, self).__init__(uri, auth, **kwargs)

    # ============================================ #
    # ============ Auxiliary Function ============ #
    # ============================================ #
    # TODO: document auxiliary functions
    @staticmethod
    def _extract_nodes_from_results_(query_result: Result) -> List[Node]:
        return [
            node for node_list in query_result.values()
            for node in node_list
        ]

    @staticmethod
    def _extract_edges_from_results_(query_result: Result) -> List[
            Relationship]:
        """
        get query results from relationship queries and extra the
        `Relationship` objects from it

        Parameters
        ----------
        query_result: Result
            :obj:`neo4j.Result` object conatining relationship query result

        Returns
        -------
        List[Relationship]
            list of relationships contained in `query_result`
        """
        return [
            edge for edge_list in query_result.values()
            for edge in edge_list
        ]

    def _extract_nodes_from_edge_(
        self, edge: Relationship, as_strings: bool = False
    ) -> Union[Tuple[Node, Node], Tuple[str, str]]:
        """
        Extract the nodes from a neo4j relationship object

        Parameters
        ----------
        edge: Relationship
            Edge as a neo4j relationship

        as_strings: bool, default False
            Whether to return nodes as strings (True) :obj:`neo4j.grpah.Node`

        Returns
        -------
        Union[Tuple[Node, Node], Tuple[str, str]]
            1-tuple of source and target node, either strings or as
            :obj:`neo4j.graph.Node`
        """
        node_ids = [node.id for node in edge.nodes]
        query_str = "MATCH (n) WHERE ID(n) = $node RETURN n"
        if as_strings:
            return (
                self.run(query_str, node=node_ids[0]).value()[0][
                    'nodeLabel'],
                self.run(query_str, node=node_ids[1]).value()[0][
                    'nodeLabel']
            )
        return (
            self.run(query_str, node=node_ids[0]).value()[0],
            self.run(query_str, node=node_ids[1]).value()[0]
        )

    @staticmethod
    def _get_node_type(node: Node):
        """
        Extract the node type from a neo4j node

        Parameters
        ----------
        node: Node
            :obj:`neo4j.graph.Node`

        Returns
        -------
        str
            Node type, specified by node.labels
        """
        # NOTE: this only works correctly if node has __exactly one__ label
        #       => ensured in `Verifier`
        return next(iter(node.labels))

    def _get_node_ids_by_type(self, node_type: str) -> Dict[int, str]:
        """
        Get all nodes of specific node type and return their ids

        Parameters
        ----------
        node_type: str
            name of the node type to query

        Returns
        -------
        Dict[int, str]
            dictionary of node IDs (:obj:`int`) as keys and node labels
            as values
        """
        return {
            node.id: node['nodeLabel'] for node_list in
            # NOTE: we can do this safely, since there is a fixed number of
            # node type options
            self.run(
                f"MATCH (n: {NODE_TYPE_NAMES[node_type]}) RETURN n").values()
            for node in node_list
        }

    @staticmethod
    def _edge_query_to_set(results: Result, src_type: str, tgt_type: str) -> \
            Set[Edge]:
        """
        Takes the results of a query returning source and target nodes as an
        iterable, where each element contains exactly two nodes

        Parameters
        ----------
        results: Result
            query results

        src_type: str
            source node type

        tgt_type: str
            target node type

        Returns
        -------
        Set[Edge]
            set of edges contained in the query
        """
        edges = set()
        for edge in results:
            # sanity checking node types
            if src_type in edge[0].labels and tgt_type in edge[1].labels:
                edges.add(Edge(edge[0]['nodeLabel'], edge[1]['nodeLabel']))
        return edges

    def _second_order_neighbours(
        self, src_nodes: Set[str], src_type: str, tgt_type: str,
        intermediate_type: str, edge_types: Tuple[str, str] = None,
        targets: Set[str] = None
    ) -> Set[Edge]:
        """
        Query all second order neighbours of a node belonging to a specific
        node type

        Parameters
        ----------
        src_nodes: str
            Source node labels to query
        src_type: str
            node type of the query node
        tgt_type: str
            node type of the second order neighbours
        intermediate_type: str
            node type of the nodes connecting enrichment nodes and their 2-hop
            neighbours
        edge_types: Tuple[str, str], Optional
            edge types connecting the query node to intermediate nodes (first
            position) and intermediate nodes to the 2-hop neighbours
            (second position)
        targets: Set[str], Optional
            target node labels to choose from. If not given, all possible
            targets will be considered

        Returns
        -------
        Set[Edge]
            Set of Edge objects representing second order edges
        """
        _check_node_types([src_type, tgt_type, intermediate_type])

        params = {"labels": list(src_nodes)}

        src_match = f"(s:{src_type})"
        intermediate = f"(g:{intermediate_type})"
        tgt_match = f"(t:{tgt_type})"
        if edge_types:
            _check_edge_types(edge_types)
            src_edge = f"rs:{edge_types[0]}"
            tgt_edge = f"rt:{edge_types[1]}"
        else:
            src_edge = ""
            tgt_edge = ""
        query = f"MATCH {src_match}-[{src_edge}]-{intermediate}-" \
                f"[{tgt_edge}]-{tgt_match} " \
                "WHERE s.nodeLabel IN $labels"
        if targets:
            params["target_labels"] = list(targets)
            query += " AND t.nodeLabel IN $target_labels" \
                     " AND s.nodeLabel<>t.nodeLabel"
        query += " RETURN t, s"
        return self._edge_query_to_set(
            self.run(query, **params), src_type, tgt_type)

    def _edges_to_higher_order_neighbours(
        self, src_nodes: Set[str], src_type: str, tgt_type: str,
        intermediate_types: Union[str, Tuple[str, str]],
        targets: Set[str] = None,
        edges: Set[Edge] = None
    ) -> Union[None, Set[Edge]]:
        """
        Wrapper for _third_order_neighbours and _second_order_neighbours,
        restricting outputs to a given set of possible higher order neighbours

        Parameters
        ----------
        src_nodes: Set[str]
            All source node labels to query

        src_type: str
            node type of the query node

        tgt_type: str
            node type of the third order neighbours

        intermediate_types: Union[str, Tuple[str, str]]
            If a string is given, second order neighbours of that node type
            will be queried Else, the 2-tuple represents the node types of the
            nodes connecting enrichment nodes and their 3-hop neighbours. First
            position indicates first order neighbour node type, second position
            second order neighbour node type

        targets: Set[str], Optional
            Set of target nodes, if not None only nodes, whose label is in this
            set will be returned

        edges: Set[Edge], Optional
            If given edges will be stored in this set as and
            :obj:`Edge`(node, neighbour) and None will be returned. Otherwise,
            edges are saved and returned in a new set.

        Returns
        -------
        Union[None, Set[Edge]]
            None, if a set of edges to fill is provided, else a set of
            :obj:`Edge` representing higher order connections
        """
        return_ = False
        if edges is None:
            edges = set()
            return_ = True
        if isinstance(intermediate_types, str):
            if src_type == NODE_TYPE_NAMES['metabolite'] and tgt_type == \
                    NODE_TYPE_NAMES['metabolite']:
                neighbours = self._second_order_neighbours(
                    src_nodes, src_type, tgt_type, intermediate_types,
                    (EDGE_TYPE_NAMES['substrate'], EDGE_TYPE_NAMES['product']),
                    targets
                )
            else:
                neighbours = self._second_order_neighbours(
                    src_nodes, src_type, tgt_type, intermediate_types,
                    targets=targets
                )
        elif len(intermediate_types) == 2:
            _check_node_types([src_type, tgt_type] + list(intermediate_types))
            # Currently, the only third oder connections that make sense are
            # organism -> metabolites
            params = {"src_nodes": list(src_nodes)}
            query = f"MATCH (s:{src_type})-[]-(io:{intermediate_types[0]})-[]"\
                    f"-(it:{intermediate_types[1]})-[]-(t:{tgt_type}) " \
                    "WHERE s.nodeLabel IN $src_nodes"
            if targets:
                params["targets"] = list(targets)
                query += " AND t.nodeLabel IN $targets"\
                         " AND s.nodeLabel<>t.nodeLabel"  # remove self-loops
            query += " RETURN t, s"
            neighbours = self._edge_query_to_set(
                self.run(query, **params), src_type, tgt_type)
        else:
            raise ValueError("Only second and third order neighbours possible")
        edges = edges.union(neighbours)
        if return_:
            return edges

    def _get_neighbours_of_type(
        self, src_type: str, tgt_type: str, edge_type: str
    ) -> Set[Tuple[Node, Node]]:
        """
        Query all direct neighbour pairs between two specific node types

        Parameters
        ----------
        src_type: str
            source node type

        tgt_type: str
            target node type

        edge_type: str
            edge type connection source and target nodes

        Returns
        -------
        Set[Tuple[Node, Node]]
            set of 2-tuples representing queried edges
        """
        _check_node_types([src_type, tgt_type])
        _check_edge_types([edge_type])
        edges = set()
        src_ = f"(s:{src_type})"
        tgt_ = f"(t:{tgt_type})"
        query = f"MATCH {src_}-[r:{edge_type}]-{tgt_} RETURN t, s"
        results = self.run(query)
        for edge in results:
            if edge[0] != edge[1]:
                edges.add(Edge(edge[0], edge[1]))
        return edges

    # ===================================================== #
    # ============ Whole Node/Edge set queries ============ #
    # ===================================================== #
    @property
    def n_nodes(self) -> Dict[str, int]:
        """
        Returning the number of nodes per node type.

        Returns
        -------
        Dict[str, int]
            Node types are keys and node counts are values
            If no nodes are found and empty :obj:`dict` will be returned.

        """
        return {
            label: self.run(
                f"MATCH (n:{label}) RETURN count(n) as count").value()[0]
            for label in NODE_TYPES
        }

    @property
    def nodes(self) -> Union[Dict[str, Set[str]], Dict[str, Set[Node]]]:
        """
        Querying all nodes in the database

        Returns
        -------
        Dict[str, Set[str]]
            nodeLabels of all nodes by node types, where nodes types are the
            keys and the values are sets of strings

            If no nodes are found and empty :obj:`dict` will be returned.
        """
        return {
            label: {
                node['nodeLabel']
                for node in self._extract_nodes_from_results_(
                    self.run(f"MATCH (n:{label}) RETURN n"))
            }
            for label in NODE_TYPES
        }

    @property
    def neo4j_nodes(self) -> Dict[str, Set[Node]]:
        """
        Returns
        -------
        Dict[str, Set[Node]]
            nodeLabels of all nodes by node types, where nodes types are the
            keys and the values are sets of  :obj:`neo4j.graph.Node`.
            If no nodes are found and empty :obj:`dict` will be returned.
        """
        return {
            label: {
                node for node in self._extract_nodes_from_results_(
                        self.run(f"MATCH (n:{label} RETURN n")
                )
            }
            for label in NODE_TYPES
        }

    @property
    def n_edges(self) -> Dict[str, int]:
        """
        Returning the number of relations per type

        Returns
        -------
        Dict[str, int]
            Relation types are keys and counts values
            If no edges are found and empty :obj:`dict` will be returned.
        """
        return {
            label: self.run(
                f"MATCH ()-[r:{label}]->() RETURN count(r)").value()[0]
            for label in EDGE_TYPES
        }

    @property
    def edges(self, as_string: bool = True) -> Union[
            Dict[str, Set[str]], Dict[str, Edge]]:
        """
        Querying all nodes in the database

        Returns
        -------
        Union[Dict[str, Set[str]], Dict[str, Edge]]
            All edges contained in the database by type.

            Edges are either represented as strings of
            "source.nodeLabel -> target.nodeLabel"
            or as :obj:`Edge` (:obj:`NamedTuple` with attributes
            `source` at position 0 and `target` at position 1)

            If no edges are found and empty :obj:`dict` will be returned.
        """
        edges = {
            label: {edge for edge in
                    self._extract_edges_from_results_(self.run(
                        f"MATCH ()-[r:{label}]->() RETURN DISTINCT r"))}
            for label in EDGE_TYPES
        }
        if as_string:
            return {
                label: {
                    Edge(
                        *self._extract_nodes_from_edge_(edge, as_strings=True)
                    )
                    for edge in edges_}
                for label, edges_ in edges.items()
            }
        else:
            return {
                label: {self._extract_nodes_from_edge_(edge, as_strings=False)
                        for edge in edges_}
                for label, edges_ in edges.items()
            }

    def get_node_by_id(self, node_id: int, as_string: bool = False) -> Union[
            str, Node]:
        """
        Query a node by its ID

        Parameters
        ----------
        node_id: int
            Node ID to query

        as_string: bool, Optional, default False
            If true the node label is returned, else the
            :obj:`neo4j.graph.Node`

        Returns
        -------
        Union[str, Node]
            Node with the respective ID as :obj:`str` (nodeLabel) if
            `as_string` is True, else the :obj:`neo4j.graph.Node`
        """
        query = "MATCH (n) WHERE ID(n) = $node_id RETURN n"
        result = self.run(query, node_id=node_id).value()[0]
        if as_string:
            return result['nodeLabel']
        return result

    @property
    def metabolite_ids(self):
        """Get all metabolites in the database by their ID and node label

        Returns
        -------
        Dict[int, str]
            Dictionary of all metabolite in the database as ID, name pairs
        """
        return self._get_node_ids_by_type(NODE_TYPE_NAMES['metabolite'])

    @property
    def organism_ids(self):
        """Get all organisms in the database by their ID and node label

        Returns
        -------
        Dict[int, str]
            Dictionary of all organisms in the database as ID, name pairs
        """
        return self._get_node_ids_by_type(NODE_TYPE_NAMES['organism'])

    @property
    def gene_ids(self):
        """Get all gene in the database by their ID and node label

        Returns
        -------
        Dict[int, str]
            Dictionary of all genes in the database as ID, name pairs
        """
        return self._get_node_ids_by_type(NODE_TYPE_NAMES['gene'])

    @property
    def reaction_ids(self):
        """Get all reactions in the database by their ID and node label

        Returns
        -------
        Dict[int, str]
            Dictionary of all reactions in the database as ID, name pairs
        """
        return self._get_node_ids_by_type(NODE_TYPE_NAMES['reaction'])

    # ================================================================ #
    # ============ Query Functions for network Generation ============ #
    # ================================================================ #
    def get_all_edges(self, edge_type: str, limit: int = None) -> Set[Edge]:
        """
        Query all relationships of a specific type

        Parameters
        ----------
        edge_type : str
            Must be one of the elements in `utils.EDGE_TYPES`
        limit : int, Optional
            If specified it represents the maximum number of edges to return,
            otherwise all edges are returned (default)

        Returns
        -------
        Set[Edge]
            All edges of the respective edge type represented as namedtuple of
            size with attributes `source` and `target`, which are both
            :obj:`str` of the respective nodeLabels.  If no edges are found and
            empty :obj:`set` will be returned.
        """
        params = {}
        if edge_type not in EDGE_TYPES:
            raise IncorrectEdgeType(
                f"Unknown edge type: {edge_type}. "
                f"Valid edge types are '{EDGE_TYPES}'")
        query = f"MATCH (s)-[r:{edge_type}]-(t) RETURN DISTINCT r"
        if limit is not None:
            params["limit"] = limit
            query += "LIMIT $limit;"
        else:
            query += ";"
        results = self._extract_edges_from_results_(
            self.run(query, **params))
        return {
            self._extract_nodes_from_edge_(edge, as_strings=False)
            for edge in results
        }

    def get_node_attributes(
        self, node: str, node_type: str = None
    ) -> Dict[str, any]:
        """Get the attributes of a specific node

        Parameters
        ----------
        node: str
            Name of the node (i.e. internal iD/species name)
        node_type: str, optional
            Node type. Specifying this will speed up the computation, since the
            number of nodes filtered by neo4j are reduced

        Returns
        -------
        Dict[str, any]
            Node attribute dictionary
        """
        # FIXME: change these names in the database itself
        node = node.replace("'", "")
        if node_type is None:
            query = "MATCH (n {nodeLabel:$node}) RETURN n"
        else:
            _check_node_types([node_type])
            query = f"MATCH (n:{node_type} {{nodeLabel:$node}}) RETURN n"
        try:
            return dict(self.run(query, node=node).value()[0])
        except IndexError:
            return {}

    def get_node_neighbours(
        self, node: str, node_type: str = None, as_strings: bool = False
    ) -> Union[Dict[str, Set[str]], Dict[str, Set[Node]]]:
        """
        Query all neighbours of a specific node by node label, irrespective of
        their node type

        Parameters
        ----------
        node : str
            Node label of the node to query
        node_type: str, Optional
            Query node type. Query results should be the same, since
            nodeLabels are supposed to be unique across all node types,
            however, speed might be different
        as_strings: bool, Optional, default False
            If True nodes will be returned as their nodeLabels, else as
            :obj:neo4j.graph.Node` objects.

        Returns
        -------
        Union[Dict[str, Set[str]], Dict[str, Set[Node]]]
            All direct neighbours by node type (`dict.keys`). If as_strings is
            True nodes are :obj:`set` of :obj:`str`, else :obj:`set` of
            :obj:`neo4j.graph.Node`.

            If no neighbours are found and empty :obj:`dict` will be returned.
        """
        # NOTE: the following queries depend on the uniqueness of node labels
        if node_type:
            _check_node_types([node_type])
            query = f"MATCH (s:{node_type} {{nodeLabel:$node}})-[]-(t) " \
                    "RETURN DISTINCT t"
        else:
            query = \
                "MATCH (s {nodeLabel:$node})-[]-(t) RETURN DISTINCT t"
        results = self.run(query, node=node).values()
        neighbours = {}
        if as_strings:
            for node_ in results:
                neighbours.setdefault(
                    self._get_node_type(node_[0]), set()).add(
                        node_[0]['nodeLabel'])
        else:
            for node_ in results:
                neighbours.setdefault(
                    self._get_node_type(node_[0]), set()).add(node_[0])
        return neighbours

    def get_node_edges(self, node, node_type: str = None, **kwargs) -> Dict[
            str, Set[Edge]]:
        """
        Query all edges of a given node (by node label), irrespective of edge
        types

        Parameters
        ----------
        node: str
            Node label of the node to query

        node_type: str, Optional
            Query results should be the same, since nodeLabels
            are supposed to be unique across all node types, however, speed
            might be different

        kwargs
            Optional keyword arguments

        Returns
        -------
        Dict[str, Set[Edge]]
            All edges going out of or to the given input node by edge type
        """
        if node_type:
            if node_type not in NODE_TYPES:
                raise IncorrectNodeType(
                    f"Invalid node type {node_type}. "
                    f"Please use only {NODE_TYPES}"
                )
            query = f"MATCH (n:{node_type} {{nodeLabel:$node}})-[r]-() " \
                    f"RETURN DISTINCT r"
        else:
            query = "MATCH (n {nodeLabel:$node})-[r]-() " \
                    "RETURN DISTINCT r"
        results = self.run(query, node=node, **kwargs).values()
        edges = {}
        for edge in results:
            edges.setdefault(edge[0].type, set()).add(
                self._extract_nodes_from_edge_(edge[0], False))
        return edges

    def get_organism_metabolite_connections(
        self, organisms: Set[str], metabolites: Set[str]
    ) -> Set[Edge]:
        """
        Query all pairwise connections between organisms and metabolites
        given as input

        Parameters
        ----------
        organisms: Set[str]
            Set of all organisms to query

        metabolites: Set[str]
            Set of all metabolites to query

        Returns
        -------
        Set[Edge]
            Set of all organism - metabolite connections. Edges are assumed
            to have **no direction**
        """
        if not isinstance(organisms, set):
            organisms = set(organisms)
        if not isinstance(metabolites, set):
            metabolites = set(metabolites)
        return self._edges_to_higher_order_neighbours(
            organisms, NODE_TYPE_NAMES['organism'],
            NODE_TYPE_NAMES['metabolite'],
            (NODE_TYPE_NAMES['gene'], NODE_TYPE_NAMES['reaction']), metabolites
        )

    def get_gene_metabolite_connections(self, genes: Set[str],
                                        metabolites: Set[str]):
        """
        Query all pairwise connections between genes and metabolites
        given as input

        Parameters
        ----------
        genes: Set[str]
            Set of all genes to query

        metabolites: Set[str]
            Set of all metabolites to query

        Returns
        -------
        Set[Edge]
            Set of all gene - metabolite connections. Edges are assumed
            to have **no direction**
        """
        if not isinstance(genes, set):
            genes = set(genes)
        if not isinstance(metabolites, set):
            metabolites = set(metabolites)
        return self._edges_to_higher_order_neighbours(
            genes, NODE_TYPE_NAMES['gene'], NODE_TYPE_NAMES['metabolite'],
            NODE_TYPE_NAMES['reaction'], metabolites
        )

    def get_metabolite_metabolite_connection(self, metabolites: Set[str]):
        """
        Query all pairwise connections between metabolites given as input

        Parameters
        ----------
        metabolites: Set[str]
            Set of all metabolites to query

        Returns
        -------
        Set[Edge]
            Set of all metabolite - metabolite connections. Edges are assumed
            to have **no direction**
        """
        return self._edges_to_higher_order_neighbours(
            metabolites, NODE_TYPE_NAMES['metabolite'],
            NODE_TYPE_NAMES['metabolite'],
            NODE_TYPE_NAMES['reaction'], metabolites
        )

    def get_reaction_subgraph(
        self, organisms: Set[str], genes: Set[str], metabolites: Set[str],
        reaction_organism: Optional[Tuple[str, str]] = None
    ) -> Dict[str, Set[Edge]]:
        """Extract the edges for a given set of entities

        Query a subgraph with all genes, organisms and metabolites given
        and retain the original graph structure with reaction nodes.

        **Important**: gene - organism, organism - reaction and gene - reaction
        edges are of opposite direction outside the database. The database
        structure is made to allow efficient queries, which do not reflect the
        'passing' directions required for quantitative metabolic-network style
        analyses.

        Parameters
        ----------
        organisms: Set[str]
            A set of all organisms to be included in the subgraph.
            The names must correspond to the nodeLabel property in
            the database.
        genes: Set[str]
            A set of all genes to be included in the subgraph.
            The names must correspond to the nodeLabel property in
            the database.
        metabolites: Set[str]
            A set of all metabolites to be included in the subgraph.
            The names must correspond to the nodeLabel property in
            the database.
        reaction_organism: Tuple[str, str], optional
            Specify an organism for which the metabolic reactions should be
            extracted as a 2-tuple of [ID type, ID], where ID type must be
            'Abbreviation_KEGG' or 'KeggID' and ID the KEGG organism code or
            T number, respectively. For human this would thus either be
            ['Abbreviation_KEGG', 'hsa'] or ['KeggID', 'T01001'].
            If `organisms` is not empty the specified organism will be added
            on top.

        Returns
        -------
        Dict[str, Set[Edge]]
            A dictionary, where keys represent edge types as specified in
            `utils.EDGE_TYPES` pointing to a set of :obj:`Edge` representing
            all edges of the respective type contained in the subgraph.

        Examples
        --------
        >>> generator = NetworkGenerator(
        ...     "bolt://localhost:7687", auth=('<user>', '<password>'))
        >>> metabos = {'FDMO3', 'h2o', 'FDMO2', 'fald', 'FDMO6', 'so3',
        >>>            'FMNRx', 'nad', 'FMNRx2', 'fmn', 'nadp'}
        >>> gs = {'1576', '1557', '1559'}
        >>> orgs = {"Streptomyces tsukubensis", "Bacillus smithii",
        >>>         "Streptomyces fulvissimus"}
        >>> edges_ = generator.get_reaction_subgraph(orgs, gs, metabos)
        """
        if reaction_organism:
            if reaction_organism[0] not in ("Abbreviation_KEGG", "KeggID"):
                raise ValueError(
                    "Unknown organism identification type "
                    f"{reaction_organism[0]}. Please use either "
                    "'Abbreviation_KEGG' or 'KeggID'"
                )
            # getting all reactions from the "main"/host organism
            reaction_query = \
                f"MATCH (r:{NODE_TYPE_NAMES['reaction']})" \
                f"-[:{EDGE_TYPE_NAMES['organism_reaction']}]-" \
                f"(o:{NODE_TYPE_NAMES['organism']} " \
                f"{{{reaction_organism[0]}:$ro_id}}) RETURN r"
            org_reaction_query = self.run(
                reaction_query, ro_id=reaction_organism[1]).values()
            org_reactions = {
                str(reaction[0].id) for reaction in org_reaction_query}
            if organisms:
                # getting all reactions from the "parasitic" organism
                reaction_query = \
                    f"MATCH (r:{NODE_TYPE_NAMES['reaction']})" \
                    f"-[:{EDGE_TYPE_NAMES['organism_reaction']}]-" \
                    f"(o:{NODE_TYPE_NAMES['organism']})" \
                    f"WHERE o.nodeLabel IN $organisms " \
                    "RETURN r"
                org_reaction_query = self.run(
                    reaction_query, organisms=list(organisms)).values()
                org_reactions.union({
                    str(reaction.id) for _, reaction, __ in org_reaction_query
                })

            reaction_restriction = "AND ID(r) in $target_reactions "
        else:
            reaction_restriction = ""

        # metabolite connections
        # substrates
        metabo_query = f"MATCH (s:{NODE_TYPE_NAMES['metabolite']})" \
                       f"-[:{EDGE_TYPE_NAMES['substrate']}]-" \
                       f"(r:{NODE_TYPE_NAMES['reaction']})" \
                       f"-[:{EDGE_TYPE_NAMES['product']}]-" \
                       f"(p:{NODE_TYPE_NAMES['metabolite']}) " \
                       "WHERE s.nodeLabel IN $metabolite_list " \
                       "AND p.nodeLabel IN $metabolite_list " \
                       "AND s.nodeLabel <> p.nodeLabel " \
                       f"{reaction_restriction}" \
                       f"RETURN s, r, p"
        if reaction_organism is None:
            metabolite_relations = self.run(
                metabo_query, metabolite_list=list(metabolites)).values()
        else:
            metabolite_relations = self.run(
                metabo_query, metabolite_list=list(metabolites),
                target_reactions=list(org_reactions)
            ).values()
        reaction_nodes = ','.join(
            [str(reaction.id) for _, reaction, __ in metabolite_relations])

        edges = {
            EDGE_TYPE_NAMES['substrate']: set(),
            EDGE_TYPE_NAMES['product']: set()
        }
        # NOTE: queried nodes come back in inverse order
        for product, reaction, substrate in metabolite_relations:
            edges[EDGE_TYPE_NAMES['substrate']].add(
                Edge(substrate['nodeLabel'], reaction['nodeLabel']))
            edges[EDGE_TYPE_NAMES['product']].add(
                Edge(reaction['nodeLabel'], product['nodeLabel']))
        # gene connections
        if genes and organisms:
            go_query = f"MATCH (r:{NODE_TYPE_NAMES['reaction']})" \
                       f"-[:{EDGE_TYPE_NAMES['gene_reaction']}]-" \
                       f"(g:{NODE_TYPE_NAMES['gene']})" \
                       f"-[:{EDGE_TYPE_NAMES['organism_gene']}]-" \
                       f"(o:{NODE_TYPE_NAMES['organism']}) " \
                       "WHERE ID(r) IN $reaction_nodes " \
                       "AND g.nodeLabel IN $genes " \
                       "AND o.nodeLabel IN $organisms " \
                       f"RETURN r, g, o"
            go_relations = self.run(
                go_query, reaction_nodes=list(reaction_nodes),
                genes=list(genes), organisms=list(organisms)
            ).values()

            edges[EDGE_TYPE_NAMES['gene_reaction']] = set()
            edges[EDGE_TYPE_NAMES['organism_gene']] = set()
            for organism, gene, reaction in go_relations:
                edges[EDGE_TYPE_NAMES['gene_reaction']].add(
                    Edge(reaction['nodeLabel'], gene['nodeLabel']))
                edges[EDGE_TYPE_NAMES['organism_gene']].add(
                    Edge(gene['nodeLabel'], organism['nodeLabel']))
        elif genes:
            g_query = f"MATCH (r:{NODE_TYPE_NAMES['reaction']})" \
                      f"-[:{EDGE_TYPE_NAMES['gene_reaction']}]-" \
                      f"(g:{NODE_TYPE_NAMES['gene']})" \
                      f"WHERE ID(r) IN $reaction_nodes " \
                      f"AND g.nodeLabel IN $genes " \
                      f"RETURN r, g"
            g_relations = self.run(
                g_query, reaction_nodes=list(reaction_nodes), genes=list(genes)
            ).values()

            edges[EDGE_TYPE_NAMES['gene_reaction']] = set()
            for gene, reaction in g_relations:
                edges[EDGE_TYPE_NAMES['gene_reaction']].add(
                    Edge(reaction['nodeLabel'], gene['nodeLabel']))
        # organism connections if genes are NOT included
        elif organisms:
            o_query = f"MATCH (r:{NODE_TYPE_NAMES['reaction']})" \
                      f"-[:{EDGE_TYPE_NAMES['organism_reaction']}]-" \
                      f"(o:{NODE_TYPE_NAMES['organism']}) " \
                      f"WHERE ID(r) IN $reaction_nodes " \
                      f"AND o.nodeLabel IN $organisms " \
                      f"RETURN r, o"
            o_relations = self.run(
                o_query, reaction_nodes=list(reaction_nodes),
                organisms=list(organisms)
            ).values()
            edges[EDGE_TYPE_NAMES['organism_reaction']] = set()
            for organism, reaction in o_relations:
                edges[EDGE_TYPE_NAMES['organism_reaction']].add(
                    Edge(reaction['nodeLabel'], organism['nodeLabel']))
        return edges

    def get_subgraph(
        self, organisms: Set[str], genes: Set[str], metabolites: Set[str]
    ) -> Dict[str, Set[Edge]]:
        """Returns a subgraph with all nodes given plus the reaction nodes
        required to connect them

        Parameters
        ----------
        organisms: Set[str]
            Set of all organisms to query
        genes: Set[str]
            Set of all genes to query
        metabolites: Set[str]
            Set of all metabolites to query

        Returns
        -------
        Dict[str, Set[Edge]]
            All connections between organisms, genes and metabolites contained
            in the database. organism - metabolite are third order connections
            (via gene and reaction nodes), all other connections are second
            order (via reaction nodes)

        Examples
        --------
        >>> generator = NetworkGenerator(
        ...     "bolt://localhost:7687", auth=('<user>', '<password>'))
        >>> metabos = {'FDMO3', 'h2o', 'FDMO2', 'fald', 'FDMO6', 'so3',
        >>>            'FMNRx', 'nad', 'FMNRx2', 'fmn', 'nadp'}
        >>> gs = {'1576', '1557', '1559'}
        >>> orgs = {"Streptomyces tsukubensis", "Bacillus smithii",
        >>>         "Streptomyces fulvissimus"}
        >>> edges_ = generator.get_subgraph(orgs, gs, metabos)
        """
        if not isinstance(organisms, set):
            organisms = set(organisms)
        if not isinstance(genes, set):
            genes = set(genes)
        if not isinstance(organisms, set):
            metabolites = set(metabolites)
        connections = {}
        if organisms and genes:
            connections[DIRECT_EDGE_TYPE_NAMES['organism_gene']] = \
                self.get_organism_metabolite_connections(
                    organisms, genes)
        elif organisms and metabolites:
            connections[DIRECT_EDGE_TYPE_NAMES['organism_metabolite']] = \
                self.get_organism_metabolite_connections(
                    organisms, metabolites)
        if genes and metabolites:
            connections[DIRECT_EDGE_TYPE_NAMES['gene_metabolite']] = \
                self.get_gene_metabolite_connections(
                    genes, metabolites)
        if metabolites:
            connections[DIRECT_EDGE_TYPE_NAMES['metabolite_metabolite']] = \
                self.get_metabolite_metabolite_connection(
                    metabolites
                )
        return connections

    def as_networkx(
        self, nodes: Dict[str, Set[str]] = None,
        edges: Dict[str, Set[Edge]] = None, include_attributes: bool = True,
        reaction_subgraph: bool = False, reduce: bool = True
    ) -> nx.DiGraph:
        """Convert a set of nodes or edges to a networkx Graph

        Parameters
        ----------
        nodes: Dict[str, Set[str]], Optional
            Nodes to include by node type. Generally optional, but either
            `nodes` or `edges` need to be given.
            Please note: if `edges` is not specified, and `reaction_subgraph`
            is True reaction nodes given in `nodes` will NOT be considered.
        edges: Dict[str, Set[Edge]], Optional
            Edges to include by edge type. Generally optional, but either
            `nodes` or `edges` need to be given. If not specified, edges will
            be queried from the database using the specified nodes using either
            `get_subgraph` or `get_reaction_subgraph` depending on
            `include_attributes`.
            The only edge attribute currently included is `edge_type`.
        include_attributes: bool, default True
            If True, the nx.Graph.nodes contain the attributes specified in the
            database. Else `node_type` will be the only node attribute in the
            output graph. Please be aware that if True and `edges` are None,
            this might make the function much less efficient.
        reaction_subgraph: bool, default False
            Only relevant if edges is None. If True subgraph edges queried
            result in a reaction subgraph (see
            :py:meth:`~NetworkGenerator.get_reaction_subgraph`) else
            the subgraph will not contain reaction nodes
            (:py:meth:`~NetworkGenerator.get_subgraph`)
        reduce: bool, False
            Whether to reduce the reaction nodes at the end of the

        Returns
        -------
        nx.DiGraph
            Subgraph as a :obj:`nx.DiGraph`

        # TODO: add sample data
        Examples
        --------
        >>> edges_ = {
        >>>     EDGE_TYPE_NAMES['substrate']: {
        >>>         # TODO: example edges
        >>>     },
        >>>     EDGE_TYPE_NAMES['product']: {
        >>>         # TODO: example edges
        >>>     }
        >>> }
        >>> generator = NetworkGenerator(
        ...     "bolt://localhost:7687", auth=('<user>', '<password>'))
        >>> generator.as_networkx(edges=edges_)
        """
        replace_dict = {'metabolites': 'metabolite', 'genes': 'gene'}
        g = nx.DiGraph()
        if nodes:
            # adding nodes and obtaining node attributes from database
            # if include_attributes
            if include_attributes:
                for node_type, nodes_ in nodes.items():
                    for node in nodes_:
                        attrs = self.get_node_attributes(node, node_type)
                        attrs["node_type"] = replace_dict.get(
                            node_type, node_type)
                        g.add_node(node, **attrs)
            else:
                for node_type, nodes_ in nodes.items():
                    g.add_nodes_from(nodes_, node_type=node_type)
            if not edges:
                # edges are calculated by get_*subgraph if not provided
                if reaction_subgraph:
                    edges = self.get_reaction_subgraph(
                        nodes.get(NODE_TYPE_NAMES['organism'], set()),
                        nodes.get(NODE_TYPE_NAMES['gene'], set()),
                        nodes.get(NODE_TYPE_NAMES['metabolite'], set()),
                    )
                else:
                    edges = self.get_subgraph(
                        nodes.get(NODE_TYPE_NAMES['organism'], set()),
                        nodes.get(NODE_TYPE_NAMES['gene'], set()),
                        nodes.get(NODE_TYPE_NAMES['metabolite'], set()),
                    )
            for edge_type, edges_ in edges.items():
                g.add_edges_from(edges_, edge_type=edge_type)
        elif edges:
            if include_attributes:
                for edge_type, edges_ in edges.items():
                    for edge in edges_:
                        if reaction_subgraph:
                            if edge_type in {EDGE_TYPE_NAMES["organism_gene"],
                                             EDGE_TYPE_NAMES["gene_reaction"]}:
                                # the reversion of node types accounts for the
                                # inversion of edge directions of
                                # organism - gene - reaction relations outside
                                # the database (see self.get_reaction_subgraph)
                                iter_tups = zip(
                                    edge,
                                    reversed(NODE_TYPES_BY_EDGE[edge_type])
                                )
                            else:
                                iter_tups = zip(
                                    edge, NODE_TYPES_BY_EDGE[edge_type])
                        else:
                            iter_tups = zip(
                                edge, NODE_TYPES_BY_EDGE[edge_type])
                        for node, node_type in iter_tups:
                            if node not in g.nodes:
                                attrs = self.get_node_attributes(
                                    node, node_type)
                                g.add_node(
                                    node, node_type=node_type, **attrs)
                        g.add_edge(*edge, edge_type=edge_type)
            else:
                for edge_type, edges_ in edges.items():
                    g.add_edges_from(edges_, edge_type=edge_type)
                    if edge_type:
                        if reaction_subgraph:
                            if edge_type in {
                                EDGE_TYPE_NAMES["organism_gene"],
                                EDGE_TYPE_NAMES["gene_reaction"],
                                EDGE_TYPE_NAMES["organism_reaction"]
                            }:
                                # the reversion of node types accounts for the
                                # inversion of edge directions of
                                # organism - gene - reaction relations outside
                                # the database (see self.get_reaction_subgraph)
                                types = reversed(NODE_TYPES_BY_EDGE[edge_type])
                            else:
                                types = NODE_TYPES_BY_EDGE[edge_type]
                        else:
                            types = DIRECT_NODE_TYPES_BY_EDGE[edge_type]
                        node_types = {
                            node: node_type for nodes_ in edges_
                            for node, node_type in zip(nodes_, types)
                        }
                        nx.set_node_attributes(
                            g, values=node_types, name='node_type')
                    else:
                        for node in nodes:
                            node_type = self.get_node_attributes(node).get(
                                'node_type')
                            if not node_type:
                                warnings.warn(
                                    f"Node type for {node} could neither be "
                                    f"inferred from edge type nor found in "
                                    "the database."
                                )
                                g.add_node(node)
                            else:
                                g.add_node(node, node_type=node_type)
        else:
            raise ValueError("Either 'nodes' or 'edges' must be given.")

        if reduce:
            return reduce_reaction_nodes(g)
        return g
