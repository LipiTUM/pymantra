from abc import ABC
from typing import Dict, Set, Tuple
import warnings
import os
import json
import requests
import networkx as nx
from html.parser import HTMLParser

from pymantra.statics import Edge


def _json_safe_dict(py_dict: Dict[str, any]):
    return {
        k: list(v) if isinstance(v, set) else v for k, v in py_dict.items()}


class _HTMLBodyTextReader(HTMLParser, ABC):
    """Helper class to read html to text

    adapted from https://stackoverflow.com/a/55825140 and
    https://stackoverflow.com/questions/16773583
    """
    def __init__(self):
        super().__init__()
        self.text = ""
        self.inbody = False

    def handle_starttag(self, tag, attrs):
        self.inbody = False
        if tag == "body":
            self.inbody = True

    def handle_endtag(self, tag):
        if tag == "body":
            self.inbody = False

    def handle_data(self, data):
        # get data but only after <body> and before </body>
        if self.inbody and data.strip():
            self.text += data


def _get_max_lines(to_print: str, n: int):
    lines = to_print.split(os.linesep)
    return os.linesep.join(lines[:min(len(lines) - 1, n)])


def _parse_html_response(response: requests.Response, n: int = 10):
    hp = _HTMLBodyTextReader()
    hp.feed(response.text)
    return _get_max_lines(hp.text, n)


class APINetworkGenerator:
    """API mirror for :class:`~NetworkGenerator`

    Querying the mantra online neo4j database containing the reference network
    generated with the Neo4jGenerator class.

    Most query functions depend on the requirements test with the `Verifier`
    class. To ensure that all functions work as expected, only databases test
    for their correctness should be used.

    Attributes
    ----------
    url: str
        Base URL to where requests go
    """
    __slots__ = ["url", "_using_local_api"]

    url: str

    def __init__(self, base_url: str = "https://exbio.wzw.tum.de/pymantradb"):
        """Initialize a new APINetworkGenerator instance

        Initialize a new instance to run queries to the neo4j mantra-db API.

        Parameters
        ----------
        base_url: str, https://exbio.wzw.tum.de/pymantradb
            Set the root URL where the server is located
        """
        # local API requires different setting to allow for connections
        self._using_local_api: bool = False
        if "127.0.0.1" in base_url:
            self._local_api("127.0.0.1")
        elif "localhost" in base_url:
            self._local_api("localhost")

        # only http/https URIs are allowed
        if not base_url.startswith("http"):
            if not self._using_local_api:
                raise ValueError(
                    "'base_url' needs to be a full web-address it no local "
                    "API is used. Did not find the required starting "
                    "'https://' or 'http://'."
                )
            else:
                # NOTE: local django API cannot handle https
                base_url = f"http://{base_url}"

        # checking whether given base is working
        self.url = base_url
        if not self.verify_connection():
            raise ConnectionError(
                f"Base URL '{base_url}' seems to be invalid. Connection could "
                "not be verified!"
            )

    def _local_api(self, local_option: str):
        """Make local APIs available"""
        warnings.warn(
            f"Setting 'NO_PROXY' to '{local_option}' to ensure that requests "
            "is able to reach API server"
        )
        os.environ['NO_PROXY'] = local_option
        self._using_local_api = True

    def verify_connection(self) -> bool:
        """Check whether the given base URL is correct

        Returns
        -------
        bool
            True if status code is 200 (connection verified)
        """
        try:
            test_req = requests.get(f"{self.url}/verify-connection")
        except requests.exceptions.RequestException:
            raise ConnectionError(
                f"Connection to {self.url} could not be established! Please "
                "select a valid URL. Aborting verification process..."
            )
        return test_req.status_code == 200

    def _query(self, data: Dict[str, any], subpath: str):
        """Run a query and return the data from json"""
        # dicts are not send properly if passed as python objects
        _prepped_data = {
            k: json.dumps(_json_safe_dict(v)) if isinstance(v, dict) else v
            for k, v in data.items()
        }
        response = requests.post(
            f"{self.url}/{subpath}", data=_prepped_data)
        if response.status_code != 200:
            if response.status_code == 404:
                raise ConnectionError(
                    f"URL '{self.url}/{subpath}' was not found!")
            if response.status_code == 414:
                if self._using_local_api:
                    raise ConnectionError(
                        "Request-URI Too Large. Please change the "
                        "configuration in the APIs nginx container."
                    )
                else:
                    raise ConnectionError(
                        "Request-URI Too Large. Please contact the "
                        "developers to increase the allowed size or find an "
                        "alternative for processing your data."
                    )
            if response.status_code == 422:
                raise ValueError(
                    f"Invalid data reached the server: {response.text}")
            if response.status_code == 500:
                raise ConnectionError(
                    "An internal server error occurred wile processing your "
                    f"request: {_parse_html_response(response)}"
                )
            raise ValueError(
                "Something went wrong while processing your request, failed "
                f"with exit code {response.status_code} ("
                f"{_parse_html_response(response)})."
            )
        return json.loads(response.text)

    def get_reaction_subgraph(
        self, organisms: Set[str], genes: Set[str], metabolites: Set[str],
        reaction_organism: Tuple[str, str] = None
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
        >>> generator = APINetworkGenerator()
        >>> metabos = {'FDMO3', 'h2o', 'FDMO2', 'fald', 'FDMO6', 'so3',
        ...            'FMNRx', 'nad', 'FMNRx2', 'fmn', 'nadp'}
        >>> gs = {'1576', '1557', '1559'}
        >>> orgs = {"Streptomyces tsukubensis", "Bacillus smithii",
        ...         "Streptomyces fulvissimus"}
        >>> edges_ = generator.get_reaction_subgraph(orgs, gs, metabos)
        """
        data = {
            "organisms": organisms, "genes": genes, "metabolites": metabolites,
            "reaction_organism": reaction_organism
        }
        return {
            edge_type: {Edge(*edge) for edge in edges}
            for edge_type, edges in
            self._query(data, "database/get-reaction-subgraph").items()
        }

    def get_subgraph(
        self, organisms: Set[str], genes: Set[str], metabolites: Set[str]
    ) -> Dict[str, Set[Edge]]:
        """Returns a subgraph with all nodes given plus the reaction nodes
        required to connect them

        Parameters
        ----------
        organisms: Set[str]
            Set of all organisms to query
            The names must correspond to the nodeLabel property in
            the database.
        genes: Set[str]
            Set of all genes to query
            The names must correspond to the nodeLabel property in
            the database.
        metabolites: Set[str]
            Set of all metabolites to query
            The names must correspond to the nodeLabel property in
            the database.

        Returns
        -------
        Dict[str, Set[Edge]]
            All connections between organisms, genes and metabolites contained
            in the database. organism - metabolite are third order connections
            (via gene and reaction nodes), all other connections are second
            order (via reaction nodes)

        Examples
        --------
        >>> generator = APINetworkGenerator("http://127.0.0.1:8084")
        >>> metabos = {'FDMO3', 'h2o', 'FDMO2', 'fald', 'FDMO6', 'so3',
        ...            'FMNRx', 'nad', 'FMNRx2', 'fmn', 'nadp'}
        >>> gs = {'1576', '1557', '1559'}
        >>> orgs = {"Streptomyces tsukubensis", "Bacillus smithii",
        ...         "Streptomyces fulvissimus"}
        >>> edges_ = generator.get_subgraph(orgs, gs, metabos)
        """
        data = {
            "organisms": organisms, "genes": genes, "metabolites": metabolites}
        return {
            edge_type: {Edge(*edge) for edge in edges}
            for edge_type, edges in
            self._query(data, "database/get-subgraph").items()
        }

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
            (:py:meth:`~APINetworkGenerator.get_subgraph`)
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
        >>> generator = APINetworkGenerator()
        >>> generator.as_networkx(edges=edges_)
        """
        if nodes is None and edges is None:
            raise ValueError(
                "Either 'nodes' or 'edges' must be given to run 'as_networkx")
        data = {
            "nodes": nodes, "edges": edges,
            "include_attributes": include_attributes,
            "reaction_subgraph": reaction_subgraph,
            "reduce": reduce
        }
        graph_data = self._query(data, "database/as-networkx")

        graph = nx.DiGraph()
        for node, node_data in graph_data["nodes"].items():
            graph.add_node(node, **node_data)
        for edge, edge_data in graph_data["edges"].items():
            src, tgt = edge.split(graph_data["split_str"])
            graph.add_edge(src, tgt, **edge_data)
        return graph
