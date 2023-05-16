# distutils: language = c++


from pymantra.network.enrichment.LSO.lso cimport (
    LocalSearch, extract_reaction_graph)
from pymantra.network.enrichment.spearman cimport HAS_PARALLEL
from libcpp.string cimport string
from libcpp.vector cimport vector
import os
import warnings
import networkx as nx
from typing import Tuple


_OBJECTIVE_FUNCTIONS = {
    "metabolic_reactions": 0,
    "ld_reactions": 1,
    "reaction_microbiome": 2,
    "reaction_transcriptome": 3,
    "precomputed_objectives": 4
}


def reaction_graph_extraction(
    graph: nx.Graph, include_attributes: bool = True
):
    """Extract the reaction-reaction graph from a metabolite-reaction graph

    Parameters
    ----------
    graph : nx.Graph | nx.DiGraph
        Bipartite metabolite-reaction from which reaction-reaction connections
        are extracted
    include_attributes : bool, True
        Whether to include no
    Returns
    -------
    Set[Tuple[str, str]]
        reaction-reaction edges as a set of 2-tuples. Reactions are always
        represented as strings
    """
    edge_props = {
        (src, tgt): data
        for src, tgt, data in graph.edges(data=True)
    }
    node_data = dict(graph.nodes(data=True))
    edge_data = {
        (src, tgt): data
        for src, tgt, data in graph.edges(data=True)
    }
    return extract_reaction_graph(node_data, edge_data, include_attributes)


cdef class LocalSearchOptimization:
    """
    Interface for running local search with (predefined) objective functions
    using pre-computed node and edge values.

    Usually these values are coming from reaction activity approximation and
    the approximation of reaction/metabolite associations to other omics
    entities approximation, as implemented in the 'reaction_estimation' and
    'reaction_associations' modules.

    If you want to use other node/edge metrics, you need to set/overwrite the
    values stored as 'data'. These values are used by the pre-defined objective
    functions **without** any further checks or corrections.

    If you are intending to use this function with a graph generated manually
    (i.e. with functions outside this module) it must contain the following
    node attributes:
      * 'data'
      * 'node_tupe'


    and the following edge attributes:
      * 'data'
      * 'edge_type'

    """
    __slots__ = ['lso', 'groups']

    cdef LocalSearch* lso
    groups: set

    def __cinit__(
        self, network: nx.Graph, temp: float, delta_min: float, l_min: int,
        l_max: int, max_iter: int, objective_function: str, p: float = 1,
        min_reactions: int = 5, *args, **kwargs
    ):
        """Initialize a LocalSearchOptimization object

        The actual documentation is found in __init__
        """
        # sanity check objective function selection
        if objective_function not in _OBJECTIVE_FUNCTIONS.keys():
            raise ValueError(
                f"Objective function '{objective_function}' is unknown. "
                f"Please use one of the following options: "
                f"{', '.join(_OBJECTIVE_FUNCTIONS.keys())}"
            )

        network_data = dict(network.nodes(data=True))
        # sanity check input graph
        if not network_data:
            raise ValueError(
                "'network' does not contain any node data!"
            )
            # TODO: further sanity checks for network
        
        # NOTE: this simply uses the groups from the first node
        # => consistency assumed!
        try:
            data = nx.get_node_attributes(network, "data")
            if data:
                dnode = list(data.keys())[0]
                self.groups = set(data[dnode].keys())
            else:
                vdata = nx.get_node_attributes(network, "vec_group_data")
                vdnode = list(vdata.keys())[0]
                self.groups = set(vdata[vdnode].keys())
        except KeyError:
            raise KeyError(
                "A 'data' attribute must be present in the network attributes"
                ", which is a dictionary with one (numeric) value per group "
                "(key)"
            )
        if len(self.groups) < 2:
            raise ValueError(
                "At least two groups are required to run the enrichment, "
                f"but only {len(self.groups)} were found. Please provide "
                "the groups as node attributes with the name 'value'."
            )
        # TODO: check per_sample data
        # TODO: check data, groups and max_iter for types
        edge_data = {
            (src, tgt): data
            for src, tgt, data in network.edges(data=True)
        }
        # build c++ object
        self.lso = new LocalSearch(
            network_data, edge_data, _OBJECTIVE_FUNCTIONS[objective_function],
            temp, delta_min, l_min, l_max, min_reactions, max_iter, p
        )

    def __init__(
        self, network: nx.Graph, temp: float, delta_min: float, l_min: int,
        l_max: int, max_iter: int, objective_function: str, p: float = 1,
        min_reactions: int = 5, *args, **kwargs
    ):
        """Initialize a LocalSearchOptimization object

        Regarding the local search objective it is possible to adapt or
        add new objective functions in general, but currently requires C++
        functions and re-compilation.

        The current implemented objective functions are

        1. 'metabolic_reactions': TODO

        2. 'ld_reactions': TODO

        3. 'reaction_microbe': currently not implemented

        4. 'reaction_transcriptome': currently not implemented

        Parameters
        ----------
        network: nx.Graph
            Reaction network on which the local search should be computed.

            The graph is assumed to have the following node attributes:
              * 'data'
              * 'node_type'

            and the following edge attributes:
              * 'data'
              * 'edge_type'

            Usually `network` will be computed using
            :py:class:`pymantra.databases.NetworkGenerator`
        temp : float
            Initial simulated annealing temperature, exponentially
            decaying every iteration. The higher `temp` the more likely
            it is to a solution with a lower score at any iteration.
            Intuitively more 'hops' will be performed at higher temperature.
        delta_min : float
            Minimum improvement per iteration
        l_min : int
            Minimal solution size
        l_max : int
            Maximal solution size
        max_iter : int
            Maximum number of iterations before local search is stopped,
            if the (sub)optimal results has not been found
        objective_function : str
            Which objective functions to use. Possible options are currently:
            "ld_reactions", "reaction_microbiome", "reaction_gene" and
            "precomputed_objectives"
        min_reactions : int
            Minimum number of reactions to be contained in the solution.
            This has no effect for metabolomics-only experiments.
        p : Optional[float], default 2.
            Which :math:`L^p`(Minkowski)-norm to use in the objective function.
            Might not be relevant for all objective functions
        """
        # just a dummy to recognize the signature of __cinit__
        pass

    @property
    def objective_values(self):
        # objective values as dictionary
        return self.lso.get_objective_values()

    def score_final_solution(self, groups):
        cdef vector[string] group_vec = [
            str.encode(groups[i]) for i in range(2)]
        return self.lso.score_solution(group_vec)

    def run_local_search(self, groups: Tuple[str, str] = None,
                         n_threads: int = 1):
        # NOTE for docs: random seed is only generated once => if you want
        #                different results you need to set the seed via
        #                'set_seed' between calls to run_local_search
        if n_threads > os.cpu_count():
            raise ValueError(
                f"Number of threads to use ({n_threads}) cannot be larger "
                f"than the number of threads available ({os.cpu_count()})"
            )
        elif n_threads > 1 and not HAS_PARALLEL():
            warnings.warn(
                "Parallel processing is only possible when compiled with "
                "OpenMP support", category=RuntimeWarning
            )
            n_threads = 1
        if groups is None:
            if len(self.groups) == 2:
                groups = tuple(self.groups)
            else:
                raise ValueError(
                    "Groups must be specified as a 2-tuple, "
                    "when the total number of groups "
                    f"in the data is >2 ({len(self.groups)} groups found)."
                )
        elif len(groups) > 2:
            raise ValueError(
                f"'groups' must be a 2-tuple, not a {len(groups)}-tuple! "
                "Please run pairwise enrichment only."
            )
        elif groups[0] not in self.groups:
            raise ValueError(
                f"Group '{groups[0]}' not found in the first node of the "
                "graph. Please make sure all group values are available for "
                "all nodes."
            )
        elif groups[1] not in self.groups:
            raise ValueError(
                f"Group '{groups[1]}' not found in the first node of the "
                "graph. Please make sure all group values are available for "
                "all nodes."
            )
        # prepare group array for passing
        cdef vector[string] group_vec = [
            str.encode(groups[i]) for i in range(2)]
        # actual local search
        self.lso.run_local_search(group_vec, n_threads)

        # extracting results
        subgraph = self.lso.get_best_solution()
        score = self.lso.get_best_score()
        converged = self.lso.get_converged()

        return subgraph, score, converged

    def set_seed(self, seed, seed_size=None):
        if seed_size is None:
            if not isinstance(seed, set) and not isinstance(seed, list):
                raise ValueError(
                    "'seed_size' must be set when 'seed' is empty"
                )
            self.lso.set_seed_py(seed)
        elif seed_size < self.lso.get_lmin():
            warnings.warn(
                "'seed_size' is smaller than l_min. The generated "
                "seed will be of size l_min"
            )
        # TODO: check seed and seed_size for types
        seed = str(seed)
        self.lso.set_seed_py(seed.encode('utf-8'), seed_size)

    def score_progression(self):
        return self.lso.get_score_progression()

    def precompute(self, groups):
        if len(groups) != 2:
            raise ValueError(
                f"Only binary comparisons supported. Please provide exactly "
                f"2 groups, not {len(groups)}."
            )
        cdef vector[string] group_vec = [
            str(group).encode('utf-8') for group in groups]
        self.lso.precompute_objectives_py(group_vec)

    def converged(self) -> bool:
        """Whether the local search converged or terminated due to reaching
        max_iter
        """
        return self.lso.get_converged()

    def l_min(self) -> int:
        """Current choice of minimum solution size"""
        return self.lso.get_lmin()

    def min_reactions(self) -> int:
        """Current choice of minimum number of reactions in the solution"""
        return self.lso.get_min_reactions()

    def l_max(self) -> int:
        """Current choice of maximum solution size"""
        return self.lso.get_lmax()

    def temp(self) -> float:
        """Current choice of simulated annealing temperature"""
        return self.lso.get_temp()

    def max_iter(self) -> int:
        """Current choice of maximum number of iteration for local search"""
        return self.lso.get_maxiter()

    def delta_min(self) -> float:
        """Current choice of minimum progress per iteration"""
        return self.lso.get_deltamin()

    def p(self) -> float:
        """Current choice of lp-norm"""
        return self.lso.get_lp_norm()

    def set_l_min(self, lmin: int):
        """Set the choice of minimum solution size"""
        self.lso.set_lmin(lmin)

    def set_min_reactions(self, min_: int):
        """Set the choice of minimum number of reactions in the solution"""
        self.lso.set_min_reactions(min_)

    def set_l_max(self, lmax: int):
        """Set the choice of maximum solution size"""
        self.lso.set_lmax(lmax)

    def set_temp(self, temp: float):
        """Set the choice of simulated annealing temperature"""
        self.lso.set_temp(temp)

    def set_max_iter(self, max_iter: int):
        """Set the choice of maximum number of iteration for local search"""
        self.lso.set_maxiter(max_iter)

    def set_delta_min(self, delta_min: float):
        """Set the choice of minimum progress per iteration"""
        self.lso.set_deltamin(delta_min)

    def set_p(self, p: float):
        """Set the choice of lp-norm to use. Only relevant for certain
        objective functions
        """
        self.lso.set_lp_norm(p)
