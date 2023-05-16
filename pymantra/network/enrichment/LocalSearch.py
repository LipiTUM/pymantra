import json
import pathlib
from warnings import warn
from typing import Tuple, List, Set, Union, Dict
from math import ceil
from abc import ABC, abstractmethod
from dataclasses import dataclass, fields
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

from pymantra.plotting import plot_undirected_graph, plot_directed_graph
from pymantra.statics import NODE_TYPE_NAMES, EDGE_BY_NODE_TYPE
from pymantra.network.exceptions import GraphComponentError

from .LSO.lso import LocalSearchOptimization, reaction_graph_extraction


_REACTION_EDGE = "REACTION_REACTION"
_FILE_TYPE = Union[str, pathlib.Path]


@dataclass(frozen=True)
class EnrichmentBaseResults(ABC):
    """Base class handling IO operations for enrichment result classes

    This base class generically enables to read from and write to json files.
    It is not meant to be used directly and is thus an abstract base class.
    """
    def as_dict(self, json_compatible: bool = False) -> dict:
        """Get the class attributes as a dictionary

        Parameters
        ----------
        json_compatible: bool, False
            Whether to make all objects json compatible. In the current
            implementation this means all sets are cast to lists
        """
        if json_compatible:
            dict_repr = {}
            for field in fields(self):
                # sets are cast to lists
                if field.type.__name__ == "set":
                    dict_repr[field.name] = list(getattr(self, field.name))
                else:
                    dict_repr[field.name] = getattr(self, field.name)
            return dict_repr
        return {
            field.name: getattr(self, field.name) for field in fields(self)}

    def to_json(self, file: _FILE_TYPE):
        """Save an enrichment result to a json file

        Parameters
        ----------
        file: Union[str, pathlib.Path]
            Path to the file in which the results will be saved
        """
        json.dump(self.as_dict(json_compatible=True), open(file, "w"))

    @staticmethod
    def _cast_field_dict(cls, field_data: Dict[str, any]):
        """Cast field data from a json file to the correct field types"""
        for field in fields(cls):
            if field.type.__name__ == "set":
                try:
                    field_data[field.name] = set(field_data[field.name])
                except KeyError:
                    raise ValueError(
                        f"Field data for '{field.name}' not found"
                    )
        return field_data

    @classmethod
    @abstractmethod
    def from_json(cls, file: _FILE_TYPE):
        """Read a previously computed enrichment result from a json file"""
        pass


@dataclass(frozen=True)
class EnrichmentResults(EnrichmentBaseResults):
    """Object holding local optimization results

    Attributes
    ----------
    subgraph: Set[int]
        Subnetwork with the best found objective function value
    score: float
        Objective function value of the solution
    converged: bool
        **currently unused** Boolean indicating whether algorithm converged
    """
    subgraph: set
    score: float
    converged: bool

    @classmethod
    def from_json(cls, file: _FILE_TYPE):
        """Read a previously computed enrichment result from a json file

        Parameters
        ----------
        file: Union[str, pathlib.Path]
            Path to the file containing the enrichment results

        Returns
        -------
        EnrichmentResults
        """
        content = json.load(open(file, "r"))
        try:
            content = EnrichmentBaseResults._cast_field_dict(cls, content)
        except KeyError as kerr:
            raise KeyError(
                f"Missing data in '{file}'. {kerr}")
        return cls(**content)


@dataclass(frozen=True)
class RepeatedEnrichmentResults(EnrichmentBaseResults):
    """Object holding local optimization results over multiple repeats

    Attributes
    ----------
    subgraph: Set[int]
        Subnetwork with the best found objective function value per repeat
    scores: List[float]
        Objective function value of the solution per repeat
    score_progressions: List[List[float]]
        Objective score progressions per repeat
    """
    subgraph: set
    scores: list
    score_progressions: list

    @classmethod
    def from_json(cls, file: _FILE_TYPE):
        """Read a previously computed enrichment result from a json file

        Parameters
        ----------
        file: Union[str, pathlib.Path]
            Path to the file containing the enrichment results

        Returns
        -------
        RepeatedEnrichmentResults
        """
        content = json.load(open(file, "r"))
        try:
            content = EnrichmentBaseResults._cast_field_dict(cls, content)
        except KeyError as kerr:
            raise KeyError(
                f"Missing data in '{file}'. {kerr}")
        return cls(**content)


# Helper decorator functions
def _has_ls_results(func):
    """Decorator to check whether local search results are available"""
    def _check_results(*args, **kwargs):
        if args[0].non_empty_results():
            return func(*args, **kwargs)
        warn("No local search results yet!")
        return None

    return _check_results


def _check_numeric_type(func):
    """Decorator to check whether the passed parameter is the correct numeric
    type
    """
    def _check_size_t(*args, **kwargs):
        type_ = list(func.__annotations__.values())[0]
        if len(args) > 1:
            num_ = args[1]
        else:
            num_ = list(kwargs.values())[0]
        if not isinstance(num_, type_):
            try:
                num_ = type_(num_)
            except ValueError:
                raise ValueError(
                    "Parameter must be a numeric value, but found "
                    f"object of type '{type(num_).__name__}'"
                )
        if num_ < 0:
            raise ValueError(
                f"Passing a negative value is not allowed for {func.__name__}")
        # first positional argument => self
        return func(args[0], num_)

    return _check_size_t


# helper functions to simplify repeated enrichment merging
def _intersection(a, b):
    return a.intersection(b)


def _union(a, b):
    return a.union(b)


class LocalSearch(ABC):
    """Local search base class

    Interface for running local search with (predefined) objective functions
    using pre-computed node and edge values.

    Usually these values are coming from reaction activity approximation and
    the approximation of reaction/metabolite associations to other omics
    entities approximation, as implemented in
    :py:func:`pymantra.network.compute_reaction_estimates` and
    :py:func:`pymantra.network.compute_multiomics_associations`.

    If you want to use other node/edge metrics, you need to set/overwrite the
    values stored as 'data'. These values are used by the pre-defined objective
    functions **without** any further checks or corrections.

    If you are intending to use this function with a graph generated manually
    (i.e. with functions outside this module) it must contain the following
    node attributes:
    * 'data' or 'vec_data'
    * 'node_type'

    and the following edge attributes:
    * 'data'
    * 'edge_type'

    Attributes
    ----------
    lso : pymantra.network.enrichment.lso.LocalSearchOptimization
        Interface to the c++ class handling the local search
    """
    __slots__ = [
        'lso', '_solution', '_score_progression',
        '_objective_values', '_base_seed_size',
        '_reaction_nodes', "_data_nodes"
    ]

    lso: LocalSearchOptimization
    _solution: Union[
        None, List[EnrichmentResults], List[RepeatedEnrichmentResults]]
    _score_progression: Union[None, List[List[float]]]
    _objective_values: Union[
        None, Dict[str, Dict[Union[str, Tuple[str, str]], float]]
    ]
    _reaction_nodes: Union[None, Set[str]]
    _data_nodes: Set[str]

    def __init__(
        self, network: nx.Graph,
        temp: float, delta_min: float,
        l_min: int, l_max: int, max_iter: int,
        objective_function: str, min_reactions: int,
        p: float = 2, seed_size=10, *args, **kwargs
    ):
        """Initialize a LocalSearchOptimization object

        Regarding the local search objective it is possible to adapt or
        add new objective functions in general, but currently requires C++
        functions and re-compilation.

        The current implemented objective functions are

        1. 'metabolic_reactions': 'ld_reactions'

        2. 'reaction_microbe': 'reaction_microbiome'

        3. 'reaction_transcriptome': 'reaction_transcriptome'

        Parameters
        ----------
        network : nx.Graph
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
        # sanity check objective function selection
        self.lso = LocalSearchOptimization(
            network, temp, delta_min, l_min,
            l_max, max_iter, objective_function,
            p, min_reactions, *args, **kwargs
        )

        self._solution = None
        self._score_progression = None
        self._objective_values = None
        self._base_seed_size = seed_size

    def score_final_solution(self, groups: Tuple[str, str]) -> float:
        """Recompute the score of the final solution

        Recomputes the score fo the final subnetwork with the given
        groups. If `groups` contains the same groups as during the
        local search, the result will be equal to the score in
        :attr:`solution`. However, this function also enables the
        calculation of the objective score with *other* group combinations
        and the computed solution.

        Parameters
        ----------
        groups : Tuple[str, str]
            2-tuple containing the group names

        Returns
        -------
        float
            The objective function score
        """
        return self.lso.score_final_solution(groups)

    def non_empty_results(self) -> bool:
        """Check whether both solution and score progression are not empty"""
        return self._solution is not None and \
            self._score_progression is not None

    @abstractmethod
    def _prepare_graph(self):
        raise NotImplementedError(
            "'_prepare_graph' is only implemented in subclasses"
        )

    def run_local_search(
        self, groups: Tuple[str, str] = None, n_threads: int = 1,
        seed_size: int = None, min_comp_size: int = 25
    ):
        """Run a local search

        Do a local search optimization with parameters given by
        the instance's attributes for the binary comparison specified
        by `groups`

        Parameters
        ----------
        groups : Tuple[str, str]
            2-tuple containing the group names
        n_threads : int, default 1
            Number of threads to use.
            NOTE: currently multi-threading can cause unexpected termination
        seed_size : int, optional
            Option to specify the seed size. If None, the seed_size at
            :obj:`LocalSearch` object initialization is used.
        min_comp_size : int, default 25
            Minimum component size to run a local search on the component

        """
        graph = self._prepare_graph()

        if seed_size is None:
            seed_size = self._base_seed_size
        self._run_local_search(
            graph, groups, n_threads, seed_size, min_comp_size)

    def run_repeated_local_search(
        self, n_repeats: int, groups: Tuple[str, str] = None,
        combine_mode: str = "union", n_threads: int = 1,
        seed_size: int = None, min_comp_size: int = 25
    ):
        """Run a local search repeatedly n-times

        Do a local search optimization with parameters given by
        the instance's attributes for the binary comparison specified
        by `groups` n times to get a more robust result.

        Before each repetition a new random seed will be set.

        Results from the different repeats will be merged either by the union
        or intersection of the nodes contained in their subgraph-solutions.

        Parameters
        ----------
        n_repeats: int
            Number of repeated local searches to perform
        groups : Optional[Tuple[str, str]]
            2-tuple containing the group names
        combine_mode: str, "union"
            How to combine the results of all iterations. Either "union" or
            "intersection"
        n_threads : int, default 1
            Number of threads to use.
            NOTE: currently multi-threading can cause unexpected termination
        seed_size : int, optional
            Option to specify the seed size. If None, the seed_size at
            :obj:`LocalSearch` object initialization is used.
        min_comp_size : int, default 25
            Minimum component size to run a local search on the component
        """
        graph = self._prepare_graph()

        if seed_size is None:
            seed_size = self._base_seed_size
        if combine_mode == "union":
            self._run_repeated_local_search(
                graph, groups, n_repeats, n_threads, seed_size, min_comp_size,
                _union
            )
        elif combine_mode == "intersection":
            self._run_repeated_local_search(
                graph, groups, n_repeats, n_threads, seed_size, min_comp_size,
                _intersection
            )
        else:
            raise ValueError(
                f"Unknown option '{combine_mode}' for 'combine_mode'. Only "
                "'union' and 'intersection' are currently supported."
            )

    @staticmethod
    def _get_components_with_min_size(network, min_comp_size):
        components = [
            comp for comp in nx.connected_components(network)
            # TODO: what is a good size here?
            if len(comp) >= min_comp_size
        ]
        if not components:
            raise ValueError(
                f"No components with minimum size {min_comp_size} found. "
                "Consider increasing 'min_comp_size' to get valid components "
                "for the enrichment step."
            )
        return components

    def _run_local_search(
        self, network, groups, n_threads, seed_size, min_comp_size
    ):
        self._solution = list()
        self._score_progression = list()
        components = self._get_components_with_min_size(network, min_comp_size)
        for i, comp in tqdm(
                enumerate(components), total=len(components),
                desc="Local search per connected component"):
            # re-setting seed - _set_seed will throw and error if there
            # is no seed available which meets the minimum node
            # requirements
            try:
                self._set_seed(seed_size, comp)
            except GraphComponentError as gce:
                warn(
                    f"{gce.message} (component {i + 1}). "
                    "Continuing with next component"
                )
                continue

            # actual enrichment call + storing in class attributes
            subgraph, score, converged = \
                self.lso.run_local_search(groups, n_threads)
            self._solution.append(
                EnrichmentResults(subgraph, score, converged)
            )
            self._score_progression.append(self.lso.score_progression())

    def _run_repeated_local_search(
        self, network, groups, n_repeats, n_threads, seed_size, min_comp_size,
        combine_function
    ):
        self._solution = list()
        self._score_progression = list()

        components = self._get_components_with_min_size(network, min_comp_size)

        for i, comp in tqdm(
                enumerate(components), total=len(components),
                desc="Repeated local search per connected component"):
            solution = set()
            scores = []
            score_progressions = []
            for j in range(n_repeats):
                # re-setting seed - _set_seed will throw and error if there
                # is no seed available which meets the minimum node
                # requirements
                try:
                    self._set_seed(seed_size, comp)
                except GraphComponentError as gce:
                    warn(
                        f"{gce.message} (component {i}). "
                        "Continuing with next component"
                    )
                    break

                # actual enrichment call + storing
                subgraph, score, _ = \
                    self.lso.run_local_search(groups, n_threads)
                solution = combine_function(solution, subgraph)
                scores.append(score)
                score_progressions.append(self.lso.score_progression())

            # adding results to class attributes
            self._solution.append(
                RepeatedEnrichmentResults(
                    solution, scores, score_progressions)
            )
            self._score_progression.append(score_progressions)

    def set_seed(self, seed: str, seed_size: int = None):
        """Set the seed at which local search is starting

        'Manually' setting the local search seed. When calling
        :meth:`run_local_search` a seed is automatically computed and
        cached for re-usage when running another local search.
        If you want multiple local search runs with independent seeds
        this method should be called in between runs.

        Parameters
        ----------
        seed : Union[List[str], Set[str], str]
            Either a node (specified by name) or an iterable collection of
            nodes. If given a single node `seed_size` must be set to specify
            the number of neighbours of the seed node to draw.
            Otherwise, the size of the collection must bei in the range of
            [l_min, l_max]
        seed_size : Optional[int]
            Size of the seed subgraph. Required if `seed` is a single node,
            otherwise ignored.
        """
        self.lso.set_seed(seed, seed_size)

    def _set_seed(
        self, seed_size: int, seed_options: Set[str]
    ):
        # setting the seed to be within the largest connected component
        comp_nodes = list(seed_options)
        np.random.shuffle(comp_nodes)

        seed_node = None
        for node in comp_nodes:
            if node in self._data_nodes:
                seed_node = node
                break
        if seed_node is None:
            raise GraphComponentError(
                "No reaction node (node_type == '"
                f"{NODE_TYPE_NAMES['reaction']}') with data found in "
                "the component"
            )
        self.lso.set_seed(seed_node, seed_size)

    @property
    @_has_ls_results
    def solution(self):
        return self._solution

    @property
    @_has_ls_results
    def score_progression(self):
        if self._score_progression is None:
            warn("No local search results yet!")
        return self._score_progression

    @property
    def objective_values(self):
        if self._objective_values is None:
            try:
                self._objective_values = self.lso.objective_values
            except ValueError as ve:
                if self._solution is not None and \
                        self._score_progression is not None:
                    raise ValueError(
                        "Objective values not computed yet. Please call "
                        "'run_local_search' first!"
                    )
                else:
                    raise ValueError(
                        f"Unknown error occurred: {str(ve)}"
                    )
        return self._objective_values

    @property
    def converged(self) -> bool:
        """Whether the local search converged or terminated due to reaching
        max_iter
        """
        return self.lso.converged()

    @property
    def l_min(self) -> int:
        """Current choice of minimum solution size"""
        return self.lso.l_min()

    @property
    def min_reactions(self) -> int:
        """Current choice of minimum number of reactions in the solution"""
        return self.lso.min_reactions()

    @property
    def l_max(self) -> int:
        """Current choice of maximum solution size"""
        return self.lso.l_max()

    @property
    def temp(self) -> float:
        """Current choice of simulated annealing temperature"""
        return self.lso.temp()

    @property
    def max_iter(self) -> int:
        """Current choice of maximum number of iteration for local search"""
        return self.lso.max_iter()

    @property
    def delta_min(self) -> float:
        """Current choice of minimum progress per iteration"""
        return self.lso.delta_min()

    @property
    def p(self) -> float:
        """Current choice of lp-norm"""
        return self.lso.p()

    @_check_numeric_type
    def set_l_min(self, lmin: int):
        """Set the choice of minimum solution size"""
        return self.lso.set_l_min(lmin)

    @_check_numeric_type
    def set_min_reactions(self, min_: int):
        """Set the choice of minimum number of reactions in the solution"""
        self.lso.set_min_reactions(min_)

    @_check_numeric_type
    def set_l_max(self, lmax: int):
        """Set the choice of maximum solution size"""
        return self.lso.set_l_max(lmax)

    @_check_numeric_type
    def set_temp(self, temp: float):
        """Set the choice of simulated annealing temperature"""
        return self.lso.set_temp(temp)

    @_check_numeric_type
    def set_max_iter(self, max_iter: int):
        """Set the choice of maximum number of iteration for local search"""
        return self.lso.set_max_iter(max_iter)

    @_check_numeric_type
    def set_p(self, p: float):
        """Set the choice of lp-norm to use. Only relevant for certain
        objective functions
        """
        return self.lso.set_p(p)

    @_check_numeric_type
    def set_delta_min(self, delta_min: float):
        """Set the choice of minimum progress per iteration"""
        return self.lso.set_delta_min(delta_min)

    @_has_ls_results
    def plot_score_progression(self, ax=None):
        """Plot the progression of the objective score function
        per iteration as a combined line- and scatter-plot

        Parameters
        ----------
        ax : Optional[Union[plt.axis, List[plt.axis]]]
            matplotlib.pyplot axis to plot on

        Returns
        -------
        plt.axis | List[plt.axis]
            Plot axis containing score progression plot or list of axes if
            repeated enrichment was called last
        """
        # this means we have a non-repeated local search
        if isinstance(self._solution[0], EnrichmentResults):
            if ax is None:
                _, ax = plt.subplots(figsize=(16, 9))
            for i, prog in enumerate(self._score_progression):
                n_iter = np.arange(len(prog))
                ax.plot(n_iter, prog, "o-", label=f"Component {i}")
            ax.set_xlabel("Iteration")
            ax.set_ylabel("Objective Function Score")
            ax.legend()
        # repeated local search
        else:
            # pre-plotting checks
            n_comps = len(self._solution)
            if ax is None:
                # TODO: add an option here to have at max n columns
                fig, ax = plt.subplots(figsize=(16, 9), ncols=n_comps)
                if n_comps == 1:
                    ax = [ax]
            else:
                if not isinstance(ax, list):
                    raise ValueError(
                        "'ax' must be a list of matplotlib.axis objects if "
                        "the last enrichment call was to a repeated "
                        "enrichment method, not a single axis object."
                    )
                fig = ax[0].get_figure()
            # plotting one subplot per component and each iteration with a
            # different colour
            for i, solution in enumerate(self._solution):
                for j, prog in enumerate(solution.score_progressions):
                    n_iter = np.arange(len(prog))
                    ax[i].plot(n_iter, prog, "o-", label=f"Repetition {j + 1}")
                ax[i].set_title(f"Component {i}")
                # legend is only generated for the last component
                if i == n_comps - 1:
                    ax[i].legend()
            # common axis labelling
            fig.text(0.5, 0.02, 'Iteration', ha='center')
            fig.text(
                0.02, 0.5, 'Objective Function Score', va='center',
                rotation='vertical'
            )

        return ax

    @staticmethod
    def _extract_subgraph(
        network: nx.DiGraph, subnodes: Set[str],
        from_reaction_graph_only: bool = False
    ) -> nx.DiGraph:
        """Extract the metabolite-reaction subgraph indicated by the
        local search solution
        """
        if from_reaction_graph_only:
            def _check_sugraph_nodes(src_node, tgt_node):
                return src_node in subnodes and tgt_node in subnodes
        else:
            def _check_sugraph_nodes(src_node, tgt_node):
                return src_node in subnodes or tgt_node in subnodes

        subgraph = nx.DiGraph()
        for src, tgt in network.edges:
            if _check_sugraph_nodes(src, tgt):
                if src not in subgraph.nodes:
                    subgraph.add_node(
                        src, **network.nodes[src]
                    )
                if tgt not in subgraph.nodes:
                    subgraph.add_node(
                        tgt, **network.nodes[tgt]
                    )
                subgraph.add_edge(
                    src, tgt, **network.edges[(src, tgt)]
                )
        return subgraph

    @abstractmethod
    def _build_from_solution(self, subgraph, *args, **kwargs):
        raise NotImplementedError(
            "Building a LocalSearch object from a solution is currently "
            "not implemented."
        )

    def _plot_subnetwork(
        self, network: nx.DiGraph = None,
        subplot_args: dict = None, **kwargs
    ) -> Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]:
        """Plot the subgraph returned by local search optimization

        Plotting the edge subgraph returned by local search optimization.
        Either plots a reaction-reaction-organism graph containing exactly the
        edges in :attr:`solution` or a metabolite-reaction network
        containing all nodes :attr:`solution` and their metabolite
        neighbours.

        Parameters
        ----------
        network : nx.DiGraph, optional
            Same metabolite-reaction graph used as input in the constructor.
            It is assumed to only contain *directed* metabolite-reaction edges.
        subplot_args : dict, optional
            Keyword arguments to pass to plt.subplots.
        kwargs
            Optional keyword arguments passed to
            :meth:`pymantra.plotting.plot_undirected_graph` or
            :meth:`pymantra.plotting.plot_directed_graph` depending on
            whether `network` is None.
            Note that `reaction_graph` cannot be passed.

        Returns
        -------
        Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]
            2-Tuple with the first element being the figure on which the
            subplots are lying and the second a list or array of plt.axis
            which are drawn in the figure.
        """
        kwargs.pop('reaction_graph', None)

        if subplot_args is None:
            subplot_args = {}

        ncomps = len(self._solution)
        max_cols = 3
        plot_on_input = False
        if ncomps == 1:
            if "ax" in kwargs.keys():
                axes = [kwargs.pop("ax")]
                fig = None
                plot_on_input = True
            else:
                fig, ax = plt.subplots(
                    figsize=subplot_args.pop('figsize', (16, 9)),
                    **subplot_args
                )
                axes = [ax]
        elif ncomps <= max_cols:
            if kwargs.pop("ax", None) is not None:
                warn(
                    "'ax' cannot be set when the input graph contains "
                    "multiple components"
                )
            fig, axes = plt.subplots(
                figsize=subplot_args.pop('figsize', (16, 9)),
                ncols=subplot_args.pop('ncols', ncomps),
            )
        else:
            if kwargs.pop("ax", None) is not None:
                warn(
                    "'ax' cannot be set when the input graph contains "
                    "multiple components"
                )
            nrows = ceil(ncomps / max_cols)
            fig, axes = plt.subplots(
                figsize=subplot_args.pop('figsize', (16, 9)),
                ncols=subplot_args.pop('ncols', max_cols),
                nrows=subplot_args.pop('nrows', nrows),
            )
        if network is not None:
            for i, (solution, ax) in enumerate(zip(self._solution, axes)):
                graph = self._extract_subgraph(network, solution.subgraph)
                ax = plot_directed_graph(
                    graph, reaction_graph=True, ax=ax, **kwargs
                )
                if len(self._solution) > 1:
                    ax.set_title(f"Component {i + 1}")
                ax.axis('off')
        else:
            # NOTE: this case should NOT happen right now => we want the node
            # attributes for plotting
            for i, (solution, ax) in enumerate(zip(self._solution, axes)):
                node_types = kwargs.pop("node_types", None)
                if node_types is None:
                    raise ValueError(
                        "The optional 'node_types' argument must be given, "
                        "when 'network' is not parsed"
                    )
                graph = self._build_from_solution(
                    solution.subgraph, node_types=node_types)
                ax = plot_undirected_graph(
                    graph, reaction_graph=False, ax=ax, **kwargs
                )
                ax.set_title(f"Component {i}")
                ax.axis('off')
        if plot_on_input:
            return axes[0]
        return fig, axes

    @abstractmethod
    def plot_subnetwork(self, *args, **kwargs):
        raise NotImplementedError


class MetaboliteLocalSearch(LocalSearch):
    """Local search class for metabolomics-data only

    Interface for running local search for metabolomics-only data with
    (predefined) objective functions using pre-computed node values.

    Usually these values are coming from reaction activity approximation as
    implemented in :py:func:`pymantra.network.compute_reaction_estimates`.

    If you want to use other node/edge metrics, you need to set/overwrite the
    values stored as 'data'. These values are used by the pre-defined objective
    functions **without** any further checks or corrections.

    If you are intending to use this function with a graph generated manually
    (i.e. with functions outside this module) it must contain the following
    node attributes:
    * 'data' or 'vec_data'
    * 'node_type'

    and the following edge attributes:
    * 'data'
    * 'edge_type'

    Attributes
    ----------
    lso : LocalSearchOptimization
        Interface to the c++ class handling the local search
    reaction_edges: Set[Tuple[str, str]]
        reaction-reaction edges either passed by or extracted from the graph
        passed to the constructor
    """
    __slots__ = ["reaction_edges"]

    reaction_edges: set

    def __init__(
        self, network: nx.Graph, temp: float, delta_min: float,
        l_min: int, l_max: int, max_iter: int, p: float = 2,
        objective_function="ld_reactions",
        is_reaction_graph: bool = False, **kwargs
    ):
        """Initialize a MetaboliteLocalSearch object

        Initializes a specialized subclass of :class:`LocalSearchOptimization`
        meant to compute a local search on metabolite-reaction graphs without
        additional node types.

        Parameters
        ----------
        network : nx.Graph
            Either a metabolite-reaction graph or a reaction-reaction
            graph extracted from a metabolite-reaction graph.
            Usually `network` will be computed using
            :class:`pymantra.database.NetworkGenerator` and/or
            :meth:`reaction_graph_extraction`
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
        p : Optional[float], default 2.
            Which :math:`L^p`(Minkowski)-norm to use in the objective
            function.
        objective_function : Optional[str], default 'metabolic_reactions'
            Currently changing this option is not supported
        is_reaction_graph : Optional[bool], default False
            Whether `network` is a reaction-reaction or a
            metabolite-reaction graph
        kwargs
            Keyword arguments to pass to :class:LocalSearchOptimization`.
            Note that "min_reactions" is automatically set to `l_min`, passing
            it will have no effect.
        """
        # automatically overwriting "min_reactions
        kwargs["min_reactions"] = l_min

        if is_reaction_graph:
            self.reaction_edges = set(network.edges)
            self._reaction_nodes = set(network.nodes)

            self._data_nodes = set()
            for node in self._reaction_nodes:
                if "data" in network.nodes[node].keys() or \
                        "vec_group_data" in network.nodes[node].keys():
                    self._data_nodes.add(node)

            super().__init__(
                network, objective_function=objective_function,
                temp=temp, delta_min=delta_min, l_min=l_min,
                l_max=l_max, max_iter=max_iter, p=p, **kwargs
            )
            if not self.reaction_edges:
                raise ValueError(
                    "No reaction-reaction edges were passed, but are required "
                    "when 'is_reaction_graph' is True."
                )
        else:
            # computing reactions-reaction edges
            self.reaction_edges = reaction_graph_extraction(
                network, include_attributes=False)
            self._reaction_nodes = {
                node for edge in self.reaction_edges
                for node in edge
            }

            reaction_graph = nx.Graph()
            self._data_nodes = set()
            for node in self._reaction_nodes:
                if "data" in network.nodes[node].keys() or \
                        "vec_group_data" in network.nodes[node].keys():
                    reaction_graph.add_node(
                        node, **network.nodes[node]
                    )
                    self._data_nodes.add(node)
            if len(self._data_nodes) <= l_min:
                raise ValueError(
                    f"Only {len(self._data_nodes)} reaction nodes with data "
                    "found, but the minimum size of the solution is set to "
                    f"{l_min}. Either update the reaction nodes with data or "
                    "decrease the 'l_min' parameter."
                )
            for edge in self.reaction_edges:
                if edge[0] in self._data_nodes and edge[1] in self._data_nodes:
                    reaction_graph.add_edge(
                        edge[0], edge[1], edge_type=_REACTION_EDGE)

            # build LSO object only on reactions nodes
            super().__init__(
                reaction_graph,
                objective_function=objective_function,
                temp=temp, delta_min=delta_min, l_min=l_min,
                l_max=l_max, max_iter=max_iter, p=p, **kwargs
            )
            # sanity checking subgraph
            if not self.reaction_edges:
                raise ValueError(
                    "No reaction-reaction and reaction-microbe edges were "
                    "found in 'network'. "
                    "Please make sure node types are assigned correctly and "
                    "graph generation worked properly."
                )

    def _prepare_graph(self):
        graph = nx.from_edgelist(self.reaction_edges)
        nx.set_node_attributes(graph, values='reaction', name='node_type')
        return graph

    def plot_subnetwork(
        self, network: nx.DiGraph = None, subplot_args: dict = None, **kwargs
    ) -> Union[plt.axis, Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]]:
        """Plot the subgraph returned by local search optimization

        Plotting the edge subgraph returned by local search optimization.
        Either plots a reaction-reaction-organism graph containing exactly the
        edges in :attr:`solution` or a metabolite-reaction network
        containing all nodes :attr:`solution` and their metabolite
        neighbours.

        Parameters
        ----------
        network : nx.DiGraph, optional
            Same metabolite-reaction graph used as input in the constructor.
            It is assumed to only contain *directed* metabolite-reaction edges.
        subplot_args : dict, optional
            Keyword arguments to pass to plt.subplots.
        kwargs
            Optional keyword arguments passed to
            :meth:`pymantra.plotting.plot_undirected_graph` or
            :meth:`pymantra.plotting.plot_directed_graph` depending on
            whether `network` is None.
            Note that `reaction_graph` cannot be passed.

        Returns
        -------
        Union[plt.axis, Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]]
            Either a single matplotlib axis object, if 'ax' is given as a
            keyword argument or a
            2-tuple with the first element being the figure on which the
            subplots are lying and the second a list or array of plt.axis
            which are drawn in the figure.
        """
        return self._plot_subnetwork(network, subplot_args, **kwargs)

    def _build_from_solution(self, subgraph, *args, **kwargs):
        graph = nx.Graph()
        for src, tgt in self.reaction_edges:
            if src in subgraph and tgt in subgraph:
                graph.add_edge(
                    src, tgt, edge_type=_REACTION_EDGE
                )
        nx.set_node_attributes(graph, NODE_TYPE_NAMES["reaction"], "node_type")
        return graph


class MultiOmicsLocalSearch(LocalSearch):
    """Local search class for metabolomics with multi-omics data

    Interface for running local search with (predefined) objective functions
    using pre-computed node and edge values.

    Usually these values are coming from reaction activity approximation and
    the approximation of reaction/metabolite associations to other omics
    entities approximation, as implemented in
    :py:func:`pymantra.network.compute_reaction_estimates` and
    :py:func:`pymantra.network.compute_multiomics_associations`.

    If you want to use other node/edge metrics, you need to set/overwrite the
    values stored as 'data'. These values are used by the pre-defined objective
    functions **without** any further checks or corrections.

    If you are intending to use this function with a graph generated manually
    (i.e. with functions outside this module) it must contain the following
    node attributes:
    * 'data' or vec_data
    * 'node_type'

    and the following edge attributes:
    * 'data'
    * 'edge_type'

    Attributes
    ----------
    lso : LocalSearchOptimization
        Interface to the c++ class handling the local search
    reaction_multiomics_edges: Set[Tuple[str, str]]
        reaction-reaction and multi omics-reactions edges either passed by or
        extracted from the graph passed to the constructor
    """
    __slots__ = ['reaction_multiomics_edges', '_node_type', '_edge_type']

    reaction_multiomics_edges: set
    _node_type: str
    _edge_type: str

    def __init__(
        self, network: nx.Graph, omics: str, temp: float, delta_min: float,
        l_min: int, l_max: int, max_iter: int, min_reactions: int,
        p: float = 2, is_reaction_graph: bool = False, **kwargs
    ):
        """Initialize a MetaboliteLocalSearch object

        Initializes a specialized subclass of :class:`LocalSearchOptimization`
        meant to compute a local search on metabolite-reaction graphs without
        additional node types.

        Parameters
        ----------
        network : nx.Graph
            Either a metabolite-reaction graph or a reaction-reaction
            graph extracted from a metabolite-reaction graph.
            Usually `network` will be computed using
            :class:`NetworkGenerator` and/or
            :meth:`reaction_graph_extraction`
        omics : str
            Type of multi-omics association. Currently, this must be either
            "organism" (for microbiome) or "gene" (for transcriptome or
            metagenome)
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
        min_reactions : int
            Minimum number of reactions to be contained in the solution.
        p : float, default 2.
            Which :math:`L^p`(Minkowski)-norm to use in the objective
            function.
        is_reaction_graph : bool, default False
            Whether `network` is a reaction-reaction or a
            metabolite-reaction graph
        kwargs
            Keyword arguments to pass to :class:LocalSearchOptimization`
        """
        if omics == "organism":
            objective_function = "reaction_microbiome"
            self._node_type = NODE_TYPE_NAMES['organism']
        elif omics == "gene":
            objective_function = "reaction_transcriptome"
            self._node_type = NODE_TYPE_NAMES['gene']
        else:
            raise ValueError(
                f"{omics} is not a supported option for 'omics'. At the "
                "moment only 'organism' and 'gene' are supported."
            )

        self._edge_type = EDGE_BY_NODE_TYPE[
            (NODE_TYPE_NAMES['reaction'], self._node_type)]

        if is_reaction_graph:
            self.reaction_multiomics_edges = set(network.edges)
            self._reaction_nodes = {
                node for node in network.nodes
                if network.nodes['node_type'] == NODE_TYPE_NAMES['reaction']
            }

            self._data_nodes = set()
            for node in network.nodes:
                if "data" in network.nodes[node].keys() or \
                        "vec_group_data" in network.nodes[node].keys():
                    self._data_nodes.add(node)

            super().__init__(
                network=network,
                objective_function=objective_function,
                temp=temp, delta_min=delta_min, l_min=l_min,
                l_max=l_max, max_iter=max_iter, p=p,
                min_reactions=min_reactions, **kwargs
            )
        else:
            # TODO: checks for type names etc.
            # computing reactions-reaction edges
            self.reaction_multiomics_edges = reaction_graph_extraction(
                network, include_attributes=False)
            self._reaction_nodes = {
                node for edge in self.reaction_multiomics_edges
                for node in edge
            }
            # adding in reaction-microbe edges
            _mo_reactions = set()
            for edge, edge_type in \
                    nx.get_edge_attributes(network, "edge_type").items():
                if edge_type == self._edge_type:
                    self.reaction_multiomics_edges.add(edge)
                    _mo_reactions.add(edge)
            # sanity checking subgraph
            if not self.reaction_multiomics_edges:
                raise ValueError(
                    "No reaction-reaction and reaction-microbe edges were "
                    "found in 'network'. "
                    "Please make sure node types are assigned correctly "
                    "and graph generation worked properly."
                )

            reaction_nodes = {
                node for edge in self.reaction_multiomics_edges
                for node in edge
            }
            reaction_graph = nx.Graph()
            self._data_nodes = set()
            for node in reaction_nodes:
                if "data" in network.nodes[node].keys() or \
                        "vec_group_data" in network.nodes[node].keys():
                    reaction_graph.add_node(
                        node, **network.nodes[node]
                    )
                    self._data_nodes.add(node)
                elif network.nodes[node]["node_type"] != \
                        NODE_TYPE_NAMES["reaction"]:
                    reaction_graph.add_node(
                        node, **network.nodes[node]
                    )

            if len(self._data_nodes) <= min_reactions:
                raise ValueError(
                    f"Only {len(self._data_nodes)} reaction nodes with "
                    "data found, but the minimum number of reaction nodes in "
                    f"the solution is set to {min_reactions}. Either update "
                    "the reaction nodes with data or decrease the "
                    "'min_reactions' parameter."
                )
            n_nodes = reaction_graph.number_of_nodes()
            if n_nodes <= l_min:
                raise ValueError(
                    f"Only {n_nodes} total nodes, but the "
                    f"minimum size of the solution is set to {l_min}. "
                    "Either update the reaction nodes with data or decrease "
                    "the 'l_min' parameter."
                )
            for edge in self.reaction_multiomics_edges:
                if edge in _mo_reactions:
                    if edge[0] in self._data_nodes:
                        reaction_graph.add_edge(
                            edge[0], edge[1], **network.edges[edge])
                else:
                    if edge[0] in self._data_nodes and \
                            edge[1] in self._data_nodes:
                        reaction_graph.add_edge(
                            edge[0], edge[1], edge_type=_REACTION_EDGE)

            # TODO: sanity check node/edge values
            super().__init__(
                network=reaction_graph,
                objective_function=objective_function,
                temp=temp, delta_min=delta_min, l_min=l_min,
                l_max=l_max, max_iter=max_iter, p=p,
                min_reactions=min_reactions, **kwargs
            )

    def _prepare_graph(self):
        graph = nx.from_edgelist(self.reaction_multiomics_edges)
        node_types = {
            node: NODE_TYPE_NAMES["reaction"] if node in self._reaction_nodes
            else self._node_type
            for node in graph.nodes
        }
        nx.set_node_attributes(graph, values=node_types, name='node_type')
        return graph

    def _get_edge_type(self, src, tgt, network):
        src_type = network.nodes[src]['node_type']
        tgt_type = network.nodes[tgt]['node_type']
        if src_type == tgt_type:
            return _REACTION_EDGE
        return self._edge_type

    def plot_subnetwork(
        self, network: nx.DiGraph = None, subplot_args: dict = None, **kwargs
    ) -> Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]:
        """Plot the subgraph returned by local search optimization

        Plotting the edge subgraph returned by local search optimization.
        Either plots a reaction-reaction-organism graph containing exactly the
        edges in :attr:`solution` or a metabolite-reaction network
        containing all nodes :attr:`solution` and their metabolite
        neighbours.

        Parameters
        ----------
        network : nx.DiGraph, optional
            Same metabolite-reaction graph used as input in the constructor.
            It is assumed to only contain *directed* metabolite-reaction edges.
        subplot_args : dict, optional
            Keyword arguments to pass to plt.subplots.
        kwargs
            Optional keyword arguments passed to
            :meth:`pymantra.plotting.plot_undirected_graph` or
            :meth:`pymantra.plotting.plot_directed_graph` depending on
            whether `network` is None.
            Note that `reaction_graph` cannot be passed.

        Returns
        -------
        Union[plt.axis, Tuple[plt.figure, Union[np.ndarray, List[plt.axis]]]]
            Either a single matplotlib axis object, if 'ax' is given as a
            keyword argument or a
            2-Tuple with the first element being the figure on which the
            subplots are lying and the second a list or array of plt.axis
            which are drawn in the figure.
        """
        return self._plot_subnetwork(
            network, subplot_args, **kwargs)

    def _build_from_solution(self, subgraph, *args, **kwargs):
        node_types = kwargs.get("node_types")
        graph = nx.Graph()
        for src, tgt in self.reaction_multiomics_edges:
            if src in subgraph and tgt in subgraph:
                if src not in graph.nodes:
                    graph.add_node(src, node_type=node_types.get(src))
                if tgt not in graph.nodes:
                    graph.add_node(tgt, node_type=node_types.get(tgt))
                if graph.nodes[tgt]["node_type"] == \
                        NODE_TYPE_NAMES["reaction"]:
                    edge_type = _REACTION_EDGE
                else:
                    edge_type = self._edge_type
                graph.add_edge(src, tgt, edge_type=edge_type)
        return graph
