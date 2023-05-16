#ifndef ENRICHMENT_LOCALSEARCH_HPP
#define ENRICHMENT_LOCALSEARCH_HPP


#define PY_SSIZE_T_CLEAN


#include <functional>
#include <list>
#include <utility>
#include <vector>
#include <string>
#include <boost/unordered_set.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <Python.h>
#include "reaction_graph.hpp"
#include "objective_functions.hpp"
#include "lso_utils.hpp"


using std::function;
using std::list;
using std::vector;
using std::vector;
using std::pair;
using boost::unordered_set;
using boost::unordered_map;
using boost::tuple;
using boost::vertex_index_t;


class LocalSearch {
public:
    BestAction solution;
    vector<double> score_progression;
    bool converged{};
    vector<vertex_subgraph> solutions;
    vector<Action> action_progression;
    vertex_subgraph seed{};

    bool precompute_objectives{};
    objective_precompute precompute{};
    string objectives_computed[2]{"", ""};
    float p_norm{};

    /**
     * Construct a `LocalSearch` object from python objects
     *
     * @param node_properties Nodes and their properties as a python dictionary where keys are the
     *                        node names and values are the node properties as a dictionary
     * @param edge_properties Edges and their properties as a python dictionary where keys are
     *                        2-tuples representing the edges and values are the edge properties as
     *                        a python dictionary
     * @param objective_function objective function specifier
     * @param temperature simulated annealing temperature
     * @param delta_min minimum improvement to continue running the local search
     * @param l_min minimum solution size
     * @param l_max maximum solution size
     * @param min_reactions minimum number of reaction nodes in the solution
     * @param max_iter maximum number of iterations to run
     * @param p which p-norm to use when computing objective function scores
     */
    LocalSearch(
        PyObject *node_properties,
        PyObject *edge_properties,
        int objective_function,
        double temperature,
        double delta_min,
        size_t l_min,
        size_t l_max,
        size_t min_reactions,
        size_t max_iter,
        float p
    );

    /**
     * Construct a `LocalSearch` object from pure c++ code
     *
     * @param node_properties Vector of node properties. Nodes appear in the graph in the same order
     *                        as provided here.
     * @param edge_properties Edges and their properties as an unordered map where keys are
     *                        pairs representing the edges (node indices) and values are the
     *                        edge properties
     * @param objective_function objective function specifier
     * @param temperature simulated annealing temperature
     * @param delta_min minimum improvement to continue running the local search
     * @param l_min minimum solution size
     * @param l_max maximum solution size
     * @param min_reactions minimum number of reaction nodes in the solution
     * @param max_iter maximum number of iterations to run
     * @param p which p-norm to use when computing objective function scores
     */
    LocalSearch(
        vector<vertex_props> &node_properties,
        unordered_map<pair<size_t, size_t>, edge_props> &edge_properties,
        int objective_function,
        double temperature,
        double delta_min,
        size_t l_min,
        size_t l_max,
        size_t min_reactions,
        size_t max_iter,
        float p
    );

    // TODO: optional constructors that allow to set the objective function freely

    bool get_converged() const;

    /**
     * Get the best solution-subgraph as a python set of nodes (strings)
     *
     * @return set of nodes forming the subgraph of the final solution
     */
    PyObject *get_best_solution();

    /**
     * Get the progression of objective scores over the local search
     *
     * @return vector of scores in order of computation throughout the optimization
     */
    vector<double> get_score_progression() const;

    /**
     * Score the final solution for a comparison of two groups. These do not
     * necessarily need to be the groups over which the local search was performed.
     * In this case, it is of course not guaranteed that the solution is optimal
     *
     * @param groups subgraph according to the objective function
     *
     * @return objective function score of the solution for the given groups
     */
    double score_solution(vector<string>& groups);

    /**
     * Score a given solution
     *
     * @param groups sample groups for which the comparison should be performed. Only the first two
     *               are considered
     * @param solution_ subgraph to be scored
     *
     * @return objective function score for the given subgraph with the specified groups
     */
    double score_solution(vector<string>& groups, vertex_subgraph &solution_);

    /**
     * Retrieve the best score found during the local search
     * @return best objective function score
     */
    double get_best_score() const;

    /**
     * Run a local search optimizing the objective function
     *
     * @param groups sample groups to compare. Only the first two are considered
     * @param n_threads number of threads to use. This is only relevant if compiled with OpenMP support
     */
    void run_local_search(vector<string> &groups, int n_threads = 1);

    /**
     * Set the initial local search optimization seed
     */
    void set_seed(vertex_subgraph seed_);

    /**
     * Set the initial local search optimization seed based
     * on a single node
     */
    void set_seed(size_t seed_size);

    /**
     * Set the initial local search optimization seed
     * from python
     */
    void set_seed_py(PyObject* seed_init);

    /**
     * Set the initial local search optimization seed based
     * on a single node from python
     */
    void set_seed_py(const string& seed_init, size_t seed_size);

    size_t get_lmax() const;
    void set_lmax(size_t lmax);

    size_t get_lmin() const;
    void set_lmin(size_t lmin);

    size_t get_min_reactions() const;
    void set_min_reactions(size_t min_reactions_);

    double get_temp() const;
    void set_temp(double temp_);

    size_t get_maxiter() const;
    void set_maxiter(size_t iter);

    double get_deltamin() const;
    void set_deltamin(double delta);

    double get_lp_norm() const;
    void set_lp_norm(float p);

    void precompute_objectives_py(vector<string>& groups);
    PyObject *get_objective_values();

    /**
     * Merge three **sorted** vectors into one sorted vector.
     * The function fail silently if the input vectors are not sorted
     * according to the comparison function!
     *
     * @tparam T type of vectors to be merged
     * @tparam F function type used for merging
     * @param ins first vector to merge
     * @param del second vector to merge
     * @param subs third vector to merge
     * @param comp_fun comparison function. Needs to implement a `greater` function between two instances of T
     *
     * @return sorted vector with a size of the sum of sizes of the three input vectors
     */
    template<class T, class F>
    static vector<T> merge_solutions(vector<T>& ins, vector<T>& del, vector<T>& subs, F comp_fun);

private:
    Graph graph;
    objective_ objective_function{};
    double delta_min{};
    double temp{};
    size_t l_min{};
    size_t l_max{};
    size_t min_reactions{};
    size_t max_iter{};
    bool multi_omics{};

    /**
     * Helper function for the constructors implementing all construction parts that are independent
     * of whether the constructor is called from pure c++ or the python bindings
     *
     * @param objective_function_ objective function specifier
     * @param temperature_ simulated annealing temperature
     * @param delta_min_ minimum improvement to continue running the local search
     * @param l_min_ minimum solution size
     * @param l_max_ maximum solution size
     * @param min_reactions_ minimum number of reaction nodes in the solution
     * @param max_iter_ maximum number of iterations to run
     * @param p which p-norm to use when computing objective function scores
     */
    void construct_params(int objective_function_, double temperature_, double delta_min_, size_t l_min_,
                          size_t l_max_, size_t min_reactions_, size_t max_iter_, float p);

    /**
     * Subset this->graph to a vertex subgraph containing
     * all vertices given. The returned graph has the same
     * type as this->graph
     *
     * @param vertices vertices to be included in the subgraph
     *
     * @return subset of this->graph containing only the vertices and all connections between them
     */
    Graph subgraph_from_vertices(vertex_subgraph& vertices);

    /**
     * Find all articulation points in a (sub)graph.
     * This function is used to exclude vertices whose deletion
     * would trigger a disconnection of the solution
     *
     * @param subgraph subgraph as a set of vertices
     *
     * @return list of vertices to be excluded from deletion
     */
    unordered_set<size_t> find_articulation_points(vertex_subgraph& vertices);

    /**
     * Generate a random subgraph of a given size (compliant with min_reactions) by doing a
     * random walk
     *
     * @param vertex initial vertex to start from
     * @param subgraph_size desired size of the random subgraph
     *
     * @return set of vertices defining the subgraph
     */
    vertex_subgraph random_walk_subgraph(vertex_index_t vertex, size_t subgraph_size);

    /**
     * Find all direct neighbours of a subset of vertices in a given (sub)graph.
     *
     * @param curr_solution current solution in the form of a vertex subgraph
     *
     * @return unordered set of first order neighbours
     */
    unordered_set<size_t> one_hop_neighbours(vertex_subgraph& curr_solution);

    /**
     * Update the currently best found solution
     *
     * @param prev previously best solution
     * @param score new objective function score
     * @param action action taken to get from the best solution in the last iteration to the
     *               new solution
     * @param solution new solution as a vertex subgraph
     * @param vertex_1 first vertex affected. If action is a deletion or insertion this is the
     *                 affected node, in case of a substitution this is the inserted node
     * @param vertex_2 second vertex affected. For deletion and insertion this is -1, for a
     *                 substitution the deleted node index
     */
    static void update_solution(BestAction& prev, double score, Action action,
                                vertex_subgraph solution_, int vertex_1, int vertex_2,
                                size_t n_reactions);

    /**
     * Score the change of the objective function for any
     * possible insertion
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param neighbours directly neighbouring nodes of the current solution
     * @param best_solution optimal solution for the current iteration
     * @param moved whether a valid move was already found in this iteration
     *
     * @return whether a valid was taken of previously found in this iteration
     */
    bool score_insertions(vector<string> &groups, vertex_subgraph curr_solution, unordered_set<size_t> &neighbours,
                          BestAction &best_solution, bool moved);

    /**
     * Score the change of the objective function for any possible deletions. This function is
     * only available when compiled with OpenMP support.
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param neighbours directly neighbouring nodes of the current solution
     * @param threads number of threads to use
     * @param n_reactions number of reactions in the current solution
     *
     * @return vector of all valid actions in decreasing objective score order
     */
    vector<BestAction>
    score_insertions(vector<string> &groups, vertex_subgraph &curr_solution, unordered_set<size_t> &neighbours,
                     int threads, unsigned n_reactions);

    /**
     * Score the change of the objective function for any
     * possible deletion
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param best_solution optimal solution for the current iteration
     * @param moved whether a valid move was already found in this iteration
     *
     * @return whether a valid was taken of previously found in this iteration
     */
    bool score_deletions(vector<string> &groups, vertex_subgraph curr_solution, BestAction &best_solution, bool moved);

    /**
     * Score the change of the objective function for any
     * possible deletion
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param neighbours directly neighbouring nodes of the current solution
     * @param threads number of threads to use
     * @param n_reactions number of reactions in the current solution
     *
     * @return vector of all valid actions in decreasing objective score order
     */
    vector<BestAction>
    score_deletions(vector<string> &groups, vertex_subgraph curr_solution, int threads, unsigned n_reactions);

    /**
     * Score the change of the objective function for any
     * possible substitution
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param neighbours directly neighbouring nodes of the current solution
     * @param best_solution optimal solution for the current iteration
     * @param moved whether a valid move was already found in this iteration
     *
     * @return whether a valid was taken of previously found in this iteration
     */
    bool score_substitutions(vector<string> &groups, vertex_subgraph curr_solution, unordered_set<size_t> &neighbours,
                             BestAction &best_solution, bool moved);

    /**
     * Score the change of the objective function for any
     * possible substitution
     *
     * @param groups group comparison. Only the first two entries are considered
     * @param curr_solution current solution as a vertex subgraph
     * @param neighbours directly neighbouring nodes of the current solution
     * @param threads number of threads to use
     * @param n_reactions number of reactions in the current solution
     *
     * @return vector of all valid actions in decreasing objective score order
     */
    vector<BestAction>
    score_substitutions(vector<string>& groups, vertex_subgraph curr_solution, vertex_subgraph neighbours,
                        int threads, unsigned n_reactions);

    /**
     * Check if a subgraph is connected
     *
     * @param vertices vertices specifying the subgraph
     *
     * @return whether the subgraph is connected or disconnected
     */
    bool subgraph_is_connected(vertex_subgraph& vertices);

    /**
     * Find the best local action based on a given objective function
     *
     * @param groups sample groups. Only the first two entries are considered
     * @param curr_solution current best solution
     * @param n_iter current iteration count
     * @param t simulated annealing temperature in the current iteration
     *
     * @return best found action in the given iteration
     */
    BestAction find_best_action(vector<string>& groups, BestAction &curr_solution,
                                int n_iter, double t);

    /**
     * Find the best local action based on a given objective function, evaluating all options in
     * a multi-threaded fashion. This is only available when compiled with OpenMP support.
     *
     * @param groups sample groups. Only the first two entries are considered
     * @param curr_solution current best solution
     * @param n_iter current iteration count
     * @param t simulated annealing temperature in the current iteration
     * @param threads number of threads to use
     *
     * @return best found action in the given iteration
     */
    BestAction find_best_action(vector<string>& groups, BestAction curr_solution, int iter,
                                double t, int threads);

    /**
     * Apply simulated annealing to a solution. Possible changes are applied inplace.
     *
     * @param best_solution best action found in the current iteration
     * @param curr_solution optimal subnetwork as determined in the previous iteration
     * @param n_iter current iteration counter
     * @param t simulated annealing temperature
     */
    void simulated_annealing(BestAction &best_solution, vertex_subgraph curr_solution,
                             int n_iter, double t);

    template<typename T>
    static bool contains(unordered_set<T> set, T key);

    /**
     * Check whether a given solution was used before
     *
     * @param curr_solution solution as a vertex subgraph
     * @return whether the solution was used
     */
    bool solution_used(const vertex_subgraph& curr_solution) const;

};


#include "LocalSearch.tpp"


#endif //ENRICHMENT_LOCALSEARCH_HPP
