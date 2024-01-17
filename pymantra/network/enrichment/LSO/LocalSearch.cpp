#include <iostream>
#include <limits>
#include <random>
#include <utility>
#include <random>
#include <algorithm>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/biconnected_components.hpp>
#include <boost/graph/subgraph.hpp>
#include <boost/graph/random.hpp>
#include <boost/unordered_map.hpp>
#include <boost/algorithm/cxx11/all_of.hpp>
#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/format.hpp>

#include "LocalSearch.hpp"
#include "../Exceptions.hpp"
#include "../pyutils.hpp"


#ifdef USE_OPENMP
    #include <omp.h>
#endif


using std::vector;
using std::numeric_limits;
using std::pow;
using std::exp;
using std::copy;
using std::back_inserter;
using std::max;
using boost::unordered_map;
using boost::add_vertex;
using boost::add_edge;
using boost::source;
using boost::target;
using boost::adjacency_iterator;
using boost::adjacent_vertices;
using boost::tie;
using boost::algorithm::all_of;
using boost::algorithm::any_of;
using boost::format;
using boost::str;


LocalSearch::LocalSearch(
    PyObject *node_properties, PyObject *edge_properties,
    int objective_function, double temperature, double delta_min,
    size_t l_min, size_t l_max, size_t min_reactions, size_t max_iter, float p
) : solution() {
    this->construct_params(objective_function, temperature, delta_min,
                           l_min, l_max, min_reactions, max_iter, p);
    this->graph = pygraph_to_boost(node_properties, edge_properties);
}


LocalSearch::LocalSearch(
    vector<vertex_props> &node_properties,
    unordered_map<pair<size_t, size_t>, edge_props> &edge_properties,
    int objective_function, double temperature, double delta_min, size_t l_min,
    size_t l_max, size_t min_reactions, size_t max_iter, float p
)  : solution() {
    this->construct_params(objective_function, temperature, delta_min,
                           l_min, l_max, min_reactions, max_iter, p);
    this->graph = Graph();
    Graph::vertex_descriptor vertex_;
    unsigned i;
    for (i = 0; i < node_properties.size(); i++) {
        vertex_ = add_vertex(graph);
        graph[vertex_].index = static_cast<vertex_index_t>(i);
        graph[vertex_].name = node_properties[i].name;
        graph[vertex_].value = node_properties[i].value;
        graph[vertex_].vec_value = node_properties[i].vec_value;
        graph[vertex_].type = node_properties[i].type;
        graph[vertex_].objective = 0.;
    }
    Graph::edge_descriptor edge_;
    for (auto& edge : edge_properties) {
        add_edge(edge.first.first, edge.first.second, edge.second, graph);
    }
}


void LocalSearch::construct_params(
        int objective_function_, double temperature_, double delta_min_, size_t l_min_,
        size_t l_max_, size_t min_reactions_, size_t max_iter_, float p
) {
    switch (objective_function_) {
        case 0:
            this->objective_function = reaction_dysregulation;
            this->precompute_objectives = true;
            this->precompute = set_reaction_kcl_objectives;
            this->multi_omics = false;
            break;
        case 1:
            this->objective_function = reaction_dysregulation;
            this->precompute_objectives = true;
            this->precompute = set_reaction_ld_objectives;
            this->multi_omics = false;
            break;
        case 2:
            this->objective_function = microbiome_reaction_dysregulation;
            this->precompute_objectives = true;
            this->precompute = set_microbiome_reaction_dysregulation_objectives;
            this->multi_omics = true;
            break;
        case 3:
            this->objective_function = transcriptome_reaction_dysregulation;
            this->precompute_objectives = true;
            this->precompute = set_transcriptome_reaction_dysregulation_objectives;
            this->multi_omics = true;
            break;
        case 4:
            // TODO: option to pass multi-omics here
            this->objective_function = reaction_dysregulation;
            this->precompute_objectives = false;
            this-> precompute = precomputed_objectives;
            this->multi_omics = false;
            break;
        default:
            throw InvalidObjectiveFunction(
                    "Only the following objective functions are allowed:\n"
                    "0: metabolic reactions only\n"
                    "1: metabolite-microbiome dysregulation\n"
                    "2: metabolite-transcriptome dysregulation (currently not implemented)"
            );
    }
    this->temp = temperature_;
    if (l_min_ < l_max_) {
        this->l_min = l_min_;
        this->l_max = l_max_;
    } else{
        // TODO: add warning
        this->l_max = l_min_;
        this->l_min = l_max_;
    }
    // min_reactions_ > l_min_ this wouldn't make any sense
    if (min_reactions_ > l_min_) {
        // TODO: add warning
        this->min_reactions = l_min_;
    } else {
        this->min_reactions = min_reactions_;
    }
    this->max_iter = max_iter_;
    this->delta_min = delta_min_;
    this->converged = false;
    this->p_norm = p;
}


template<typename T>
bool LocalSearch::contains(unordered_set<T> set, T key) {
    return set.find(key) != set.end();
}


bool LocalSearch::solution_used(const vertex_subgraph& curr_solution) const {
    return any_of(
        solutions, [curr_solution](const vertex_subgraph& solution_) {
            // return subgraphs_equal(solution_, curr_solution);
            return solution_ == curr_solution;
        }
    );
}


Graph LocalSearch::subgraph_from_vertices(vertex_subgraph &vertices) {
    Graph g;
    unordered_map<size_t, size_t> v_map;
    size_t i = 0;
    for (auto v : vertices) {
        v_map.insert(pair<size_t, size_t>(v, i));
        i++;
    }
    // get all edges of the subgraph
    size_t vi, ni;
    for (auto v : vertices) {
        vi = v_map[v];
        for (auto neighbour : boost::make_iterator_range(boost::adjacent_vertices(v, graph))) {
            ni = v_map[neighbour];
            if (contains(vertices, neighbour)) {
                if (vi < ni) {
                    // NOTE: we only work with an undirected graph, hence
                    //       should avoid parallel edges => produces un
                    //       expected results on directed graphs!
                    boost::add_edge(vi, ni, g);
                }
            }
        }
    }
    // Graph g{edge_arr.begin(), edge_arr.end(), vertices.size()};
    return g;
}


unordered_set<size_t> LocalSearch::find_articulation_points(vertex_subgraph& vertices) {
    unordered_set<size_t> articulation_points;
    // TODO: is it more efficient to have an own implementation e.g.:
    // based on https://www.geeksforgeeks.org/articulation-points-or-cut-vertices-in-a-graph/
    Graph subgraph = subgraph_from_vertices(vertices);
    vector<size_t> art_point_vec;
    boost::articulation_points(subgraph, back_inserter(art_point_vec));
    unordered_set<size_t> art_points;
    for (auto art_point : art_point_vec) {
        art_points.insert(art_point);
    }
    return art_points;
}


void LocalSearch::set_seed(vertex_subgraph seed_) {
    seed = std::move(seed_);
}


vertex_subgraph LocalSearch::random_walk_subgraph(vertex_index_t vertex, size_t subgraph_size) {
    // random neighbour setup
    AdjacencyIter neighbour, neighbour_end;
    // auxiliary variables
    vertex_index_t current = vertex;
    size_t n_reactions = 0, n_non_reactions = 0, max_non_reactions, target_reactions;

    if (subgraph_size < min_reactions) {
        throw std::runtime_error(
            "seed size cannot be smaller than minimum number of reactions in the solution"
        );
    }

    max_non_reactions = subgraph_size - min_reactions;
    if (multi_omics) target_reactions = min_reactions;
    else target_reactions = subgraph_size;

    // subgraph vertex set
    vertex_subgraph subgraph;
	subgraph.insert(vertex);
    if (graph[vertex].type == "reaction") ++n_reactions;
    else ++n_non_reactions;

    bool iter_success;
    while (subgraph.size() < subgraph_size) {
        iter_success = false;

        // actually draw the neighbours
        for (tie(neighbour, neighbour_end) = adjacent_vertices(current, graph); neighbour != neighbour_end; ++neighbour) {
            if (subgraph.find(graph[*neighbour].index) != subgraph.end())
                continue;

            if (graph[*neighbour].type == "reaction") {
                if (n_reactions < target_reactions) {
                    current = graph[*neighbour].index;
                    subgraph.insert(current);
                    iter_success = true;
                    ++n_reactions;
                    break;
                }
            } else if (n_non_reactions < max_non_reactions) {
                current = graph[*neighbour].index;
                subgraph.insert(current);
                iter_success = true;
                ++n_non_reactions;
                break;
            }
        }
        if (!iter_success) {
            return subgraph;
        }
    }
    return subgraph;
}


void LocalSearch::set_seed(size_t seed_size) {
    vector<size_t> seeds(boost::num_vertices(graph));
    for (unsigned i = 0; i < boost::num_vertices(graph); ++i) {
        seeds[i] = i;
    }
    std::shuffle(seeds.begin(), seeds.end(), std::mt19937(std::random_device()()));

    // TODO: choose the biggest connected component
    for (auto seed_ : seeds) {
        seed = random_walk_subgraph((vertex_index_t) seed_, seed_size);
        if (seed.size() >= seed_size)
            break;
    }
    if (seed.size() < seed_size) {
        throw std::runtime_error(str(format("No connected subgraph of size %1% found!") % seed_size));
    }
}


void LocalSearch::set_seed_py(PyObject* seed_init) {
    vertex_subgraph seed_;
    size_t vertex;
    if (PyList_Check(seed_init) | PySet_Check(seed_init)) {
        PyObject *iter = PyObject_GetIter(seed_init);
        if (iter == nullptr) {
            throw PyObjectType("empty seed not allowed");
        }
        PyObject *pyvertex = PyIter_Next(iter);
        while (pyvertex != nullptr) {
            vertex = PyLong_AsLong(pyvertex);
            if (PyErr_Occurred()) {
                PyErr_Clear();
                throw PyObjectType("All objects in the seed must be python integers");
            }
            seed_.insert(vertex);
        }
    } else {
        throw PyObjectType("`seed_init` must a list or a set");
    }
    seed = seed_;
}


void LocalSearch::set_seed_py(const string& seed_init, size_t seed_size) {
    // string seed_node = PyBytes_AsString(seed_init);
    // if (PyErr_Occurred()) {
    //     PyErr_Clear();
    //     throw PyObjectType("seed_init must be a python string in utf-8 encoding");
    // }
    // we initialize a seed in case the given seed is not found in the graph
    static std::mt19937 rng { std::random_device{}() };
    boost::graph_traits<Graph>::vertex_descriptor seed_ = boost::random_vertex(graph, rng);
    vertex_index_t seed_idx = graph[seed_].index;
    for (auto vertex : boost::make_iterator_range(boost::vertices(graph))) {
        if (graph[vertex].name == seed_init) seed_idx = graph[vertex].index;
    }
    if (seed_idx == graph[seed_].index) {
        std::cout << "Starting seed node '" << seed_init << "' not found. ";
        std::cout << "Continuing with a random node" << std::endl;
    }
    // we initialize every node with an index property => warning can be ignored
    seed = random_walk_subgraph(seed_idx, seed_size);
}


bool LocalSearch::subgraph_is_connected(vertex_subgraph& vertices) {
    Graph subgraph = subgraph_from_vertices(vertices);
    // get all vertices of the subgraph
    vector<int> vertex_vector(vertices.size());
    int i = 0;
    for (auto vertex : vertices) {
        vertex_vector[i] = (int) vertex;
        i++;
    }
    auto iter_map = boost::make_iterator_property_map(
        vertex_vector.begin(), boost::get(boost::vertex_index, subgraph),
        vertex_vector[0]
    );
    // find the number of connected components in the subgraph
    int n_components = boost::connected_components(subgraph, iter_map);
    return n_components == 1;
}


unordered_set<size_t> LocalSearch::one_hop_neighbours(vertex_subgraph& curr_solution) {
    unordered_set<size_t> neighbours;
    std::pair<AdjacencyIter, AdjacencyIter> v_neighbours;
    AdjacencyIter neighbour, neighbour_end;
    for(auto vertex : curr_solution) {
        v_neighbours = adjacent_vertices(vertex, graph);
        for (tie(neighbour, neighbour_end) = v_neighbours; neighbour != neighbour_end; ++neighbour) {
            if (curr_solution.find(*neighbour) == curr_solution.end()) {
                neighbours.insert(*neighbour);
            }
        }
    }
    return neighbours;
}


void LocalSearch::simulated_annealing(
        BestAction &best_solution, vertex_subgraph curr_solution, int n_iter, double t
) {
    // In the first iteration we ignore SA
    if (n_iter > 0) {
        // difference between current and last score
        double delta = best_solution.score - score_progression[n_iter];
        // if score did not improve apply simulated annealing
        if (delta < 0) {
            // compute test value
            double test_value;
            if (t != 0 && delta != 0)
                // NOTE: originally this is -delta, but we compute delta as current - previous
                // instead of previous - current
                test_value = exp(delta / t);
            else
                test_value = 0.;

            // draw at random from a uniform distribution over [0, 1)
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_real_distribution<double> dist(0., 1.);
            double rand_val = dist(gen);

            if (test_value < rand_val) {
                // TODO: why does overriding not work properly here?
                // random test failed => discarding the un-improved solution
                best_solution.action = Action::no_move;
                best_solution.score = score_progression[n_iter];
                // can this assignment potentially blow up?
                best_solution.solution = std::move(curr_solution);
            } else {
                // accept the solution, but mark it as a simulated annealing solution
                best_solution.action = Action::simulated_annealing;
            }
        }
    }
}


/*
 * ============================
 * single-thread implementation
 * ============================
 */
void LocalSearch::update_solution(
    BestAction& prev, double score, Action action,
    vertex_subgraph solution_, int vertex_1, int vertex_2,
    size_t n_reactions
) {
    prev.score = score;
    prev.action = action;
    // TODO: add node i to curr_solution to yield solution
    //       make sure the reference is not changed!
    prev.solution = std::move(solution_);
    prev.vertex_affected[0] = vertex_1;
    prev.vertex_affected[1] = vertex_2;
    prev.n_reactions = n_reactions;
}


bool LocalSearch::score_insertions(
    vector<string> &groups, vertex_subgraph curr_solution, unordered_set<size_t> &neighbours,
    BestAction &best_solution, bool moved
) {
    double tmp_score;
    for (auto neighbour : neighbours) {
        // do insertion
        curr_solution.insert(neighbour);
        // score action and update if increased
        if (!solution_used(curr_solution)) {
            tmp_score = objective_function(groups, graph, curr_solution, p_norm);
            if (tmp_score > best_solution.score) {
                update_solution(
                    best_solution, tmp_score, Action::insertion,
                    curr_solution, (int) neighbour, -1,
                    best_solution.n_reactions + (unsigned)(graph[neighbour].type == "reaction")
                );
                moved = true;
            }
        }
        // restore initial solution
        curr_solution.erase(neighbour);
    }
    return moved;
}


bool LocalSearch::score_deletions(
    vector<string> &groups, vertex_subgraph curr_solution, BestAction &best_solution,
    bool moved
) {
    unordered_set<size_t> articulation_points = find_articulation_points(curr_solution);
    double tmp_score;
    vector<size_t> vertices;
    unsigned n_reactions;
    vertices.insert(vertices.end(), curr_solution.begin(), curr_solution.end());
    for (auto v : vertices) {
        n_reactions = best_solution.n_reactions - (unsigned)(graph[v].type == "reaction");
        // skip to next iteration if vertex is an articulation point
        if (!contains(articulation_points, (size_t)v) && (n_reactions >= min_reactions)) {
            // do deletion
            curr_solution.erase(v);
            if (!solution_used(curr_solution)) {
                // score action and update if increased
                tmp_score = objective_function(groups, graph, curr_solution, p_norm);
                if (tmp_score > best_solution.score) {
                    update_solution(
                        best_solution, tmp_score, Action::deletion,
                        curr_solution, (int) v, -1, n_reactions
                    );
                    moved = true;
                }
            }
            // restore initial solution
            curr_solution.insert(v);
        }
    }
    return moved;
}


bool LocalSearch::score_substitutions(
    vector<string> &groups, vertex_subgraph curr_solution,
    unordered_set<size_t> &neighbours, BestAction &best_solution, bool moved
) {
    double tmp_score;
    vector<size_t> vertices;
    unsigned add_int, rem_int, n_reactions;
    vertices.insert(vertices.end(), curr_solution.begin(), curr_solution.end());
    // scored are all pairwise combinations of insertion and deletion
    for (auto v : vertices) {
        for (auto neighbour : neighbours) {
            add_int = (unsigned)(graph[neighbour].type == "reaction");
            rem_int = (unsigned)(graph[v].type == "reaction");
            n_reactions = best_solution.n_reactions + add_int - rem_int;
            // do substitution
            curr_solution.insert(neighbour);
            curr_solution.erase(v);
            // score action if resulting graph is connected
            // => this cannot be checked easily beforehand
            //    in contrast to deletion only action
            bool connected, seen;
            connected = subgraph_is_connected(curr_solution);
            seen = solution_used(curr_solution);

            // if (subgraph_is_connected(curr_solution) && !(solution_used(curr_solution)) &&
            //         ((unsigned)n_reactions >= min_reactions)) {
            if (connected && !seen && n_reactions >= min_reactions) {
                tmp_score = objective_function(groups, graph, curr_solution, p_norm);
                if (tmp_score > best_solution.score) {
                    update_solution(
                        best_solution, tmp_score, Action::substitution,
                        curr_solution, (int) neighbour, (int) v, n_reactions
                    );
                    moved = true;
                }
            }
            // restore initial solution
            curr_solution.erase(neighbour);
            curr_solution.insert(v);
        }
    }
    return moved;
}


BestAction LocalSearch::find_best_action(vector<string>& groups, BestAction &curr_solution,
                                         int n_iter, double t) {
    BestAction best_solution;
    if (n_iter == 0)
        best_solution = BestAction(
            0, Action::initial_seed, curr_solution.n_reactions);
    else
        best_solution = BestAction(
            0, Action::no_move, curr_solution.n_reactions);

    unordered_set<size_t> neighbours = one_hop_neighbours(curr_solution.solution);
    bool moved = false;
    // compute all possible actions and update best_solution if
    // a score better than the previous one is found
    // TODO: add in min_reaction criterion for multi-omics
    if (curr_solution.solution.size() < l_max)
        moved = score_insertions(groups, curr_solution.solution, neighbours, best_solution, moved);
    if (curr_solution.solution.size() > l_min)
        moved = score_deletions(groups, curr_solution.solution, best_solution, moved);
    // substitutions do not change the size of the solution
    moved = score_substitutions(groups, curr_solution.solution, neighbours, best_solution, moved);

    if (!moved) {
        // this means no new solution found
        return best_solution;
    }

    // apply simulated annealing
    simulated_annealing(best_solution, curr_solution.solution, n_iter, t);

    return best_solution;
}


/*
 * ===========================
 * parallelized implementation
 * ===========================
 */
#ifdef USE_OPENMP

    vector<BestAction>
    LocalSearch::score_insertions(vector<string> &groups, vertex_subgraph &curr_solution,
                                  unordered_set<size_t> &neighbours, int threads,
                                  unsigned n_reactions) {
        // no neighbours => return
        if (neighbours.empty()) return{};

        vector<BestAction> actions;
        // shared variables
        vector<size_t> neighbour_vec(neighbours.begin(), neighbours.end());
        // private variables
        vertex_subgraph priv_solution = curr_solution;
        double score;
        BestAction action;

        #pragma omp parallel num_threads(threads) private(priv_solution, score, action)
        {

            #pragma omp for nowait
            for (size_t &i: neighbour_vec) {
                // copy shared object to private
                priv_solution = curr_solution;
                priv_solution.insert(i);

                if (!solution_used(priv_solution)) {

                    score = objective_function(groups, graph, priv_solution, p_norm);
                    action = BestAction(
                        score, priv_solution, Action::insertion, (int) i, -1,
                        n_reactions + (int)(graph[i].type == "reaction")
                    );

                    #pragma omp critical
                    actions.push_back(action);
                }
            }
        }

        std::sort(actions.begin(), actions.end(), action_greater);
        return actions;
    }


    vector<BestAction>
    LocalSearch::score_deletions(vector<string> &groups, vertex_subgraph curr_solution, int threads,
                                 unsigned n_reactions) {
        vector<BestAction> actions;
        // shared variables
        vector<size_t> solution_vec;
        solution_vec.insert(solution_vec.end(), curr_solution.begin(), curr_solution.end());
        unordered_set<size_t> art_points = find_articulation_points(curr_solution);
        // private variables
        vertex_subgraph priv_solution;
        double score;
        BestAction action;
        unsigned is_reaction;

        #pragma omp parallel num_threads(threads) private(priv_solution, score, action, is_reaction)
        {

            #pragma omp for nowait
            for (auto &i: solution_vec) {
                // copy shared object to private
                priv_solution = curr_solution;

                if (!contains(art_points, i)) {
                    is_reaction = (unsigned)(graph[i].type == "reaction");

                    if ((n_reactions - is_reaction) >= min_reactions) {
                        priv_solution.erase(i);

                        if (!solution_used(priv_solution)) {
                            score = objective_function(groups, graph, priv_solution, p_norm);
                            action = BestAction(
                                score, priv_solution, Action::deletion, (int) i, -1,
                                n_reactions - is_reaction
                            );

                            #pragma omp critical
                            actions.push_back(action);
                        }
                    }
                }
            }
        }
        std::sort(actions.begin(), actions.end(), action_greater);
        return actions;
    }


    vector<BestAction>
    LocalSearch::score_substitutions(
        vector<string>& groups, vertex_subgraph curr_solution, vertex_subgraph neighbours,
        int threads, unsigned n_reactions
    ) {
        // no neighbours => no substitutions
        if (neighbours.empty()) return {};

        vector<BestAction> actions;

        // initializing shared variables
        vector<size_t> neighbour_vec(neighbours.begin(), neighbours.end());
        vector<size_t> vertices(curr_solution.begin(), curr_solution.end());

        // initializing private variables
        vertex_subgraph priv_solution;
        double score;
        BestAction action;
        int n_reaction_change, add_int, rem_int;

        #pragma omp parallel num_threads(threads) private(priv_solution, score, action, n_reaction_change, add_int, rem_int)
        {
            #pragma omp for nowait
            for(size_t &i : vertices) {
                for (size_t &j : neighbour_vec) {
                    // copy shared solution to private
                    priv_solution = curr_solution;

                    add_int = (int)(graph[j].type == "reaction");
                    rem_int = (int)(graph[i].type == "reaction");
                    n_reaction_change = add_int - rem_int;

                    priv_solution.insert(j);
                    priv_solution.erase(i);

                    if ((n_reactions + n_reaction_change) >= min_reactions) {
                        if (subgraph_is_connected(priv_solution)) {
                            if (!solution_used(priv_solution)) {

                                score = objective_function(groups, graph, priv_solution, p_norm);
                                action = BestAction(
                                    score, priv_solution, Action::substitution,
                                    (int) i, (int) j,
                                    n_reactions + n_reaction_change
                                );

                                #pragma omp critical
                                actions.push_back(action);
                            }
                        }
                    }
                }
            }
        }

        std::sort(actions.begin(), actions.end(), action_greater);
        return actions;
    }


    BestAction
    LocalSearch::find_best_action(vector<string>& groups, BestAction curr_solution,
                                  int iter, double t, int threads) {
        vertex_subgraph neighbours = one_hop_neighbours(curr_solution.solution);

        // compute all possible actions and update best_solution if
        // a score better than the previous one is found
        vector<BestAction> insertions, deletions, substitutions;

        if (curr_solution.solution.size() < l_max)
            insertions = score_insertions(
                    groups, curr_solution.solution, neighbours, threads, curr_solution.n_reactions);
        if (curr_solution.solution.size() > l_min)
            deletions = score_deletions(
                    groups, curr_solution.solution, threads, curr_solution.n_reactions);
        // substitutions do not change the size of the solution
        substitutions = score_substitutions(
                groups, curr_solution.solution, neighbours, threads, curr_solution.n_reactions);
        // merging all actions sorted by score (decreasing)
        vector<BestAction> actions = merge_solutions(
                insertions, deletions, substitutions, action_greater);

        // collecting best action information
        BestAction best_solution(-1., Action::no_move, curr_solution.n_reactions);

        bool moved = false;
        unsigned i;
        for (i = 0; i < actions.size(); i++) {
            // this means we have found the best unused solution
            // apply_action(curr_solution.solution, best_solution, actions[i]);
            best_solution = std::move(actions[i]);
            // NOTE: this shouldn't be necessary more of a double check
            if (!solution_used(best_solution.solution)) {
                moved = true;
                break;
            }
        }

        if (!moved)
            // no unused solution found
            return best_solution;

        simulated_annealing(best_solution, curr_solution.solution, iter, t);

        return best_solution;
    }
#endif


void LocalSearch::run_local_search(vector<string> &groups, int n_threads) {
    // computing objectives and adding to the graph
    if (precompute_objectives && ((objectives_computed[0] != groups[0]) && (objectives_computed[1] != groups[1]))) {
        precompute(groups, graph);
        objectives_computed[0] = groups[0];
        objectives_computed[1] = groups[1];
    }


    // re-setting logs
    score_progression.clear();
    solutions.clear();
    converged = false;

    // generating a random seed if not generated before
    // NOTE: seed is NOT overwritten! If you want to have a different seed
    //       for each run you need to set the seed manually in between runs
    if (seed.empty()) set_seed((int) l_min);

    // computing seed values
    double seed_score = this->score_solution(groups, seed);
    this->score_progression.push_back(seed_score);
    BestAction current_solution(seed_score, Action::initial_seed, 0);
    current_solution.solution = seed;
    for (auto v : current_solution.solution) {
        if (graph[boost::vertex(vertex_index_t(v), graph)].type == "reaction") {
            ++current_solution.n_reactions;
        }
    }
    solution = current_solution;
    solutions.push_back(current_solution.solution);

    // starting local search
    converged = true;
    double t0 = temp;
    double t = t0;
    int i;
    for (i = 0; i < (int) max_iter; ++i) {
#ifdef USE_OPENMP
        // TODO: is there a way to avoid if-else in every iteration? => lambda outside of loop
        // find update
        if (n_threads == 1)
            current_solution = find_best_action(groups, current_solution, i, t);
        else
            current_solution = find_best_action(groups, current_solution, i, t, n_threads);
#else
        current_solution = find_best_action(groups, current_solution, i, t);
#endif
        if (current_solution.action == Action::no_move &&
                action_progression[max(0, i - 1)] == Action::no_move)
        {
            std::cout << "No more local solutions to test in iteration " << i;
            std::cout << " or too cold for simulated annealing. Returning..." << std::endl;
            break;
        }

        // log progression of optimal score at each iteration
        score_progression.push_back(current_solution.score);
        // log solution subnetwork
        solutions.push_back(current_solution.solution);
        // log solution action
        action_progression.push_back(current_solution.action);
        // update temperature
        t = t0 * (pow(.9, (float) i));
    }

    if ((size_t)i == max_iter)
        converged = false;

    if (score_progression.empty() || solutions.empty()) {
        throw std::runtime_error(
                "No solutions matching the size criteria found!");
    }
    if (score_progression.size() != solutions.size()) {
        throw std::runtime_error(
                "Internal error: number of scores does not match number of solutions");
    }
    vector<double>::iterator best_score;
    best_score = std::max_element(score_progression.begin(), score_progression.end());
    size_t best_idx = std::distance(score_progression.begin(), best_score);

    if (best_idx >= score_progression.size()) {
        throw std::runtime_error("'best_idx' search failed");
    }
    // Setting the best score found as the final solution
    solution = BestAction(score_progression[best_idx], solutions[best_idx], Action::initial_seed);
}


double LocalSearch::score_solution(vector<string>& groups, vertex_subgraph& solution_) {
    return objective_function(groups, graph, solution_, p_norm);
}


double LocalSearch::score_solution(vector<string>& groups) {
    if (solution.score == -1.) {
        throw MissingComputation(
                "No enrichment computed yet. Please call 'run_local_search'."
        );
    }
    // if (!converged) {
    //     throw ConvergenceError(
    //             "Local search did not converge! Please re-run with different seed and/or parameters."
    //     );
    // }
    return objective_function(groups, graph, solution.solution, p_norm);
}

/*  ==================
 *  getters for python
 *  ==================
 */
bool LocalSearch::get_converged() const {
    return converged;
}

PyObject *LocalSearch::get_best_solution() {
    PyObject *vertex_set = PySet_New(nullptr);
    Graph::vertex_descriptor vertex;
    for (size_t vertex_idx : solution.solution) {
        vertex = boost::vertex(vertex_idx, graph);
        PySet_Add(vertex_set, PyUnicode_FromString(graph[vertex].name.c_str()));
    }
    return vertex_set;
}

vector<double> LocalSearch::get_score_progression() const {
    return score_progression;
}

double LocalSearch::get_best_score() const {
    return solution.score;
}

size_t LocalSearch::get_lmax() const {
    return l_max;
}

size_t LocalSearch::get_lmin() const {
    return l_min;
}

size_t LocalSearch::get_min_reactions() const {
    return this->min_reactions;
}

double LocalSearch::get_temp() const {
    return temp;
}

size_t LocalSearch::get_maxiter() const {
    return max_iter;
}

double LocalSearch::get_deltamin() const {
    return delta_min;
}

double LocalSearch::get_lp_norm() const {
    return this->p_norm;
}


/*  ==================
 *  setters for python
 *  ==================
 */
void LocalSearch::set_lmax(size_t lmax) {
    this->l_max = lmax;
}

void LocalSearch::set_lmin(size_t lmin) {
    this->l_min = lmin;
}

void LocalSearch::set_min_reactions(size_t min_reactions_) {
    this->min_reactions = min_reactions_;
}

void LocalSearch::set_temp(double temp_) {
    this->temp = temp_;
}

void LocalSearch::set_maxiter(size_t iter) {
    this->max_iter = iter;
}

void LocalSearch::set_deltamin(double delta) {
    this->delta_min = delta;
}

void LocalSearch::set_lp_norm(float p) {
    this->p_norm = p;
}


void LocalSearch::precompute_objectives_py(vector<string>& groups) {
    precompute(groups, graph);
}

PyObject *LocalSearch::get_objective_values() {
    if (precompute_objectives && objectives_computed[0].empty() && objectives_computed[1].empty()) {
        throw std::invalid_argument(
            "No values computed yet, please call 'precompute' or 'run_local_search' first!"
        );
    }
    PyObject *objectives = PyDict_New(), *vertices = PyDict_New();
    Graph::vertex_iterator vertex, v_end;
    for (tie(vertex, v_end) = boost::vertices(graph); vertex != v_end; ++vertex) {
        if (graph[*vertex].type == "reaction") {
            PyDict_SetItem(vertices,
                           PyUnicode_FromString(graph[*vertex].name.c_str()),
                           PyFloat_FromDouble(graph[*vertex].objective));
        }
    }
    PyDict_SetItem(objectives, PyUnicode_FromString("nodes"), vertices);
    // only in this case edge_values are required
    if (this->objective_function == microbiome_reaction_dysregulation) {
        PyObject *edges = PyDict_New();
        Graph::edge_iterator edge, e_end;
        string src, tgt;
        for (tie(edge, e_end) = boost::edges(graph); edge != e_end; ++edge) {
            if (graph[*edge].type == "REACTION_ORGANISM") {
                PyObject *pyedge = PyTuple_New(2);
                src = graph[source(*edge, graph)].name;
                tgt = graph[target(*edge, graph)].name;
                PyTuple_SetItem(pyedge, 0, PyUnicode_FromString(src.c_str()));
                PyTuple_SetItem(pyedge, 1, PyUnicode_FromString(tgt.c_str()));
                PyDict_SetItem(edges, pyedge, PyFloat_FromDouble(graph[*edge].objective));
            }
        }
        PyDict_SetItem(objectives, PyUnicode_FromString("edges"), edges);
    }
    else if (this->objective_function != reaction_dysregulation) {
        throw std::runtime_error("Currently only single-omics and metabolome-microbiome analyses are supported");
    }
    return objectives;
}
