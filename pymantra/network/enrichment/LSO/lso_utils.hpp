#ifndef ENRICHMENT_LSO_UTILS_HPP
#define ENRICHMENT_LSO_UTILS_HPP


#include <vector>
#include "reaction_graph.hpp"
#include "objective_functions.hpp"


using std::vector;


// objective function template
// NOTE: the objective functions can be redefined as long as they only the
//       following vertex properties
//          * node type ()
//          * node value (vector<double> => one value per group or sample)
//       and edge properties
//          * edge type ()
//          * edge value (vector<double> => one value per group or sample)
typedef double (* objective_)(vector<string>&, Graph&, vertex_subgraph&, float p);

typedef void (* objective_precompute)(vector<string>&, Graph&);

enum Action {
    insertion,
    deletion,
    substitution,
    initial_seed,
    no_move,
    simulated_annealing
};


struct BestAction {
    double score;
    vertex_subgraph solution{};
    Action action;
    size_t n_reactions{};
    int vertex_affected[2] = {-1, -1};

    BestAction();
    BestAction(double score, Action action, unsigned n_reactions);
    BestAction(double score, vertex_subgraph solution, Action action);
    BestAction(double score, Action action, int v1, int v2, size_t n_reactions);
    BestAction(double score, vertex_subgraph solution, Action action, int v1, int v2, size_t n_reactions);

    bool operator==(const BestAction& a) const;
};



bool action_greater(BestAction &a, BestAction &b);


AdjacencyIter draw_random_neighbour();


// inline bool operator==(const BestAction& a, BestAction &b) { return a.score == b.score; }

typedef bool (* comp_double)(double, double);
typedef bool (* comp_action)(BestAction, BestAction);


#endif //ENRICHMENT_LSO_UTILS_HPP
