#include "lso_utils.hpp"


BestAction::BestAction() {
    this->score = -1.;
    this->action = initial_seed;
}

BestAction::BestAction(double score, Action action, unsigned n_reactions) {
    this->score = score;
    this->action = action;
    this->n_reactions = n_reactions;
}

BestAction::BestAction(double score, vertex_subgraph solution, Action action) {
    this->score = score;
    this->solution = std::move(solution);
    this->action = action;
}

BestAction::BestAction(double score, Action action, int v1, int v2, size_t n_reactions) {
    this->score = score;
    this->action = action;
    this->vertex_affected[0] = v1;
    this->vertex_affected[1] = v2;
    this->n_reactions = n_reactions;
}

BestAction::BestAction(double score, vertex_subgraph solution, Action action, int v1, int v2, size_t n_reactions) {
    this->score = score;
    this->solution = std::move(solution);
    this->action = action;
    this->vertex_affected[0] = v1;
    this->vertex_affected[1] = v2;
    this->n_reactions = n_reactions;
}

bool BestAction::operator==(const BestAction &a) const {
    return this->score == a.score;
}

bool action_greater(BestAction &a, BestAction &b) {
    return a.score > b.score;
}
