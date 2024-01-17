#include "lso_utils.hpp"

#include <utility>


BestAction::BestAction()
	: score(-1), action(initial_seed)
{}

BestAction::BestAction(double score, Action action, unsigned n_reactions)
	: score(score), action(action), n_reactions(n_reactions)
{}

BestAction::BestAction(double score, vertex_subgraph solution, Action action)
	: score(score), solution(std::move(solution)), action(action)
{}

BestAction::BestAction(double score, Action action, int v1, int v2, size_t n_reactions)
	: score(score), action(action), n_reactions(n_reactions)
{
    this->vertex_affected[0] = v1;
    this->vertex_affected[1] = v2;
}

BestAction::BestAction(double score, vertex_subgraph solution, Action action, int v1, int v2, size_t n_reactions)
	: score(score), solution(std::move(solution)), action(action), n_reactions(n_reactions)
{
    this->vertex_affected[0] = v1;
    this->vertex_affected[1] = v2;
}

bool BestAction::operator==(const BestAction &a) const {
    return this->score == a.score;
}

bool action_greater(BestAction &a, BestAction &b) {
    return a.score > b.score;
}
