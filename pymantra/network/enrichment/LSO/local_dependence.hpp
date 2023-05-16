#ifndef ENRICHMENT_LOCAL_DEPENDENCE_HPP
#define ENRICHMENT_LOCAL_DEPENDENCE_HPP

#include <vector>
#include <boost/graph/adjacency_list.hpp>

using std::vector;
namespace boost {

typedef std::pair<size_t, size_t> Edge;
// => vertex property is a uint describing node type
typedef adjacency_list<mapS, vecS, directedS, size_t> DirectedGraph;

// directed subgraph (without reaction nodes?)
adjacency_list<> get_k_neighbourhood();

double conditional_mutual_information();
double conditional_dependence();

vector<double> local_dependence();

}

#endif //ENRICHMENT_LOCAL_DEPENDENCE_HPP