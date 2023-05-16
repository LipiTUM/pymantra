#include <vector>
#include <tuple>
#include <string>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include "Exceptions.hpp"

#ifndef SRC_UTILS_HPP
#define SRC_UTILS_HPP


using std::vector;
using std::string;


vector<vector<double>> edgelist_to_incidence_matrix(
    const vector<std::tuple<int, int>> &edgeList,
    unsigned long nNodes,
    bool directed = false,
    const vector<double> &edgeWeights = vector<double>(0)
);


double mutualInformation(const vector<double> &srcData, const vector<double> &tgtData);


void extract_enrichment(
    const vector<vector<double>> &microbe_data,
    const vector<int> &microbe_map,
    unsigned long nNodes, float threshold,
    vector<double>& enrichmentScores,
    boost::unordered_set<std::tuple<int, int>>& enrichedEdges,
    boost::unordered_map<int, std::tuple<vector<int>, vector<int>>>& metaboliteNeighbours
);


// Helper to use enums as hashes in maps
struct EnumHash {
    template <typename T>
    size_t operator()(T t) const
    {
        return static_cast<size_t>(t);
    }
};


enum class EdgeType {
    microbe_microbe,
    metabolite_metabolite,
    metabolite_microbe,
    // should always stay at the end
    // adapt when the number of edge types are changed
    NUMBER_OF_EDGE_TYPES=3
};
EdgeType NAME_TO_EDGETYPE(string&);
string EDGETYPE_TO_NAME(EdgeType);
int EDGETYPE_TO_INT(EdgeType);


enum NodeType {
    microbe,
    metabolite,
    // should always stay at the end
    // adapt when the number of edge types are changed
    NUMBER_OF_NODE_TYPES=2
};
NodeType NAME_TO_NODETYPE(const string&);
string NODETYPE_TO_NAME(NodeType);
int NODETYPE_TO_INT(NodeType);

vector<NodeType> EDGE_NODE_TYPES(EdgeType);

void check_edge_nodes(EdgeType, NodeType, NodeType);


#endif // SRC_UTILS_HPP
