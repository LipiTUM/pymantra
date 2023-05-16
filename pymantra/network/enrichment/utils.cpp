#include <vector>
#include <tuple>
#include <set>
#include <map>
#include "utils.hpp"
#include "Exceptions.hpp"


std::vector<std::vector<double>> edgelist_to_incidence_matrix(
    const std::vector<std::tuple<int, int>> &edgeList,
    unsigned long nNodes,
    bool directed,
    const std::vector<double> &edgeWeights
) {
    std::vector<std::vector<double>> incidenceMatrix(edgeList.size(), std::vector<double>(nNodes));
    bool weights = !edgeWeights.empty();
    for (unsigned i = 0; i < edgeList.size(); i++) {
        if (weights) {
            if (directed) {
                incidenceMatrix[i][std::get<0>(edgeList[i])] = edgeWeights[i];
                incidenceMatrix[i][std::get<1>(edgeList[i])] = -edgeWeights[i];
            } else {
                incidenceMatrix[i][std::get<0>(edgeList[i])] = edgeWeights[i];
                incidenceMatrix[i][std::get<1>(edgeList[i])] = edgeWeights[i];
            }
        }
        else {
            if (directed) {
                incidenceMatrix[i][std::get<0>(edgeList[i])] = 1;
                incidenceMatrix[i][std::get<1>(edgeList[i])] = -1;
            } else {
                incidenceMatrix[i][std::get<0>(edgeList[i])] = 1;
                incidenceMatrix[i][std::get<1>(edgeList[i])] = 1;
            }
        }
    }
    return incidenceMatrix;
}


double mutualInformation(const std::vector<double> &srcData, const std::vector<double> &tgtData) {
    // TODO: https://github.com/Craigacp/MIToolbox/blob/master/src/MutualInformation.c
    return 0.0;
}


void extract_enrichment(
    const std::vector<std::vector<double>> &microbe_data,
    const std::vector<int> &microbe_map,
    unsigned long nNodes, float threshold,
    std::vector<double>& enrichmentScores,
    std::set<std::tuple<int, int>>& enrichedEdges,
    std::map<int, std::tuple<std::vector<int>, std::vector<int>>>& metaboliteNeighbours
) {
    // extracting enriched reactions
    std::vector<int> enrichedReactions;
    for (unsigned i = 0; i < nNodes; i++) {
        if (enrichmentScores[i] > threshold) {
            enrichedReactions.push_back(i);
        }
    }
    // extracting enriched reaction-metabolite
    for (int reaction : enrichedReactions) {
        for (int substrate : std::get<0>(metaboliteNeighbours[reaction])) {
            enrichedEdges.insert(std::make_tuple(reaction, substrate));
        }
        for (int product : std::get<1>(metaboliteNeighbours[reaction])) {
            enrichedEdges.insert(std::make_tuple(product, reaction));
        }
    }
}


void check_edge_nodes(EdgeType edge_type, NodeType src_type, NodeType tgt_type) {
    switch (edge_type) {
        case EdgeType::microbe_microbe:
            if ((src_type != NodeType::microbe) | (tgt_type != NodeType::microbe)) {
                throw InvalidEdgeNode("metabolite", "microbe_microbe");
            }
        case EdgeType::metabolite_metabolite:
            if ((src_type != NodeType::metabolite) | (tgt_type != NodeType::metabolite)) {
                throw InvalidEdgeNode("microbe", "metabolite_metabolite");
            }
        case EdgeType::metabolite_microbe:
            if (src_type == tgt_type) {
                string message = "Only '" + NODETYPE_TO_NAME(src_type) + "' nodes found for a metabolite_microbe edge";
                throw InvalidEdgeNode(message);
            }
		default: {
                string message = "Unknown edge type'" + NODETYPE_TO_NAME(src_type) + "'";
                throw InvalidEdgeNode(message);
			}
    }
}


EdgeType NAME_TO_EDGETYPE(string &edge_type) {
    if (edge_type == "microbe_microbe") return EdgeType::microbe_microbe;
    else if (edge_type == "metabolite_microbe") return EdgeType::metabolite_microbe;
    else if (edge_type == "metabolite_metabolite") return EdgeType::metabolite_metabolite;
    throw InvalidEdgeType(edge_type);
}

string EDGETYPE_TO_NAME(EdgeType edge_type) {
    switch (edge_type) {
        case EdgeType::microbe_microbe:
            return "microbe_microbe";
        case EdgeType::metabolite_microbe:
            return "metabolite_microbe";
        case EdgeType::metabolite_metabolite:
            return "metabolite_metabolite";
        default:
            throw InvalidEdgeType("Unknown");
    }
}

int EDGETYPE_TO_INT(EdgeType edge_type) {
    switch (edge_type) {
        case EdgeType::microbe_microbe:
            return 0;
        case EdgeType::metabolite_microbe:
            return 1;
        case EdgeType::metabolite_metabolite:
            return 2;
        default:
            throw InvalidEdgeType("Unknown");
    }
}


NodeType NAME_TO_NODETYPE(const string &node_type) {
    if (node_type == "microbe") return NodeType::microbe;
    else if(node_type == "metabolite") return NodeType::metabolite;
    else throw InvalidNodeType(node_type);
}

string NODETYPE_TO_NAME(NodeType node_type) {
    switch(node_type) {
        case NodeType::microbe:
            return "microbe";
        case NodeType::metabolite:
            return "metabolite";
        default:
            throw InvalidNodeType("Unknown");
    }
}

int NODETYPE_TO_INT(NodeType node_type) {
    switch (node_type) {
        case NodeType::microbe:
            return 0;
        case NodeType::metabolite:
            return 1;
        default:
            throw InvalidNodeType("Unknown");
    }
}


vector<NodeType> EDGE_NODE_TYPES(EdgeType edge_type) {
    switch (edge_type) {
        case EdgeType::microbe_microbe:
            return vector<NodeType>({NodeType::microbe});
        case EdgeType::metabolite_metabolite:
            return vector<NodeType>({NodeType::metabolite});
        case EdgeType::metabolite_microbe:
            return vector<NodeType>({NodeType::metabolite, NodeType::microbe});
        default:
            throw InvalidEdgeType("Unknown");
    }
}
