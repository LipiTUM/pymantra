#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN


#include "../doctest.h"
#include "../RFA/ReferenceNetwork.hpp"
#include "../RFA/ReferenceNetwork.cpp"
#include "../Exceptions.hpp"
#include "../Exceptions.cpp"
#include "../utils.hpp"
#include "../utils.cpp"
#include "../statsutils.cpp"
#include "../pyutils.hpp"
#include "../pyutils.cpp"
#include "boost/unordered_set.hpp"
#include <iostream>


using boost::unordered_set;
using std::cout;
using std::endl;


bool networks_identical(unordered_set<Edge> net1, vector<Edge> &net2) {
    for (const auto& edge : net1) {
        for (const auto& vec_edge : net2) {
            if (edge.first == vec_edge.first & edge.second == vec_edge.second) {
                return false;
            }
        }
    }
    for (const auto& edge : net2) {
        if (net1.find(edge) != net1.end()) {
            return false;
        }
    }
    return true;
}


TEST_SUITE_BEGIN("ReferenceNetwork");


TEST_CASE("Network Initialization") {
    GIVEN("a vector of edges and node types") {
        vector<strEdge> edges = {
            {"1", "2"},
            {"0", "3"},
            {"1", "3"},
            {"0", "2"},
            {"0", "4"}
        };

        vector<Edge> int_edges = {
                {1, 2},
                {0, 3},
                {1, 3},
                {0, 2},
                {0, 4}
        };

        string node_types[5] = {
            "microbe", "metabolite", "metabolite",
            "microbe", "metabolite"
        };
        string edge_types[5] {
           "metabolite_metabolite", "microbe_microbe",
           "metabolite_microbe", "metabolite_microbe",
           "metabolite_microbe"
        };
        vector<vector<double>> correlations = {
            {1, .4, -.9, .1, -.3},
            {.4, 1, .6, -.4, -.7},
            {-.9, .6, 1, .85, .78},
            {.1, -.4, .85, 1, -.65},
            {-.3, -.7, .78, -.65, 1}
        };
        unordered_map<string, size_t> node_order = {
                {"0", 0},
                {"1", 1},
                {"2", 2},
                {"3", 3},
                {"4", 4},
        };
        try {
            ReferenceNetwork refnet = ReferenceNetwork(&edges[0], 5, edge_types, 5,
                                                       node_order, node_types);
            WHEN("Initializing") {
                THEN("Generate Network") {
                    CHECK_FALSE(networks_identical(refnet.network.edges, int_edges));
                }
            }
            WHEN("Computing cutoffs") {
                vector<double> best_cutoffs(2);
                refnet.compute_optimal_cutoff(correlations, 5, best_cutoffs);
                for (int i = 0; i < 2; i++) {
                    cout << "Cutoff " << i << ": " << best_cutoffs[i] << endl;
                }
                // TODO: compute cutoffs manually and check if they are computed correctly
                CHECK_EQ(best_cutoffs[0], 0);
                CHECK_EQ(best_cutoffs[1], 1);
            }
        } catch (const std::exception& e) {
            FAIL("Unexpected exception: ", e.what());
        }
    }
}

// TODO: test enrichment
