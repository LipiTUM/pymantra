#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN


#include "../doctest.h"
#include <iostream>
#include <string>
#include <vector>
#include "../LSO/objective_functions.hpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;


TEST_SUITE_BEGIN("objective functions");

TEST_CASE("objective functions from c++") {
    GIVEN("metabolite-metabolite") {
        vector<string> groups = {"x", "y"};
        unordered_set<pair<size_t, size_t>> edges = {
            pair(0, 1), pair(0, 2), pair(0, 5),
            pair(1, 2), pair(2, 3), pair(2, 5),
            pair(3, 4), pair(5, 6)
        };

        unordered_map<string, vector<double>> vec_dummy;

        vector<vertex_props> node_properties = {
            vertex_props{vertex_index_t(0), unordered_map<string, double>({pair("x", -0.42), pair("y", .15)}), vec_dummy, "a", "reaction"},
            vertex_props{vertex_index_t(1), unordered_map<string, double>({pair("x", .76), pair("y", -3.)}), vec_dummy, "b", "reaction"},
            vertex_props{vertex_index_t(2), unordered_map<string, double>({pair("x", 3.34), pair("y", -3.7)}), vec_dummy, "c", "reaction"},
            vertex_props{vertex_index_t(3), unordered_map<string, double>({pair("x", .78), pair("y", -.36)}), vec_dummy, "d", "reaction"},
            vertex_props{vertex_index_t(4), unordered_map<string, double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "e", "reaction"},
            vertex_props{vertex_index_t(5), unordered_map<string, double>({pair("x", .16), pair("y", .37)}), vec_dummy, "f", "reaction"},
            vertex_props{vertex_index_t(6), unordered_map<string, double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "g", "reaction"}
        };

        unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
            pair(pair(0, 1), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(0, 2), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(0, 5), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(1, 2), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(2, 3), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(2, 5), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(3, 4), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(5, 6), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"})
        };

        Graph graph;
        Graph::vertex_descriptor vertex_;
        int i;
        for (i = 0; i < node_properties.size(); i++) {
            vertex_ = add_vertex(graph);
            graph[vertex_].index = static_cast<vertex_index_t>(i);
            graph[vertex_].name = node_properties[i].name;
            graph[vertex_].value = node_properties[i].value;
            graph[vertex_].type = node_properties[i].type;
        }
        Graph::edge_descriptor edge_;
        for (auto& edge : edge_properties) {
            add_edge(edge.first.first, edge.first.second, edge.second, graph);
        }
        vertex_subgraph vsg = {1, 2, 3};
        double dysreg = reaction_dysregulation(groups, graph, vsg, 2);
        CHECK(dysreg - 2.8388971 < 1e-6);
    }


    GIVEN("metabolome-microbiome") {
        vector<string> groups = {"x", "y"};
        Graph graph;
        vertex_subgraph vsg;
        double dysreg = microbiome_reaction_dysregulation(groups, graph, vsg, 2);
        CHECK_EQ(dysreg, 0.0);
    }

    GIVEN("metabolome-transcriptome") {
        vector<string> groups = {"x", "y"};
        Graph graph;
        vertex_subgraph vsg;
        double dysreg = transcriptome_reaction_dysregulation(groups, graph, vsg, 2);
        CHECK_EQ(dysreg, 0.0);
    }
}
