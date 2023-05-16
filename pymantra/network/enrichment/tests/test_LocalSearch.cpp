#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN


#include "../doctest.h"
#include <iostream>
#include <string>
#include "../LSO/LocalSearch.hpp"


using std::cout;
using std::endl;
using std::string;


template<typename T>
bool contains(unordered_set<T> &set, T key) {
    return set.find(key) != set.end();
}


bool solution_equals(BestAction& computed, vertex_subgraph& expected) {
    vector<size_t> fp{};
    vector<size_t> fn{};
    for (auto vertex : computed.solution) {
        cout << vertex << " -- ";
        if (!contains(expected, vertex)) {
            fp.push_back(vertex);
        }
    }
    cout << endl;
    for (auto vertex : expected) {
        if (!contains(computed.solution, vertex)) {
            fn.push_back(vertex);
        }
        cout << vertex << " -- ";
    }
    cout << endl;
    if (fp.empty() &fn.empty()) return true;
    if (!fn.empty()) {
        cout << "Expected vertices missing" << endl;
    }
    if (!fn.empty()) {
        cout << "Unexpected vertices found" << endl;
    }
    return false;
}


void report_score_progression(LocalSearch& lso) {
    if (lso.score_progression.empty()) {
        cout << R"(\033[1;31m Oops, empty score progression! \033[0m)" << endl;
    }
    cout << "Score progression with maximum score " << lso.get_best_score() << ":" << endl << "-- ";
    for (auto score : lso.get_score_progression()) {
        cout << score << " --";
    }
    cout << endl;
}


void report_solution(LocalSearch& lso) {
    for (auto vertex : lso.solution.solution) {
        cout << vertex << " -- ";
    }
    cout << endl;
}


/*
 *  NOTE: This suite only tests the C++ constructor.
 *        The python C API constructor is tested with the bindings.
 */
TEST_SUITE_BEGIN("LocalSearch class");


TEST_CASE("metabolites only") {

    Py_Initialize();

    unordered_set<pair<size_t, size_t>> edges = {
            pair(0, 1), pair(0, 2), pair(0, 5),
            pair(1, 2), pair(2, 3), pair(2, 5),
            pair(3, 4), pair(5, 6)
    };

    unordered_map<string, vector<double>> vec_dummy;
    vector<vertex_props> node_properties = {
            vertex_props{vertex_index_t(0), unordered_map<string,
                    double>({pair("x", -0.42), pair("y", .15)}), vec_dummy, "a", "reaction"},
            vertex_props{vertex_index_t(1), unordered_map<string,
                    double>({pair("x", .76), pair("y", -3.)}), vec_dummy, "b", "reaction"},
            vertex_props{vertex_index_t(2), unordered_map<string,
                    double>({pair("x", 3.34), pair("y", -3.7)}), vec_dummy, "c", "reaction"},
            vertex_props{vertex_index_t(3), unordered_map<string,
                    double>({pair("x", .78), pair("y", -.36)}), vec_dummy, "d", "reaction"},
            vertex_props{vertex_index_t(4), unordered_map<string,
                    double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "e", "reaction"},
            vertex_props{vertex_index_t(5), unordered_map<string,
                    double>({pair("x", .16), pair("y", .37)}), vec_dummy, "f", "reaction"},
            vertex_props{vertex_index_t(6), unordered_map<string,
                    double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "g", "reaction"}
    };

    unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
            pair(pair(0, 1), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(0, 2), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(0, 5), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(1, 2), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(2, 3), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(2, 5), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(3, 4), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
            pair(pair(5, 6), edge_props{unordered_map<string,
                    double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"})
    };

    int objective_function = 0;
    double temperature = 1.;
    double delta_min = 1e-10;
    size_t l_min = 3;
    size_t l_max = 5;
    size_t max_iter = 25;

    LocalSearch lso = LocalSearch(
        node_properties, edge_properties, objective_function, temperature, delta_min,
        l_min, l_max, 3, max_iter, 1.
    );

    vertex_subgraph expected_solution = {1, 2, 3};
    vector<double> progression;

    GIVEN("single-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 1);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }
    GIVEN("multi-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 4);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }

    Py_Finalize();
}


TEST_CASE("microbiome") {

    Py_Initialize();

    unordered_set<pair<size_t, size_t>> edges = {
        pair(0, 1), pair(0, 2), pair(0, 5),
        pair(1, 2), pair(2, 3), pair(2, 5),
        pair(3, 4), pair(5, 6),
        pair(2, 7), pair(2, 9), pair(3, 7), pair(5, 8)
    };

    unordered_map<string, vector<double>> vec_dummy;
    vector<vertex_props> node_properties = {
        vertex_props{vertex_index_t(0), unordered_map<string,
                double>({pair("x", -0.42), pair("y", .15)}), vec_dummy, "a", "reaction"},
        vertex_props{vertex_index_t(1), unordered_map<string,
                double>({pair("x", .76), pair("y", -3.)}), vec_dummy, "b", "reaction"},
        vertex_props{vertex_index_t(2), unordered_map<string,
                double>({pair("x", 3.34), pair("y", -3.7)}), vec_dummy, "c", "reaction"},
        vertex_props{vertex_index_t(3), unordered_map<string,
                double>({pair("x", .78), pair("y", -.36)}), vec_dummy, "d", "reaction"},
        vertex_props{vertex_index_t(4), unordered_map<string,
                double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "e", "reaction"},
        vertex_props{vertex_index_t(5), unordered_map<string,
                double>({pair("x", .16), pair("y", .37)}), vec_dummy, "f", "reaction"},
        vertex_props{vertex_index_t(6), unordered_map<string,
                double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "g", "reaction"},
        vertex_props{vertex_index_t(7), unordered_map<string,
                double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "x", "organism"},
        vertex_props{vertex_index_t(8), unordered_map<string,
                double>({pair("x", .16), pair("y", .37)}), vec_dummy, "y", "organism"},
        vertex_props{vertex_index_t(9), unordered_map<string,
                double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "z", "organism"}
    };

    unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
        pair(pair(0, 1), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 2), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 5), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(1, 2), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 3), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 5), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(3, 4), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(5, 6), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 7), edge_props{unordered_map<string,
                double>({pair("x", .4), pair("y", -.3)}), "reaction_organism"}),
        pair(pair(2, 9), edge_props{unordered_map<string,
                double>({pair("x", .85), pair("y", .9)}), "reaction_organism"}),
        pair(pair(3, 7), edge_props{unordered_map<string,
                double>({pair("x", -.9), pair("y", 8)}), "reaction_organism"}),
        pair(pair(5, 8), edge_props{unordered_map<string,
                double>({pair("x", .7), pair("y", -.5)}), "reaction_organism"})
    };

    int objective_function = 2;
    double temperature = 1.;
    double delta_min = 1e-10;
    size_t l_min = 5;
    size_t l_max = 7;
    size_t max_iter = 25;

    vertex_subgraph expected_solution = {1, 2, 7, 8, 9};

    LocalSearch lso = LocalSearch(
        node_properties, edge_properties, objective_function, temperature, delta_min,
        l_min, l_max, 3, max_iter, 1.
    );

    GIVEN("single-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 1);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }
    GIVEN("multi-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 4);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }

    Py_Finalize();

}


TEST_CASE("gene") {

    Py_Initialize();

    unordered_set<pair<size_t, size_t>> edges = {
        pair(0, 1), pair(0, 2), pair(0, 5),
        pair(1, 2), pair(2, 3), pair(2, 5),
        pair(3, 4), pair(5, 6),
        pair(2, 7), pair(2, 9), pair(3, 7), pair(5, 8)
    };

    unordered_map<string, vector<double>> vec_dummy;
    vector<vertex_props> node_properties = {
        vertex_props{vertex_index_t(0), unordered_map<string,
                double>({pair("x", -0.42), pair("y", .15)}), vec_dummy, "a", "reaction"},
        vertex_props{vertex_index_t(1), unordered_map<string,
                double>({pair("x", .76), pair("y", -3.)}), vec_dummy, "b", "reaction"},
        vertex_props{vertex_index_t(2), unordered_map<string,
                double>({pair("x", 3.34), pair("y", -3.7)}), vec_dummy, "c", "reaction"},
        vertex_props{vertex_index_t(3), unordered_map<string,
                double>({pair("x", .78), pair("y", -.36)}), vec_dummy, "d", "reaction"},
        vertex_props{vertex_index_t(4), unordered_map<string,
                double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "e", "reaction"},
        vertex_props{vertex_index_t(5), unordered_map<string,
                double>({pair("x", .16), pair("y", .37)}), vec_dummy, "f", "reaction"},
        vertex_props{vertex_index_t(6), unordered_map<string,
                double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "g", "reaction"},
        vertex_props{vertex_index_t(7), unordered_map<string,
                double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "x", "gene"},
        vertex_props{vertex_index_t(8), unordered_map<string,
                double>({pair("x", .16), pair("y", .37)}), vec_dummy, "y", "gene"},
        vertex_props{vertex_index_t(9), unordered_map<string,
                double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "z", "gene"}
    };

    unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
        pair(pair(0, 1), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 2), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 5), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(1, 2), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 3), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 5), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(3, 4), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(5, 6), edge_props{unordered_map<string,
                double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 7), edge_props{unordered_map<string,
                double>({pair("x", .4), pair("y", .5)}), "reaction_gene"}),
        pair(pair(2, 9), edge_props{unordered_map<string,
                double>({pair("x", -.7), pair("y", .9)}), "reaction_gene"}),
        pair(pair(3, 7), edge_props{unordered_map<string,
                double>({pair("x", -.7), pair("y", .5)}), "reaction_gene"}),
        pair(pair(5, 8), edge_props{unordered_map<string,
                double>({pair("x", -1), pair("y", 1)}), "reaction_gene"})
    };

    int objective_function = 2;
    double temperature = 1.;
    double delta_min = 1e-10;
    size_t l_min = 5;
    size_t l_max = 7;
    size_t max_iter = 25;

    vertex_subgraph expected_solution = {1, 2, 7, 8, 9};

    LocalSearch lso = LocalSearch(
        node_properties, edge_properties, objective_function, temperature, delta_min,
        l_min, l_max, 3, max_iter, 1.
    );

    GIVEN("single-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 1);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }
    GIVEN("multi-thread") {
        vector<string> groups = {"x", "y"};
        lso.run_local_search(groups, 4);
        report_score_progression(lso);

        CHECK(solution_equals(lso.solution, expected_solution));
    }

    Py_Finalize();

}


TEST_CASE("Property getters and setters") {

    Py_Initialize();

    unordered_set<pair<size_t, size_t>> edges = {
        pair(0, 1), pair(0, 2), pair(0, 5),
        pair(1, 2), pair(2, 3), pair(2, 5),
        pair(3, 4), pair(5, 6)
    };

    unordered_map<string, vector<double>> vec_dummy;
    vector<vertex_props> node_properties = {
        vertex_props{vertex_index_t(0), unordered_map<string,
            double>({pair("x", -0.42), pair("y", .15)}), vec_dummy, "a", "reaction"},
        vertex_props{vertex_index_t(1), unordered_map<string,
            double>({pair("x", .76), pair("y", -3.)}), vec_dummy, "b", "reaction"},
        vertex_props{vertex_index_t(2), unordered_map<string,
            double>({pair("x", 3.34), pair("y", -3.7)}), vec_dummy, "c", "reaction"},
        vertex_props{vertex_index_t(3), unordered_map<string,
            double>({pair("x", .78), pair("y", -.36)}), vec_dummy, "d", "reaction"},
        vertex_props{vertex_index_t(4), unordered_map<string,
            double>({pair("x", -.32), pair("y", .02)}), vec_dummy, "e", "reaction"},
        vertex_props{vertex_index_t(5), unordered_map<string,
            double>({pair("x", .16), pair("y", .37)}), vec_dummy, "f", "reaction"},
        vertex_props{vertex_index_t(6), unordered_map<string,
            double>({pair("x", -.3), pair("y", .04)}), vec_dummy, "g", "reaction"}
    };

    unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
        pair(pair(0, 1), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 2), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(0, 5), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(1, 2), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 3), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(2, 5), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(3, 4), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"}),
        pair(pair(5, 6), edge_props{unordered_map<string,
            double>({pair("x", 0), pair("y", 0)}), "reaction_reaction"})
    };

    int objective_function = 0;
    double temperature = 1.;
    double delta_min = 1e-10;
    size_t l_min = 3;
    size_t l_max = 5;
    size_t max_iter = 25;
    float p_norm = 2.;

    LocalSearch lso = LocalSearch(
        node_properties, edge_properties, objective_function, temperature,
        delta_min, l_min, l_max, 3, max_iter, p_norm
    );

    CHECK_EQ(lso.get_lmin(), l_min);
    lso.set_lmin(4);
    CHECK_EQ(lso.get_lmin(), 4);

    CHECK_EQ(lso.get_lmax(), l_max);
    lso.set_lmax(7);
    CHECK_EQ(lso.get_lmax(), 7);

    CHECK_EQ(lso.get_maxiter(), max_iter);
    lso.set_maxiter(55);
    CHECK_EQ(lso.get_maxiter(), 55);

    CHECK_EQ(lso.get_temp(), temperature);
    lso.set_temp(25.);
    CHECK_EQ(lso.get_temp(), 25.);

    CHECK_EQ(lso.get_min_reactions(), 3);
    lso.set_min_reactions(2);
    CHECK_EQ(lso.get_min_reactions(), 2);

    CHECK_EQ(lso.get_deltamin(), delta_min);
    lso.set_deltamin(.001);
    CHECK_EQ(lso.get_deltamin(), .001);

    CHECK_EQ(lso.get_lp_norm(), p_norm);
    lso.set_lp_norm(2.);
    CHECK_EQ(lso.get_lp_norm(), 2.);

    DOCTEST_CHECK_THROWS_AS(
        lso.get_objective_values(), typename std::invalid_argument);

    Py_Finalize();
}
