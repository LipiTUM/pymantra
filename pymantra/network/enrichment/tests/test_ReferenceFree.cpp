#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define PY_SSIZE_T_CLEAN


#include "../doctest.h"
#include "../RFA/ReferenceFree.hpp"
#include "../RFA/ReferenceNetwork.hpp"
#include "../Exceptions.hpp"
#include "../utils.hpp"
#include <Python.h>
#include <numpy/ndarrayobject.h>
#include <boost/unordered_set.hpp>
#include <iostream>


using std::cout;
using std::endl;


bool networks_identical(unordered_set<Edge> net1, vector<Edge> &net2) {
    for (auto edge : net1) {
        for (auto vec_edge : net2) {
            if (edge.first == vec_edge.first & edge.second == vec_edge.second) {
                return false;
            }
        }
    }
    for (auto edge : net2) {
        if (net1.find(edge) != net1.end()) {
            return false;
        }
    }
    return true;
}


TEST_CASE("") {
    vector<strEdge> edges = {
            {"1", "2"},
            {"0", "3"},
            {"1", "3"},
            {"0", "2"},
            {"0", "4"}
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
    vector<string> node_order = vector<string>({"0", "1", "2", "3", "4"});
    GIVEN("a vector of edges and node types") {
        try {
            ReferenceFree reffree = ReferenceFree(
                correlations, node_order, &edges[0], 5, edge_types, node_types
            );
            reffree.set_n_threads(2);

            WHEN("Computing cutoffs") {
                reffree.compute_optimal_network();
                vector<double> best_cutoffs = reffree.get_optimal_cutoff();
                // TODO: compute cutoffs manually and check if they are computed correctly
                CHECK_EQ(best_cutoffs[0], 0.);
                CHECK_EQ(best_cutoffs[1], 1.);
                // fisher's exact progression
                int i, j;
                vector<vector<double>> progressions = reffree.get_overlap_progression();
                printf("Overlap score progressions:");
                for (i = 0; i < progressions.size(); i++) {
                    for (j = 0; j < progressions[i].size(); j++) {
                        cout << progressions[i][j] << ", ";
                    }
                    cout << endl << endl;
                }
            }

            WHEN("Computing enrichment") {
                // TODO
            }
        } catch (const std::exception& e) {
            FAIL("Unexpected exception: ", e.what());
        }
    }
    GIVEN("Python data structures") {
        Py_Initialize();
        _import_array();

        int i;
        npy_intp dims[2] = {5, 5};
        double corr_arr[25] = {
                1, .4, -.9, .1, -.3,
                .4, 1, .6, -.4, -.7,
                -.9, .6, 1, .85, .78,
                .1, -.4, .85, 1, -.65,
                -.3, -.7, .78, -.65, 1
        };
        PyObject *pycorrs = PyArray_SimpleNewFromData(2, &dims[0], NPY_DOUBLE, &corr_arr);
        PyArrayObject *pycorrarr = PyArray_GETCONTIGUOUS(reinterpret_cast<PyArrayObject *>(pycorrs));
        PyObject *pyorder = PyList_New(5);
        for (i = 0; i < 5; i++) {
            PyList_SetItem(pyorder, i, PyUnicode_FromString(node_order[0].c_str()));
        }
        PyObject *edge_tuple = PyTuple_New(2);
        PyObject *pyedges = PyList_New((long)edges.size());
        for (i = 0; i < 5; i++) {
            PyTuple_SetItem(edge_tuple, 0, PyUnicode_FromString(edges[i].first.c_str()));
            PyTuple_SetItem(edge_tuple, 1, PyUnicode_FromString(edges[i].second.c_str()));
            PyList_SetItem(pyedges, i, edge_tuple);
        }
        PyObject *pyetypes = PyList_New(5);
        for (i = 0; i < 5; i++) {
            PyList_SetItem(pyetypes, i, PyUnicode_FromString(edge_types[0].c_str()));
        }
        PyObject *pyntypes = PyList_New(5);
        for (i = 0; i < 5; i++) {
            PyList_SetItem(pyntypes, i, PyUnicode_FromString(node_types[0].c_str()));
        }
        try {
            ReferenceFree reffree = ReferenceFree(
                pycorrarr, pyorder, pyedges,
                pyetypes, pyntypes
            );
            reffree.set_n_threads(2);

            WHEN("Computing cutoffs") {
                reffree.compute_optimal_network();
                vector<double> best_cutoffs = reffree.get_optimal_cutoff();
                // TODO: compute cutoffs manually and check if they are computed correctly
                CHECK_EQ(best_cutoffs[0], 0.);
                CHECK_EQ(best_cutoffs[1], 1.);
                // fisher's exact progression
                int j;
                vector<vector<double>> progressions = reffree.get_overlap_progression();
                printf("Overlap score progressions:");
                for (i = 0; i < progressions.size(); i++) {
                    for (j = 0; j < progressions[i].size(); j++) {
                        cout << progressions[i][j] << ", ";
                    }
                    cout << endl << endl;
                }
            }

            WHEN("Computing enrichment") {
                // TODO
            }
        } catch (const std::exception& e) {
            FAIL("Unexpected exception during ReferenceFree test: ", e.what());
        }
    }
}