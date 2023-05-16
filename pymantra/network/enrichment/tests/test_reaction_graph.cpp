#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#define PY_SSIZE_T_CLEAN


#include "../doctest.h"
#include <iostream>
#include <string>
#include <vector>
#include <Python.h>
#include "../LSO/reaction_graph.hpp"


using std::cout;
using std::endl;
using std::string;
using std::vector;
using boost::vertex;


typedef unordered_map<string, double>::iterator data_iter;
typedef unordered_map<pair<size_t, size_t>, edge_props>::iterator edge_iter;


bool py_subgraphs_equal(PyObject *computed, PyObject *expected) {
    Py_ssize_t c_size = PySet_Size(computed);
    Py_ssize_t e_size = PySet_Size(expected);
    if (c_size == 0 | c_size != e_size) {
        cout << "Incorrect graph sizes: ";
        cout << "expected: " << e_size;
        cout << " computed: " << c_size;
        return false;
    }
    bool eq = true;
    Py_ssize_t i;
    PyObject *edge, *iter;
    iter = PyObject_GetIter(computed);
    for (i = 0; i < c_size; i++) {
       edge = PyIter_Next(iter);
       if (!PySet_Contains(computed, edge)) {
           eq = false;
           break;
       }
    }
    iter = PyObject_GetIter(expected);
    for (i = 0; i < c_size; i++) {
        edge = PyIter_Next(iter);
        if (!PySet_Contains(expected, edge)) {
            eq = false;
            break;
        }
    }
    if (!eq) {
        cout << "Computed subgraph is incorrect!" << endl;
        cout << "Expected:" << endl;
        iter = PyObject_GetIter(expected);
        for (i = 0; i < c_size; i++) {
            edge = PyIter_Next(iter);
            cout << "(" << PyTuple_GetItem(edge, 0) << ",";
            cout << PyTuple_GetItem(edge, 1) << ") ";
        }
        cout << "Computed:" << endl;
        iter = PyObject_GetIter(computed);
        for (i = 0; i < c_size; i++) {
            edge = PyIter_Next(iter);
            cout << "(" << PyTuple_GetItem(edge, 0) << ",";
            cout << PyTuple_GetItem(edge, 1) << ") ";
        }
    }
    return eq;
}


void data_to_pydict(unordered_map<string, double> data, PyObject *dict) {
    data_iter iter;
    for (iter = data.begin(); iter != data.end(); iter++) {
        PyDict_SetItem(
            dict, PyUnicode_FromString(iter->first.c_str()),
            PyFloat_FromDouble(iter->second)
        );
    }
}


TEST_SUITE_BEGIN("reaction graph");

TEST_CASE("reaction graph from multi-omics network") {
    Py_Initialize();

    unordered_map<string, vector<double>> vec_dummy;
    vector<vertex_props> node_properties = {
        vertex_props{static_cast<vertex_index_t>(0), unordered_map<string, double>({pair("x", .42)}), vec_dummy, "1", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(1), unordered_map<string, double>({pair("x", .76)}), vec_dummy, "2", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(2), unordered_map<string, double>({pair("x", .34)}), vec_dummy, "3", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(3), unordered_map<string, double>({pair("x", .78)}), vec_dummy, "4", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(4), unordered_map<string, double>({pair("x", .32)}), vec_dummy, "5", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(5), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "6", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(6), unordered_map<string, double>({pair("x", .32)}), vec_dummy, "7", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(7), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "8", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(8), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "9", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(9), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "10", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(10), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "11", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(11), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "12", "metabolite"},
        vertex_props{static_cast<vertex_index_t>(12), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "a", "reaction"},
        vertex_props{static_cast<vertex_index_t>(13), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "b", "reaction"},
        vertex_props{static_cast<vertex_index_t>(14), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "c", "reaction"},
        vertex_props{static_cast<vertex_index_t>(15), unordered_map<string, double>({pair("x", .16)}), vec_dummy, "d", "reaction"}
    };

    unordered_map<pair<size_t, size_t>, edge_props> edge_properties = {
        pair(pair(0, 14), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(0, 15), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(1, 13), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(1, 14), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(2, 12), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(2, 14), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(3, 15), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(4, 15), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(5, 12), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(5, 15), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(6, 14), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(6, 15), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(7, 13), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(7, 14), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(8, 13), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(9, 13), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(10, 12), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
        pair(pair(11, 12), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"})
    };

    Graph graph;

    // expected edges: pair("a", "c"), pair("a", "d"), pair("b", "c"), pair("c", "d")
    vector<string> exp_src = {"a", "a", "b", "c"};
    vector<string> exp_tgt = {"c", "d", "c", "d"};

    // build the boost graph
    int i;
    Graph::vertex_descriptor v;
    for (auto &v_prop : node_properties) {
        v = boost::add_vertex(graph);
        graph[v].index = v_prop.index;
        graph[v].name = v_prop.name;
        graph[v].type = v_prop.type;
    }
    for (auto &e_prop : edge_properties) {
        boost::add_edge(
            vertex(e_prop.first.first, graph), vertex(e_prop.first.second, graph),
            e_prop.second, graph
        );
    }

    PyObject *computed_graph, *expected_graph;

    // convert the expected edges to a python object representing the subgraph
    expected_graph = PySet_New(nullptr);
    for (i = 0; i < exp_src.size(); i++) {
        PyObject *pyedge = PyTuple_New(2);
        PyTuple_SetItem(pyedge, 0, PyUnicode_FromString(exp_src[i].c_str()));
        PyTuple_SetItem(pyedge, 1, PyUnicode_FromString(exp_tgt[i].c_str()));
        PySet_Add(expected_graph, pyedge);
    }

    GIVEN("c++ objects") {
        // computing the reaction subgraph
        subgraph reaction_graph;
        extract_edges(graph, reaction_graph);

        // converting subgraph to python object
        computed_graph = subgraph_to_py(reaction_graph);

        // actual comparison
        bool eq = py_subgraphs_equal(computed_graph, expected_graph);
        CHECK(eq);
    }

    GIVEN("python objects") {

        unordered_map<pair<string, string>, edge_props> edge_ps = {
            pair(pair("1", "c"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("1", "d"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("2", "b"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("2", "c"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("3", "a"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("3", "c"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("4", "d"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("5", "d"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("6", "a"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("6", "d"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("7", "c"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("7", "d"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("8", "b"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("8", "c"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("9", "b"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("10", "b"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("11", "a"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"}),
            pair(pair("12", "a"), edge_props{unordered_map<string, double>({pair("x", 0), pair("y", 0)}), "metabolite_reaction"})
        };
        PyObject *pynodes, *pyedges;
        pynodes = PyDict_New();
        for (i = 0; i < node_properties.size(); i++) {
            PyObject *pynode = PyDict_New(), *pydata = PyDict_New();
            data_to_pydict(node_properties[i].value, pydata);
            PyDict_SetItem(pynode, PyUnicode_FromString("index"),
                           PyLong_FromSize_t(i));
            PyDict_SetItem(pynode, PyUnicode_FromString("data"), pydata);
            PyDict_SetItem(pynode, PyUnicode_FromString("node_type"),
                           PyUnicode_FromString(node_properties[i].type.c_str()));
            PyDict_SetItem(pynodes, PyUnicode_FromString(node_properties[i].name.c_str()),
                           pynode);
        }
        pyedges = PyDict_New();
        unordered_map<pair<string, string>, edge_props>::iterator edge_;
        for (edge_ = edge_ps.begin(); edge_ != edge_ps.end(); edge_++) {
            PyObject *pyedge = PyDict_New(), *pydata = PyDict_New();
            PyObject *nodes = PyTuple_New(2);

            data_to_pydict(edge_->second.value, pydata);
            PyDict_SetItem(pyedge, PyUnicode_FromString("data"), pydata);
            PyDict_SetItem(pyedge, PyUnicode_FromString("edge_type"),
                           PyUnicode_FromString(edge_->second.type.c_str()));

            string src = edge_->first.first;
            string tgt = edge_->first.second;
            PyTuple_SetItem(nodes, 0, PyUnicode_FromString(edge_->first.first.c_str()));
            PyTuple_SetItem(nodes, 1, PyUnicode_FromString(edge_->first.second.c_str()));

            PyDict_SetItem(pyedges, nodes, pyedge);
        }
        PyObject *comp_graph = extract_reaction_graph(pynodes, pyedges);
        bool eq = py_subgraphs_equal(comp_graph, expected_graph);
        CHECK(eq);
    }

    Py_Finalize();
}
