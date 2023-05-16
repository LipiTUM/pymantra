#ifndef ENRICHMENT_REACTION_GRAPH_HPP
#define ENRICHMENT_REACTION_GRAPH_HPP


#define PY_SSIZE_T_CLEAN


#include <string>
#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/unordered_set.hpp>
#include <boost/unordered_map.hpp>
#include <Python.h>


using boost::vertex_index_t;
using std::string;
using std::pair;
using std::vector;
using boost::unordered_set;
using boost::unordered_map;


// boost graph types
struct vertex_props {
    vertex_index_t index;
    unordered_map<string, double> value;
    unordered_map<string, vector<double>> vec_value;
    string name;
    string type;
    double objective;
};


struct edge_props {
    unordered_map<string, double> value;
    string type;
    double objective;
};


/* NOTE: since we do not do any add/remove operations on the main graph
 *       the time complexity of listS does not come into play. Hence,
 *       using vecS to reduce per-vertex space overhead and have a lower
 *       constant time factor on iterator operators.
 */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
                              vertex_props, edge_props> Graph;
typedef typename boost::graph_traits<Graph> GraphTrait;

typedef typename GraphTrait::vertex_iterator VertexIter;
typedef typename GraphTrait::adjacency_iterator AdjacencyIter;


typedef pair<string, string> edge_str;
typedef pair<vertex_index_t, vertex_index_t> edge;
typedef unordered_set<edge_str> subgraph;


string pystring_to_string(PyObject *pystring);

unordered_map<string, double> get_group_values(PyObject *values);
unordered_map<string, vector<double>> get_vec_group_values(PyObject *values);

Graph pygraph_to_boost(PyObject *node_properties, PyObject *edge_properties,
                       bool include_attributes = true);


void step(Graph& graph, subgraph& edges, vertex_index_t reaction_vertex);

void extract_edges(Graph& graph, subgraph& pyedges);


PyObject *subgraph_to_py(subgraph& edges);



/**
 * Main method. Do a BFS and collect reaction - reaction edges on the way
 * @param node_properties
 * @param edges
 * @param include_attributes
 * @return
 */
PyObject *extract_reaction_graph(PyObject *node_properties, PyObject *edge_properties,
                                 bool include_attributes = true);


#endif //ENRICHMENT_REACTION_GRAPH_HPP
