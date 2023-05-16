#include <boost/format.hpp>
#include "reaction_graph.hpp"
#include "../Exceptions.hpp"
#include "../pyutils.hpp"


using std::vector;
using boost::unordered_map;
using boost::format;
using boost::str;
using boost::tie;
using boost::make_iterator_range;
using boost::adjacent_vertices;
using boost::vertex;


string pystring_to_string(PyObject *pystring) {
    PyObject *pyrepr = PyObject_Str(pystring);
    return PyUnicode_AsUTF8(pyrepr);
}


unordered_map<string, double> get_group_values(PyObject *values) {
    if(!PyDict_Check(values)) {
        throw IncorrectPyObjectElementType("Values must be provided as a dictionary");
    }
    unordered_map<string, double> vals;
    PyObject *group, *val;
    Py_ssize_t pos = 0;
    while(PyDict_Next(values, &pos, &group, &val)) {
        vals.insert(pair<string, double>(pystring_to_string(group), PyFloat_AsDouble(val)));
    }
    return vals;
}

unordered_map<string, vector<double>> get_vec_group_values(PyObject *values) {
    if(!PyDict_Check(values)) {
        throw IncorrectPyObjectElementType("Values must be provided as a dictionary");
    }
    unordered_map<string, vector<double>> vals;
    PyObject *group, *val;
    Py_ssize_t pos = 0, i;
    while(PyDict_Next(values, &pos, &group, &val)) {
        vector<double> tmp_vec(PyList_Size(val));
        for (i = 0; i < PyList_Size(val); i++) {
            tmp_vec[i] = PyFloat_AsDouble(PyList_GetItem(val, i));
        }
        vals.insert(pair<string, vector<double>>(pystring_to_string(group), tmp_vec));
    }
    return vals;
}

Graph pygraph_to_boost(PyObject *node_properties, PyObject *edge_properties, bool include_attributes) {
    string DATA_KEY = "data";
    string VEC_DATA_KEY = "vec_group_data";
    string NODE_TYPE_KEY = "node_type";
    string EDGE_TYPE_KEY = "edge_type";

    Graph graph;
    if(!PyDict_Check(node_properties)) {
        throw IncorrectPyObjectElementType("'node_properties' must be a python dictionary");
    }
    if(!PyDict_Check(edge_properties)) {
        throw IncorrectPyObjectElementType("'edge_properties' must be a python dictionary");
    }
    PyObject *node, *props, *edge, *group_data, *type_, *vec_group_data;
    string node_name, node_type;
    Py_ssize_t pos = 0;
    int i = 0;
    unordered_map<string, int> vertices;
    Graph::vertex_descriptor vertex_;
    while(PyDict_Next(node_properties, &pos, &node, &props)) {
        node_name = pystring_to_string(node);
        vertices.insert(pair<string, int>(node_name, i));
        i++;
        vertex_ = boost::add_vertex(graph);
        graph[vertex_].index = static_cast<vertex_index_t>(vertex_);
        graph[vertex_].name = node_name;
        graph[vertex_].objective = 0.;

        type_ = dict_get(props, &NODE_TYPE_KEY);
        node_type = pystring_to_string(type_);
        if (type_ == nullptr) {
            throw MissingAttribute(
                str(
                    format(
                        "Node type key '%1%' not found for node '%2%'. Please make sure all nodes "
                        "contain an attribute named '%1%'!"
                    ) % NODE_TYPE_KEY % node_name
                )
            );
        }
        graph[vertex_].type = node_type;

        if (include_attributes) {
            // NOTE: we only the data for reaction nodes, if we do this kind of extraction
            // in case we want genes or organisms this would be very simple to add in, however
            // they are most likely going to have edge values instead (i.e. associations)
            group_data = dict_get(props, &DATA_KEY);
            vec_group_data = dict_get(props, &VEC_DATA_KEY);
            if (node_type == "reaction") {
                if ((group_data == nullptr) & (vec_group_data == nullptr)) {
                    throw MissingAttribute(
                        str(
                            format(
                                "Group data key '%1%' or '%2%' not found for reaction node '%3%'. "
                                "Please make sure all nodes contain an attribute named '%1%' or '%2%'!"
                            ) % DATA_KEY % VEC_DATA_KEY % node_name
                        )
                    );
                }
                if (group_data != nullptr) {
                    unordered_map<string, double> gd = get_group_values(group_data);
                    graph[vertex_].value = get_group_values(group_data);
                }
                if (vec_group_data != nullptr) {
                    graph[vertex_].vec_value = get_vec_group_values(vec_group_data);
                }
            }
        }
    }
    vector<edge_props> edges;
    string src, tgt, edge_type;
    int src_id, tgt_id;
    pos = 0;
    while(PyDict_Next(edge_properties, &pos, &edge, &props)) {
        src = pystring_to_string(PyTuple_GetItem(edge, 0));
        src_id = vertices[src];
        tgt = pystring_to_string(PyTuple_GetItem(edge, 1));
        tgt_id = vertices[tgt];

        type_ = dict_get(props, &EDGE_TYPE_KEY);
        edge_type = pystring_to_string(type_);
        if (type_ == nullptr) {
            throw MissingAttribute(
                str(
                    format(
                        "Edge type key '%1%' not found for edge ('%2%', '%3%'). Please make sure all edges contain"
                        " an attribute named '%1%'!"
                    ) % EDGE_TYPE_KEY % src % tgt
                )
            );
        }
        if (include_attributes) {

            group_data = dict_get(props, &DATA_KEY);
            // TODO: do we need the data attributes for other edge types too?
            if ((edge_type == "REACTION_ORGANISM" || (edge_type == "REACTION_CATALYST"))) {
                if (group_data == nullptr) {
                    throw MissingAttribute(
                        str(
                            format(
                                "Group data key '%1%' not found for edge ('%2%', '%3%'). Please make sure all edges "
                                "contain an attribute named '%1%'!"
                            ) % DATA_KEY % src % tgt
                        )
                    );
                }
                boost::add_edge(
                    vertex(src_id, graph), vertex(tgt_id, graph),
                    {get_group_values(group_data), pystring_to_string(type_), 0.},
                    graph
                );
            } else {
                boost::add_edge(
                    vertex(src_id, graph), vertex(tgt_id, graph),
                    {{}, pystring_to_string(type_), 0.},
                    graph
                );
            }
        }
        else {
            boost::add_edge(
                    vertex(src_id, graph), vertex(tgt_id, graph),
                    {{}, "", 0.},
                    graph
            );
        }
    }
    // TODO: is there a more efficient way than adding nodes/edges from a loop?
    // e.g.:
    // edge_pair edge_arr[n_edges];
    // for (i = 0; i < n_edges; i++) {
    //     src_node = get_list_pair_int(edges, i, 0);
    //     tgt_node = get_list_pair_int(edges, i, 1);
    //     edge_arr[i] = edge_pair(src_node, tgt_node);
    // }
    // graph = Graph(edge_arr, edge_arr + sizeof(edge_arr) / sizeof(edge_pair),
    //               PyDict_GET_SIZE(node_properties));

    return graph;
}


void step(Graph& graph, subgraph& edges, vertex_index_t metabolite_vertex) {
    unordered_set<vertex_index_t> adjacent_reactions;
    AdjacencyIter vertex, v_end;
    for (tie(vertex, v_end) = adjacent_vertices(metabolite_vertex, graph); vertex != v_end; vertex++) {
        if (graph[*vertex].type == "reaction") adjacent_reactions.insert(graph[*vertex].index);
    }
    // TODO: get all pairwise combinations => add edge
    for (auto src_reaction : adjacent_reactions) {
        for (auto tgt_reaction : adjacent_reactions) {
            if (src_reaction != tgt_reaction) {
                // this guarantees ordering of the nodes, since we work with undirected
                // graphs without parallel edges
                if (src_reaction < tgt_reaction) {
                    edges.insert(edge_str(graph[src_reaction].name, graph[tgt_reaction].name));
                } else {
                    edges.insert(edge_str(graph[tgt_reaction].name, graph[src_reaction].name));
                }
            }
        }
    }
}


void extract_edges(Graph& graph, subgraph& edges) {
    vector<bool> seen(boost::num_vertices(graph), false);
    VertexIter start_node, s_end;
    AdjacencyIter vertex, v_end;
    string node_type;
    for (tie(start_node, s_end) = vertices(graph); start_node != s_end; start_node++) {
        node_type = graph[*start_node].type;
        if (graph[*start_node].type != "reaction") continue;
        for (tie(vertex, v_end) = adjacent_vertices(*start_node, graph); vertex != v_end; vertex++) {
            if ((graph[*vertex].type != "metabolite") || (seen[graph[*vertex].index])) continue;
            // => unseen metabolite node
            step(graph, edges, graph[*vertex].index);
            seen[graph[*vertex].index] = true;
        }
    }
}


PyObject *subgraph_to_py(subgraph& edges) {
    PyObject *pyedges = PySet_New(nullptr);
    for (const auto& edge : edges) {
        // TODO: does this object survive the loop?
        PyObject *pyedge = PyTuple_New(2);
        PyTuple_SetItem(pyedge, 0, PyUnicode_FromString(edge.first.c_str()));
        PyTuple_SetItem(pyedge, 1, PyUnicode_FromString(edge.second.c_str()));
        PySet_Add(pyedges, pyedge);
    }
    return pyedges;
}


PyObject *extract_reaction_graph(PyObject *node_properties, PyObject *edge_properties, bool include_attributes) {
    // generate boost graph
    Graph graph = pygraph_to_boost(node_properties, edge_properties, include_attributes);
    subgraph reaction_graph;
    // edges for the final subgraph
    extract_edges(graph, reaction_graph);
    // convert to python object
    return subgraph_to_py(reaction_graph);
}
