#ifndef ENRICHMENT_SRC_PYUTILS_HPP
#define ENRICHMENT_SRC_PYUTILS_HPP


#define PY_SSIZE_T_CLEAN


#include <string>
#include <vector>
#include <Python.h>
#include <boost/unordered_set.hpp>

using boost::unordered_set;
using std::string;
using std::vector;

typedef std::pair<size_t, size_t> Edge;

void edgelist_to_pyobject(const unordered_set<Edge>&, PyObject*);

string get_list_element(PyObject *, unsigned);

string get_list_pair_string(PyObject *list, unsigned list_idx, unsigned tuple_idx);

size_t get_list_pair_int(PyObject *list, unsigned list_idx, unsigned tuple_idx);

Edge tuple_to_edge(PyObject *);

PyObject *dict_get(PyObject *dict, string *key);

string dict_get_str(PyObject *dict, string *key);

vector<vector<double>> nested_list_to_vector(PyObject *list);
vector<double> list_to_vector(PyObject *list);

PyObject *vector_to_list(const vector<vector<double>>&);

#endif //ENRICHMENT_SRC_PYUTILS_HPP
