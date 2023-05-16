#include "pyutils.hpp"


void edgelist_to_pyobject(const unordered_set<Edge>& edges, PyObject *pyedges) {
    Py_ssize_t i = 0;
    for (auto edge: edges) {
        PyObject *pyedge = PyTuple_New(2);
        PyTuple_SetItem(pyedge, 0, PyLong_FromUnsignedLong(edge.first));
        PyTuple_SetItem(pyedge, 0, PyLong_FromUnsignedLong(edge.second));
        PyList_SetItem(pyedges, i, pyedge);
        i++;
    }
}

string get_list_element(PyObject *list, unsigned position) {
    PyObject *byte_item = PyList_GetItem(list, position);
    if (!byte_item) {
        return nullptr;
    }
    string item = PyBytes_AsString(PyUnicode_AsASCIIString(byte_item));
    if (item.empty()) {
        return nullptr;
    }
    return item;
}

string get_list_pair_string(PyObject *list, unsigned list_idx, unsigned tuple_idx) {
    PyObject *tuple = PyList_GetItem(list, list_idx);
    return PyBytes_AsString(PyUnicode_AsASCIIString(PyTuple_GET_ITEM(tuple, tuple_idx)));
}

size_t get_list_pair_int(PyObject *list, unsigned list_idx, unsigned tuple_idx) {
    PyObject *tuple = PyList_GetItem(list, list_idx);
    return PyLong_AsLong(PyTuple_GET_ITEM(tuple, tuple_idx));
}

Edge tuple_to_edge(PyObject *tuple) {
     return std::make_pair(PyLong_AS_LONG(PyTuple_GET_ITEM(tuple, 0)),
                           PyLong_AS_LONG(PyTuple_GET_ITEM(tuple, 1)));
}

PyObject *dict_get(PyObject *dict, string* key) {
    PyObject *elem = PyDict_GetItem(dict, PyUnicode_FromString(key->c_str()));
    return elem;
}

string dict_get_str(PyObject *dict, string* key) {
    return PyBytes_AsString(PyDict_GetItem(dict, PyUnicode_FromString(key->c_str())));
}


vector<vector<double>> nested_list_to_vector(PyObject *list) {
    PyObject *row, *elem;
    row = PyList_GetItem(list, 0);
    if (!PyList_Check(row)) {
        throw std::runtime_error(
            "Expected a nested list, but got a flat list."
        );
    }
    elem = PyList_GetItem(row, 0);
    if (!PyFloat_Check(elem)) {
        throw std::runtime_error(
            "Expected a nested list of floats, but got a different element type."
        );
    }

    long i, j, nrow, ncol;
    nrow = PyList_Size(list);
    ncol = PyList_Size(row);
    vector<vector<double>> vec(nrow, vector<double>(ncol));
    for (i = 0; i < nrow; i++) {
        row = PyList_GetItem(list, i);
        for (j = 0; j < ncol; j++) {
            elem = PyList_GetItem(row, j);
            vec[i][j] = PyFloat_AsDouble(elem);
        }
    }
    return vec;
}

vector<double> list_to_vector(PyObject *list) {
    PyObject *elem = PyList_GetItem(list, 0);
    if (!PyFloat_Check(elem)) {
        throw std::runtime_error(
            "Expected a list of floats, but got a different element type."
        );
    }

    long n = PyList_Size(list), i;
    vector<double> vec(n);
    for (i = 0; i < n; i++) {
        elem = PyList_GetItem(list, i);
        vec[i] = PyFloat_AsDouble(elem);
    }
    return vec;
}

PyObject *vector_to_list(const vector<vector<double>>& x) {
    PyObject *list;
    long i, j, nrow, ncol;
    nrow = static_cast<long>(x.size());
    ncol = static_cast<long>(x[0].size());
    list = PyList_New(nrow);
    for (i = 0; i < nrow; i++) {
        PyObject *row = PyList_New(ncol);
        for (j = 0; j < ncol; j++) {
            PyList_SetItem(row, j, PyFloat_FromDouble(x[i][j]));
        }
        PyList_SetItem(list, i, row);
    }
    return list;
}