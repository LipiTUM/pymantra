#include "Exceptions.hpp"
#include "utils.hpp"
#include <iostream>
#include <boost/format.hpp>


using std::string;


const char *BaseException::what() const noexcept {
    return error.c_str();
}

MissingComputation::MissingComputation(const string &function_call) {
    BaseException::error = "Function '" + function_call + "' was not called yet!";
}

InvalidNodeType::InvalidNodeType(const string &node_type) {
    BaseException::error = "'" + node_type + "' is not a valid node type";
}

InvalidEdgeNode::InvalidEdgeNode(const string &message) {
    BaseException::error = message;
}

InvalidEdgeNode::InvalidEdgeNode(const string &node_type, const string &edge_type) {
    BaseException::error = "'" + node_type + "' is not a valid end point for an edge of type '" + edge_type + "'!";
}

InvalidEdgeType::InvalidEdgeType(const string &edge_type) {
    BaseException::error = "'" + edge_type + "' is not a valid edge type";
}

PyObjectType::PyObjectType(const string &message) {
    BaseException::error = message;
}

IncorrectPyObjectElementType::IncorrectPyObjectElementType(const string &message) {
    BaseException::error = message;
}

IncorrectPyArrayDimensions::IncorrectPyArrayDimensions(const string &message) {
    BaseException::error = message;
}

ConvergenceError::ConvergenceError(const string &message) {
    BaseException::error = message;
}

InvalidObjectiveFunction::InvalidObjectiveFunction(const string &message) {
    InvalidObjectiveFunction::error = message;
}

MissingAttribute::MissingAttribute(const string &message) {
    MissingAttribute::error = message;
}


void LOG_CERR(const string& m) {
    std::cerr << boost::str(boost::format("warning: %1%") % m) << std::endl;
}
