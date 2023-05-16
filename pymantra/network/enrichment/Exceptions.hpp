#include <exception>
#include <string>


#ifndef ENRICHMENT_SRC_EXCEPTIONS_HPP
#define ENRICHMENT_SRC_EXCEPTIONS_HPP


using std::string;

struct BaseException : public std::exception
{
    string error;
    const char * what () const noexcept override;
};

struct MissingComputation : BaseException
{
    explicit MissingComputation(const string&);
};

struct InvalidNodeType : BaseException
{
    explicit InvalidNodeType(const string&);
};

struct InvalidEdgeNode : BaseException {
    explicit InvalidEdgeNode(const string&);
    InvalidEdgeNode(const string&, const string&);
};

struct InvalidEdgeType : BaseException
{
    explicit InvalidEdgeType(const string&);
};

struct PyObjectType : BaseException
{
    explicit PyObjectType(const string&);
};

struct IncorrectPyObjectElementType : BaseException
{
    explicit IncorrectPyObjectElementType(const string&);
};

struct IncorrectPyArrayDimensions : BaseException
{
    explicit IncorrectPyArrayDimensions(const string&);
};


struct ConvergenceError : BaseException
{
    explicit ConvergenceError(const string&);
};

struct InvalidObjectiveFunction : BaseException
{
    explicit InvalidObjectiveFunction(const string&);
};

struct MissingAttribute : BaseException
{
    explicit MissingAttribute(const string&);
};


void LOG_CERR(const string&);

#endif //ENRICHMENT_SRC_EXCEPTIONS_HPP
