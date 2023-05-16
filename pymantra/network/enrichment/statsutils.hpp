
#ifndef ENRICHMENT_SRC_STATSUTILS_HPP
#define ENRICHMENT_SRC_STATSUTILS_HPP


bool HAS_PARALLEL();


#include <vector>
#include <Python.h>


using std::vector;
using std::pair;


typedef pair<vector<vector<double>>, vector<vector<double>>> SpearmansResults;


bool isnan_(double);
vector<double> check_nan_safe(vector<double>&);

double fishers_exact(const unsigned int *);

double mean(vector<double>&);
double var(vector<double>&);
double var(vector<double>&, double);

vector<double> norm_cdf(vector<double>&);
double norm_cdf(double);

double welch_test(vector<double>&, vector<double>&, bool return_statistic=false);

// auxiliaries
pair<bool, vector<vector<bool>>> get_nan_masks(vector<vector<double>>&);
vector<double> to_ranks(vector<double>&);
// actual coefficient calculation
double spearman_coefficient(vector<double>&, vector<double>&);
// p-value from coefficient + df
double spearman_pvalue(double, unsigned);
// wrapper for 1D x 1D
pair<double, double> spearmans_rank(vector<double>&, vector<double>&);
// wrapper with precomputed nans
pair<double, double> spearmans_rank(vector<double>&, vector<double>&, vector<bool>&, vector<bool>&);
// wrapper for 2D x 2D
SpearmansResults spearmans_rank(vector<vector<double>>&, vector<vector<double>>&,
                                      int n_threads=1);
// wrapper with pre-computed nans
SpearmansResults spearmans_rank(vector<vector<double>>&, vector<vector<double>>&,
                                      vector<vector<bool>>&, vector<vector<bool>>&, int n_threads=1);

// for python bindings
PyObject* spearmans_rank_1by1(PyObject*, PyObject*);
PyObject* spearmans_rank_2by2(PyObject*, PyObject*, int n_threads=1);


#endif //ENRICHMENT_SRC_STATSUTILS_HPP
