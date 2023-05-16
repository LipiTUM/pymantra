#include "statsutils.hpp"
#include "Exceptions.hpp"
#include "pyutils.hpp"
#include <cmath>
#include <boost/math/distributions/hypergeometric.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/statistics/t_test.hpp>
#include <boost/math/special_functions/beta.hpp>


#ifdef USE_OPENMP
    #include <omp.h>

    bool HAS_PARALLEL() {
        return true;
    }
#else
    bool HAS_PARALLEL() {
        return false;
    }
#endif


using boost::math::hypergeometric_distribution;
using boost::math::pdf;
using boost::math::normal;
using boost::math::cdf;
using boost::math::ibeta;
using boost::math::students_t;
using boost::math::statistics::two_sample_t_test;
using boost::math::isnan;


bool isnan_(double x) {
   return isnan(x);
}


vector<double> check_nan_safe(vector<double> &x) {
    vector<double> y;
    unsigned i;
    for (i = 0; i < x.size(); i++) {
        if (!isnan(x[i]))
            y.push_back(x[i]);
    }
    if (y.empty()) {
        throw std::invalid_argument(
            "Found 'nan' values only, but non-nan values required to compute statistics!"
        );
    } else if (y.size() < 3) {
        throw std::invalid_argument(
            "Less than 3 non-nan values found, but at least 3 are required to compute statistics!"
        );
    }
    return y;
}


/**
 * One-sided Fisher's exact test
 * @param cont array of doubles of size four indicating the contingency table.
 *             In the context of this library the values are (in order):
 *             0: True Positives
 *             1: False Positives
 *             2: False Negatives
 *             3: True Negatives
 *
 * @return double, p-value
 */
double fishers_exact(const unsigned *cont) {
    unsigned n_total = cont[0] + cont[1] + cont[2] + cont[3];
    unsigned n1 = cont[0] + cont[2];
    unsigned n2 = cont[2] + cont[3];
    hypergeometric_distribution<> hg_dist(n1, n2, n_total);

    unsigned zero = 0;
    unsigned max_k = std::min(n1, n2);
    unsigned min_k = std::max(zero, n1 + n2 - n_total);
    double thresh = pdf(hg_dist, cont[2]);

    double pval = 0.;
    double pi;
    unsigned k;
    for (k = min_k; k <= max_k; k++) {
        pi = pdf(hg_dist, k);
        if (pi <= thresh) pval += pi;
    }
    return pval;
}


/**
 * Computing the mean of a vector
 * @param x numeric vector to compute mean
 * @return mean value of `x`
 */
double mean(vector<double> &x) {
    double sum = .0;
    unsigned i;
    for (i = 0; i < x.size(); i++) {
        sum += x[i];
    }
    return sum / (double)i;
}


/**
 *
 * @param x
 * @param mu
 * @return
 */
double var(vector<double> &x, double mu) {
    double se = .0;
    unsigned i;
    for (i = 0; i < x.size(); i++) {
        se += std::pow(x[i] - mu, 2);
    }
    return se / (double) i;
}


double var(vector<double> &x) {
    double mu = mean(x);
    return var(x, mu);
}

vector<double> norm_cdf(vector<double>& x) {
    vector<double> cdf_vals(x.size());
    unsigned i;
    for (i = 0; i < x.size(); i++) {
        cdf_vals[i] = norm_cdf(x[i]);
    }
    return cdf_vals;
}

double norm_cdf(double z) {
    normal norm_dist;
    if (isnan(z))
        return 1.;
    return cdf(norm_dist, z);
}


double welch_test(vector<double> &x1, vector<double> &x2, bool return_statistic) {
    // TODO: check for nan => remove or similar
    vector<double> xs1, xs2;
    xs1 = check_nan_safe(x1);
    xs2 = check_nan_safe(x2);
    std::pair<double, double> welch = two_sample_t_test(xs1, xs2);
    if (return_statistic) return welch.first;
    return welch.second;
}


pair<bool, vector<vector<bool>>> get_nan_masks(vector<vector<double>>& x) {
    bool contains_nan = false;
    vector<vector<bool>> nans(x.size(), vector<bool>(x[0].size()));
    bool isna;
    unsigned i, j;
    for (i = 0; i < x.size(); i++) {
        for (j = 0; j < x[i].size(); j++) {
            isna = isnan(x[i][j]);
            nans[i][j] = isna;
            if (isna) contains_nan = true;
        }
    }
    return {contains_nan, nans};
}


vector<double> to_ranks(vector<double> &x) {
    // TODO: use merge sort to compute ranks
    unsigned n = x.size(), i, r, s, j;
    vector<double> ranked(n);
    for (i = 0; i < n; i++) {
        r = 1;
        s = 1;
        for (j = 0; j < i; j++) {
            if (x[j] < x[i]) r++;
            if (x[j] == x[i]) s++;
        }
        for (j = i + 1; j < n; j++) {
            if (x[j] < x[i]) r++;
            if (x[j] == x[i]) s++;
        }
        ranked[i] = r + (s - 1.) / 2.;
    }
    return ranked;
}

double spearman_coefficient(vector<double>& x, vector<double>& y) {
    unsigned n = x.size(), i;
    if (n != y.size()) {
        LOG_CERR(
            "Both arrays must have the same size to compute spearman's correlation coefficient"
        );
        return NAN;
    } else if (n < 3) {
        LOG_CERR(
            "At least three non-nan observations required to compute spearman's correlation coefficient"
        );
        return NAN;
    }
    double sx = 0, sy = 0, sxy = 0, ssx = 0, ssy = 0;
    for (i = 0; i < n; i++) {
        // sum of elements
        sx += x[i];
        sy += y[i];
        sxy += x[i] * y[i];
        // square sums
        ssx += pow(x[i], 2);
        ssy += pow(y[i], 2);
    }
    double nom = (double)n * sxy - sx * sy;
    double denom = (n * ssx - pow(sx, 2)) * (n * ssy - pow(sy, 2));
    return nom / sqrt(denom);
}


double spearman_pvalue(double coeff, unsigned n) {
    unsigned df = n - 2;
    double t = coeff * sqrt(df / (1 - pow(coeff, 2)));
    if (df < 1)
        return NAN;
    students_t dist(df);
    // two-sided p-value
    if (isnan(t))
        return NAN;
    return cdf(dist, -abs(t)) * 2;
}


pair<double, double> spearmans_rank(vector<double> &x, vector<double> &y) {
    vector<double> xclean, yclean, xrank, yrank;
    double coeff, pval;
    unsigned i;

    for (i = 0; i < x.size(); i++) {
        if (!(isnan(x[i]) | isnan(y[i]))) {
            xclean.push_back(x[i]);
            yclean.push_back(y[i]);
        }
    }

    xrank = to_ranks(xclean);
    yrank = to_ranks(yclean);
    coeff = spearman_coefficient(xrank, yrank);
    pval = spearman_pvalue(coeff, xrank.size());
    return {coeff, pval};
}


SpearmansResults spearmans_rank(vector<vector<double>>& x, vector<vector<double>>& y, int n_threads) {
    // NOTE: correlations are computed between the ROWS of x and y!
    vector<vector<double>> xranks(x.size(), vector<double>(x[0].size()));
    vector<vector<double>> yranks(y.size(), vector<double>(y[0].size()));

    pair<bool, vector<vector<bool>>> x_nan_masks = get_nan_masks(x);
    pair<bool, vector<vector<bool>>> y_nan_masks = get_nan_masks(y);

    if (x_nan_masks.first | y_nan_masks.first) {
        return spearmans_rank(
            x, y, x_nan_masks.second,
            y_nan_masks.second, n_threads
        );
    }

    unsigned i, j;
    for (i = 0; i < x.size(); i++) {
        xranks[i] = to_ranks(x[i]);
    }
    for (j = 0; j < y.size(); j++) {
        yranks[j] = to_ranks(y[j]);
    }

    size_t n = x[0].size();
    vector<vector<double>> coeffs(x.size(), vector<double>(y.size())), pvals(x.size(), vector<double>(y.size()));

#ifdef USE_OPENMP
    if (n_threads == 1) {
        double coeff;
        for (i = 0; i < x.size(); i++) {
            for (j = 0; j < y.size(); j++) {
                coeff = spearman_coefficient(xranks[i], yranks[j]);
                coeffs[i][j] = coeff;
                pvals[i][j] = spearman_pvalue(coeff, n);
            }
        }
    } else {
        #pragma omp parallel num_threads(n_threads)
        {
            int io, jo;
            double coeff;
            #pragma omp for nowait
            for (io = 0; io < (int)x.size(); io++) {
                for (jo = 0; jo < (int)y.size(); jo++) {
                    coeff = spearman_coefficient(xranks[io], yranks[jo]);
                    coeffs[io][jo] = coeff;
                    pvals[io][jo] = spearman_pvalue(coeff, n);
                }
            }
        }
    }
#else
    double coeff;
    for (i = 0; i < x.size(); i++) {
        for (j = 0; j < y.size(); j++) {
            coeff = spearman_coefficient(xranks[i], yranks[j]);
            coeffs[i][j] = coeff;
            pvals[i][j] = spearman_pvalue(coeff, n);
        }
    }
#endif

    return {coeffs, pvals};
}


// with nan computations
pair<double, double> spearmans_rank(
    vector<double>& x, vector<double>& y, vector<bool>& xnan, vector<bool>& ynan
) {
    vector<double> xclean, yclean, xrank, yrank;
    double coeff, pval;
    unsigned i;

    for (i = 0; i < x.size(); i++) {
        if (!(xnan[i] | ynan[i])) {
            xclean.push_back(x[i]);
            yclean.push_back(y[i]);
        }
    }

    xrank = to_ranks(xclean);
    yrank = to_ranks(yclean);
    coeff = spearman_coefficient(xrank, yrank);
    pval = spearman_pvalue(coeff, xrank.size());
    return {coeff, pval};
}


SpearmansResults spearmans_rank(
    vector<vector<double>>& x, vector<vector<double>>& y,
    vector<vector<bool>> &xnans, vector<vector<bool>> &ynans, int n_threads
) {
    // NOTE: correlations are computed between the ROWS of x and y!
    vector<vector<double>> coeffs(x.size(), vector<double>(y.size())), pvals(x.size(), vector<double>(y.size()));

#ifdef USE_OPENMP
    if (n_threads == 1) {
        unsigned i, j;
        pair<double, double> res;
        vector<double> xranks, yranks;
        for (i = 0; i < x.size(); i++) {
            for (j = 0; j < y.size(); j++) {
                res = spearmans_rank(x[i], y[j], xnans[i], ynans[j]);
                coeffs[i][j] = res.first;
                pvals[i][j] = res.second;
            }
        }
    } else {
        #pragma omp parallel num_threads(n_threads)
        {
            pair<double, double> res;
            int io, jo;
            #pragma omp for nowait
            for (io = 0; io < (int)x.size(); io++) {
                for (jo = 0; jo < (int)y.size(); jo++) {
                    res = spearmans_rank(x[io], y[jo], xnans[io], ynans[jo]);
                    coeffs[io][jo] = res.first;
                    pvals[io][jo] = res.second;
                }
            }
        }
    }
#else
    unsigned i, j;
    pair<double, double> res;
    vector<double> xranks, yranks;
    for (i = 0; i < x.size(); i++) {
        for (j = 0; j < y.size(); j++) {
            res = spearmans_rank(x[i], y[j], xnans[i], ynans[j]);
            coeffs[i][j] = res.first;
            pvals[i][j] = res.second;
        }
    }
#endif

    return {coeffs, pvals};
}


/* ===================
 * wrappers for python
 * ===================
 */
PyObject* spearmans_rank_1by1(PyObject* x, PyObject* y) {
    vector<double> xvec = list_to_vector(x);
    vector<double> yvec = list_to_vector(y);
    pair<double, double> res = spearmans_rank(xvec, yvec);

    PyObject* tup = PyTuple_New(2);
    PyTuple_SetItem(tup, 0, PyFloat_FromDouble(res.first));
    PyTuple_SetItem(tup, 1, PyFloat_FromDouble(res.second));

    return tup;
}

PyObject* spearmans_rank_2by2(PyObject* x, PyObject* y, int n_threads) {
    vector<vector<double>> xvec = nested_list_to_vector(x);
    vector<vector<double>> yvec = nested_list_to_vector(y);
    SpearmansResults cor = spearmans_rank(xvec, yvec, n_threads);

    PyObject* tup = PyTuple_New(2);
    PyTuple_SetItem(tup, 0, vector_to_list(cor.first));
    PyTuple_SetItem(tup, 1, vector_to_list(cor.second));

    return tup;
}
