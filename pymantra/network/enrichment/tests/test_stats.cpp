#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN

#include "../doctest.h"
#include "../statsutils.hpp"
#include <iostream>
#include <limits>


using std::cout;
using std::endl;

#define TOLERANCE 0.000009


template <typename T>
vector<vector<T>> transpose(vector<vector<T>>& x) {
    unsigned xi = x.size(), xj = x[0].size(), i, j;
    vector<vector<T>> xt (xj, vector<T>(xi));
    for (j = 0; j < xj; j++) {
        for (i = 0; i < xi; i++) {
            xt[j][i] = x[i][j];
        }
    }
    return xt;
}


template <typename T>
bool all_equal(vector<T> a, vector<T> b) {
    if (a.size() != b.size()) return false;
    int i;
    for (i = 0; i < a.size(); i++) {
        if (std::abs(a[i] - b[i]) > TOLERANCE)  {
            // cout << "mismatch: " << a[i] << " -- " << b[i] << endl;
            return false;
        }
    }
    return true;
}


TEST_SUITE_BEGIN("Statistical Functions");

TEST_CASE("Fisher's Exact Test") {
    GIVEN("Contingency Table") {
        const unsigned cont_table[4] = {1, 9, 11, 3};
        // this is the result from both the R as well as the scipy implementation
        double expected_result = 0.002759;
        double diff = fabs(expected_result - fishers_exact(cont_table));
        THEN("Test works correctly") {
            CHECK(diff < 1e-5);
        }
    }
}

TEST_CASE("Cumulative density function on gaussians") {
    vector<double> res = {
        -0.88088039,  0.33381562,  0.27940721,  1.07386361, -0.46342849,
        -1.88026051,  1.33370453, -2.04006419,  0.80943612, -0.50718711
    };
    // expected values are computed using scipy.norm.cdf(x)
    vector<double> expected_cdf_vals = {
        0.18919128, 0.63074065, 0.61003383, 0.85855809, 0.32152863,
        0.03003629, 0.90884964, 0.02067197, 0.79086783, 0.30601177
    };
    vector<double> cdf_vals = norm_cdf(res);

    CHECK(all_equal(expected_cdf_vals, cdf_vals));

}

TEST_CASE("Welch test") {
    vector<double> res_a = {
        -0.10444687,  0.20528081,  0.71117937, -0.22149706, -0.02231591,
        1.2024399 , -0.97986801, -0.58016235, -0.2708519 ,  0.06024202
    };
    vector<double> res_b = {
        -1.41016552, -0.08713225, -0.18770675,  0.30592808, -0.57377686,
        0.02058034,  0.9230752 ,  0.69495933, -1.16316011, -0.71472771
    };

    // expected values are computed using scipy.stats.ttest_ind(res_a, res_b)
    double expected_statistic = 0.709113;
    double statistic = welch_test(res_a, res_b, true);
    cout << statistic << endl;
    CHECK(std::abs(expected_statistic - statistic) < TOLERANCE);

    double expected_pval = 0.487338;
    double pval = welch_test(res_a, res_b);
    cout << pval << endl;
    CHECK(std::abs(expected_pval - pval) < TOLERANCE);
}


TEST_CASE("Welch test with nan values") {
    vector<double> res_a = {
        -0.10444687,  0.20528081,  0.71117937, -0.22149706, -0.02231591,
        1.2024399 , -0.97986801, -0.58016235, NAN,
        -0.2708519 ,  0.06024202
    };
    vector<double> res_b = {
        -1.41016552, -0.08713225, -0.18770675,  0.30592808, -0.57377686,
        0.02058034,  0.9230752 ,  0.69495933, -1.16316011, -0.71472771,
        NAN
    };

    // expected values are computed using scipy.stats.ttest_ind(res_a, res_b)
    double statistic = welch_test(res_a, res_b, true);
    cout << statistic << endl;

    double pval = welch_test(res_a, res_b);
    cout << pval << endl;
}


// TODO: test for spearman's functions
TEST_CASE("Spearmans's Correlation") {
    vector<double> x1, y1, x1n, y1n;
    SpearmansResults spear_res;
    vector<vector<double>> coeffs;
    double coeff, exp_coeff, p, exp_p;
    unsigned i;

    GIVEN("Two 1D vectors") {
        x1 = {
            -0.18256708,  1.81042663, -0.26167468, -1.69506019,  1.2356669 ,
            -2.42780505, -0.32866334, -0.03021432, -1.78031191, -0.70281443
        };
        y1 = {
            -0.60513949, -0.59049128,  0.38576032, -1.15587999, -1.09633867,
            -0.44189462,  0.28603914, -0.27911719, -0.07665268,  0.45849339
        };

        exp_coeff = -0.236363;
        coeff = spearman_coefficient(x1, y1);
        CHECK(std::abs(exp_coeff - coeff) < TOLERANCE);

        exp_p =  0.510885;
        p = spearman_pvalue(coeff, x1.size());
        CHECK(std::abs(exp_p - p) < TOLERANCE);
    }

    GIVEN("Two 1D vectors with nans") {
        x1 = {
            NAN,  1.81042663, -0.26167468, -1.69506019,  1.2356669,
            -2.42780505, NAN, -0.03021432, -1.78031191, -0.70281443
        };
        y1 = {
            NAN, -0.59049128,  0.38576032, -1.15587999, -1.09633867,
            -0.44189462,  0.28603914, -0.27911719, -0.07665268, NAN
        };

        // from scipy.stats.spearmanr
        exp_coeff = -0.25;
        coeff = spearman_coefficient(x1, y1);

        cout << coeff << endl;
        CHECK(std::abs(exp_coeff - coeff) < TOLERANCE);

        exp_p = 0.58872;
        p = spearman_pvalue(coeff, 7);
        CHECK(std::abs(exp_p - p) < TOLERANCE);
    }

    GIVEN("Two 2D vectors") {
        vector<vector<double>> x2{
            {1.37905698,  2.10038556,  0.16040252},
            {-0.75498863,  1.77267099, -1.53857276},
            {-0.10407954,  0.58776359,  0.37617236},
            {-0.41908092,  1.17346653,  1.16457266}
        };
        x2 = transpose(x2);
        vector<vector<double>> y2{
            { 0.35664017,  0.17592594, -2.23915764},
            { 0.31112835, -0.90583765, -0.67430889},
            {-1.66821649, -0.00911379, -0.10622856},
            {-0.77695438,  0.78963106,  1.35621344}
        };
        y2 = transpose(y2);

        vector<vector<double>> exp_coeffs{
            {0.2,  1.,  -0.6},
            {0.4,  0.,   0.8},
            {-0.4, -0.8,  0.8}
        };
        exp_coeffs = transpose(exp_coeffs);

        spear_res = spearmans_rank(x2, y2);
        coeffs = spear_res.first;

        for (i = 0; i < coeffs.size(); i++) {
            for (int j = 0; j < coeffs[1].size(); j++) {
                cout << "[" << coeffs[i][j] << "," << exp_coeffs[i][j] << "] ";
            }
            cout << endl;
        }
        cout << endl;

        for (i = 0; i < coeffs.size(); i++) {
            CHECK(all_equal(coeffs[i], exp_coeffs[i]));
        }

        // TODO: test
        // vector<vector<double>> exp_pvals{
        //    {0.8, 0.,  0.4,},
        //    {0.6, 1.,  0.2,},
        //    {0.6, 0.2, 0.2}
        // };
        // for (i = 0; i < coeffs.size(); i++) {
        //     CHECK(all_equal(pvals[i], pvals[i]));
        // }
    }

    GIVEN("Two 2D vectors with nans") {
        vector<vector<double>> x2n{
            {NAN,  2.10038556,  0.16040252},
            {-0.75498863,  NAN, -1.53857276},
            {-0.10407954,  0.58776359,  0.37617236},
            {-0.41908092,  1.17346653,  1.16457266}
        };
        x2n = transpose(x2n);

        vector<vector<double>> y2n{
            {NAN,  0.17592594, NAN},
            {0.31112835, -0.90583765, -0.67430889},
            {-1.66821649, -0.00911379, -0.10622856},
            {-0.77695438,  0.78963106,  1.35621344}
        };
        y2n = transpose(y2n);

        vector<vector<double>> exp_coeffs{
            {-1.,   1.,  -0.5},
            {0.5,  0.5,  0.8},
            {0.5,  1.,   1.}
        };
        exp_coeffs = transpose(exp_coeffs);

        spear_res = spearmans_rank(x2n, y2n);

        for (i = 0; i < spear_res.first.size(); i++) {
            for (int j = 0; j < spear_res.first[1].size(); j++) {
                cout << "[" << spear_res.first[i][j] << "," << exp_coeffs[i][j] << "] ";
            }
            cout << endl;
        }
        cout << endl;

        // TODO: add in checks for p-values
        for (i = 0; i < spear_res.first.size(); i++) {
            CHECK(all_equal(spear_res.first[i], exp_coeffs[i]));
        }
        // {0.,                nan 0.66666,}
        // {0.66666667, 0.66666667 0.2,       }
        // {0.66666667,        nan 0.,        }
    }
}

TEST_CASE("Spearmans's Correlation with OpenMP") {
    vector<vector<double>> coeffs;
    SpearmansResults spear_res;
    unsigned i;

    GIVEN("Two 2D vectors") {
        vector<vector<double>> x2{
                {1.37905698,  2.10038556,  0.16040252},
                {-0.75498863,  1.77267099, -1.53857276},
                {-0.10407954,  0.58776359,  0.37617236},
                {-0.41908092,  1.17346653,  1.16457266}
        };
        x2 = transpose(x2);
        vector<vector<double>> y2{
                { 0.35664017,  0.17592594, -2.23915764},
                { 0.31112835, -0.90583765, -0.67430889},
                {-1.66821649, -0.00911379, -0.10622856},
                {-0.77695438,  0.78963106,  1.35621344}
        };
        y2 = transpose(y2);

        vector<vector<double>> exp_coeffs{
                {0.2,  1.,  -0.6},
                {0.4,  0.,   0.8},
                {-0.4, -0.8,  0.8}
        };
        exp_coeffs = transpose(exp_coeffs);

        spear_res = spearmans_rank(x2, y2, 4);
        coeffs = spear_res.first;

        for (i = 0; i < coeffs.size(); i++) {
            for (int j = 0; j < coeffs[1].size(); j++) {
                cout << "[" << coeffs[i][j] << "," << exp_coeffs[i][j] << "] ";
            }
            cout << endl;
        }
        cout << endl;

        for (i = 0; i < coeffs.size(); i++) {
            CHECK(all_equal(coeffs[i], exp_coeffs[i]));
        }
    }

    GIVEN("Two 2D vectors with nans") {
        vector<vector<double>> x2n{
                {NAN,         2.10038556, 0.16040252},
                {-0.75498863, NAN,        -1.53857276},
                {-0.10407954, 0.58776359, 0.37617236},
                {-0.41908092, 1.17346653, 1.16457266}
        };
        x2n = transpose(x2n);

        vector<vector<double>> y2n{
                {NAN,         0.17592594, NAN},
                {0.31112835,  -0.90583765, -0.67430889},
                {-1.66821649, -0.00911379, -0.10622856},
                {-0.77695438, 0.78963106,  1.35621344}
        };
        y2n = transpose(y2n);

        vector<vector<double>> exp_coeffs{
                {-1., 1.,  -0.5},
                {0.5, 0.5, 0.8},
                {0.5, 1.,  1.}
        };
        exp_coeffs = transpose(exp_coeffs);

        spear_res = spearmans_rank(x2n, y2n, 4);
        coeffs = spear_res.first;

        for (i = 0; i < coeffs.size(); i++) {
            for (int j = 0; j < coeffs[1].size(); j++) {
                cout << "[" << coeffs[i][j] << "," << exp_coeffs[i][j] << "] ";
            }
            cout << endl;
        }
        cout << endl;

        for (i = 0; i < coeffs.size(); i++) {
            CHECK(all_equal(coeffs[i], exp_coeffs[i]));
        }
    }
}

