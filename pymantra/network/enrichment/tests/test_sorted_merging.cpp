#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN


#include "../doctest.h"
#include <iostream>
#include "../LSO/LocalSearch.hpp"


using std::cout;
using std::endl;


void print_vector(vector<double> vec) {
    int i;
    for (i = 0; i < vec.size(); i++) {
        cout << vec[i] << ", " << endl;
    }
}


void print_vector(vector<BestAction> vec) {
    int i;
    for (i = 0; i < vec.size(); i++) {
        cout << vec[i].score << ", " << endl;
    }
}


template<typename T>
void report_false_sorting(vector<T> computed, vector<T> expected) {
    cout << "Expected:" << endl;
    print_vector(expected);
    cout << "Computed:" << endl;
    print_vector(computed);
}


bool equals(double a, double b) {
    return a == b;
}


bool equals(BestAction& a, BestAction& b) {
    return a.score == b.score;
}


template<typename T>
bool merge_results_sorted(vector<T> computed, vector<T> expected) {
    bool eq = true;
    if (computed.size() != expected.size()) {
        cout << R"(\033[1;31m Something went completely wrong, resulting vectors have different sizes!\033[0m)" << endl;
        report_false_sorting(computed, expected);
        eq = false;
    } else {
        int i;
        for (i = 0; i < computed.size(); i++) {
            if (!equals(computed[i], expected[i])) {
                eq = false;
                break;
            }
        }
    }
    if (!eq) {
        report_false_sorting(computed, expected);
    }
    return eq;
}

TEST_SUITE_BEGIN("sorted merging");

TEST_CASE("merging") {
    GIVEN("double vectors") {
        vector<double> a = {1., 4., 7.};
        vector<double> b = {3., 5., 6.};
        vector<double> c = {2., 8., 9.};
        vector<double> merged = {1., 2., 3., 4., 5., 6., 7., 8., 9.};

        vector<double> computed_merge = LocalSearch::merge_solutions(a, b, c, std::less<>());
        CHECK(merge_results_sorted(computed_merge, merged));
    }
    GIVEN("BestAction vectors") {
        vector<BestAction> a = {
            BestAction(2., Action::no_move),
            BestAction(1., Action::no_move)
        };
        vector<BestAction> b = {
            BestAction(4., Action::no_move),
            BestAction(3., Action::no_move)
        };
        vector<BestAction> c = {
            BestAction(6., Action::no_move),
            BestAction(5., Action::no_move)
        };
        vector<BestAction> merged = {
            BestAction(6., Action::no_move),
            BestAction(5., Action::no_move),
            BestAction(4., Action::no_move),
            BestAction(3., Action::no_move),
            BestAction(2., Action::no_move),
            BestAction(1., Action::no_move)
        };

        vector<BestAction> computed_merge = LocalSearch::merge_solutions(a, b, c, action_greater);
        CHECK(merge_results_sorted(computed_merge, merged));
    }
}
