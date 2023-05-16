#include <vector>
#include <iostream>

#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/students_t.hpp>
#include <boost/math/statistics/t_test.hpp>

using std::vector;
using std::cout;
using std::endl;

using boost::math::cdf;
using boost::math::students_t;
using boost::math::statistics::two_sample_t_test;


int main() {
	vector<double> x = {.1, .2, .3, .4, .5};
	vector<double> y = {.3, .4, .5, .6, .7};
	std::pair<double, double> test_res = two_sample_t_test(x, y);
	cout << "t-stats: " << test_res.first << " ";
	cout << "pval: " << test_res.second << endl;


	boost::math::normal norm_dist;
	double z = 1.9;
	double pn = boost::math::cdf(norm_dist, z);
	cout << "normal cdf: " << pn << endl;

	boost::math::students_t t_dist(x.size() - 1);
	double tn = boost::math::cdf(t_dist, z);
	cout << "student's cdf: " << tn << endl;
}

