#include <omp.h>
#include <vector>


using std::vector;


int main() {

	vector<vector<int>> x(3, vector<int>(4));
#pragma omp parallel num_threads(4)
	{
		int i, j;
		int xij;
		#pragma omp for nowait
		for (i = 0; i < x.size(); i++) {
			for (j = 0; j < x[0].size(); j++) {
				xij = i * (j + 1);
				x[i][j] = xij;
			}
		}
	}
}
