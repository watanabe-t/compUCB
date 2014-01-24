
// -------------- compUCB.cpp --------------

#include <iostream>
#include "execute.h"
using namespace std;

int main() {
	const char* dataname = "data.csv";
	int M = 100;
	int N = 500;

	int times = 1;
	double eps[4] = {0.99, 0.95, 0.90, 0.85};
	double K = 1;
	const char* output_name = "c:/users/mistpc/aaa.csv";

	exe::data_setup(dataname, M, N);
	exe::getResults(times, eps, K, output_name);

	return 0;
}
