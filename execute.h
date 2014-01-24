
#include <iostream>
#include <vector>
#include <string>
#include "hoge.h"
#include "ucbs.h"
using namespace std;

namespace exe{

vector<vector<int> > Y; // input data
int len;

void data_setup(const char* input_name, int M, int N){
	vector<vector<int> > tempD = hoge::importData(input_name);
	vector<vector<int> > tempY = hoge::Data2UIMatrix(tempD);
	Y = hoge::top_N_item(hoge::top_N_user(tempY, M), N);
	len = N;
}

void getResults(int times, double* eps, double K, const char* output_name){
	vector<int> arr;
	for(int i = 0; i < len; ++i){
		arr.push_back(i);
	}

	int n = sizeof(eps);
	double Acc1[n];
	double Acc2[n];
	double Acc3[n];
	double Chosen1[n];
	double Chosen2[n];
	double Chosen3[n];
	double Num1[n];
	double Num2[n];
	double Num3[n];
	for(int i = 0; i < n; ++i){
		Acc1[i] = 0;
		Acc2[i] = 0;
		Acc3[i] = 0;
		Chosen1[i] = 0;
		Chosen2[i] = 0;
		Chosen3[i] = 0;
		Num1[i] = 0;
		Num2[i] = 0;
		Num3[i] = 0;
	}
	for(int t = 0; t < times; ++t){
		// using built-in random generator:
		random_shuffle(arr.begin(), arr.end());
		// using myrandom:
		random_shuffle(arr.begin(), arr.end(), hoge::myrandom);

		// shuffle columns of Y
		vector<vector<int> > Z = hoge::submatrix_column(Y, arr);

		for(int i = 0; i < n; ++i){
			ucbs::construct(Z, eps[i], K, len);
			ucbs::init();
			ucbs::setResults_nocomp();
			Acc1[i] += ucbs::getAcc();
			Chosen1[i] += ucbs::getChosen();
			Num1[i] += ucbs::getNum();
			ucbs::all_clear();

			ucbs::construct(Z, eps[i], K, len);
			ucbs::init();
			ucbs::setResults_withGL_ver1();
			Acc2[i] += ucbs::getAcc();
			Chosen2[i] += ucbs::getChosen();
			Num2[i] += ucbs::getNum();
			ucbs::all_clear();

			ucbs::construct(Z, eps[i], K, len);
			ucbs::init();
			ucbs::setResults_withGL_ver2();
			Acc3[i] += ucbs::getAcc();
			Chosen3[i] += ucbs::getChosen();
			Num3[i] += ucbs::getNum();
			ucbs::all_clear();

			cout << t << "	" << i << endl;
		}
	}

	for(int i = 0; i < n; ++i){
		Acc1[i] /= times;
		Acc2[i] /= times;
		Acc3[i] /= times;
		Chosen1[i] /= times;
		Chosen2[i] /= times;
		Chosen3[i] /= times;
		Num1[i] /= times;
		Num2[i] /= times;
		Num3[i] /= times;
	}

	vector<string> Strs1, Strs2, Strs3;
	for(int i = 0; i < n; ++i){
		string Acc1_str = hoge::double2string(Acc1[i]);
		string Acc2_str = hoge::double2string(Acc2[i]);
		string Acc3_str = hoge::double2string(Acc3[i]);
		string Chosen1_str = hoge::double2string(Chosen1[i]);
		string Chosen2_str = hoge::double2string(Chosen2[i]);
		string Chosen3_str = hoge::double2string(Chosen3[i]);
		string Num1_str = hoge::double2string(Num1[i]);
		string Num2_str = hoge::double2string(Num2[i]);
		string Num3_str = hoge::double2string(Num3[i]);

		string str1 = hoge::double2string(eps[i])+","+Acc1_str+","+Chosen1_str+","+Num1_str;
		string str2 = hoge::double2string(eps[i])+","+Acc2_str+","+Chosen2_str+","+Num2_str;
		string str3 = hoge::double2string(eps[i])+","+Acc3_str+","+Chosen3_str+","+Num3_str;
		Strs1.push_back(str1);
		Strs2.push_back(str2);
		Strs3.push_back(str3);
	}

	const char* file1 = "c:/users/mistpc/results/result1_1.csv";
	const char* file2 = "c:/users/mistpc/results/result2_2.csv";
	const char* file3 = "c:/users/mistpc/results/result3_3.csv";
	hoge::exportData(file1, Strs1);
	hoge::exportData(file2, Strs2);
	hoge::exportData(file3, Strs3);
}

}
