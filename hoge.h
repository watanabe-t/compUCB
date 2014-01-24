
// Useful functions for UCB1-Thresholding

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
using namespace std;

namespace hoge{

string int2string(int a){
	stringstream ss;
	ss << a;
	return ss.str();
}

string double2string(double a){
	stringstream ss;
	ss << a;
	return ss.str();
}

// compute the inner production of two vectors u and v
double innerprod(vector<double> u, vector<double> v){
	if(u.size() == 0 || v.size() == 0 || u.size() != v.size()){
		cout << "invalid operation" << endl;
		return 0;
	}else{
		double sum = 0;
		for(unsigned int i = 0; i < u.size(); ++i){
			sum += u[i]*v[i];
		}
		return sum;
	}
}

// compute the 2-norm of a vector u
double getNorm(vector<double> u){
	double res = 0;
	for(unsigned int i = 0; i < u.size(); ++i){
		res += u[i]*u[i];
	}
	return sqrt(res);
}

// compute the mean value of a vector u
double getMean(vector<double> u){
	double res = 0;
	for(unsigned int i = 0; i < u.size(); ++i){
		res += u[i];
	}
	return res/u.size();
}

// compute the max value of a vector u
double getMax(vector<double> u){
	double res = -infinity();
	for(unsigned int i = 0; i < u.size(); ++i){
		if(res < u[i]){
			res = u[i];
		}
	}
	return res;
}

// compute intersection of A and B (A\cap B)
vector<int> intersect(vector<int> a, vector<int> b){
	vector<int> v(a.size()+b.size());
	vector<int>::iterator it;

	sort(a.begin(), a.end());
	sort(b.begin(), b.end());

	it = set_intersection(a.begin(), a.end(), b.begin(), b.end(), v.begin());
	v.resize(it-v.begin());
	return v;
}

// compute set difference of A and B (A\B)
vector<int> setdiff(vector<int> a, vector<int> b){
	vector<int> v(a.size()+b.size());
	vector<int>::iterator it;

	sort(a.begin(), a.end());
	sort(b.begin(), b.end());

	it = set_difference(a.begin(), a.end(), b.begin(), b.end(), v.begin());
	v.resize(it-v.begin());

	return v;
}

// read csv file to make a data list
vector<vector<int> > importData(const char* filename){
	ifstream ifs(filename, ifstream::in);
	vector<vector<string> > values;
	string str;
	int p;

	if(!ifs){
		cout << "Error:Input data file not found" << endl;
		vector<vector<int> > emp(0);
		return emp;
	}

	while(getline(ifs, str)){
		if((p = str.find("//")) != str.npos){
			continue;
		}
		vector<string> inner;

		while((p = str.find(",")) != str.npos){
			inner.push_back(str.substr(0, p));
			str = str.substr(p+1);
		}

		inner.push_back(str);
		values.push_back(inner);
	}

	vector<vector<int> > Y(values.size());
	for(unsigned int i = 0; i < values.size(); ++i){
		Y[i].reserve(values[i].size());
		for(unsigned int j = 0; j < values[i].size(); ++j){
			istringstream iss(values[i][j]);
			int num = 0;
			iss >> num;
			Y[i].push_back(num);
		}
	}

	return Y;
}

// make user-item matrix from the data
vector<vector<int> > Data2UIMatrix(vector<vector<int> > D){
	int maxuser = -1;
	int maxitem = -1;
	for(unsigned int i = 0; i < D.size(); ++i){
		if(D[i][0] > maxuser){
			maxuser = D[i][0];
		}
		if(D[i][1] > maxitem){
			maxitem = D[i][1];
		}
	}

	vector<vector<int> > Y(maxuser);
	for(int i = 0; i < maxuser; ++i){
		Y[i].reserve(maxitem);
	}

	for(int i = 0; i < maxuser; ++i){
		for(int j = 0; j < maxitem; ++j){
			Y[i].push_back(0);;
		}
	}

	for(unsigned int i = 0; i < D.size(); ++i){
		if(D[i][2] > 2){
			Y[D[i][0]-1][D[i][1]-1] = 1;
		}else{
			Y[D[i][0]-1][D[i][1]-1] = -1;
		}
	}

	return Y;
}

// compute the sub-matrix of X (sub-row)
vector<vector<int> > submatrix_row(vector<vector<int> > X, vector<int> arr){
	if(arr.size() > X.size()){
		cout << "Array size over row size!!" << endl;
		return X;
	}
	else{
		sort(arr.begin(), arr.end());
		vector<vector<int> > Y(arr.size());
		for(unsigned int i = 0; i < arr.size(); ++i){
			for(unsigned int j = 0; j < X[arr[i]].size(); ++j){
				Y[i].push_back(X[arr[i]][j]);
			}
		}

		return Y;
	}
}

// compute the sub-matrix of X (sub-column)
vector<vector<int> > submatrix_column(vector<vector<int> > X, vector<int> arr){
	if(arr.size() > X[0].size()){
		cout << "Array size over column size!!" << endl;
		return X;
	}
	else{
		vector<vector<int> > Y(X.size());
		for(unsigned int i = 0; i < X.size(); ++i){
			for(unsigned int j = 0; j < arr.size(); ++j){
				Y[i].push_back(X[i][arr[j]]);
			}
		}

		return Y;
	}
}

// make a text file
void exportData(const char* filename, vector<string> Strs){
	ofstream ofs(filename);
	vector<string>::iterator it;
	for(it = Strs.begin(); it != Strs.end(); ++it){
		ofs << (*it) << endl;
	}
	ofs.close();

	return;
}

// generate a random integer from 0 to a-1
int myrandom(int a){
	return rand() % a;
}

// cut off workers who did not label a lot (top-N workers are chosen)
vector<vector<int> > top_N_user(vector<vector<int> > Y, unsigned int N){
	if(N > Y.size()){
		cout << "N is larger than the row size of Y." << endl;
		return Y;
	}

	vector<int> Sums;
	for(unsigned int i = 0; i < Y.size(); ++i){
		int sum = 0;
		for(unsigned int t = 0; t < Y[i].size(); ++t){
			if(Y[i][t] != 0){
				sum++;
			}
		}
		Sums.push_back(sum);
	}

	multimap<int, int> mMap;
	for(unsigned int i = 0; i < Sums.size(); ++i){
		mMap.insert(multimap<int, int>::value_type(Sums[i], i));
	}

	vector<int> array;
	multimap<int, int>::iterator it = --mMap.end();
	for(unsigned int i = 0; i < N; ++i){
		array.push_back((*it).second);
		--it;
	}

	sort(array.begin(), array.end());
	return submatrix_row(Y, array);
}

// cut off items with few workers' label (top-N items are chosen)
vector<vector<int> > top_N_item(vector<vector<int> > Y, unsigned int N){
	if(N > Y[0].size()){
		cout << "N is larger than the row size of Y." << endl;
		return Y;
	}

	vector<int> Sums;
	for(unsigned int t = 0; t < Y[0].size(); ++t){
		int sum = 0;
		for(unsigned int i = 0; i < Y.size(); ++i){
			if(Y[i][t] != 0){
				sum++;
			}
		}
		Sums.push_back(sum);
	}

	multimap<int, int> mMap;
	for(unsigned int t = 0; t < Sums.size(); ++t){
		mMap.insert(multimap<int, int>::value_type(Sums[t], t));
	}

	vector<int> array;
	multimap<int, int>::iterator it = --mMap.end();
	for(unsigned int i = 0; i < N; ++i){
		array.push_back((*it).second);
		--it;
	}

	sort(array.begin(), array.end());
	return submatrix_column(Y, array);

}

}
