
#include <iostream>
#include <vector>
#include <cmath>
using namespace std;

namespace ucbs{

vector<vector<int> > UserList, ItemList;
vector<vector<double> > Y;
vector<vector<int> >isObserved;

// outputs
vector<double> Accs;
vector<int> Chosens;
int Num = 0;

// algorithm parameters
double eps, K;
int len;

void construct(vector<vector<int> > Y0, double eps0, double K0, int len0){
	eps = eps0;
	K = K0;
	len = len0;

	vector<vector<double> > tempX(Y0.size());
	for(unsigned int i = 0; i < Y0.size(); ++i){
		for(unsigned int j = 0; j < Y0[0].size(); ++j){
			tempX[i].push_back((double)Y0[i][j]);
		}
	}
	Y = tempX;
}

void init(){
	vector<vector<int> > tempU(Y[0].size());
	for(unsigned int t = 0; t < Y[0].size(); ++t){
		for(unsigned int i = 0; i < Y.size(); ++i){
			if(Y[i][t] != 0){
				tempU[t].push_back(i);
			}
		}
	}
	UserList = tempU;

	vector<vector<int> > tempI(Y.size());
	for(unsigned int i = 0; i < Y.size(); ++i){
		for(unsigned int t = 0; t < Y[i].size(); ++t){
			if(Y[i][t] != 0){
				tempI[i].push_back(t);
			}
		}
	}
	ItemList = tempI;

	vector<vector<int> > tempZ(Y.size());
	for(unsigned int i = 0; i < Y.size(); ++i){
		for(unsigned int t = 0; t < Y[i].size(); ++t){
			tempZ[i].push_back(0);
		}
	}
	isObserved = tempZ;
}

double getAcc(){
	return Accs[Accs.size()-1];
}

int getChosen(){
	int sum = 0;
	for(unsigned int i = 0; i < Chosens.size(); ++i){
		sum += Chosens[i];
	}
	return sum;
}

int getNum(){
	return Num;
}

// set results (without completing missing labels)
void setResults_nocomp(){
	int m = Y.size();
	vector<double> tempAccs;
	vector<int> tempChosens;

	vector<vector<double> > Rew(m); // reward of each worker
	for(int i = 0; i < m; ++i){
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(0);
		}
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(1);
		}
	}

	vector<double> UCB; // upper confidence bound
	for(int i = 0; i < m; ++i){
		UCB.push_back(0);
	}

	int labeled = 0;
	int count = 0;

	// main loop
	for(int t = 0; t < len; ++t){
		// compute the upper confidence bound
		for(int i = 0; i < m; ++i){
			if(Rew[i].empty()){
				UCB[i] = 0;
			}else{
				UCB[i] = hoge::getMean(Rew[i]) + sqrt(2*log(1+t)/Rew[i].size());
			}
		}

		// choose workers who have large UCB
		double M = hoge::getMax(UCB);
		vector<int> Sel;
		for(int i = 0; i < m; ++i){
			if(UCB[i] >= eps*M){
				Sel.push_back(i);
			}
		}

		vector<int> L = hoge::intersect(UserList[t], Sel); // exclude workers who did not label
		if(L.empty()){
			continue;
		}

		labeled++;
		// compute the majority vote of L
		double score = 0;
		for(unsigned int i = 0; i < L.size(); ++i){
			score += Y[L[i]][t];
		}

		int vote;
		if(score > 0){
			vote = 1;
		}else if(score < 0){
			vote = -1;
		}else{
			vote = hoge::myrandom(2)*2 - 1;
		}


		// compute the majority vote of U_t
		score = 0;
		for(unsigned int i = 0; i < UserList[t].size(); ++i){
			score += Y[UserList[t][i]][t];
		}
		int entire;
		if(score > 0){
			entire = 1;
		}else if(score < 0){
			entire = -1;
		}else{
			entire = hoge::myrandom(2)*2 - 1;
		}

		if(vote == entire){
			count++;
		}

		tempAccs.push_back((double)count/labeled);
		tempChosens.push_back(L.size());

		// update rewards
		for(unsigned int i = 0; i < L.size(); ++i){
			double rew_i = 1 - 0.5*fabs(Y[L[i]][t] - vote);
			Rew[L[i]].push_back(rew_i);
		}
	}

	Accs = tempAccs;
	Chosens = tempChosens;
	Num = labeled;
}

// set results (with GroupLens and [observed labels + predicted labels])
void setResults_withGL_ver1(){
	int m = Y.size();
	vector<double> tempAccs;
	vector<int> tempChosens;

	vector<vector<double> > Rew(m); // reward of each worker
	for(int i = 0; i < m; ++i){
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(0);
		}
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(1);
		}
	}

	vector<double> UCB; // upper confidence bound
	for(int i = 0; i < m; ++i){
		UCB.push_back(0);
	}

	vector<vector<double> > R(m);
	for(int i = 0; i < m; ++i){
		for(int j = 0; j < m; ++j){
			R[i].push_back(0);
		}
	}

	int labeled = 0;
	int count = 0;

	// main loop
	for(int t = 0; t < len; ++t){
		// complete missing labels
		if(t != 0){
			// compute the similarity
			for(int i = 0; i < m; ++i){
				for(int j = i; j < m; ++j){
					vector<int> I_i, I_j;
					for(int s = 0; s < t; ++s){
						if(isObserved[i][s] != 0){
							I_i.push_back(s);
						}
						if(isObserved[j][s] != 0){
							I_j.push_back(s);
						}
					}
					vector<int> I_set = hoge::intersect(I_i, I_j);
					if(I_set.empty()){
						continue;
					}

					double m_i = 0;
					double m_j = 0;
					int l_i = 0;
					int l_j = 0;
					for(int s = 0; s < t; ++s){
						if(isObserved[i][s] != 0){
							m_i += Y[i][s];
							l_i++;
						}
						if(isObserved[j][s] != 0){
							m_j += Y[j][s];
							l_j++;
						}
					}
					if(l_i > 0){
						m_i /= l_i;
					}
					if(l_j > 0){
						m_j /= l_j;
					}

					vector<double> v_i, v_j;
					for(unsigned int s = 0; s < I_set.size(); ++s){
						v_i.push_back(Y[i][I_set[s]] - m_i);
						v_j.push_back(Y[j][I_set[s]] - m_j);
					}

					R[i][j] = hoge::innerprod(v_i, v_j)/hoge::getNorm(v_i)/hoge::getNorm(v_j);
					if(hoge::getNorm(v_i) == 0 || hoge::getNorm(v_j) == 0){
						R[i][j] = 0;
					}
					R[j][i] = R[i][j];
				}
			}

			// predict labels
			for(int s = 0; s < t; ++s){
				vector<int> W, W_s;
				for(int i = 0; i < m; ++i){
					if(isObserved[i][s] != 0){
						W_s.push_back(i);
					}
					W.push_back(i);
				}

				vector<int> compW_s = hoge::setdiff(W, W_s);
				for(unsigned int i = 0; i < compW_s.size(); ++i){
					double sum_i = 0;
					double div_i = 0;
					for(unsigned int j = 0; j < W_s.size(); ++j){
						sum_i += R[compW_s[i]][W_s[j]]*Y[W_s[j]][s];
						div_i += fabs(R[compW_s[i]][W_s[j]]);
					}
					if(div_i == 0){
						sum_i = 0;
					}else{
						sum_i /= div_i;
					}

					double m_i = 0;
					int l_i = 0;
					for(int u = 0; u < t; ++u){
						if(isObserved[compW_s[i]][u] != 0){
							m_i += Y[compW_s[i]][u];
							l_i++;
						}
					}
					if(l_i == 0){
						m_i = 0;
					}else{
						m_i /= l_i;
					}

					Y[compW_s[i]][s] = m_i + sum_i;
				}
			}
		}

		// compute each worker's reward & compute upper confidence bound
		vector<double> Votes;
		for(int s = 0; s < t; ++s){
			double m_s = 0;
			for(int i = 0; i < m; ++i){
				m_s += Y[i][s];
			}
			int val;
			if(m_s > 0){
				val = 1;
			}else if(m_s < 0){
				val = -1;
			}else{
				val = hoge::myrandom(2)*2 - 1;
			}
			Votes.push_back(val);
		}
		for(int i = 0; i < m; ++i){
			if(t == 0){
				UCB[i] = hoge::getMean(Rew[i]);
			}else{
				vector<double> temp = Rew[i];
				for(int s = 0; s < t; ++s){
					temp.push_back(1 - fabs(Y[i][s] - Votes[s]));
				}
				UCB[i] = hoge::getMean(temp) + sqrt(2*log(1+t)/t);
			}
		}

		// choose workers
		double M = hoge::getMax(UCB);
		vector<int> Sel;
		for(int i = 0; i < m; ++i){
			if(UCB[i] >= eps*M){
				Sel.push_back(i);
			}
		}
		vector<int> L = hoge::intersect(Sel, UserList[t]);
		if(L.empty()){
			continue;
		}
		for(unsigned int i = 0; i < L.size(); ++i){
			isObserved[L[i]][t] = 1;
		}
		labeled++;

		// compute the majority vote of L
		double score = 0;
		for(unsigned int i = 0; i < L.size(); ++i){
			score += Y[L[i]][t];
		}

		int vote;
		if(score > 0){
			vote = 1;
		}else if(score < 0){
			vote = -1;
		}else{
			vote = hoge::myrandom(2)*2 - 1;
		}

		// compute the majority vote of U_t
		score = 0;
		for(unsigned int i = 0; i < UserList[t].size(); ++i){
			score += Y[UserList[t][i]][t];
		}
		int entire;
		if(score > 0){
			entire = 1;
		}else if(score < 0){
			entire = -1;
		}else{
			entire = hoge::myrandom(2)*2 - 1;
		}

		if(vote == entire){
			count++;
		}

		tempAccs.push_back((double)count/labeled);
		tempChosens.push_back(L.size());
	}

	Accs = tempAccs;
	Chosens = tempChosens;
	Num = labeled;
}

// set results (with GroupLens and observed labels only)
void setResults_withGL_ver2(){
	int m = Y.size();
	vector<double> tempAccs;
	vector<int> tempChosens;

	vector<vector<double> > Rew(m); // reward of each worker
	for(int i = 0; i < m; ++i){
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(0);
		}
		for(int times = 0; times < K; ++times){
			Rew[i].push_back(1);
		}
	}

	vector<double> UCB; // upper confidence bound
	for(int i = 0; i < m; ++i){
		UCB.push_back(0);
	}

	vector<vector<double> > R(m);
	for(int i = 0; i < m; ++i){
		for(int j = 0; j < m; ++j){
			R[i].push_back(0);
		}
	}

	int labeled = 0;
	int count = 0;

	// main loop
	for(int t = 0; t < len; ++t){
		// complete missing labels
		if(t != 0){
			// compute the similarity
			for(int i = 0; i < m; ++i){
				for(int j = i; j < m; ++j){
					vector<int> I_i, I_j;
					for(int s = 0; s < t; ++s){
						if(isObserved[i][s] != 0){
							I_i.push_back(s);
						}
						if(isObserved[j][s] != 0){
							I_j.push_back(s);
						}
					}
					vector<int> I_set = hoge::intersect(I_i, I_j);
					if(I_set.empty()){
						continue;
					}

					double m_i = 0;
					double m_j = 0;
					int l_i = 0;
					int l_j = 0;
					for(int s = 0; s < t; ++s){
						if(isObserved[i][s] != 0){
							m_i += Y[i][s];
							l_i++;
						}
						if(isObserved[j][s] != 0){
							m_j += Y[j][s];
							l_j++;
						}
					}
					if(l_i > 0){
						m_i /= l_i;
					}
					if(l_j > 0){
						m_j /= l_j;
					}

					vector<double> v_i, v_j;
					for(unsigned int s = 0; s < I_set.size(); ++s){
						v_i.push_back(Y[i][I_set[s]] - m_i);
						v_j.push_back(Y[j][I_set[s]] - m_j);
					}

					R[i][j] = hoge::innerprod(v_i, v_j)/hoge::getNorm(v_i)/hoge::getNorm(v_j);
					if(hoge::getNorm(v_i) == 0 || hoge::getNorm(v_j) == 0){
						R[i][j] = 0;
					}
					R[j][i] = R[i][j];
				}
			}

			// predict labels
			for(int s = 0; s < t; ++s){
				vector<int> W, W_s;
				for(int i = 0; i < m; ++i){
					if(isObserved[i][s] != 0){
						W_s.push_back(i);
					}
					W.push_back(i);
				}

				vector<int> compW_s = hoge::setdiff(W, W_s);
				for(unsigned int i = 0; i < compW_s.size(); ++i){
					double sum_i = 0;
					double div_i = 0;
					for(unsigned int j = 0; j < W_s.size(); ++j){
						sum_i += R[compW_s[i]][W_s[j]]*Y[W_s[j]][s];
						div_i += fabs(R[compW_s[i]][W_s[j]]);
					}
					if(div_i == 0){
						sum_i = 0;
					}else{
						sum_i /= div_i;
					}

					double m_i = 0;
					int l_i = 0;
					for(int u = 0; u < t; ++u){
						if(isObserved[compW_s[i]][u] != 0){
							m_i += Y[compW_s[i]][u];
							l_i++;
						}
					}
					if(l_i == 0){
						m_i = 0;
					}else{
						m_i /= l_i;
					}

					Y[compW_s[i]][s] = m_i + sum_i;
				}
			}
		}

		// compute each worker's reward & compute upper confidence bound
		vector<double> Votes;
		for(int s = 0; s < t; ++s){
			double m_s = 0;
			for(unsigned int i = 0; i < UserList[s].size(); ++i){
				m_s += Y[UserList[s][i]][s];
			}
			int val;
			if(m_s > 0){
				val = 1;
			}else if(m_s < 0){
				val = -1;
			}else{
				val = hoge::myrandom(2)*2 - 1;
			}
			Votes.push_back(val);
		}
		for(int i = 0; i < m; ++i){
			if(t == 0){
				UCB[i] = hoge::getMean(Rew[i]);
			}else{
				vector<double> temp = Rew[i];
				for(int s = 0; s < t; ++s){
					if(isObserved[i][s] != 0){
						temp.push_back(1 - fabs(Y[i][s] - Votes[s]));
					}
				}
				UCB[i] = hoge::getMean(temp) + sqrt(2*log(1+t)/temp.size());
			}
		}

		// choose workers
		double M = hoge::getMax(UCB);
		vector<int> Sel;
		for(int i = 0; i < m; ++i){
			if(UCB[i] >= eps*M){
				Sel.push_back(i);
			}
		}
		vector<int> L = hoge::intersect(Sel, UserList[t]);
		if(L.empty()){
			continue;
		}
		for(unsigned int i = 0; i < L.size(); ++i){
			isObserved[L[i]][t] = 1;
		}
		labeled++;

		// compute the majority vote of L
		double score = 0;
		for(unsigned int i = 0; i < L.size(); ++i){
			score += Y[L[i]][t];
		}

		int vote;
		if(score > 0){
			vote = 1;
		}else if(score < 0){
			vote = -1;
		}else{
			vote = hoge::myrandom(2)*2 - 1;
		}

		// compute the majority vote of U_t
		score = 0;
		for(unsigned int i = 0; i < UserList[t].size(); ++i){
			score += Y[UserList[t][i]][t];
		}
		int entire;
		if(score > 0){
			entire = 1;
		}else if(score < 0){
			entire = -1;
		}else{
			entire = hoge::myrandom(2)*2 - 1;
		}

		if(vote == entire){
			count++;
		}

		tempAccs.push_back((double)count/labeled);
		tempChosens.push_back(L.size());
	}

	Accs = tempAccs;
	Chosens = tempChosens;
	Num = labeled;
}

// clear method
void all_clear(){
	UserList.clear();
	ItemList.clear();
	Y.clear();
	isObserved.clear();
	Accs.clear();
	Chosens.clear();
	Num = 0;
}

}
