#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>
using namespace std;


int scale = 2;
int N = pow(scale, 3);
vector<int> bitState(N, 0);
vector<int> spin = {1, -1};
int S = spin.size();
int n = pow(S, N);

vector<int> state(N);
void bitToState();
void showState();

vector<vector<int>> stateSet;
void addState();

vector<vector<vector<int>>> lattice(scale, vector<vector<int>>(scale, vector<int>(scale)));
vector<int> energySet(n);
vector<int> magSet(n);
void calEnergy();
void calMag();

double Z(double k);
double u(double k);
double c(double k);
double m(double k);


int main() {
	bitToState();
	showState();
	addState();
	int idx = 0;
	while (idx < N) {
		if (bitState.at(idx) < S-1) {
			bitState.at(idx) += 1;
			bitToState();
			showState();
			addState();
			if (idx > 0) {
				idx = 0;
			}
		}
		else if (bitState.at(idx) == S-1) {
			bitState.at(idx) = 0;
			idx += 1;
		}
	}

	calEnergy();
	calMag();

	ofstream file("file.csv");
	for (int d = 0; d < 3000; d++) {
		double k = d / 1000.0;
		double U = u(k);
		double C = c(k);
		double M = m(k);
		file << k << " , " << U << " , " << C << " , " << M << endl;
	}
	file.close();

	return 0;
}


void bitToState() {
	for (int i = 0; i < N; i++) {
		state.at(i) = spin.at(bitState.at(i));
	}
	return;
}


void showState() {
	for (int i = 0; i < N; i++) {
		int j = N - 1 - i;
		cout << state.at(j);
		if (i < N-1) {
			cout << "  "; 
		}
		else if (i == N-1) {
			cout << endl;
		}
	}
	return;
}


void addState() {
	stateSet.push_back(state);
	return;
}


void calEnergy() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < N; j++) {
			int x = j / (scale*scale);
			int y = (j - x*scale*scale) / scale;
			int z = j - x*scale*scale - y*scale;
			lattice.at(x).at(y).at(z) = stateSet.at(i).at(j);
		}
		int energy = 0;
		for (int x = 0; x < scale; x++) {
			for (int y = 0; y < scale; y++) {
				for (int z = 0; z < scale; z++) {
					if (x < scale-1) {
						energy -= lattice.at(x).at(y).at(z)*lattice.at(x+1).at(y).at(z);
					}
					if (y < scale-1) {
						energy -= lattice.at(x).at(y).at(z)*lattice.at(x).at(y+1).at(z);
					}
					if (z < scale-1) {
						energy -= lattice.at(x).at(y).at(z)*lattice.at(x).at(y).at(z+1);
					}
				}
			}
		}
		energySet.at(i) = energy;
	}
	return;
}


void calMag() {
	for (int i = 0; i < n; i++) {
		int mag = 0;
		for (int j = 0; j < N; j++) {
			mag += stateSet.at(i).at(j);
		}
		magSet.at(i) = mag;
	}
	return;
}


double Z(double k) {
	double sum = 0;
	for (int i = 0; i < n; i++) {
		double e = energySet.at(i);
		sum += exp(-k*e);
	}
	return sum;
}


double u(double k) {
	double sum1 = 0;
	for (int j = 0; j < n; j++) {
		double e1 = energySet.at(j);
		sum1 += e1*exp(-k*e1);
	}
	double x = sum1 / Z(k);
	return x;
}


double c(double k) {
	double sum2 = 0;
	for (int l = 0; l < n; l++) {
		double e2 = energySet.at(l);
		sum2 += pow(e2, 2)*exp(-k*e2);
	}
	double y = pow(k, 2)*(sum2/Z(k) - pow(u(k), 2));
	return y;
}


double m(double k) {
	double sum1 = 0;
	for (int j = 0; j < n; j++) {
		double e1 = energySet.at(j);
		double m1 = pow(magSet.at(j), 2);
		sum1 += m1*exp(-k*e1);
	}
	double x = sum1 / Z(k);
	return x;
}
