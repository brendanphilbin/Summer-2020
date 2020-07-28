// J_{ij} matrix class header file
// Brendan Philbin

// Include files
#include <stdlib.h>
#include <time.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

// Jij Matrix class
class JijMatrix {
    public:
        vector<vector<double>> j_values;
        int N;
        bool ferromagnetic;
	    mt19937 j_rng;

        JijMatrix();
        JijMatrix(int num_spins, bool ferro, int j_seed);
};

// Default constructor
JijMatrix::JijMatrix() {
    N = 0;
    ferromagnetic = false;
    vector<double> temp;
    temp.push_back(0);
    j_values.push_back(temp);
}

// Constructor w/ parameters
JijMatrix::JijMatrix(int num_spins, bool ferro, int j_seed) {
    N = num_spins;
    ferromagnetic = ferro;
    j_values.resize(N, vector<double>(N));
    int random_j;
    j_rng.seed(j_seed);

    for(int k = 0; k < N; k++)
        j_values[k][k] = 0;

    for(int i = 0; i < N; i++) {
        for(int j = i+1; j < N; j++) {
            if(ferromagnetic)
                j_values[i][j] = (double)1 / (double)N;
            else {
                random_j = (int)(j_rng() % 2);
                if(random_j == 0)
                    j_values[i][j] = (double)1 / sqrt(N);
                else
                    j_values[i][j] = (double)-1 / sqrt(N);
            }
            j_values[j][i] = j_values[i][j];
        }
    }
}
