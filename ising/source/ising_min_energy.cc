// Brendan Philbin
// 29 June 2020
// Simple Ising Model MC Simulation
// min_energy analysis

#include <stdlib.h>
#include <time.h>
#include <random>

using namespace std;

class SpinConfiguration {
    public:
        vector<double> spins;
        int N;

        SpinConfiguration(int num_spins, int spins_seed) {
            N = num_spins;
            spins.resize(N);
            srand(spins_seed);
            for(int i = 0; i < N; i++) {
                int random_spin = rand() % 2;
                if(random_spin == 0)
                    spins[i] = 1;
                else
                    spins[i] = -1;
            }
        }

        void flip(int index) {
            spins[index] *= -1;
        }

        double get(int index) {
            return spins[index];
        }

        void printSpins() {
            printf("Spins:\n");
            for(int i = 0; i < N; i++)
                printf("%f, ", spins[i]);
            printf("\n");
        }
};

class JijMatrix {
    public:
        vector<vector<double>> j_values;
        int N;
        bool ferromagnetic;

        JijMatrix() {}

        void initialize(int num_spins, int j_seed, bool ferro) {
            N = num_spins;
            ferromagnetic = ferro;
            j_values.resize(N, vector<double>(N));
            srand(j_seed);

            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++) {
                    if(i == j)
                        j_values[i][j] = 0;
                    else {
                        if(ferromagnetic)
                            j_values[i][j] = 1 / (double) N;
                        else {
                            int random_j = rand() % 2;
                            if(random_j == 0)
                                j_values[i][j] = 1 / sqrt(N);
                            else
                                j_values[i][j] = -1 / sqrt(N);
                        }
                    }
                }
            }
        }

        double get(int i, int j) {
            return j_values[i][j];
        }

        void printMatrix() {
            printf("Jij Matrix:\n");
            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++)
                    printf("%f,", j_values[i][j]);
                printf("\n");
            }
            printf("\n");
        }

};

double computeEnergy(SpinConfiguration spin_config, JijMatrix jij);
bool attemptSwap(SpinConfiguration& spin_config, JijMatrix jij, double beta_value);
double randfrom(double min, double max);

double computeEnergy(SpinConfiguration spin_config, JijMatrix jij) {

    vector<double> spin_values = spin_config.spins;
    vector<vector<double>> jij_values = jij.j_values;
    int N = spin_config.N;

    double energy = 0;

    vector<double> temp_vector;
    temp_vector.assign(N, 0);

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++)
            temp_vector[i] += spin_values[j] * jij_values[i][j];
    }

    for(int i = 0; i < N; i++)
        energy += spin_values[i] * temp_vector[i];

    energy *= -1;

    return energy;
}

bool attemptSwap(SpinConfiguration& spin_config, JijMatrix jij, double beta_value) {

    double current_energy = computeEnergy(spin_config, jij);
    int random_particle = rand() % spin_config.N;
    spin_config.flip(random_particle);
    double new_energy = computeEnergy(spin_config, jij);
    double dE = new_energy - current_energy;

    if(dE < 0)
        return true;
    else {
        double comparison = exp(-1 * dE * beta_value);
        double random_d = randfrom(0,1);

        if(random_d >= comparison) {
            spin_config.flip(random_particle);
            return false;
        }
        else
            return true;
    }
}

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char** argv) {

    // Parse parameters
    
    int num_spins, iterations, j_seed, mc_seed, spins_seed, ferro;
    double beta;

    if(argc != 8) {
        printf("ERROR: This program requires 7 command line arguments.\n");
        printf("USAGE: ./min_energy SPINS ITERATIONS BETA J_SEED MC_SEED SPINS_SEED FERRO\n");
        printf("SPINS, ITERATIONS, J_SEED, MC_SEED, SPINS_SEED must be integers. BETA must be a double.\n");
        printf("For ferromagnetic model (J = 1), set FERRO = 1. For spin glass (J = +/- 1), set FERRO = 0)\n");
        return 0;
    }

    else {
        num_spins = atoi(argv[1]);
        iterations = atoi(argv[2]);
        beta = atof(argv[3]);
        j_seed = atoi(argv[4]);
        mc_seed = atoi(argv[5]);
        spins_seed = atoi(argv[6]);
        ferro = atoi(argv[7]);
    }

    if(num_spins == 0 || iterations == 0 || beta == 0 || j_seed == 0 || mc_seed == 0 || spins_seed == 0 || !(ferro == 0 || ferro == 1)) {
        printf("ERROR: Invalid argument provided.\n");
        printf("USAGE: ./min_energy SPINS ITERATIONS BETA J_SEED MC_SEED SPINS_SEED FERRO\n");
        printf("SPINS, ITERATIONS, J_SEED, MC_SEED, SPINS_SEED must be integers. BETA must be a double.\n");
        printf("For ferromagnetic model (J = 1), set FERRO = 1. For spin glass (J = +/- 1), set FERRO = 0)\n");
        return 0;
    }

    // Create initial SpinConfiguration
    SpinConfiguration spins(num_spins, spins_seed);

    // Create Jij matrix
    JijMatrix j_values;
    if(ferro == 0)
        j_values.initialize(num_spins, j_seed, false);
    else
        j_values.initialize(num_spins, j_seed, true);

    // Create energies array and define first value, before any sweeps are performed
    double energies[iterations + 1];
    energies[0] = computeEnergy(spins, j_values);

    // Monte Carlo loop
    srand(mc_seed);
    for(int i = 1; i < iterations + 1; i++) {
        for(int j = 0; j < num_spins; j++)
            attemptSwap(spins, j_values, beta);
        energies[i] = computeEnergy(spins, j_values);
    }

    // Print final energies in CSV-readable format
    for(int i = 0; i < iterations + 1; i++) {
        if(i != iterations)
            printf("%f,", energies[i]);
        else
            printf("%f", energies[i]);
    }
}
