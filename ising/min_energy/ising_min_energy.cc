// Brendan Philbin
// 29 June 2020
// Simple Ising Model MC Simulation
// min_energy analysis

#include <stdlib.h>
#include <time.h>
#include <random>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;

// SpinConfiguration Class
class SpinConfiguration {
    public:
        vector<double> spins;
        int N;

        // Constructor
        // "num_spins" = Number of particles in desired configuration
        // "spins_seed" = RAND seed for random initialization of spins
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

        // Flip spin at index specified by "index"
        void flip(int index) {
            spins[index] *= -1;
        }

        // Return spin value at index specified by "index"
        double get(int index) {
            return spins[index];
        }

        // Return number of particles present in the configuration
        int size() {
            return N;
        }

        // Print the spin configuration in CSV form
        void printSpins() {
            printf("Spins:\n");
            for(int i = 0; i < N; i++)
                printf("%f, ", spins[i]);
            printf("\n");
        }

        // Returns an int value representing the unique configuration
        // Binary value constructed by letting -1 -> 0, 1 -> 1
        unsigned long long int toInt() {
            unsigned long long int configuration = 0;
            for(int i = 0; i < N; i++) {
                if(spins[i] == 1)
                    configuration += (int) pow(2, N - i - 1);
            }
            return configuration;
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

        int size() {
            return N;
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
    int N = spin_config.size();

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
    int random_particle = rand() % spin_config.size();
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

double randfrom(double min, double max) {

    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char** argv) {

    // Parse parameters
    
    int num_spins, iterations, j_seed, mc_seed, spins_seed, ferro, trial;
    double beta;

    if(argc != 9) {
        printf("ERROR: This program requires 8 command line arguments.\n");
        printf("USAGE: ./min_energy SPINS ITERATIONS BETA J_SEED MC_SEED SPINS_SEED FERRO TRIAL\n");
        printf("SPINS, ITERATIONS, J_SEED, MC_SEED, SPINS_SEED TRIAL must be integers. BETA must be a double.\n");
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
        trial = atoi(argv[8]);
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

    // TO DO:
    //  Instead of recording energies in an array to be printed later, write the energy to disk after each sweep (i iteration)
    //  Keep a running record of the minimum energy reached during each sweep (within the j loop).
    //      To store this, record the energy and its corresponding configuration
    //      Need to find best way to convert the configuration into an integer (like binary) - HELPER FUNCTION
    //      If multiple states have same minimum energy, record their configurations as a vector

    // Declare minimum energy as max double and create int vector of min energy configurations
    double min_energy = computeEnergy(spins, j_values);
    vector<unsigned long long int> min_configs;
    min_configs.push_back(spins.toInt());
    bool swapped;
    double current_energy;
    srand(mc_seed);

    // Prepare for file writing
    ofstream testfile;

    // Open file for sweep energy CSV
    testfile.open("min_energy_output" + to_string(trial) + ".csv");
    testfile << computeEnergy(spins, j_values) << ",";

    for(int i = 1; i < iterations + 1; i++) {
        for(int j = 0; j < num_spins; j++) {
            attemptSwap(spins, j_values, beta);
            current_energy = computeEnergy(spins, j_values);
            unsigned long long int conf_int = spins.toInt();

            if(current_energy < min_energy) {
                min_energy = current_energy;
                min_configs.clear();
                min_configs.push_back(conf_int);
            }
            else if(current_energy == min_energy) {
                if( find(min_configs.begin(), min_configs.end(), conf_int) == min_configs.end() )
                    min_configs.push_back(conf_int);
            }
        }
        
        // Write energy after each sweep to CSV file
        if(i == iterations)
            testfile << computeEnergy(spins, j_values);
        else
            testfile << computeEnergy(spins, j_values) << ","; 
    }

    // Close CSV file
    testfile.close();

    // Write results file w/ minimum energy information
    testfile.open("results.txt", ios_base::app);
    testfile << "Trial #" + to_string(trial) + " : min_energy_output" + to_string(trial) + ".csv : N = " + to_string(num_spins) + ", M = " + to_string(iterations) + ", B = " + to_string(beta) + ", J_SEED = " + to_string(j_seed) + ", MC_SEED = " + to_string(mc_seed) + ", SPIN_SEED = " + to_string(spins_seed) + ", FERRO = " + to_string(ferro) + "\n\n";
    testfile << "Minimum energy reached: " + to_string(min_energy) + "\n";
    testfile << "Configurations: ";
    for(int i = 0; i < min_configs.size(); i++)
        testfile << to_string(min_configs[i]) + ", ";
    testfile << "\n\n------------------------------------\n\n";
    testfile.close();

    printf("Trial #%d completed and stored to disk.\n", trial);

}
