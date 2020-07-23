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
#include <cxxopts.hpp>

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

// Jij Matrix class
class JijMatrix {
    public:
        vector<vector<double>> j_values;
        int N;
        bool ferromagnetic;

        // Empty constructor
        JijMatrix(int num_spins, int j_seed, bool ferro) {
            N = num_spins;
            ferromagnetic = ferro;
            j_values.resize(N, vector<double>(N));
            srand(j_seed);

            for(int i = 0; i < N; i++) {
                for(int j = i+1; j < N; j++) {
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
                        j_values[j][i] = j_values[i][j];
                    }
                }
            }
        }

        // Returns a specific J values at index {i,j}
        double get(int i, int j) {
            return j_values[i][j];
        }

        // Returns side length of the Matrix, the number of particles present
        int size() {
            return N;
        }

        // Print the matrix
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

// Implicit method declarations
double computeEnergy(SpinConfiguration spin_config, JijMatrix jij);
double attemptFlip(SpinConfiguration& spin_config, JijMatrix jij, double beta_value);
double randfrom(double min, double max);

// Returns the energy as a double given a spin configuration and Jij matrix
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

// Attempt to flip a random spin, seeded in the MC loop, depending on the acceptance algorithm
// Returns change in energy, dE, if particle was flipped. Returns DOUBLE_MAX if not.
// TO DO: CHANGE TO ONLY TAKE IN ONE ROW OF jij

double attemptFlip(SpinConfiguration& spin_config, JijMatrix jij, double beta_value) {

    double current_energy = computeEnergy(spin_config, jij);
    int random_particle = rand() % spin_config.size();
    spin_config.flip(random_particle);
    double new_energy = computeEnergy(spin_config, jij);
    double dE = new_energy - current_energy;

    if(dE < 0)
        return dE;
    else {
        double comparison = exp(-1 * dE * beta_value);
        double random_d = randfrom(0,1);

        if(random_d >= comparison) {
            spin_config.flip(random_particle);
            return numeric_limits<double>::max();
        }
        else
            return dE;
    }
}

// Returns a random double in range [min,max]
double randfrom(double min, double max) {

    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char** argv) {

    // Declare parameters
    
    int num_spins, iterations, j_seed, mc_seed, spins_seed, trial;
    double beta;
    bool ferro = false;

    // Parse parameters

    cxxopts::Options options("IsingMinEnergy", "Performs single replica minimum energy analysis");

    options.add_options()
        ("n", "number of spins", cxxopts::value<int>()->default_value("-1"))
        ("m", "number of sweeps", cxxopts::value<int>()->default_value("-1"))
        ("t", "trial number", cxxopts::value<int>()->default_value("-1"))
        ("j", "Jij matrix seed", cxxopts::value<int>()->default_value("-1"))
        ("c", "MC sweep seed", cxxopts::value<int>()->default_value("-1"))
        ("s", "initial spin seed", cxxopts::value<int>()->default_value("-1"))
        ("b", "beta value", cxxopts::value<double>()->default_value("-1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    iterations = parameters["m"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seed = parameters["c"].as<int>();
    spins_seed = parameters["s"].as<int>();
    beta = parameters["b"].as<double>();
    ferro = parameters["f"].as<bool>();

    if(num_spins == -1 || iterations == -1 || trial == -1 || j_seed == -1 || mc_seed == -1 || spins_seed == -1 || beta == -1) {
        printf("\nERROR: missing one or more required options\n");
        printf("Required options:\n");
        printf("    -n : number of spins\n");
        printf("    -m : number of sweeps\n");
        printf("    -t : trial number\n");
        printf("    -j : Jij matrix seed\n");
        printf("    -mc : MC sweep seed\n");
        printf("    -s : initial spin seed\n");
        printf("    -b : beta value\n");
        printf("    -f : add if ferromagnetic\n\n");
        return 0;
    }

    // Create initial SpinConfiguration
    SpinConfiguration spins(num_spins, spins_seed);

    // Create Jij matrix
    JijMatrix j_values(num_spins, j_seed, ferro);

    // Monte Carlo loop

    // Declare minimum energy as current energy and create long long int vector of min energy configurations
    // NOTE: energy configurations can only encompass up to N = 64 particles
    double min_energy = computeEnergy(spins, j_values);
    vector<unsigned long long int> min_configs;
    min_configs.push_back(spins.toInt());

    double current_energy;
    srand(mc_seed);

    // Prepare for file writing
    ofstream testfile;

    // Open file for sweep energy CSV
    testfile.open("../min_energy/min_energy_output" + to_string(trial) + ".csv");
    testfile << computeEnergy(spins, j_values) << ",";

    for(int i = 1; i < iterations + 1; i++) {
        for(int j = 0; j < num_spins; j++) {

            // Attempt to flip one spin
            // Calculate energy after attempted flip
            // Convert spin configuration to a long long int
            attemptFlip(spins, j_values, beta);
            current_energy = computeEnergy(spins, j_values);
            unsigned long long int conf_int = spins.toInt();

            // If current energy is less than minimum energy to this point, update min_energy,
            // clear configuration vector, and push configuration to the vector
            if(current_energy < min_energy) {
                min_energy = current_energy;
                min_configs.clear();
                min_configs.push_back(conf_int);
            }

            // If current energy is same as minimum energy, push this configuration to the vector
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
    testfile.open("../min_energy/results.txt", ios_base::app);
    testfile << "Trial #" + to_string(trial) + " : min_energy_output" + to_string(trial) + ".csv : N = " + to_string(num_spins) + ", M = " + to_string(iterations) + ", B = " + to_string(beta) + ", J_SEED = " + to_string(j_seed) + ", MC_SEED = " + to_string(mc_seed) + ", SPIN_SEED = " + to_string(spins_seed) + ", FERRO = " + to_string(ferro) + "\n\n";
    testfile << "Minimum energy reached: " + to_string(min_energy) + "\n";
    testfile << "Configurations: ";
    for(int i = 0; i < min_configs.size(); i++)
        testfile << to_string(min_configs[i]) + ", ";
    testfile << "\n\n------------------------------------\n\n";
    testfile.close();

    // Print trial(s) status to the command line
    printf("Trial #%d completed and stored to disk.\n", trial);

    return 1;

}
