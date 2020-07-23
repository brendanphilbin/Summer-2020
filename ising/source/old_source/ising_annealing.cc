// Brendan Philbin
// 29 June 2020
// Simple Ising Model MC Simulation
// min_energy analysis with several seeds and thermal annealing

#include <stdlib.h>
#include <time.h>
#include <random>
#include <math.h>
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
        SpinConfiguration(int num_spins, mt19937& rng) {
            N = num_spins;
            spins.resize(N);
            int random_spin;
            for(int i = 0; i < N; i++) {
                random_spin = (int)(rng() % 2);
                if(random_spin == 0)
                    spins[i] = 1;
                else
                    spins[i] = -1;
            }
        }

        // Flip spin at index specified by "index"
        void flip(int index) { spins[index] *= -1; }

        // Return spin value at index specified by "index"
        double get(int index) { return spins[index]; }

        // Return number of particles present in the configuration
        int size() { return N; }

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

        // Constructor
        JijMatrix(int num_spins, mt19937& rng, bool ferro) {
            N = num_spins;
            ferromagnetic = ferro;
            j_values.resize(N, vector<double>(N));
            int random_j;

            for(int i = 0; i < N; i++) {
                for(int j = i+1; j < N; j++) {
                    if(i == j)
                        j_values[i][j] = 0;
                    else {
                        if(ferromagnetic)
                            j_values[i][j] = (double)1 / (double)N;
                        else {
                            random_j = (int)(rng() % 2);
                            if(random_j == 0)
                                j_values[i][j] = (double)1 / sqrt(N);
                            else
                                j_values[i][j] = (double)-1 / sqrt(N);
                        }
                        j_values[j][i] = j_values[i][j];
                    }
                }
            }
        }

        // Return row "i" of J values 
        vector<double> row(int i) { return j_values[i]; }

        // Returns a specific J values at index {i,j}
        double get(int i, int j) { return j_values[i][j]; }

        // Returns side length of the Matrix, the number of particles present
        int size() { return N; }

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
double attemptFlip(SpinConfiguration& spin_config, JijMatrix jij, double beta_value, double energy, mt19937& rng);
double randfrom(double min, double max);
bool areEqual(double a, double b);

// Returns the energy as a double given a spin configuration and Jij matrix
double computeEnergy(SpinConfiguration spin_config, JijMatrix jij) {

    vector<double> spin_values = spin_config.spins;
    vector<vector<double>> jij_values = jij.j_values;
    int N = spin_config.size();
    double energy = 0;
    vector<double> temp_vector;
    temp_vector.assign(N, 0);
    for(int i = 0; i < N; i++) {
        for(int j = i+1; j < N; j++)
            temp_vector[i] += spin_values[j] * jij_values[i][j];
    }
    for(int i = 0; i < N; i++)
        energy += spin_values[i] * temp_vector[i];
    energy *= -2;
    return energy;
}

// Attempt to flip a random spin depending on the acceptance algorithm
// Returns change in energy, dE, if particle was flipped. Returns DOUBLE_MAX if not.
double attemptFlip(SpinConfiguration& spin_config, JijMatrix jij, double beta_value, double energy, mt19937& rng) {

    double current_energy = energy;
    int random_particle = (int)(rng() % spin_config.size());
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
    srand(time(NULL));
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

bool areEqual(double a, double b) {
    double epsilon = 0.0000000001;
    return fabs(a - b) < epsilon;
}

int main(int argc, char** argv) {

    // Declare parameters
    int num_spins, iterations, replicas, j_seed, trial, mc_seed, spin_seed;
    double beta;
    bool ferro = false;
    bool annealing = false;

    // Parse parameters

    cxxopts::Options options("IsingAnnealing", "Performs multiple replica minimum energy analysis with annealing");

    options.add_options()
        ("n", "number of spins", cxxopts::value<int>()->default_value("-1"))
        ("m", "number of sweeps", cxxopts::value<int>()->default_value("-1"))
        ("r", "number of replicas", cxxopts::value<int>()->default_value("1"))
        ("t", "trial number", cxxopts::value<int>()->default_value("-1"))
        ("j", "Jij matrix seed", cxxopts::value<int>()->default_value("-1"))
        ("c", "MC sweep seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("s", "initial spin seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("a", "toggle temperature annealing")
        ("b", "beta value", cxxopts::value<double>()->default_value("-1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    iterations = parameters["m"].as<int>();
    replicas = parameters["r"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seed = parameters["c"].as<int>();
    spin_seed = parameters["s"].as<int>();
    annealing = parameters["a"].as<bool>();
    beta = parameters["b"].as<double>();
    ferro = parameters["f"].as<bool>();

    if(num_spins == -1 || iterations == -1 || trial == -1 || j_seed == -1 || mc_seed == -1 || spin_seed == -1 || beta == -1) {
        printf("ERROR: missing one or more required arguments\n");
        printf("Required arguments:\n");
        printf("    -n : number of spins\n");
        printf("    -m : number of sweeps\n");
        printf("    -r : number of replicas\n");
        printf("    -t : trial number\n");
        printf("    -j : Jij matrix seed\n");
        printf("    -c : MC sweep seed\n");
        printf("    -s : initial spin seed\n");
        printf("    -a : toggle annealing\n");
        printf("    -b : beta value\n");
        printf("    -f : toggle ferromagnetic\n");
        return 0;
    }

    // Handle the annealing setup
    double increment = beta / iterations;
    if(annealing) {
        beta = 0;
    }

    // Create vector of RNGs
    // Create initial SpinConfiguration vector of replicas
    // Seed first config with provided spin_seed and increment for subsequent configs
    // Re-seed RNGs with mc_seeds after creating each replicas
    vector<mt19937> rngs;
    vector<SpinConfiguration> spin_configs;
    for(int i = 0; i < replicas; i++) {
        rngs.push_back( mt19937(spin_seed + i) );
        spin_configs.push_back( SpinConfiguration(num_spins, rngs[i]) );
        rngs[i].seed(mc_seed + i);
    }

    // Create Jij matrix
    mt19937 j_rng(j_seed);
    JijMatrix j_values(num_spins, j_rng, ferro);

    // For each replica, declare minimum energy as current energy and create long long int vector of min energy configurations
    // NOTE: energy configurations can only encompass up to N = 64 particles
    vector<double> min_energies;
    vector<double> current_energies;
    vector<vector<unsigned long long int>> min_configs;

    // Create a vector of ofstream objects for each replica
    vector<ofstream> streams;
    for(int i = 0; i < replicas; i++) {
        streams.push_back(ofstream());
        streams[i].open("../annealing/annealing_output" + to_string(trial) + "_r" + to_string(i) + ".csv", ios_base::app);
    }

    // Create ofstreams for the results file and the histogram file
    ofstream results;
    results.open("../annealing/results.txt", ios_base::app);
    ofstream histogram;
    histogram.open("../annealing/histogram" + to_string(trial) + ".csv");

    for(int i = 0; i < replicas; i++) {
        double current_energy = computeEnergy(spin_configs[i], j_values);
        min_energies.push_back(current_energy);
        current_energies.push_back(current_energy);
        vector<unsigned long long int> temp_vector(1, spin_configs[i].toInt());
        min_configs.push_back(temp_vector);

        // Write initial energy information to CSV
        // NOTE: RESULTS MUST BE CLEANED BEFORE EACH RUN: use "make cleanResults"
        streams[i] << current_energies[i] << ",";
    }

    // Monte Carlo loop
    double dE;

    // Iteration over "-m", number of sweeps
    for(int i = 1; i < iterations + 1; i++) {

        if(annealing)
            beta += increment;

        // Iteration over "-r", number of replicas
        for(int r = 0; r < replicas; r++) {

            // Iteration over "-n", number of spins
            for(int j = 0; j < num_spins; j++) {

                dE = attemptFlip(spin_configs[r], j_values, beta, current_energies[r], rngs[r]);
                if(dE != numeric_limits<double>::max())
                    current_energies[r] += dE;
                unsigned long long int conf_int = spin_configs[r].toInt();

                if(current_energies[r] < min_energies[r]) {
                    min_energies[r] = current_energies[r];
                    min_configs[r].clear();
                    min_configs[r].push_back(conf_int);
                }

                else if( areEqual(current_energies[r], min_energies[r]) ) {
                    if( find(min_configs[r].begin(), min_configs[r].end(), conf_int) == min_configs[r].end() )
                        min_configs[r].push_back(conf_int);
                }
            }

           if(i == iterations)
              streams[r] << current_energies[r];
           else
              streams[r] << current_energies[r] << ","; 

        }

        // This is where I put population annealing 

    }

    // Write results file w/ minimum energy information
    results << "Trial " + to_string(trial) + ": N = " + to_string(num_spins) + ", M = " + to_string(iterations) + ", J_SEED = " + to_string(j_seed) + ", B = " + to_string(beta) + "\n\n";
    for(int i = 0; i < replicas; i++) {
        results << "   Replica " + to_string(i) + ": SPIN_SEED = " + to_string(spin_seed+i) + ", MC_SEED = " + to_string(mc_seed+i) + "\n";
        results << "   Minimum energy reached = " + to_string(min_energies[i]) + "\n";
        results << "   Minimum energy configurations: ";
        for(int j = 0; j < min_configs[i].size(); j++)
            results << to_string(min_configs[i][j]) + ", ";
        results << "\n\n";
    }

    // Write minimum energy values to histogram file
    for(int i = 0; i < replicas; i++) {
        if(i == replicas - 1)
            histogram << to_string(min_energies[i]);
        else
            histogram << to_string(min_energies[i]) + ",";
    }

    // Close all file streams
    for(int i = 0; i < replicas; i++)
        streams[i].close();
    results.close();
    histogram.close();

    // Print trial status to the command line
    printf("Trial #%d completed and stored to disk.\n", trial);
    
    return 1;
}
