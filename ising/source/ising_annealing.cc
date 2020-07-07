// Brendan Philbin
// 29 June 2020
// Simple Ising Model MC Simulation
// min_energy analysis with several seeds and simulated annealing

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

        // Constructor
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
void anneal(vector<SpinConfiguration>& configs, vector<double> energies, double beta);
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

void anneal(vector<SpinConfiguration>& configs, vector<double> energies, double beta) {

    /*
    int min_index = min_element(energies.begin(), energies.end()) - v.begin();
    int min_energy = *min_element(energies.begin(), energies.end());
    */
}

// Returns a random double in range [min,max]
double randfrom(double min, double max) {

    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

int main(int argc, char** argv) {

    // Declare parameters
    
    int num_spins, iterations, replicas, j_seed, trial, anneal;
    vector<int> mc_seeds;
    vector<int> spin_seeds;
    double beta;
    bool ferro = false;

    // Parse parameters

    cxxopts::Options options("IsingMinEnergy", "Performs single replica minimum energy analysis");

    options.add_options()
        ("n", "number of spins", cxxopts::value<int>()->default_value("-1"))
        ("m", "number of sweeps", cxxopts::value<int>()->default_value("-1"))
        ("r", "number of replicas", cxxopts::value<int>()->default_value("1"))
        ("t", "trial number", cxxopts::value<int>()->default_value("-1"))
        ("j", "Jij matrix seed", cxxopts::value<int>()->default_value("-1"))
        ("c", "MC sweep seed(s)", cxxopts::value<vector<int>>()->default_value("-1"))
        ("s", "initial spin seed(s)", cxxopts::value<vector<int>>()->default_value("-1"))
        ("a", "iterations per annealing", cxxopts::value<int>()->default_value("-1"))
        ("b", "beta value", cxxopts::value<double>()->default_value("-1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    iterations = parameters["m"].as<int>();
    replicas = parameters["r"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seeds = parameters["c"].as<vector<int>>();
    spin_seeds = parameters["s"].as<vector<int>>();
    anneal = parameters["a"].as<int>();
    beta = parameters["b"].as<double>();
    ferro = parameters["f"].as<bool>();

    if(num_spins == -1 || iterations == -1 || trial == -1 || j_seed == -1 || mc_seeds.empty() || spin_seeds.empty() || beta == -1 || anneal == -1) {
        printf("ERROR: missing one or more required arguments\n");
        printf("Required options:\n");
        printf("    -n : number of spins\n");
        printf("    -m : number of sweeps\n");
        printf("    -r : number of replicas\n");
        printf("    -t : trial number\n");
        printf("    -j : Jij matrix seed\n");
        printf("    -c : MC sweep seed(s)\n");
        printf("    -s : initial spin seed(s)\n");
        printf("    -a : number of sweeps per annealing (set to a number > M for no annealing)\n");
        printf("    -b : beta value\n");
        printf("    -f : add if ferromagnetic\n");
    }

    if(mc_seeds.size() != replicas || spin_seeds.size() != replicas)
        printf("ERROR: number of SpinConfig seeds not equal to number of replicas\n");

    // Create initial SpinConfiguration vector of replicas
    vector<SpinConfiguration> spin_configs;
    for(int seed : spin_seeds)
        spin_configs.push_back(SpinConfiguration(num_spins, seed));

    // Create Jij matrix
    JijMatrix j_values(num_spins, j_seed, ferro);

    // For each replica, declare minimum energy as current energy and create long long int vector of min energy configurations
    // NOTE: energy configurations can only encompass up to N = 64 particles
    vector<double> min_energies;
    vector<double> current_energies;
    vector<vector<unsigned long long int>> min_configs;

    ofstream testfile;

    for(int i = 0; i < spin_configs.size(); i++) {
        double current_energy = computeEnergy(spin_configs[i], j_values);
        min_energies.push_back(current_energy);
        current_energies.push_back(current_energy);
        vector<unsigned long long int> temp_vector(1, spin_configs[i].toInt());
        min_configs.push_back(temp_vector);

        // Write initial energy information to CSV
        // NOTE: RESULTS MUST BE CLEANED BEFORE EACH RUN: use "make cleanResults"
        testfile.open("../several_seeds/several_seeds_output" + to_string(trial) + "_r" + to_string(i) + ".csv", ios_base::app);
        testfile << current_energies[i] << ",";
        testfile.close();
    }

    // Monte Carlo loop

    // Iteration over "-m", number of sweeps
    for(int i = 1; i < iterations + 1; i++) {

        // Iteration over "-r", number of replicas
        for(int r = 0; r < spin_configs.size(); r++) {
            
            // Re-seed each replica on every iteration to avoid same flips being chosen
            // Re-seeding with same seed results in same string of random numbers
            srand(mc_seeds[r]+i);

            testfile.open("../several_seeds/several_seeds_output" + to_string(trial) + "_r" + to_string(r) + ".csv", ios_base::app);

            // Iteration over "-n", number of spins
            for(int j = 0; j < num_spins; j++) {

                // Attempt to flip one spin
                // Calculate energy after attempted flip
                // Convert spin configuration to a long long int
                attemptFlip(spin_configs[r], j_values, beta);
                current_energies[r] = computeEnergy(spin_configs[r], j_values);
                unsigned long long int conf_int = spin_configs[r].toInt();

                // If current energy is less than minimum energy to this point, update min_energy,
                // clear configuration vector, and push configuration to the vector
                if(current_energies[r] < min_energies[r]) {
                    min_energies[r] = current_energies[r];
                    min_configs[r].clear();
                    min_configs[r].push_back(conf_int);
                }

                // If current energy is same as minimum energy, push this configuration to the vector
                else if(current_energies[r] == min_energies[r]) {
                    if( find(min_configs[r].begin(), min_configs[r].end(), conf_int) == min_configs[r].end() )
                        min_configs[r].push_back(conf_int);
                }
            }

           if(i == iterations)
              testfile << current_energies[r];
           else
              testfile << current_energies[r] << ","; 
           testfile.close();

        }

        // THIS IS WHERE REPLICA COMPARISON WILL GO
        // REACHED AFTER EACH ITERATION COMPLETES ON ALL REPLICAS
        // if( (i != 0) && (i % anneal == 0) )
        //     anneal(spin_configs, current_energies, beta);

    }

    // Write results file w/ minimum energy information
    
    testfile.open("../several_seeds/results.txt", ios_base::app);
    testfile << "Trial " + to_string(trial) + ": N = " + to_string(num_spins) + ", M = " + to_string(iterations) + ", J_SEED = " + to_string(j_seed) + ", B = " + to_string(beta) + "\n\n";
    for(int i = 0; i < spin_configs.size(); i++) {
        testfile << "   Replica " + to_string(i) + ": SPIN_SEED = " + to_string(spin_seeds[i]) + ", MC_SEED = " + to_string(mc_seeds[i]) + "\n";
        testfile << "   Minimum energy reached = " + to_string(min_energies[i]) + "\n";
        testfile << "   Minimum energy configurations: ";
        for(int j = 0; j < min_configs[i].size(); j++)
            testfile << to_string(min_configs[i][j]) + ", ";
        testfile << "\n\n";
    }
    testfile.close();

    // Print trial status to the command line
    printf("Trial #%d completed and stored to disk.\n", trial);
    
    return 1;

}
