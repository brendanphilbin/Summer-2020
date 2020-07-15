// Brendan Philbin
// Ising Model w/ population annealing 
// 11 July 2020

// Include files
#include <stdlib.h>
#include <time.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <cxxopts.hpp>

using namespace std;

// Implicit method declarations
double randfrom(double min, double max);
bool areEqual(double a, double b);

// Returns a random double in range [min,max]
double randfrom(double min, double max) {
    srand(time(NULL));
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Returns whether a and b are equal within an error range
bool areEqual(double a, double b) {
    double epsilon = 0.0000000001;
    return fabs(a - b) < epsilon;
}

// Jij Matrix class
class JijMatrix {

    public:
        vector<vector<double>> j_values;
        int N;
        bool ferromagnetic;
	    mt19937 j_rng;

        // Default constructor
        JijMatrix() {
            N = 0;
            ferromagnetic = false;
            vector<double> temp;
            temp.push_back(0);
            j_values.push_back(temp);
        }

        // Constructor
        JijMatrix(int num_spins, bool ferro, int j_seed) {
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
};


// Replica class
class Replica {
    
    // Instance variables
    public:
        int N, spin_seed, mc_seed;
        mt19937 spin_rng;
        mt19937 mc_rng;
        double beta, previous_beta, energy, min_energy;
        vector<double> spins;
        vector<unsigned long long int> min_configs;
        JijMatrix jij;

        // Constructor for a Replica object
        Replica(int num_spins, int spin_seed_param, int mc_seed_param, double beta_param, JijMatrix jij_param) {
            N = num_spins;
            beta = beta_param;
            jij = jij_param;
            spin_seed = spin_seed_param;
            mc_seed = mc_seed_param;
            spin_rng.seed(spin_seed);
            mc_rng.seed(mc_seed);
            int random_spin;
            for(int i = 0; i < N; i++) {
                random_spin = (int)(this->spin_rng() % 2);
                if(random_spin == 0)
                    spins.push_back(1);
                else
                    spins.push_back(-1);
            }
            energy = computeEnergy();
            min_energy = energy;
            min_configs.push_back(toInt());
        }

        // Returns # of particles in system
        int size() { return N; }

        int getSpinSeed() { return spin_seed; }

        // Returns Monte Carlo sweep seed
        int getMCSeed() { return mc_seed; }

        //Sets the Monte Carlo sweep seed and re-seeds RNG
        void setMCSeed(int seed) {
            mc_seed = seed;
            mc_rng.seed(mc_seed);
        }

        // Returns current beta value
        double getBeta() { return beta; }

        // Set beta value
        void setBeta(double beta_param) { beta = beta_param; }

        // Returns previous beta value
        double getPreviousBeta() { return previous_beta; }

        // Increments beta by a defined amount
        void incrementBeta(double dBeta) {
            previous_beta = beta;
            beta += dBeta;
        }

        // Returns total system energy
        double getEnergy() { return energy; }

        // Returns minimum energy reached
        double getMinEnergy() { return min_energy; }

        double attemptFlip(int index) {
		    double current_energy = energy;
		    spins[index] *= -1;
		    double new_energy = computeEnergy();
		    double dE = new_energy - current_energy;

		    if(dE < 0)
			    return dE;
		    else {
			    double comparison = exp(-1 * dE * beta);
                double random_d = randfrom(0, 1);
                
                if(random_d >= comparison) {
                    spins[index] *= -1;
                    return numeric_limits<double>::max();
                }
                else
                    return dE;
            }
	    }

        // Computes system energy 
        double computeEnergy() {
           vector<vector<double>> jij_values = this->jij.j_values;
           double energy_value = 0;
           vector<double> temp_vector;
           temp_vector.assign(N, 0);
           for(int i = 0; i < N; i++) {
               for(int j = i+1; j < N; j++)
                   temp_vector[i] += spins[j] * jij_values[i][j];
           }
           for(int i = 0; i < N; i++)
               energy_value += spins[i] * temp_vector[i];
           energy_value *= -2;
           return energy_value;
        }

        // Performs a single Monte Carlo sweep on the system
        void performSweep() {
            // Loop over N attempted flips
                // Attempt flip on a random particle
                // Use acceptance algorithm to retain or undo flip
                // Update energy values
            int random_particle;
            double dE;
            for(int i = 0; i < N; i++) {
                random_particle = (int)(this->mc_rng() % N);
                dE = attemptFlip(random_particle);
		        if(dE != numeric_limits<double>::max()) {
			        energy += dE;
                    if(energy < min_energy) {
                        min_energy = energy;
                        min_configs.clear();
                        min_configs.push_back(toInt());
                    }
                    else if(areEqual(energy, min_energy))
                        min_configs.push_back(toInt());
                }
            }
        }

        // Returns int value representing binary version of configuration
        // -1 -> 0 , 1 -> 1
        unsigned long long int toInt() {
	        unsigned long long int configuration = 0;
            for(int i = 0; i < N; i++) {
                if(spins[i] == 1)
                    configuration += (int) pow(2, N - i - 1);
            }
            return configuration;
        }

};


// Main function
int main(int argc, char** argv) {

    // Declare parameters
    int num_spins, iterations, num_replicas, j_seed, trial, mc_seed, spin_seed;
    double beta;
    bool ferro = false;
    bool temp_annealing = false;
    bool pop_annealing = false;

    // Parse parameters

    cxxopts::Options options("IsingPop", "Performs multiple replica simulation w/ temperature and population annealing");

    options.add_options()
        ("n", "number of spins", cxxopts::value<int>()->default_value("-1"))
        ("m", "number of sweeps", cxxopts::value<int>()->default_value("-1"))
        ("r", "number of replicas", cxxopts::value<int>()->default_value("1"))
        ("t", "trial number", cxxopts::value<int>()->default_value("-1"))
        ("j", "Jij matrix seed", cxxopts::value<int>()->default_value("-1"))
        ("c", "MC sweep seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("s", "initial spin seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("a", "toggle temperature annealing")
	    ("p", "toggle population annealing")
        ("b", "beta value", cxxopts::value<double>()->default_value("-1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    iterations = parameters["m"].as<int>();
    num_replicas = parameters["r"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seed = parameters["c"].as<int>();
    spin_seed = parameters["s"].as<int>();
    temp_annealing = parameters["a"].as<bool>();
    pop_annealing = parameters["p"].as<bool>();
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
        printf("    -a : toggle temperature annealing\n");
	    printf("    -p : toggle population annealing\n");
        printf("    -b : beta value\n");
        printf("    -f : toggle ferromagnetic\n");
        return 0;
    }

    // Create Jij Matrix
    JijMatrix jij(num_spins, ferro, j_seed); 

    // Create vector of replicas
    vector<Replica> replicas;
    for(int i = 0; i < num_replicas; i++) {
        replicas.push_back( Replica(num_spins, spin_seed + i, mc_seed + i, beta, jij) );
    }

    // Create replica output streams
    vector<ofstream> streams;
    for(int i = 0; i < replicas.size(); i++) {
        streams.push_back(ofstream());
        streams[i].open("../pop/pop_output" + to_string(trial) + "_r" + to_string(i) + ".csv");
    }

    // Handle thermal annealing incrementation
    double increment;
    if(temp_annealing) {
        increment = beta / (double) iterations;
        for(int r = 0; r < replicas.size(); r++)
            replicas[r].setBeta(0);
    }
    
    // Monte Carlo loop
    for(int i = 0; i < iterations; i++) {
        for(int r = 0; r < replicas.size(); r++) {
            if(temp_annealing)
                replicas[r].incrementBeta(increment);
            replicas[r].performSweep();
            if(i != iterations - 1)
                streams[r] << to_string(replicas[r].getEnergy()) + ",";
            else
                streams[r] << to_string(replicas[r].getEnergy());
        }
    }

    // Create streams for results and histogram files
    ofstream results, histogram;
    results.open("../pop/results" + to_string(trial) + ".txt");
    histogram.open("../pop/histogram" + to_string(trial) + ".csv");

    // Results and histogram file output
    results << "Trial " + to_string(trial) + " results\n";
    results << "N = " + to_string(num_spins) + ", M = " + to_string(iterations) + ", J_SEED = " + to_string(j_seed) + ", B = " + to_string(beta) + "\n\n";
    for(int r = 0; r < replicas.size(); r++) {
        results << "Replica " + to_string(r) + "\n";
        results << "SPIN_SEED = " + to_string(replicas[r].getSpinSeed()) + ", MC_SEED = " + to_string(replicas[r].getMCSeed()) + "\n";
        results << "Min energy reached = " + to_string(replicas[r].getMinEnergy()) + "\n\n";
        if(r != replicas.size() - 1)
            histogram << to_string(replicas[r].getMinEnergy()) + ",";
        else
            histogram << to_string(replicas[r].getMinEnergy());
    }

    // Close all streams
    for(int i = 0; i < streams.size(); i++)
        streams[i].close();
    results.close();
    histogram.close();

    return 1;
}
