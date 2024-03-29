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

// Replica class
class Replica {  
    public:

        // Instance Variables
        int N, spin_seed, mc_seed, poisson_seed;
        mt19937 spin_rng;
        mt19937 mc_rng;
        mt19937 poisson_rng;
        double beta, beta_increment, energy, min_energy, tau;
        vector<double> spins;
        vector<unsigned long long int> min_configs;
        JijMatrix jij;

        // Member functions
        Replica(int num_spins, int spin_seed_param, int mc_seed_param, int poisson_seed_param, double beta_param, double beta_increment_param, JijMatrix jij_param);
        int size();
        void setMCSeed(int seed);
        void setPoissonSeed(int seed);
        double getBeta();
        void setBeta(double beta_param);
        double getBetaIncrement();
        void incrementBeta();
        double getEnergy();
        double getMinEnergy();
        double attemptFlip(int index);
        double computeEnergy();
        void performSweep();
        double computeTau(double q);
        double getTau();
        int numCopies(int targetPop, int size);
        unsigned long long int toInt();

};

// Constructor for a Replica object
Replica::Replica(int num_spins, int spin_seed_param, int mc_seed_param, int poisson_seed_param, double beta_param, double beta_increment_param, JijMatrix jij_param) {
    N = num_spins;
    beta = beta_param;
    beta_increment = beta_increment_param;
    jij = jij_param;
    spin_seed = spin_seed_param;
    mc_seed = mc_seed_param;
    poisson_seed = poisson_seed_param;
    spin_rng.seed(spin_seed);
    mc_rng.seed(mc_seed);
    poisson_rng.seed(poisson_seed);
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
int Replica::size() { return N; }

// Sets the Monte Carlo sweep seed and re-seeds RNG
void Replica::setMCSeed(int seed) {
    mc_seed = seed;
    mc_rng.seed(mc_seed);
}

// Sets the Poisson distribution's RNG seed
void Replica::setPoissonSeed(int seed) {
    poisson_seed = seed;
    poisson_rng.seed(poisson_seed);
}

// Returns current beta value
double Replica::getBeta() { return beta; }

// Set beta value
void Replica::setBeta(double beta_param) { beta = beta_param; }

// Returns previous beta value
double Replica::getBetaIncrement() { return beta_increment; }

// Increments beta by a defined amount
void Replica::incrementBeta() { beta += beta_increment; }

// Returns total system energy
double Replica::getEnergy() { return energy; }

// Returns minimum energy reached
double Replica::getMinEnergy() { return min_energy; }

double Replica::attemptFlip(int index) {
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
double Replica::computeEnergy() {
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
void Replica::performSweep() {
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

// Compute normalized weight, tau
// Assigns to member variable and returns value
double Replica::computeTau(double q) {
    tau = ( exp( -1 * beta_increment * energy) / q );
    return tau;
}

double Replica::getTau() { return tau; }

int Replica::numCopies(int targetPop, int size) {
    double mean = tau * targetPop / size;
    poisson_distribution<int> poisson(mean);
    return poisson(poisson_rng);
}

// Returns int value representing binary version of configuration
// -1 -> 0 , 1 -> 1
unsigned long long int Replica::toInt() {
    unsigned long long int configuration = 0;
    for(int i = 0; i < N; i++) {
        if(spins[i] == 1)
            configuration += (int) pow(2, N - i - 1);
    }
    return configuration;
}

// Implicit method declarations
double computeQ(vector<Replica> replicas);
void anneal(vector<Replica>& replicas);

// Compute the normalization function
double computeQ(vector<Replica> replicas) {
    double numerator = 0;
    double dBeta = replicas[0].getBetaIncrement();
    double denominator = (double) replicas.size();
    for(int i = 0; i < replicas.size(); i++)
        numerator += exp( -1 * dBeta * replicas[i].getEnergy() );
    return numerator / denominator;
} 

// Population annealing function
void anneal(vector<Replica>& replicas, int targetPop, int& mc_seed, int& poisson_seed) {
    vector<Replica> replicaCopies;
    double q = computeQ(replicas);
    for(int r = 0; r < replicas.size(); r++) {
        replicas[r].computeTau(q);
        int copies = replicas[r].numCopies(targetPop, replicas.size());
        if(copies == 0) {
            replicas.erase(replicas.begin() + r);
            r--;
        }
        else {
            for(int i = 0; i < copies; i++) {
                replicaCopies.push_back(replicas[r]);
                replicaCopies[replicaCopies.size() - 1].setMCSeed(mc_seed++);
                replicaCopies[replicaCopies.size() - 1].setPoissonSeed(poisson_seed++);
            }
        }
    }
    replicas.insert(replicas.end(), replicaCopies.begin(), replicaCopies.end());
}

// Main function
int main(int argc, char** argv) {

    // Declare parameters
    int num_spins, sweeps, num_replicas, j_seed, trial, mc_seed, spin_seed, poisson_seed, steps;
    double beta;
    bool ferro = false;

    // Parse parameters

    cxxopts::Options options("IsingPop", "Performs multiple replica simulation w/ temperature and population annealing");

    options.add_options()
        ("n", "number of spins", cxxopts::value<int>()->default_value("-1"))
        ("m", "number of MC sweeps per temperature step", cxxopts::value<int>()->default_value("-1"))
        ("r", "number of replicas", cxxopts::value<int>()->default_value("1"))
        ("t", "trial number", cxxopts::value<int>()->default_value("-1"))
        ("j", "Jij matrix seed", cxxopts::value<int>()->default_value("-1"))
        ("c", "MC sweep seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("s", "initial spin seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("p", "poisson seed(s)", cxxopts::value<int>()->default_value("-1"))
        ("b", "target beta value", cxxopts::value<double>()->default_value("-1"))
        ("k", "number of temperature steps", cxxopts::value<int>()->default_value("-1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    sweeps = parameters["m"].as<int>();
    num_replicas = parameters["r"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seed = parameters["c"].as<int>();
    spin_seed = parameters["s"].as<int>();
    poisson_seed = parameters["p"].as<int>();
    beta = parameters["b"].as<double>();
    steps = parameters["k"].as<int>();
    ferro = parameters["f"].as<bool>();

    if(num_spins == -1 || sweeps == -1 || trial == -1 || j_seed == -1 || mc_seed == -1 || spin_seed == -1 || poisson_seed == -1 || beta == -1 || steps == -1) {
        printf("ERROR: missing one or more required arguments\n");
        printf("Required arguments:\n");
        printf("    -n : number of spins\n");
        printf("    -m : number of MC sweeps per temperature step\n");
        printf("    -r : number of replicas\n");
        printf("    -t : trial number\n");
        printf("    -j : Jij matrix seed\n");
        printf("    -c : MC sweep seed\n");
        printf("    -s : initial spin seed\n");
        printf("    -p : poisson distribution seed\n");
        printf("    -b : beta value\n");
        printf("    -k : number of temperature steps\n");
        printf("    -f : toggle ferromagnetic interaction\n");
        return 0;
    }

    // Create Jij Matrix
    JijMatrix jij(num_spins, ferro, j_seed); 

    // Set beta incrementation
    double increment = beta / (double) steps;
    beta = 0;

    // Create vector of replicas
    vector<Replica> replicas;
    for(int i = 0; i < num_replicas; i++) {
        replicas.push_back( Replica(num_spins, spin_seed++, mc_seed++, poisson_seed++, beta, increment, jij) );
    }

    // Create stream for each histogram
    vector<ofstream> histograms;
    for(int i = 0; i < sweeps * steps; i++) {
        histograms.push_back(ofstream());
        histograms[i].open("../pop/hist_t" + to_string(trial) + "_s" + to_string(i) + ".csv");
    }

    // Monte Carlo loop
    for(int k = 0; k < steps; k++) {
        anneal(replicas, num_replicas, mc_seed, poisson_seed);
        for(int i = 0; i < sweeps; i++) {
            for(int r = 0; r < replicas.size(); r++) {
                replicas[r].incrementBeta();
                replicas[r].performSweep();
                if(r != replicas.size() - 1)
                    histograms[k * sweeps + i] << to_string(replicas[r].getEnergy()) + ",";
                else
                    histograms[k * sweeps + i] << to_string(replicas[r].getEnergy());
            }
        }
        printf("After %d temperature steps, population size = %lu\n", k+1, replicas.size());
    }

    // Close all streams
    for(int i = 0; i < sweeps * steps; i++)
        histograms[i].close();

    printf("Trial #%d has completed and stored to disk.\n", trial);

    return 1;
}
