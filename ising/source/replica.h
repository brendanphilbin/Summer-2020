// Replica class header file
// Brendan Philbin

#include "aux.h"


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

// Sets the Poisson distribution seed and re-seeds RNG
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
    // See Matcha (2010) page 2 - under equation (2)
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
