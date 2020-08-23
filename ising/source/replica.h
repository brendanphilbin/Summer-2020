// Replica class header file
// Brendan Philbin

#include "aux.h"

class Replica {  
    public:

	    // Instance Variables
        int N, spin_seed, flip_seed, mc_seed;
        mt19937 spin_rng, flip_rng, mc_rng;
        double beta, beta_increment, energy, min_energy, tau;
        vector<double> spins;
        JijMatrix jij;

        // Member functions
        Replica(int num_spins, int spin_seed_param, int flip_seed_param, int mc_seed_param, double beta_param, double beta_increment_param, JijMatrix jij_param);
        int size();
        void setMCSeed(int seed);
        void setFlipSeed(int seed);
        double getBeta();
        void setBeta(double beta_param);
        double getBetaIncrement();
        void incrementBeta();
        double getEnergy();
        double getMinEnergy();
        double attemptFlip(int index);
        double computeEnergy();

        // void performSweep();
        void performSweep(int thread, int replica);
        int getFlipSeed();
        int getMCSeed();

        double computeTau(double q, int targetPop, int size);
        double getTau();
        int numCopies();

        // Get and set spins
        vector<double> getSpins();
        void setSpins(vector<double> spins_param);

        // Fixed population size
        double weight;
        double computeWeight();
        double getWeight();
        void setWeight(double weight_param);
        double probability;
        double computeProb(double normalization);
};

// Constructor for a Replica object
Replica::Replica(int num_spins, int spin_seed_param, int flip_seed_param, int mc_seed_param, double beta_param, double beta_increment_param, JijMatrix jij_param) {
    N = num_spins;
    beta = beta_param;
    beta_increment = beta_increment_param;
    jij = jij_param;
    spin_seed = spin_seed_param;
    flip_seed = flip_seed_param;
    mc_seed = mc_seed_param;
    spin_rng.seed(spin_seed);
    flip_rng.seed(flip_seed);
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
    weight = 1;
    min_energy = energy;
}

int Replica::getFlipSeed() { return flip_seed; }
int Replica::getMCSeed() { return mc_seed; }

// Returns # of particles in system
int Replica::size() { return N; }

// Returns spins
vector<double> Replica::getSpins() { return spins; }

// Set spins
void Replica::setSpins(vector<double> spins_param) { spins = spins_param; }

// Sets the Monte Carlo sweep seed and re-seeds RNG
void Replica::setMCSeed(int seed) {
    mc_seed = seed;
    mc_rng.seed(mc_seed);
}

// Sets the [0,1] seed for acceptance algorithm
void Replica::setFlipSeed(int seed) {
    flip_seed = seed;
    flip_rng.seed(flip_seed);
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
        double random_d = randfrom(0, 1, flip_rng);
        
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
void Replica::performSweep(int thread, int replica) {
    int random_particle;
    double dE;
    for(int i = 0; i < N; i++) {
        random_particle = (int)(this->mc_rng() % N);
//        printf("Thread %d, replica %d chose particle %d, ", thread, replica, random_particle);
        dE = attemptFlip(random_particle);
        if(dE != numeric_limits<double>::max()) {
 //           printf("flipped with dE = %f\n", dE);
            energy += dE;
            if(energy < min_energy)
                min_energy = energy;
        }
  //      else { printf("not flipped\n"); }
    }
}

// Compute normalized weight, tau
// Assigns to member variable and returns value
double Replica::computeTau(double q, int targetPop, int size) {
    tau = (targetPop / size) * ( exp( -1 * beta_increment * energy) / q );
    return tau;
}

double Replica::getTau() { return tau; }

int Replica::numCopies() {
    double prob_floor = ceil(tau) - tau;
    double dist = randfrom(0, 1, flip_rng);
    if(dist < prob_floor)
        return floor(tau);
    else
        return ceil(tau);
}

// Computes weight for fixed population size
double Replica::computeWeight() {
    // edit: add weight *
    double new_weight = exp(-1 * beta_increment * energy);
    weight = new_weight;
    return weight;
}

double Replica::getWeight() { return weight; }

void Replica::setWeight(double weight_param) {
    weight = weight_param;
}

double Replica::computeProb(double normalization) {
   double prob = weight / normalization;
   probability = prob;
   return probability; 
}
