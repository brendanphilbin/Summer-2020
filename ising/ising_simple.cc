// Brendan Philbin
// 23 June 2020
// Simple Ising Model MC Simulation

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

            for(int i = 0; i < N; i++) {
                for(int j = 0; j < N; j++) {
                    if(i == j)
                        j_values[i][j] = 0;
                    else {
                        if(ferromagnetic)
                            j_values[i][j] = 1 / (double) N;
                        else {
                            srand(j_seed);
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
    
    int num_spins, iterations, j_seed, mc_seed, ferro;
    double beta;

    if(argc != 7) {
        printf("ERROR: This program requires 6 command line arguments.\n");
        printf("USAGE: ./ising_simple SPINS ITERATIONS BETA J_SEED MC_SEED FERRO\n");
        printf("SPINS, ITERATIONS, J_SEED, MC_SEED must be integers. BETA must be a double.\n");
        printf("For ferromagnetic model (J = 1), set FERRO = 1. For spin glass (J = +/- 1), set FERRO = 0)\n");
        return 0;
    }

    else {
        num_spins = atoi(argv[1]);
        iterations = atoi(argv[2]);
        beta = atof(argv[3]);
        j_seed = atoi(argv[4]);
        mc_seed = atoi(argv[5]);
        ferro = atoi(argv[6]);
    }

    if(num_spins == 0 || iterations == 0 || beta == 0 || j_seed == 0 || mc_seed == 0 || !(ferro == 0 || ferro == 1)) {
        printf("ERROR: Invalid argument provided.\n");
        printf("USAGE: ./ising_simple SPINS ITERATIONS BETA J_SEED MC_SEED FERRO\n");
        printf("SPINS, ITERATIONS, J_SEED, MC_SEED must be integers. BETA must be a double.\n");
        printf("For ferromagnetic model (J = 1), set FERRO = 1. For spin glass (J = +/- 1), set FERRO = 0)\n");
        return 0;
    }

    // Create initial SpinConfiguration
    SpinConfiguration spins(num_spins, mc_seed);

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
//        for(int j = 0; j < num_spins; j++)
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

/*
int main(int argc, char** argv) {

    // Parse parameters

    int num_spins, iterations, j_seed, mc_seed;
    double beta;

    if(argc != 6) {
        printf("ERROR: This program requires 4 command line arguments.\n");
        printf("Usage: ./ising_simple SPINS ITERATIONS BETA J_SEED MC_SEED\n");
	printf("Set J_SEED to 0 for ferromagnetic spin glass.\n");
        printf("Ex: ./ising_simple 100 10000 1 5 7\n");
        return 0;
    }

    else {
        num_spins = atoi(argv[1]);
        iterations = atoi(argv[2]);
        beta = atof(argv[3]);
        j_seed = atoi(argv[4]);
        mc_seed = atoi(argv[5]);
    }

    if(num_spins == 0 || iterations == 0 || beta == 0) {
        printf("ERROR: Invalid argument provided.\n");
        printf("Usage: ./ising_simple SPINS ITERATIONS BETA J_SEED MC_SEED\n");
        printf("Ex: ./ising_simple 100 10000 1 5 7\n");
        return 0;
    }
    
    // Print parameters
    // printf("Spins: %d, Iterations: %d, Beta: %f\n\n", num_spins, iterations, beta);

    double spins[num_spins];

    srand(mc_seed);
    
    // Initialize 1D array of random spins, +1 or -1

    for(int i = 0; i < num_spins; i++) {
        int random1 = rand() % 2;
        if(random1 == 0)
            spins[i] = 1;
        else
            spins[i] = -1;
    }
    
    // Initialize 2D array of J (bond strength) values
    // If all_ones is true, all J values expect main diagonal are +1
    // If all_ones if false, J values except main diagonal are randomly +1 or -1
    
    bool all_ones = true;
    
    double j_values[num_spins][num_spins];
    srand(j_seed);

    for(int i = 0; i < num_spins; i++) {
        for(int j = 0; j < num_spins; j++) {
            if(i == j)
                j_values[i][j] = 0;
            else {
                if(all_ones)
                    j_values[i][j] = 1 / (double) num_spins;
                else {
                    int random2 = rand() % 2;
                    if(random == 0)
                        j_values[i][j] = 1 / std::sqrt((double) num_spins);
                    else
                        j_values[i][j] = -1 / std::sqrt((double) num_spins);
                }
            }
        }
    }

    // Create array of energy values
    // First index is energy before first iteration 
    // Last index is energy after final iteration
    
    double energies[iterations+1];
    double current_energy, new_energy;

    // Monte Carlo loop

    for(int p = 0; p < iterations; p++) {

        // Compute energy of initial configuration

        if(p == 0) {
            current_energy = 0;
            double temp_array_current[num_spins] = { 0 };

            for(int i = 0; i < num_spins; i++) {
                for(int j = 0; j < num_spins; j++) {

                    temp_array_current[i] += spins[j] * j_values[i][j];

                }
            }

            for(int i = 0; i < num_spins; i++) {
                current_energy += spins[i] * temp_array_current[i];
            }
            current_energy *= -1;
            energies[0] = current_energy;
        }

        // Print current energy
        // printf("Ei = %d", current_energy);

        // Choose random particle 

        int random_particle = rand() % num_spins;

        // Flip spin

        spins[random_particle] *= -1;

        // Compute change in energy with new configuration

        new_energy = 0;
        double temp_array_new[num_spins] = { 0 };

        for(int i = 0; i < num_spins; i++) {
            for(int j = 0; j < num_spins; j++) {

                temp_array_new[i] += spins[j] * j_values[i][j];

            }
        }

        for(int i = 0; i < num_spins; i++) {
            new_energy += spins[i] * temp_array_new[i];
        }
        new_energy *= -1;

        // Print proposed new energy
        // printf(" || Ef = %d", new_energy);

        // Use acceptance algorithm to determine whether to retain change

        double energy_change = (new_energy - current_energy);
        
        if(energy_change > 0) {
            double random3 = randfrom(0, 1);
            
            // Print random number
            // printf(" || R: %f", random3);

            double comparison = std::exp(-1 * energy_change * beta);

            // Print comparison values
            // printf(" || C: %f", comparison);

            if(random3 >= comparison) {
                spins[random_particle] *= -1;
                energies[p+1] = current_energy;
                
                // Print if not flipped
                // printf(" || No flip\n");

            }

            // Print if flipped
            else {
                energies[p+1] = new_energy;
                current_energy = new_energy;
                // printf(" || Flip\n");
            }
        }
       
        // Print if flipped 
        else {
            energies[p+1] = new_energy;
            current_energy = new_energy;
            // printf(" || Flip\n");
        }
    }

    //printf("\nFinal energies:\n");
    for(int i = 0; i < iterations + 1; i++) {
        if(i != iterations)
            printf("%f,", energies[i]);
        else
            printf("%f", energies[i]);
    }
    //printf("\n");
}
*/
