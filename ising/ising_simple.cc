// Brendan Philbin
// 23 June 2020
// Simple Ising Model MC Simulation

#include <stdlib.h>
#include <time.h>
#include <random>

double randfrom(double min, double max);

int main(int argc, char** argv) {

    // Parse parameters

    int num_spins, iterations;
    double beta;

    if(argc != 4) {
        printf("ERROR: This program requires 3 command line arguments.\n");
        printf("Usage: ./ising_simple SPINS ITERATIONS BETA\n");
        printf("Ex: ./ising_simple 100 10000 1\n");
        return 0;
    }

    else {
        num_spins = atoi(argv[1]);
        iterations = atoi(argv[2]);
        beta = atof(argv[3]);
    }

    if(num_spins == 0 || iterations == 0 || beta == 0) {
        printf("ERROR: Invalid argument provided.\nUsage: ./ising_simple SPINS ITERATIONS BETA\nEx: ./ising_simple 100 10000 50\n");
        return 0;
    }
    
    // Print parameters
    // printf("Spins: %d, Iterations: %d, Beta: %f\n\n", num_spins, iterations, beta);

    double spins[num_spins];

    srand(time(NULL));
    
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
    
    bool all_ones = false;
    
    double j_values[num_spins][num_spins];

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

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
