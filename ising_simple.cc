// Brendan Philbin
// 23 June 2020
// Simple Ising Model MC Simulation

#include <stdlib.h>
#include <time.h>
#include <random>

double randfrom(double min, double max);

int main() {

    // Parameters

    int num_spins = 500;
    int iterations = 500;

    int spins[num_spins];

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
    
    int j_values[num_spins][num_spins];

    for(int i = 0; i < num_spins; i++) {
        for(int j = 0; j < num_spins; j++) {
            if(i == j)
                j_values[i][j] = 0;
            else {
                int random2 = rand() % 2;
                if(random == 0)
                    j_values[i][j] = 1;
                else
                    j_values[i][j] = -1;
            }
        }
    }

    // Monte Carlo loop

    for(int p = 0; p < iterations; p++) {

        // Compute energy of current configuration

        int current_energy = 0;
        int temp_array_current[num_spins] = { 0 };

        for(int i = 0; i < num_spins; i++) {
            for(int j = 0; j < num_spins; j++) {

                temp_array_current[i] += spins[j] * j_values[i][j];

            }
        }

        for(int i = 0; i < num_spins; i++) {
            current_energy += spins[i] * temp_array_current[i];
        }
        current_energy *= -1;

        // Print current energy
        printf("Ei = %d", current_energy);

        // Choose random particle 

        int random_particle = rand() % num_spins;
       
        // Print which particle is being flipped
        printf(" || Particle %d", random_particle);

        // Flip spin

        spins[random_particle] *= -1;

        // Compute change in energy with new configuration

        int new_energy = 0;
        int temp_array_new[num_spins] = { 0 };

        for(int i = 0; i < num_spins; i++) {
            for(int j = 0; j < num_spins; j++) {

                temp_array_new[i] += spins[j] * j_values[i][j];

            }
        }

        for(int i = 0; i < num_spins; i++) {
            new_energy += spins[i] * temp_array_new[i];
        }
        new_energy *= -1;

        // Print new energy
        printf(" || Ef = %d", new_energy);
        
        // Use acceptance algorithm to determine whether to retain change

        int energy_change = (new_energy - current_energy);
        double energy_change_double = (double) energy_change;
        
        // Print energy change
        printf(" || dE: %d", energy_change);

        double temperature = 100;
        double boltzmann = 0.00000000000000000000001380649;

        if(energy_change > 0) {
            double random3 = randfrom(0, 1);
            
            // Print random number
            printf(" || R: %f", random3);

            // double comparison = std::exp((-1 * energy_change_double) / (temperature * boltzmann));
            double comparison = std::exp((-1 * energy_change_double) / temperature);

            // Print comparison values
            printf(" || C: %f", comparison);

            if(random3 >= comparison) {
                spins[random_particle] *= -1;
                
                // Print if not flipped
                printf(" || No flip\n");

            }

            // Print if flipped
            else { printf(" || Flip\n"); }

        }
       
        // Print if flipped 
        else { printf(" || Flip\n"); }

    }

}

double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
