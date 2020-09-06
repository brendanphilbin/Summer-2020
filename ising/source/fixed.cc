// Ising Model simulation
// Main file w/ OpenMP

// Recursively includes all header files
// #include "anneal.h"
#include "replica.h"
#include "cxxopts.hpp"
#include <omp.h>
#include <thread>
#include <chrono>

int main(int argc, char** argv) {
    
    // Start wall clock
    auto time_begin = chrono::high_resolution_clock::now();

    // Declare parameters
    int num_spins, sweeps, num_replicas, j_seed, trial, mc_seed, spin_seed, flip_seed, steps, threads;
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
        ("l", "flip [0,1] seed", cxxopts::value<int>()->default_value("-1"))
        ("b", "target beta value", cxxopts::value<double>()->default_value("-1"))
        ("k", "number of temperature steps", cxxopts::value<int>()->default_value("-1"))
        ("g", "number of threads", cxxopts::value<int>()->default_value("1"))
        ("f", "ferromagnetic or not");

    auto parameters = options.parse(argc, argv);

    num_spins = parameters["n"].as<int>();
    sweeps = parameters["m"].as<int>();
    num_replicas = parameters["r"].as<int>();
    trial = parameters["t"].as<int>();
    j_seed = parameters["j"].as<int>();
    mc_seed = parameters["c"].as<int>();
    spin_seed = parameters["s"].as<int>();
    flip_seed = parameters["l"].as<int>();
    beta = parameters["b"].as<double>();
    steps = parameters["k"].as<int>();
    threads = parameters["g"].as<int>();

    // Check for all required parameters being present
    if(num_spins == -1 || sweeps == -1 || trial == -1 || j_seed == -1 || mc_seed == -1 || flip_seed == -1 ||spin_seed == -1 || beta == -1 || steps == -1) {
        printf("\n");
        printf("ERROR: missing one or more required arguments\n");
        printf("\n");
        printf("Required arguments:\n");
        printf("    -n : number of spins\n");
        printf("    -r : initial / target population size\n");
        printf("    -m : number of MC sweeps per temperature step\n");
        printf("    -k : number of temperature steps\n");
        printf("    -b : target beta value\n");
        printf("    -t : trial number\n");
        printf("Required seeds:\n");
        printf("    -j : Jij matrix seed\n");
        printf("    -c : MC sweep seed\n");
        printf("    -s : initial spin seed\n");
        printf("    -l : flip [0,1] distribution seed\n");
        printf("Optional arguments:\n");
        printf("    -f : toggle ferromagnetic interaction\n");
        printf("    -g : number of threads (default = 1)\n");
        printf("\n");
        return 0;
    }

    // Check that number of threads doesn't exceed available processors
    int processor_count = std::thread::hardware_concurrency();
    if(threads > processor_count) {
        printf("\nERROR: number of threads requested exceeds available processors\n");
        printf("Number of available processors = %d\n\n", processor_count);
        return 0;
    }

    double increment = beta / (double) steps;
    beta = 0;

    JijMatrix jij(num_spins, ferro, j_seed);
    vector<Replica> population;
    ofstream energies;
    energies.open("../results/energies_t" + to_string(trial) + ".csv");
    for(int r = 0; r < num_replicas; r++) {
        population.push_back( Replica(num_spins, spin_seed++, flip_seed++, mc_seed++, beta, increment, jij) );
        energies << to_string(population[r].getEnergy()) + ",";
    }
    energies << "\n";

    // Define shared variables
    int min_replicas_per_thread = num_replicas / threads;
    int leftovers = num_replicas % threads;
    vector<double> weights(num_replicas, 0);
    vector<vector<double>> old_spins(num_replicas, vector<double>(num_spins, 0));
    // for(int i = 0; i < population.size(); i++)
    //     old_spins.push_back(population[i].getSpins());
    vector<vector<double>> new_spins(num_replicas, vector<double>(num_spins, 0));

    #pragma omp parallel num_threads(threads)
    {
       int thread_id = omp_get_thread_num();
       vector<Replica> replicas;

       // Populate private Replica vectors
       if(thread_id < leftovers) {
           for(int i = 0; i < min_replicas_per_thread + 1; i++) {
               int offset = thread_id * (min_replicas_per_thread + 1) + i;
               int private_spin_seed = spin_seed + offset;
               int private_flip_seed = flip_seed + offset;
               int private_mc_seed = mc_seed + offset;
               replicas.push_back( Replica(num_spins, private_spin_seed, private_flip_seed, private_mc_seed, beta, increment, jij) );
           }
       }
       else {
           for(int i = 0; i < min_replicas_per_thread; i++) {
               int offset = (min_replicas_per_thread + 1) * leftovers + (thread_id - leftovers) * min_replicas_per_thread + i;
               int private_spin_seed = spin_seed + offset;
               int private_flip_seed = flip_seed + offset;
               int private_mc_seed = mc_seed + offset;
               replicas.push_back( Replica(num_spins, private_spin_seed, private_flip_seed, private_mc_seed, beta, increment, jij) ); 
           }
       }
       
       // Temperature / annealing step loop
       for(int k = 0; k < steps; k++) {

           // Perform fixed population annealing

           for(int r = 0; r < replicas.size(); r++) {
               if(thread_id < leftovers) {
                   weights[ thread_id * (min_replicas_per_thread + 1) + r ] = replicas[r].computeWeight();
                   old_spins[ thread_id * (min_replicas_per_thread + 1) + r] = replicas[r].getSpins();
               }
               else {
                   weights[ (min_replicas_per_thread + 1) * leftovers + (thread_id - leftovers) * min_replicas_per_thread + r ] = replicas[r].computeWeight();
                   old_spins[ (min_replicas_per_thread + 1) * leftovers + (thread_id - leftovers) * min_replicas_per_thread + r] = replicas[r].getSpins();
               }
           }
           #pragma omp barrier

           #pragma omp master
           {
               double normalization = 0;
               for(int i = 0; i < weights.size(); i++)
                   normalization += weights[i];
               vector<double> prob_dist(num_replicas, 0);
               prob_dist[0] = weights[0] / normalization;
               for(int i = 1; i < weights.size(); i++)
                   prob_dist[i] = prob_dist[i-1] + weights[i] / normalization;
               for(int i = 0; i < weights.size(); i++)
                   new_spins[i] = old_spins[prob_sample(prob_dist, replicas[0].getFlipRng())];
           }
           #pragma omp barrier
           
           for(int r = 0; r < replicas.size(); r++) {
               if(thread_id < leftovers)
                   replicas[r].setSpins( new_spins[ thread_id * (min_replicas_per_thread + 1) + r ]);
               else
                   replicas[r].setSpins( new_spins[ (min_replicas_per_thread + 1) * leftovers + (thread_id - leftovers) * min_replicas_per_thread + r]);
           }
           #pragma omp barrier

           // Monte Carlo Loop
          
           for(int i = 0; i < sweeps; i++) {
               for(int r = 0; r < replicas.size(); r++) {

                   // Reduce temperature and perform Monte Carlo sweep
                   replicas[r].incrementBeta();
                   replicas[r].performSweep();

                   // Write energy information to file
                   #pragma omp critical
                   {
                       energies << to_string(replicas[r].getEnergy()) + ",";
                   }
               }
               #pragma omp barrier
               #pragma omp master
               {
                   energies << "\n";
               }
           }

       }
    }

    // Close all streams
    energies.close();

    // Stop wall clock
    auto time_end = chrono::high_resolution_clock::now();
    auto elapsed_time = chrono::duration_cast<chrono::microseconds>(time_end - time_begin);

    printf("Trial #%d has completed and stored to disk. Execution time = %f milliseconds\n", trial, elapsed_time.count() * 1e-3);

    return 1;
}
