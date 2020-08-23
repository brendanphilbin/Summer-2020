// Ising Model simulation
// Main file w/ OpenMP

// Recursively includes all header files
#include "anneal.h"
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
    bool fixed = false;

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
        ("f", "ferromagnetic or not")
        ("p", "toggle fixed population");

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
    ferro = parameters["f"].as<bool>();
    fixed = parameters["p"].as<bool>();

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
        // printf("    -p : toggle fixed population size\n");
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

    int min_replicas_per_thread = num_replicas / threads;
    int leftovers = num_replicas % threads;

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

    #pragma omp parallel num_threads(threads)
    {
       int thread_id = omp_get_thread_num();
       vector<Replica> replicas;
       
       // Populate private replica vectors on each thread
       /*
       if(thread_id < leftovers) {
           for(int i = 0; i < min_replicas_per_thread + 1; i++)
               replicas.push_back( population[ thread_id * (min_replicas_per_thread + 1) + i ] );
       }
       else {
           for(int i = 0; i < min_replicas_per_thread; i++)
               replicas.push_back( population[ (min_replicas_per_thread + 1) * (leftovers) + (thread_id - leftovers) * (min_replicas_per_thread) + i ] );
       }
       */

       // Monte Carlo loop
       for(int k = 0; k < steps; k++) {

           // RECALCULATE REPLICA DISTRIBUTION

           // Perform population annealing
           #pragma omp master
           {
              if(fixed)
                  anneal(population);
              else
                  anneal(population, num_replicas, population[population.size() - 1].mc_seed, population[population.size() - 1].flip_seed);
              min_replicas_per_thread = population.size() / threads;
              leftovers = population.size() % threads;
           }

           // Re-assign private replica vectors
           replicas.clear();
           if(thread_id < leftovers) {
             for(int i = 0; i < min_replicas_per_thread + 1; i++)
                replicas.push_back(population[ thread_id * (min_replicas_per_thread + 1) + i ]);
           }
           else {
             for(int i = 0; i < min_replicas_per_thread; i++)
                replicas.push_back(population[ (min_replicas_per_thread + 1) * (leftovers) + (thread_id - leftovers) * (min_replicas_per_thread) + i ]);
           }  

           // Monte Carlo loop
           for(int i = 0; i < sweeps; i++) {
               for(int r = 0; r < replicas.size(); r++) {
                   replicas[r].incrementBeta();
                   #pragma omp critical
                   {
                       replicas[r].performSweep(thread_id, r);
                   }

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

           printf("Pop size = %ld\n", population.size());

           // Update shared population replica vector
           if(thread_id < leftovers) {
               for(int i = 0; i < replicas.size(); i++)
                   population[ thread_id * (min_replicas_per_thread + 1) + i ] = replicas[i];
           }
           else {
               for(int i = 0; i < replicas.size(); i++)
                   population[ (min_replicas_per_thread + 1) * (leftovers) + (thread_id - leftovers) * (min_replicas_per_thread) + i ] = replicas[i];
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
