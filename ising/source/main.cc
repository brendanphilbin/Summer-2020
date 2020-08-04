// Ising Model simulation
// Main file
//
// Recursively includes all header files
//#include "anneal.h"
#include "fixed_anneal.h"
#include "cxxopts.hpp"

int main(int argc, char** argv) {

    // Declare parameters
    int num_spins, sweeps, num_replicas, j_seed, trial, mc_seed, spin_seed, flip_seed, steps;
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
    ferro = parameters["f"].as<bool>();

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
        printf("\n");
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
        replicas.push_back( Replica(num_spins, spin_seed++, flip_seed++, mc_seed++, beta, increment, jij) );
    }

    // Create stream for energies CSV
    ofstream energies;
    energies.open("../results/energies_t" + to_string(trial) + ".csv");

    // Create stream for population size analysis
    ofstream pop_size;
    pop_size.open("../results/pop_size_t" + to_string(trial) + ".csv");

    // TEST
    mt19937 test_rng;
    test_rng.seed(flip_seed++);

    // Monte Carlo loop
    for(int k = 0; k < steps; k++) {
        // anneal(replicas, num_replicas, mc_seed, flip_seed);
        anneal(replicas, test_rng, mc_seed, flip_seed);
        for(int i = 0; i < sweeps; i++) {
            for(int r = 0; r < replicas.size(); r++) {
                replicas[r].incrementBeta();
                replicas[r].performSweep();
                if(r != replicas.size() - 1)
                    energies << to_string(replicas[r].getEnergy()) + ",";
                else
                    energies << to_string(replicas[r].getEnergy());
            }
            energies << "\n";
        }
        if(k != steps - 1)
            pop_size << to_string(replicas.size()) + ",";
        else
            pop_size << to_string(replicas.size());
    }

    // Close all streams
    energies.close();
    pop_size.close();

    printf("Trial #%d has completed and stored to disk.\n", trial);

    return 1;
}
