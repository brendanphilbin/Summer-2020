// Population annealing header file
// Fixed population size
// Brendan Philbin

#include "replica.h"

// Compute the normalization function
double computeNormalization(vector<Replica> replicas) {
    double sum = 0;
    for(int r = 0; r < replicas.size(); r++)
        sum += replicas[r].getWeight();
    return sum;
}

// Sample the probability distribution
int sample(vector<double> prob_dist, double random_double) {
    for(int i = 0; i < prob_dist.size(); i++) {
        if(prob_dist[i] > random_double)
            return i;
    }
    return prob_dist.size() - 1;
}

void anneal(vector<Replica>& replicas) {
    vector<double> prob_dist;
    double normalization = 0;
    for(int r = 0; r < replicas.size(); r++)
        normalization += replicas[r].computeWeight();
    for(int r = 0; r < replicas.size(); r++) {
        double prob = replicas[r].computeProb(normalization);
        if(r == 0)
            prob_dist.push_back(prob);
        else
            prob_dist.push_back(prob_dist[r-1] + prob);
    }
    vector<vector<double>> new_spins;
    for(int i = 0; i < replicas.size(); i++) {
        double random_double = randfrom(0, 1, replicas[i].flip_rng);
        int chosen_replica = sample(prob_dist, random_double);
        new_spins.push_back(replicas[chosen_replica].getSpins());
    }
    for(int r = 0; r < replicas.size(); r++) {
        replicas[r].setSpins(new_spins[r]); 
        replicas[r].computeEnergy();
        replicas[r].setWeight(1);
    }
}

int threadStart(int thread_id, int threads, int num_replicas) {
    int min_replicas_per_thread = num_replicas / threads;
    int leftovers = num_replicas % threads;
    if(thread_id < leftovers)
        return (thread_id * (min_replicas_per_thread + 1));
    else
        return ((leftovers * (min_replicas_per_thread + 1)) + ((thread_id - leftovers) * min_replicas_per_thread));
    printf("ERROR\n");
    return -1;
}

int threadEnd(int thread_id, int threads, int num_replicas) {
    int start = threadStart(thread_id, threads, num_replicas);
    int min_replicas_per_thread = num_replicas / threads;
    int leftovers = num_replicas % threads;
    if(thread_id < leftovers)
        return start + min_replicas_per_thread + 1;
    else
        return start + min_replicas_per_thread;
    printf("ERROR\n");
    return -1;
}
