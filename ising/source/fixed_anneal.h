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
    return -1;
}

// Population annealing function
void anneal(vector<Replica>& replicas, int& mc_seed, int& flip_seed) {
    vector<double> prob_dist;
    for(int r = 0; r < replicas.size(); r++)
        replicas[r].computeWeight();
    double normalization = computeNormalization(replicas);
    for(int r = 0; r < replicas.size(); r++) {
        double prob = replicas[r].computeProb(normalization);
        if(r == 0)
            prob_dist.push_back(prob);
        else
            prob_dist.push_back(prob_dist[r-1] + prob);
    }
    int limit = replicas.size();
    for(int r = 0; r < limit; r++) {
        double random_double = randfrom(0, 1, replicas[r].flip_rng);
        int chosen_replica = sample(prob_dist, random_double);
        replicas.push_back(replicas[chosen_replica]);
        replicas[replicas.size() - 1].setMCSeed(mc_seed++);
        replicas[replicas.size() - 1].setFlipSeed(flip_seed++);
    }
    replicas.erase(replicas.begin(), replicas.begin() + limit);
    for(int r = 0; r < replicas.size(); r++)
        replicas[r].setWeight(1);
}
