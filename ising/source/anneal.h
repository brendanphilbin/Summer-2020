// Population annealing header file
// Brendan Philbin

#include "replica.h"

// Compute the normalization function
double computeQ(vector<Replica> replicas) {
    double numerator = 0;
    double dBeta = replicas[0].getBetaIncrement();
    double denominator = (double) replicas.size();
    for(int i = 0; i < replicas.size(); i++)
        numerator += exp( -1 * dBeta * replicas[i].getEnergy() );
    return numerator / denominator;
} 

// Population annealing function
void anneal(vector<Replica>& replicas, int targetPop, int& mc_seed, int& flip_seed) {
    vector<Replica> replicaCopies;
    double q = computeQ(replicas);
    for(int r = 0; r < replicas.size(); r++) {
        replicas[r].computeTau(q, targetPop, replicas.size());
        int copies = replicas[r].numCopies();
        if(copies == 0) {
            replicas.erase(replicas.begin() + r);
            r--;
        }
        else {
            for(int i = 1; i < copies; i++) {
                replicaCopies.push_back(replicas[r]);
                replicaCopies[replicaCopies.size() - 1].setMCSeed(mc_seed++);
                replicaCopies[replicaCopies.size() - 1].setFlipSeed(flip_seed++);
            }
        }
    }
    replicas.insert(replicas.end(), replicaCopies.begin(), replicaCopies.end());
}
