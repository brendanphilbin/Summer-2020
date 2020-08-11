// Brendan Philbin
// Parallel population annealing auxillary functions

#include "replica.h"
#include "fixed_anneal.h"

double computeWeights(vector<Replica>& replicas) {
    double normalization = 0;
    for(Replica r : replicas)
        normalization += r.computeWeight();
    return normalization;
}

void computeProbabilites(vector<Replica>& replicas, double normalization) {
    for(Replica r : replicas)
        r.computeProb(normalization);
}
