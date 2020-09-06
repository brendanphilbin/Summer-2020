// Auxillary functions header file
// Brendan Philbin

#include "jij.h"

// Implicit method declarations
double randfrom(double min, double max);
bool areEqual(double a, double b);
int prob_sample(vector<double> prob_dist, mt19937& rng);

// Returns a random double in range [min,max]
double randfrom(double min, double max, mt19937& rng) {
    double range = (max - min); 
    double div = rng.max() / range;
    return min + (rng() / div);
}

// Returns whether a and b are equal within an error range
bool areEqual(double a, double b) {
    double epsilon = 0.0000000001;
    return fabs(a - b) < epsilon;
}

// Samples the probability distribution
int prob_sample(vector<double> prob_dist, mt19937& rng) {
    uniform_real_distribution<double> unif(0, 1);
    double random_double = unif(rng);
    for(int i = 0; i < prob_dist.size(); i++) {
        if(prob_dist[i] > random_double)
            return i;
    }
    return prob_dist.size() - 1;
}