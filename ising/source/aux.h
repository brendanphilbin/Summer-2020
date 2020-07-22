// Auxillary functions header file

#include "jij.h"

using namespace std;

// Implicit method declarations
double randfrom(double min, double max);
bool areEqual(double a, double b);
double computeQ(vector<Replica> replicas);
void anneal(vector<Replica>& replicas);

// Returns a random double in range [min,max]
double randfrom(double min, double max) {
    srand(time(NULL));
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

// Returns whether a and b are equal within an error range
bool areEqual(double a, double b) {
    double epsilon = 0.0000000001;
    return fabs(a - b) < epsilon;
}