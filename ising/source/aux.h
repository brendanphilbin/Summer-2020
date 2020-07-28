// Auxillary functions header file
// Brendan Philbin

#include "jij.h"

// Implicit method declarations
double randfrom(double min, double max);
bool areEqual(double a, double b);

// Returns a random double in range [min,max]
double randfrom(double min, double max, mt19937& rng) {
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rng() / div);
}

// Returns whether a and b are equal within an error range
bool areEqual(double a, double b) {
    double epsilon = 0.0000000001;
    return fabs(a - b) < epsilon;
}
