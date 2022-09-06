// File containing small auxiliary functions

#include <stdlib.h>

// Returns a random number in interval [0.0, 1.0]
double Random() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

// Returns a random number in interval [x0, x1]
double Random(double x0, double x1) {
    double l = Random();
    return x0 + (x1 - x0) * l;
}