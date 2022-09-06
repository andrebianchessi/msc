// File containing small auxiliary functions

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <string>

// Returns a random number in interval [0.0, 1.0]
double Random() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

// Returns a random number in interval [x0, x1]
double Random(double x0, double x1) {
    double l = Random();
    return x0 + (x1 - x0) * l;
}

// Function to easily debug by printing on output
// Example:
//  p.SetInitialDisp(d.massId, d.val);
//  print("SetInitialDisp", d.massId, d.val);
template <typename T>
void print(T t) {
    std::cout << t << " ";
}
template <typename T, typename... Args>
void print(T t, Args... args) {
    std::cout << t << " ";
    print(args...);
    std::cout << std::endl;
}