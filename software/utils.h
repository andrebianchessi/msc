// File containing small auxiliary functions
#pragma once

#include <stdio.h>
#include <stdlib.h>

#include <iostream>
#include <memory>
#include <random>
#include <string>

#include "bounded.h"
#include "maybe.h"

// Returns a random number in interval [0.0, 1.0]
double Random();
Bounded RandomB();

// Returns a random number in interval [x0, x1]
double Random(double x0, double x1);

// Returns a random integer in interval [x0, x1]
int RandomInt(int x0, int x1);

double RelativeAbsError(double x, double y);

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

// Returns value in interval [0,1] by normalizing between min and max
Maybe<Bounded> Normalize(double val, double min, double max);
Maybe<double> NormalizeToDouble(double val, double min, double max);

// "Unnormalizes" value to unbounded interval
double Unnormalize(Bounded b, double min, double max);

// Convenience function to get current time
std::chrono::_V2::system_clock::time_point Now();

// Convenient function to get readable time difference from t0 to now
std::string TimeSince(std::chrono::_V2::system_clock::time_point t0);