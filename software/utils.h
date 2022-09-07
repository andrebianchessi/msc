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
double Random() {
    return static_cast<float>(rand()) / static_cast<float>(RAND_MAX);
}

// Returns a random number in interval [x0, x1]
double Random(double x0, double x1) {
    double l = Random();
    return x0 + (x1 - x0) * l;
}

int RandomInt(int x0, int x1) {
    std::random_device rd;   // obtain a random number from hardware
    std::mt19937 gen(rd());  // seed the generator
    std::uniform_int_distribution<> distr(x0, x1);  // define the range
    return distr(gen);
};

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
Maybe<std::shared_ptr<Bounded>> Normalize(double val, double min, double max) {
    Maybe<std::shared_ptr<Bounded>> r;
    if (val < min) {
        r.isError = true;
        r.errMsg = "val must be >=min";
        return r;
    }
    if (val > max) {
        r.isError = true;
        r.errMsg = "val must be <=max";
        return r;
    }
    r.val = std::make_shared<Bounded>();
    r.val->Set((val - min) / (max - min));
    return r;
}

// "Unnormalizes" value to unbounded interval
double Unnormalize(Bounded b, double min, double max) {
    return min + b.Get() * (max - min);
};