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
Bounded RandomB() {
    Bounded b;
    b.Set(Random());
    return b;
}

// Returns a random number in interval [x0, x1]
double Random(double x0, double x1) {
    double l = Random();
    return x0 + (x1 - x0) * l;
}

// Returns a random integer in interval [x0, x1]
int RandomInt(int x0, int x1) {
    std::random_device rd;   // obtain a random number from hardware
    std::mt19937 gen(rd());  // seed the generator
    std::uniform_int_distribution<> distr(x0, x1);  // define the range
    return distr(gen);
};

double RelativeAbsError(double x, double y) {
    if (x == y) {
        return 0;
    }
    if (y != 0) {
        return abs((x - y) / y);
    }
    return abs((y - x) / x);
}

// Returns value in interval [0,1] by normalizing between min and max
Maybe<Bounded> Normalize(double val, double min, double max) {
    Maybe<Bounded> r;
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
    r.val.Set((val - min) / (max - min));
    return r;
}

Maybe<double> NormalizeToDouble(double val, double min, double max) {
    Maybe<double> r;
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
    r.val = (val - min) / (max - min);
    return r;
}

// "Unnormalizes" value to unbounded interval
double Unnormalize(Bounded b, double min, double max) {
    return min + b.Get() * (max - min);
};