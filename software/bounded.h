// Represents a double in the interval [0.0, 1.0]

#pragma once
#include "maybe.h"

class Bounded {
   public:
    const double min = 0.0;
    const double max = 1.0;

    // Only has empty constructor. Value must be set using Set method.
    Bounded();

    double Get();
    Maybe<Void> Set(double val);

   private:
    double val;
};