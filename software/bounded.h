// Represents a double in the interval [0.0, 1.0]

#pragma once
#include <gtest/gtest.h>

#include "maybe.h"

class Bounded {
   public:
    static constexpr double min = 0.0;
    static constexpr double max = 1.0;

    // Empty constructor
    Bounded();

    // Helper constructor
    static Maybe<Bounded> CreateBounded(double x);

    double Get();
    static std::vector<double> Get(std::vector<Bounded>& v);
    Maybe<Void> Set(double val);

   private:
    double val;
    FRIEND_TEST(EvolutionTest, MutateTest);
};