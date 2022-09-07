#pragma once

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bounded.h"
#include "creature.h"
#include "utils.h"

using namespace std;
// Represents candidate solution that minimizes f:
// f(x,y) = x^2 + y^2 + 2x + y
// Global optimal solution:
// (x,y) = (-1,-0.5)
// f(-1,-0.5) = 0
class EqSol : public Creature {
   public:
    static constexpr double xMax = 2.0;
    static constexpr double xMin = -2.0;
    static constexpr double yMax = 2.0;
    static constexpr double yMin = -2.0;

    EqSol(Bounded x, Bounded y) {
        this->dna = vector<Bounded>(2);
        this->dna[0].Set(x.Get());
        this->dna[1].Set(y.Get());
    }

    // Random constructor
    EqSol() {
        Bounded x, y;
        x.Set(Random());
        y.Set(Random());
        this->dna = vector<Bounded>(2);
        this->dna[0].Set(x.Get());
        this->dna[1].Set(y.Get());
    }

    double GetCost() {
        double xVal = Unnormalize(this->dna[0], EqSol::xMin, EqSol::xMax);
        double yVal = Unnormalize(this->dna[1], EqSol::yMin, EqSol::yMax);
        return xVal * xVal + yVal * yVal + 2 * xVal + yVal;
    }

   private:
    shared_ptr<Bounded> x;
    shared_ptr<Bounded> y;
    FRIEND_TEST(CreatureTest, SimpleTest);
    FRIEND_TEST(CreatureTest, MateAndMutationTest);
    FRIEND_TEST(EvolutionTest, MutateTest);
    FRIEND_TEST(EvolutionTest, fitnessTest);
};