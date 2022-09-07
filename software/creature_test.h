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
class C : public Creature {
   public:
    static constexpr double xMax = 2.0;
    static constexpr double xMin = -2.0;
    static constexpr double yMax = 2.0;
    static constexpr double yMin = -2.0;

    C(Bounded x, Bounded y) {
        this->dna = vector<Bounded>(2);
        this->dna[0].Set(x.Get());
        this->dna[1].Set(y.Get());
    }

    // Random constructor
    C() {
        Bounded x, y;
        x.Set(Random());
        y.Set(Random());
        this->dna = vector<Bounded>(2);
        this->dna[0].Set(x.Get());
        this->dna[1].Set(y.Get());
    }

    double GetCost() {
        double xVal = this->GetX();
        double yVal = this->GetY();
        return abs(xVal * xVal + yVal * yVal + 2 * xVal + yVal);
    }

    double GetX() { return Unnormalize(this->dna[0], C::xMin, C::xMax); }

    double GetY() { return Unnormalize(this->dna[1], C::xMin, C::xMax); }

   private:
    shared_ptr<Bounded> x;
    shared_ptr<Bounded> y;
    FRIEND_TEST(CreatureTest, SimpleTest);
    FRIEND_TEST(CreatureTest, MateAndMutationTest);
    FRIEND_TEST(EvolutionTest, MutateTest);
    FRIEND_TEST(EvolutionTest, fitnessTest);
};