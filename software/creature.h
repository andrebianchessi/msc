// Abstract class that represents one candidate solution for the Evolution
//('real-encoded'/'continuous' genetic algorithm)
#pragma once

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bounded.h"

class Creature {
   protected:
    // For expensive GetCost operations, remember to cache the GetCost values.
    // The cache is deleted (hasGetCostCache set to false) after Mutations, and
    // and children are created without cache.
    double getCostCache;
    bool hasGetCostCache;
    FRIEND_TEST(CreatureTest, GetCostCacheTest);

   public:
    int DnaSize();

    // Returns value that we want to minimize.
    virtual double GetCost() = 0;

    // Applies mutation to the DNA of this creature at index i of the dna
    void Mutate(int i);

    // Combines this with c1 creature to create two children
    void Mate(Creature& c1, Creature* child0, Creature* child1);

    // Values that describe this creature.
    // A.K.A. Chromosome, in the literature.
    std::vector<Bounded> dna;
};