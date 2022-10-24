// Abstract class (interface) that represents one candidate solution for the
// Evolution ('real-encoded'/'continuous' genetic algorithm)
#pragma once

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bounded.h"

class Creature {
   public:
    // Values that describe this creature.
    // A.K.A. Chromosome, in the literature.
    std::vector<Bounded> dna;

    // Returns value that we want to minimize.
    virtual double GetCost() = 0;

    // For expensive GetCost operations, remember to cache the GetCost values.
    // The cache is deleted (hasGetCostCache set to false) after Mutations, and
    // and children are created without cache.
    double getCostCache;
    bool hasGetCostCache;
};

int DnaSize(Creature& c);

// Applies mutation to the DNA of this creature at index i of the dna
void Mutate(Creature& c, int i);

// Combines this with c1 creature to create two children
void Mate(Creature& c0, Creature& c1, Creature* child0, Creature* child1);