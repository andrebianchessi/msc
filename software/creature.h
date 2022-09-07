// Abstract class that represents one candidate solution for the Evolution
//('real-encoded'/'continuous' genetic algorithm)
#pragma once

#include <memory>
#include <vector>

#include "bounded.h"

using namespace std;
class Creature {
   protected:
    // Values that describe this creature.
    // A.K.A. Chromosome, in the literature.
    vector<Bounded> dna;

   public:
    int DnaSize();

    // Returns value that we want to minimize.
    virtual double GetCost() = 0;

    // Applies mutation to the DNA of this creature at index i of the dna
    void Mutate(int i);

    // Combines this with c1 creature to create two children
    void Mate(Creature& c1, Creature* child0, Creature* child1);
};