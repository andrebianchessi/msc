#include "creature.h"

#include <memory>
#include <vector>

#include "utils.h"

int DnaSize(Creature& c) { return c.dna.size(); }

void Mutate(Creature& c, int i) {
    // Uniform random mutation
    // Haupt, Randy L., Sue Ellen Haupt, and Sue Ellen Autor Haupt. 2004.
    // Practical Genetic Algorithms. Wiley.
    c.dna[i].Set(Random());

    c.hasGetCostCache = false;
}

void Mate(Creature& c0, Creature& c1, Creature* child0, Creature* child1) {
    // Radcliff blending method
    // Haupt, Randy L., Sue Ellen Haupt, and Sue Ellen Autor Haupt. 2004.
    // Practical Genetic Algorithms. Wiley.
    double beta;
    for (int i = 0; i < DnaSize(c0); i++) {
        beta = Random();

        child0->dna[i].Set(beta * c0.dna[i].Get() +
                           (1 - beta) * c1.dna[i].Get());
        child1->dna[i].Set((1 - beta) * c0.dna[i].Get() +
                           beta * c1.dna[i].Get());
    }
    child0->hasGetCostCache = false;
    child1->hasGetCostCache = false;
}