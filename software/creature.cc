#include "creature.h"

#include <memory>
#include <vector>

#include "utils.h"

int Creature::DnaSize() { return this->dna.size(); }

void Creature::Mutate(int i) {
    // Uniform random mutation
    // Haupt, Randy L., Sue Ellen Haupt, and Sue Ellen Autor Haupt. 2004.
    // Practical Genetic Algorithms. Wiley.
    this->dna[i].Set(Random());

    this->hasGetCostCache = false;
}

void Creature::Mate(Creature& c1, Creature* child0, Creature* child1) {
    // Randcliff blending method
    // Haupt, Randy L., Sue Ellen Haupt, and Sue Ellen Autor Haupt. 2004.
    // Practical Genetic Algorithms. Wiley.
    double beta;
    for (int i = 0; i < this->DnaSize(); i++) {
        beta = Random();

        child0->dna[i].Set(beta * this->dna[i].Get() +
                           (1 - beta) * c1.dna[i].Get());
        child1->dna[i].Set((1 - beta) * this->dna[i].Get() +
                           beta * c1.dna[i].Get());
    }
    child0->hasGetCostCache = false;
    child1->hasGetCostCache = false;
}