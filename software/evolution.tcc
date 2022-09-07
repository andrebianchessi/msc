#pragma once
#include <algorithm>
#include <memory>

#include "creature.h"
#include "evolution.h"
#include "utils.h"

using namespace std;

template <typename creature>
Evolution<creature>::Evolution(shared_ptr<vector<creature>> population) {
    this->population = population;
}

template <typename creature>
void Evolution<creature>::sortPopulation() {
    auto comparator = [](creature& left, creature& right) {
        return left.GetCost() < right.GetCost();
    };
    sort((this->population)->begin(), (this->population)->end(), comparator);
}

template <typename creature>
int Evolution<creature>::endFittest() {
    return floor(this->nKeep * this->population->size()) - 1;
}

template <typename creature>
void Evolution<creature>::mutate() {
    if (this->population->size() == 1) {
        // the first creature is excluded from mutation,
        // so we require at least 2 creatures for mutations
        return;
    }

    // Total number of slots that can be mutated
    // (the fittest creature is not mutated)
    int nDnaSlots =
        (this->population->size() - 1) * ((*this->population)[0].DnaSize());
    // target num of mutations
    int nMutations = int(nDnaSlots * this->mutationRate);
    // current mutation count
    int mutations = 0;

    // Set containing [creatureId, dnaPosition] pairs which were mutated
    set<tuple<int, int>> mutationsDone = set<tuple<int, int>>();
    while (mutations < nMutations) {
        // pick random creature (outside best one) and mutate
        int creatureI = RandomInt(1, this->population->size() - 1);
        int dnaI = RandomInt(0, (*this->population)[creatureI].DnaSize() - 1);
        auto pair = std::tuple<int, int>(creatureI, dnaI);
        if (mutationsDone.find(pair) == mutationsDone.end()) {
            (*this->population)[creatureI].Mutate(dnaI);
            mutations += 1;
            mutationsDone.insert(pair);
        }
    }
}