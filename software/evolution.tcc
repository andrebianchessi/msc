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
int Evolution<creature>::popSize() {
    return this->population->size();
}

template <typename creature>
creature* Evolution<creature>::GetCreature(int i) {
    return &(this->population->at(i));
}

template <typename creature>
vector<double> Evolution<creature>::fitness() {
    vector<double> costs = vector<double>(this->popSize());
    double totalCost = 0;
    for (int i = 0; i < this->popSize(); i++) {
        totalCost += (*this->population)[i].GetCost();
        costs[i] = (*this->population)[i].GetCost();
    }

    // Calculate the fitness vector, which indicates how fit a creature is,
    // higher fitness value -> better

    // Start by filing with cost values
    vector<double> fitness = vector<double>(this->popSize());
    for (int i = 0; i < this->popSize(); i++) {
        fitness[i] = costs[i];
    }
    // Remove negative values if needed
    // fitness: [-4,-1,0,2] -> [1,4,5,7]
    // fitness: [0,1,2] -> [1,2,3]
    // fitness: [4,6,7] -> [4,6,7]
    double minCost = costs[0];
    if (minCost <= 0) {
        for (int i = 0; i < this->popSize(); i++) {
            fitness[i] += -minCost + 1;
        }
    }
    // Invert values
    // fitness: [1,2,3] -> [1, 0.5, 0.33]
    // fitness: [4,6,7] -> [0.25, 0.16, 0.14]
    for (int i = 0; i < this->popSize(); i++) {
        fitness[i] = 1 / fitness[i];
    }

    // Normalize values
    double totalFitness = 0;
    for (int i = 0; i < this->popSize(); i++) {
        totalFitness += fitness[i];
    }
    for (int i = 0; i < this->popSize(); i++) {
        fitness[i] = fitness[i] / totalFitness;
    }
    return fitness;
}

template <typename creature>
tuple<creature*, creature*> Evolution<creature>::getParents() {
    // Biased roulette wheel method
    vector<double> f = this->fitness();
    vector<double> roulette = vector<double>(this->popSize());
    double accumulated = 0.0;
    for (int i = 0; i < this->popSize(); i++) {
        roulette[i] = f[i] + accumulated;
        accumulated += f[i];
    }

    auto it = lower_bound(roulette.begin(), roulette.end(), Random());
    int p1 = it - roulette.begin();
    int p2 = p1;
    while (p2 == p1) {
        it = lower_bound(roulette.begin(), roulette.end(), Random());
        p2 = it - roulette.begin();
    }

    return make_tuple(&((*this->population)[p1]), &((*this->population)[p2]));
};

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