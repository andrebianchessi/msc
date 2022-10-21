#pragma once
#include <algorithm>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "creature.h"
#include "evolution.h"
#include "maybe.h"
#include "utils.h"

template <typename creature>
Evolution<creature>::Evolution(std::vector<creature>* population) {
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
    return floor(this->survival * this->population->size()) - 1;
}

template <typename creature>
int Evolution<creature>::nFittest() {
    return this->endFittest() + 1;
}

template <typename creature>
int Evolution<creature>::PopSize() {
    return this->population->size();
}

template <typename creature>
double Evolution<creature>::FittestCost() {
    double cost = 0;
    for (int i = 0; i <= this->endFittest(); i++) {
        cost += this->GetCreature(i)->GetCost();
    }
    return cost;
}

template <typename creature>
double Evolution<creature>::TotalCost() {
    double cost = 0;
    for (int i = 0; i < this->PopSize(); i++) {
        cost += this->GetCreature(i)->GetCost();
    }
    return cost;
}

template <typename creature>
creature* Evolution<creature>::GetCreature(int i) {
    return &(this->population->at(i));
}

template <typename creature>
std::vector<double> Evolution<creature>::fitness() {
    std::vector<double> costs = std::vector<double>(this->nFittest());
    double totalCost = 0;
    for (int i = 0; i < this->nFittest(); i++) {
        totalCost += (*this->population)[i].GetCost();
        costs[i] = (*this->population)[i].GetCost();
    }

    // Calculate the fitness vector, which indicates how fit a creature is,
    // higher fitness value -> better

    // Start by filing with cost values
    std::vector<double> fitness = std::vector<double>(this->nFittest());
    for (int i = 0; i < this->nFittest(); i++) {
        fitness[i] = costs[i];
    }
    // Remove negative values if needed
    // fitness: [-4,-1,0,2] -> [1,4,5,7]
    // fitness: [0,1,2] -> [1,2,3]
    // fitness: [4,6,7] -> [4,6,7]
    double minCost = costs[0];
    if (minCost <= 0) {
        for (int i = 0; i < this->nFittest(); i++) {
            fitness[i] += -minCost + 1;
        }
    }
    // Invert values
    // fitness: [1,2,3] -> [1, 0.5, 0.33]
    // fitness: [4,6,7] -> [0.25, 0.16, 0.14]
    for (int i = 0; i < this->nFittest(); i++) {
        fitness[i] = 1 / fitness[i];
    }

    // Normalize values
    double totalFitness = 0;
    for (int i = 0; i < this->nFittest(); i++) {
        totalFitness += fitness[i];
    }
    for (int i = 0; i < this->nFittest(); i++) {
        fitness[i] = fitness[i] / totalFitness;
    }
    return fitness;
}

template <typename creature>
std::tuple<creature*, creature*> Evolution<creature>::getParents() {
    // Biased roulette wheel method
    std::vector<double> f = this->fitness();
    std::vector<double> roulette = std::vector<double>(this->nFittest());
    double accumulated = 0.0;
    for (int i = 0; i < this->nFittest(); i++) {
        roulette[i] = f[i] + accumulated;
        accumulated += f[i];
    }

    auto it = lower_bound(roulette.begin(), roulette.end(), Random());
    int p1 = it - roulette.begin();
    it = lower_bound(roulette.begin(), roulette.end(), Random());
    int p2 = it - roulette.begin();

    // If got the same parent twice from roulette,
    // get the next for the second parent.
    if (p2 == p1) {
        if (p2 == this->endFittest()) {
            p2 = 0;
        } else {
            p2 += 1;
        }
    }

    return std::make_tuple(&((*this->population)[p1]),
                           &((*this->population)[p2]));
};

template <typename creature>
void Evolution<creature>::mutate() {
    if (this->population->size() == 1) {
        // the first creature is excluded from mutation,
        // so we require at least 2 creatures for mutations
        return;
    }

    // Total number of slots that can be mutated
    // (the fittest creatures are not mutated)
    int nDnaSlots =
        this->population->size() * (1-this->survival) * ((*this->population)[0].DnaSize());
    // target num of mutations
    int nMutations = int(nDnaSlots * this->mutationRate);
    // current mutation count
    int mutations = 0;

    // Set containing [creatureId, dnaPosition] pairs which were mutated
    std::set<std::tuple<int, int>> mutationsDone =
        std::set<std::tuple<int, int>>();
    while (mutations < nMutations) {
        // pick random creature from the less fit ones
        // The fittest ones do not suffer mutation
        int creatureI =
            RandomInt(this->endFittest() + 1, this->population->size() - 1);
        int dnaI = RandomInt(0, (*this->population)[creatureI].DnaSize() - 1);
        auto pair = std::tuple<int, int>(creatureI, dnaI);
        if (mutationsDone.find(pair) == mutationsDone.end()) {
            (*this->population)[creatureI].Mutate(dnaI);
            mutations += 1;
            mutationsDone.insert(pair);
        }
    }
}

template <typename creature>
void Evolution<creature>::step() {
    // Replace less fit of the population (assuming sortPopulation was already
    // called)
    for (int i = this->endFittest() + 1; i < this->PopSize(); i++) {
        auto parents = this->getParents();
        creature* p0 = std::get<0>(parents);
        creature* p1 = std::get<1>(parents);

        creature* child0;
        creature* child1;
        child0 = this->GetCreature(i);
        // If we still have 2 or more creatures to replace, get two creatures.
        // Else, set both childs to same child.
        if (i + 1 < this->PopSize()) {
            child1 = this->GetCreature(i + 1);
        } else {
            child1 = child1;
        }

        p0->Mate(*p1, child0, child1);
    }

    this->mutate();

    this->sortPopulation();
}

template <typename creature>
Maybe<creature*> Evolution<creature>::Evolve(double stop,
                                             bool printProgression) {
    Maybe<creature*> r;
    if (stop > 1 || stop < 0) {
        r.isError = true;
        r.errMsg = "stop must be 0<=stop<=1";
        return r;
    }

    // step() requires sortPopulation() to be called first
    this->sortPopulation();
    double pastCost = -1.0;
    double currentCost = this->FittestCost();

    int i = 0;
    const int maxI = 1000000;
    if (printProgression) {
        std::cout << "generation,cost" << std::endl;
    }

    while (i == 0 || RelativeAbsError(currentCost, pastCost) > stop) {
        if (printProgression) {
            std::cout << i << "," << currentCost << std::endl;
        }
        pastCost = currentCost;
        this->step();
        currentCost = this->FittestCost();

        if (i == maxI) {
            r.isError = true;
            r.errMsg = "Evolution did not converge";
            break;
        }
        i += 1;
    }
    if (printProgression) {
        std::cout << i << "," << currentCost << std::endl;
    }

    r.val = this->GetCreature(0);
    return r;
}

template <typename creature>
Maybe<creature*> Evolution<creature>::Evolve(int generations,
                                             bool printProgression) {
    Maybe<creature*> r;
    if (generations < 1) {
        r.isError = true;
        r.errMsg = "generations must be >=1";
        return r;
    }

    // step() requires sortPopulation() to be called first
    this->sortPopulation();
    int i = 0;
    if (printProgression) {
        std::cout << "generation,cost" << std::endl;
    }
    while (i < generations) {
        if (printProgression) {
            std::cout << i << "," << this->FittestCost() << std::endl;
        }
        this->step();
        i += 1;
    }
    if (printProgression) {
        std::cout << i << "," << this->FittestCost() << std::endl;
    }

    r.val = this->GetCreature(0);
    return r;
}