#pragma once
#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "maybe.h"

using namespace std;
// Template class must be child of Creature
template <typename creature>
class Evolution {
   public:
    Evolution(vector<creature>* population);

    // Pointer to vector of creatures which this instance is optimizing
    vector<creature>* population;

    // Get creature by index
    creature* GetCreature(int i);

    // Performs Genetic Algorithm optimization until the FittestCost
    // difference is smaller than 'stop', i.e. until
    // abs((FittestCost_1-FittestCost_0)/FittestCost_0)<stop
    // Returns error if solution doesn't converge.
    // If printProgression, the progression of costs will be printed to stdout.
    Maybe<creature*> Evolve(double stop, bool printProgression);

    // Helper function which returns population size;
    int PopSize();

    // Returns the sum of the cost of the fittest creatures.
    // The number of creatures that we consider for this is determined
    // by Evolve.survival
    double FittestCost();
    // Returns the sum of the cost of all creatures
    double TotalCost();

   private:
    // Percentage of the population that will survive each generation
    double survival = 0.5;
    // Percentage of the dna slots that will be mutated;
    // Example:
    // creatures with DnaSize = 10
    // population of size 30
    // total DnaSlots = 10*30 = 300
    // 0.05*300 will be mutated
    double mutationRate = 0.05;
    FRIEND_TEST(EvolutionTest, TotalCostCount);

    // Sort population vector by increasing cost;
    // population[0] has the minimum cost after calling this method
    void sortPopulation();
    FRIEND_TEST(EvolutionTest, SimpleTest);

    // Returns the fitness vector;
    // WARNING: Assumes population vector is sorted.
    // This vector is made up of strictly positive numbers, and the sum
    // of all elements is 1.0
    // The numbers are in ascending order, since a higher fitness
    // means a better solution.
    vector<double> fitness();
    FRIEND_TEST(EvolutionTest, fitnessTest);

    // Get creatures from the population that will mate
    // A.K.A. "Select Mates" in the literature.
    // This function assumes the population is sorted, i.e. sortPopulation
    // must be called before calling this method;
    tuple<creature*, creature*> getParents();
    FRIEND_TEST(EvolutionTest, getParentsTest);

    // Returns the position (inclusive) in which we separate the creatures
    // who will survive to the next generation;
    // population[0] up to population[endFittest()] will survive to the next
    // generation.
    int endFittest();
    FRIEND_TEST(EvolutionTest, EndFittestTest);

    // Returns number of creatures that will survive to next generation
    int nFittest();

    // Randomly mutates some of the less fit creatures
    void mutate();
    FRIEND_TEST(EvolutionTest, MutateTest);

    // Performs one evolution step;
    // i.e. replace less fit with offsprings from the most fit part of the
    // population
    void step();
    FRIEND_TEST(EvolutionTest, StepTest);
};

#include "evolution.tcc"