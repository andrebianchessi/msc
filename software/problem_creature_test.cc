#include "problem_creature.h"

#include <gtest/gtest.h>

#include <memory>
#include <vector>

#include "bounded.h"
#include "evolution.h"
#include "utils.h"

TEST(ProblemCreatureTest, SimpleTest) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(1.0, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    pd.SetFixedMass(0);
    pd.AddInitialVel(1, 10.0);
    ASSERT_TRUE(pd.IsOk());

    // Since the first mass is fixed, its accelerations will be zero
    ProblemCreature c0 = ProblemCreature(&pd, 0, 0.5);
    ASSERT_FALSE(c0.hasGetCostCache);
    ASSERT_DOUBLE_EQ(c0.GetCost(), 0.0);
    ASSERT_TRUE(c0.hasGetCostCache);

    ProblemCreature c1 = ProblemCreature(&pd, 1, 0.5);
    ASSERT_FALSE(c1.hasGetCostCache);
    ASSERT_TRUE(c1.GetCost() > 0.0);
    ASSERT_TRUE(c1.hasGetCostCache);

    std::vector<ProblemCreature> pop = std::vector<ProblemCreature>();
    pop.push_back(ProblemCreature(&pd, 1, 0.5));
    pop.push_back(ProblemCreature(&pd, 1, 0.5));
    Evolution<ProblemCreature> evolution = Evolution<ProblemCreature>(&pop);
    ASSERT_NO_THROW(evolution.GetCreature(0)->GetCost());
    ASSERT_NO_THROW(evolution.GetCreature(1)->GetCost());
}

TEST(ProblemCreatureTest, EvolutionTest) {
    // Optimize problem similar to MultiBodyBibliographyDataTest2
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 1.0, 0.0);  // m2
    pd.AddMass(150, 1.0, 3.0);  // m3
    pd.AddMass(700, 2.0, 0.0);  // m4
    pd.AddMass(80, 3.0, 0.0);   // m5

    double min = 100.0;
    double max = 100000;
    pd.AddSpring(0, 1, min, max);
    pd.AddSpring(1, 2, min, max);
    pd.AddSpring(1, 3, min, max);
    pd.AddSpring(1, 4, min, max);
    pd.AddDamper(1, 4, min, max);
    pd.AddSpring(0, 2, min, max);
    pd.AddDamper(0, 2, min, max);
    pd.AddSpring(2, 4, min, max);
    pd.AddDamper(2, 4, min, max);
    pd.AddSpring(0, 3, min, max);
    pd.AddDamper(0, 3, min, max);
    pd.AddSpring(3, 4, min, max);
    pd.AddDamper(3, 4, min, max);
    pd.AddSpring(4, 5, min, max);
    pd.AddDamper(4, 5, min, max);

    pd.SetFixedMass(0);
    pd.AddInitialVel(200.0);
    ASSERT_TRUE(pd.IsOk());

    // Create population of 50 creatures
    std::vector<ProblemCreature> pop = std::vector<ProblemCreature>();
    for (int i = 0; i < 50; i++) {
        pop.push_back(ProblemCreature(&pd, 5, 0.15));
    }

    Evolution<ProblemCreature> evolution = Evolution<ProblemCreature>(&pop);
    // Evolve multiple times and check that cost is always <=
    for (int i = 0; i < 5; i++) {
        double cost0 = evolution.FittestCost();
        auto p = evolution.Evolve(1, false);
        ASSERT_FALSE(p.isError);
        ASSERT_TRUE(evolution.FittestCost() <= cost0);
    }
}

TEST(ProblemCreatureTest, EvolutionUntilConvergenceTest) {
    // Similar to EvolutionTest but using Evolve(float, true)
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 1.0, 0.0);  // m2
    pd.AddMass(150, 1.0, 3.0);  // m3
    pd.AddMass(700, 2.0, 0.0);  // m4
    pd.AddMass(80, 3.0, 0.0);   // m5

    double min = 100.0;
    double max = 100000;
    pd.AddSpring(0, 1, min, max);
    pd.AddSpring(1, 2, min, max);
    pd.AddSpring(1, 3, min, max);
    pd.AddSpring(1, 4, min, max);
    pd.AddDamper(1, 4, min, max);
    pd.AddSpring(0, 2, min, max);
    pd.AddDamper(0, 2, min, max);
    pd.AddSpring(2, 4, min, max);
    pd.AddDamper(2, 4, min, max);
    pd.AddSpring(0, 3, min, max);
    pd.AddDamper(0, 3, min, max);
    pd.AddSpring(3, 4, min, max);
    pd.AddDamper(3, 4, min, max);
    pd.AddSpring(4, 5, min, max);
    pd.AddDamper(4, 5, min, max);

    pd.SetFixedMass(0);
    pd.AddInitialVel(200.0);
    ASSERT_TRUE(pd.IsOk());

    // Create population of 20 creatures
    std::vector<ProblemCreature> pop = std::vector<ProblemCreature>();
    for (int i = 0; i < 20; i++) {
        pop.push_back(ProblemCreature(&pd, 5, 0.15));
    }

    Evolution<ProblemCreature> evolution = Evolution<ProblemCreature>(&pop);
    double cost0 = evolution.FittestCost();
    auto p = evolution.Evolve(0.1, true);
    ASSERT_FALSE(p.isError);
    ASSERT_TRUE(evolution.FittestCost() <= cost0);

    Problem best = pd.BuildFromDNA(evolution.GetCreature(0)->dna).val;
    print("k1: ", best.springs[0].Get_k());
    print("k2: ", best.springs[1].Get_k());
    print("k3: ", best.springs[2].Get_k());
    print("k4: ", best.springs[3].Get_k());
    print("c4: ", best.dampers[0].Get_c());
    print("k5: ", best.springs[4].Get_k());
    print("c5: ", best.dampers[1].Get_c());
}