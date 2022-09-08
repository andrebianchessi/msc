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
    ProblemCreature c0 = ProblemCreature(&pd, 0, 0.5, 0.01);
    ASSERT_DOUBLE_EQ(c0.GetCost(), 0.0);

    ProblemCreature c1 = ProblemCreature(&pd, 1, 0.5, 0.01);
    ASSERT_TRUE(c1.GetCost() > 0.0);
}

TEST(ProblemCreatureTest, EvolutionTest) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(1.0, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    pd.SetFixedMass(0);
    pd.AddInitialVel(1, 10.0);
    ASSERT_TRUE(pd.IsOk());

    // Create population of 5 creatures
    std::vector<ProblemCreature> pop = std::vector<ProblemCreature>();
    for (int i = 0; i < 5; i++) {
        pop.push_back(ProblemCreature(&pd, 0, 0.5, 0.01));
    }

    Evolution<ProblemCreature> evolution = Evolution<ProblemCreature>(&pop);
    auto p = evolution.Evolve(0.01, true);
}