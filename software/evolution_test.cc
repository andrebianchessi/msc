#include "evolution.h"

#include <gtest/gtest.h>

#include "creature_test.h"
using namespace std;
TEST(EvolutionTest, SimpleTest) {
    // Create 10 random creatures
    vector<C> pop = vector<C>();
    for (int i = 0; i < 10; i++) {
        pop.push_back(C());
    }

    // Instantiate Evolution
    Evolution<C> ev = Evolution<C>(&pop);
    ASSERT_EQ(ev.PopSize(), 10);

    // Sort
    ev.sortPopulation();
    for (int i = 0; i < int(ev.PopSize()) - 1; i++) {
        ASSERT_TRUE(ev.population->at(i).GetCost() <=
                    ev.population->at(i + 1).GetCost());
    }
}

TEST(EvolutionTest, EndFittestTest) {
    // Create 10 random creatures
    vector<C> pop = vector<C>();
    for (int i = 0; i < 10; i++) {
        pop.push_back(C());
    }

    // Instantiate Evolution
    Evolution<C> ev = Evolution<C>(&pop);

    ASSERT_EQ(ev.endFittest(), 4);
    ASSERT_EQ(ev.nFittest(), 5);

    // Test other values of survival
    ev.survival = 0.33;
    ASSERT_EQ(ev.endFittest(), 2);
    ASSERT_EQ(ev.nFittest(), 3);
    ev.survival = 0.1;
    ASSERT_EQ(ev.endFittest(), 0);
    ASSERT_EQ(ev.nFittest(), 1);
}

TEST(EvolutionTest, fitnessTest) {
    // Create population of size 4 in which all creatures have
    // DNA = [0.5, 0.5]
    vector<C> pop = vector<C>();
    Bounded xB = Bounded();
    xB.Set(0.5);
    // xB = 0.5 -> x = 0.0;
    for (int i = 0; i < 4; i++) {
        pop.push_back(C(xB, xB));
    }
    Evolution<C> ev = Evolution<C>(&pop);
    ASSERT_EQ(ev.PopSize(), 4);

    // 3 creatures will survive to next generation
    ev.survival = 3 / 4.0;

    // costs = [0,0,0,0]
    // fitness = [1/3,1/3,1/3]
    auto f = ev.fitness();
    ASSERT_EQ(f.size(), 3);
    ASSERT_DOUBLE_EQ(f[0], 1 / 3.0);
    ASSERT_DOUBLE_EQ(f[1], 1 / 3.0);
    ASSERT_DOUBLE_EQ(f[2], 1 / 3.0);

    // Test other values of dna
    // Set xB0 = 0.25 -> x = -1
    ev.GetCreature(0)->dna[0].Set(0.25);
    // Set xB1 = 0.75 -> x = 1
    // Set yB1 = 0.75 -> y = 1
    ev.GetCreature(1)->dna[0].Set(0.75);
    ev.GetCreature(1)->dna[1].Set(0.75);
    // f(x,y) = x^2 + y^2 + 2x + y
    // cost = [|1-2|, 5, 0, 0]
    // cost = [1, 5, 0, 0]
    ev.sortPopulation();
    // cost = [0, 0, 1, 5]
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->GetCost(), 0.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(1)->GetCost(), 0.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(2)->GetCost(), 1.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(3)->GetCost(), 5.0);
    // fitness = [1,1,2] -> [1,1,1/2] -> [1/(1+1+1/2), ...]
    f = ev.fitness();
    ASSERT_EQ(f.size(), 3);
    ASSERT_DOUBLE_EQ(f[0], 1 / (1 + 1 + 1 / 2.0));
    ASSERT_DOUBLE_EQ(f[1], 1 / (1 + 1 + 1 / 2.0));
    ASSERT_DOUBLE_EQ(f[2], 1 / 2.0 / (1 + 1 + 1 / 2.0));
}

TEST(EvolutionTest, getParentsTest) {
    vector<C> pop = vector<C>();
    Bounded x = Bounded();
    Bounded y = Bounded();

    // cost = 1
    x.Set(0.25);
    y.Set(0.25);
    pop.push_back(C(x, y));

    // cost = 0
    x.Set(0.5);
    y.Set(0.5);
    pop.push_back(C(x, y));

    // cost = 2
    x.Set(0.0);
    y.Set(0.0);
    pop.push_back(C(x, y));

    // cost = 14
    x.Set(1.0);
    y.Set(1.0);
    pop.push_back(C(x, y));

    Evolution<C> ev = Evolution<C>(&pop);
    ev.survival = 3 / 4.0;
    ev.sortPopulation();

    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->GetCost(), 0.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(1)->GetCost(), 1.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(2)->GetCost(), 2.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(3)->GetCost(), 14.0);

    // Number of parents with cost X selected each time by getParents
    int cost0ParentCount = 0;
    int cost1ParentCount = 0;
    int cost2ParentCount = 0;
    int cost14ParentCount = 0;

    for (int i = 0; i < 100000; i++) {
        auto e = ev.getParents();
        C* p1 = get<0>(e);
        if (abs(p1->GetCost() - 0.0) < 0.005f) {
            cost0ParentCount += 1;
        }
        if (abs(p1->GetCost() - 1.0) < 0.005f) {
            cost1ParentCount += 1;
        }
        if (abs(p1->GetCost() - 2.0) < 0.005f) {
            cost2ParentCount += 1;
        }
        if (abs(p1->GetCost() - 14.0) < 0.005f) {
            cost14ParentCount += 1;
        }
    }

    ASSERT_TRUE(cost0ParentCount > cost1ParentCount);
    ASSERT_TRUE(cost1ParentCount > cost2ParentCount);
    ASSERT_EQ(cost14ParentCount, 0);

    // fitness = [1/k,(1/2)/k,(1/3)/k]
    // Frequency of choosing each parent is proportional to fitness:
    // freqCost0 = alpha*1
    // freqCost1 = alpha*1/2
    // freqCost2 = alpha*1/3
    // -> cost0Count/cost1Count = 2.0
    // -> cost0Count/cost2Count = 3.0
    // -> cost1Count/cost2Count = 1.5
    ASSERT_TRUE(abs((float(cost0ParentCount) / cost1ParentCount - 2.0) / 2.0) <
                0.05);
    ASSERT_TRUE(abs((float(cost0ParentCount) / cost2ParentCount - 3.0) / 3.0) <
                0.05);
    ASSERT_TRUE(abs((float(cost1ParentCount) / cost2ParentCount - 1.5) / 1.5) <
                0.05);
}

TEST(EvolutionTest, MutateTest) {
    // Create 4 creatures with DNA = [0.0, 0.0]
    vector<C> pop = vector<C>();
    Bounded x = Bounded();
    for (int i = 0; i < 4; i++) {
        pop.push_back(C(x, x));
    }
    Evolution<C> ev = Evolution<C>(&pop);
    ASSERT_EQ(ev.PopSize(), 4);

    // Set invalid values for DNAs
    // This is done just to make sure which values were edited by the mutate
    // function, and should never be done outside tests.
    for (int i = 0; i < int(ev.PopSize()); i++) {
        ev.GetCreature(i)->dna[0].val = -99.0;
        ev.GetCreature(i)->dna[1].val = -99.0;
    }

    ev.mutationRate = 0.5;
    ev.mutate();

    // Best is not changed
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[1].val, -99.0);

    int changeCount = 0;
    for (int i = 0; i < ev.PopSize(); i++) {
        if (ev.GetCreature(i)->dna[0].val != -99.0) {
            changeCount += 1;
        }
        if (ev.GetCreature(i)->dna[1].val != -99.0) {
            changeCount += 1;
        }
    }

    // 4 creatures with 2 dna slots each
    // 1 is not mutated (elitism)
    // 3 remain -> 6 slots
    // 0.5 mutation rate -> 3 changes
    ASSERT_EQ(changeCount, 3);

    // Repeat for different mutation rates: 2/3
    for (int i = 0; i < ev.PopSize(); i++) {
        ev.GetCreature(i)->dna[0].val = -99.0;
        ev.GetCreature(i)->dna[1].val = -99.0;
    }
    ev.mutationRate = 2 / 3.0;
    ev.mutate();
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[1].val, -99.0);
    changeCount = 0;
    for (int i = 0; i < ev.PopSize(); i++) {
        if (ev.GetCreature(i)->dna[0].val != -99.0) {
            changeCount += 1;
        }
        if (ev.GetCreature(i)->dna[1].val != -99.0) {
            changeCount += 1;
        }
    }
    ASSERT_EQ(changeCount, 4);

    // Repeat for different mutation rates: 0%
    for (int i = 0; i < ev.PopSize(); i++) {
        ev.GetCreature(i)->dna[0].val = -99.0;
        ev.GetCreature(i)->dna[1].val = -99.0;
    }
    ev.mutationRate = 0;
    ev.mutate();
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ(ev.GetCreature(0)->dna[1].val, -99.0);
    changeCount = 0;
    for (int i = 0; i < ev.PopSize(); i++) {
        if (ev.GetCreature(i)->dna[0].val != -99.0) {
            changeCount += 1;
        }
        if (ev.GetCreature(i)->dna[1].val != -99.0) {
            changeCount += 1;
        }
    }
    ASSERT_EQ(changeCount, 0);
}

TEST(EvolutionTest, TotalCostCount) {
    // Create 4 creatures
    vector<C> pop = vector<C>();
    Bounded x = Bounded();
    Bounded y = Bounded();

    // cost = 1
    x.Set(0.25);
    y.Set(0.25);
    pop.push_back(C(x, y));

    // cost = 0
    x.Set(0.5);
    y.Set(0.5);
    pop.push_back(C(x, y));

    // cost = 2
    x.Set(0.0);
    y.Set(0.0);
    pop.push_back(C(x, y));

    // cost = 14
    x.Set(1.0);
    y.Set(1.0);
    pop.push_back(C(x, y));

    Evolution<C> ev = Evolution<C>(&pop);

    ASSERT_DOUBLE_EQ(ev.TotalCost(), 1 + 0 + 2 + 14);

    ev.survival = 0.5;
    ASSERT_DOUBLE_EQ(ev.FittestCost(), 1 + 0);
}

TEST(EvolutionTest, StepTest) {
    // Create random population of size 21
    // (Using odd number to check if it all works)
    vector<C> pop = vector<C>();
    for (int i = 0; i < 21; i++) {
        pop.push_back(C());
    }
    Evolution<C> ev = Evolution<C>(&pop);

    // Perform multiple evolution steps and check that total cost and
    // fittest cost are reduced
    double tc0 = ev.TotalCost();
    double fc0 = ev.FittestCost();
    for (int i = 0; i < 1000; i++) {
        ev.step();
    }
    double tc1 = ev.TotalCost();
    double fc1 = ev.FittestCost();
    ASSERT_TRUE(tc1 < tc0);
    ASSERT_TRUE(fc1 < fc0);

    C* bestCreature = ev.GetCreature(0);
    double bestCost = bestCreature->GetCost();
    // f(x,y) = x^2 + y^2 + 2x + y
    // Theoretical global min = 0.0
    double bestX = bestCreature->GetX();
    double bestY = bestCreature->GetY();
    ASSERT_DOUBLE_EQ(bestCost,
                     abs(bestX * bestX + bestY * bestY + 2 * bestX + bestY));
    EXPECT_TRUE(abs(bestCost) < 0.01);
    // print("bestCost, bestX, bestY: ", bestCost, bestCreature->GetX(),
    //   bestCreature->GetY());
}

TEST(EvolutionTest, EvolveTest) {
    // Create random population of size 30
    vector<C> pop = vector<C>();
    for (int i = 0; i < 30; i++) {
        pop.push_back(C());
    }
    Evolution<C> ev = Evolution<C>(&pop);

    // Perform multiple evolution steps and check that total cost and
    // fittest cost are reduced
    double tc0 = ev.TotalCost();
    double fc0 = ev.FittestCost();

    auto e = ev.Evolve(0.001, true);
    ASSERT_FALSE(e.isError);

    C* best = e.val;
    double tc1 = ev.TotalCost();
    double fc1 = ev.FittestCost();

    ASSERT_TRUE(tc1 < tc0);
    ASSERT_TRUE(fc1 < fc0);

    EXPECT_TRUE(abs(best->GetCost()) < 0.1);
}