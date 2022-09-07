#include "evolution.h"

#include <gtest/gtest.h>

#include "creature_test.h"
using namespace std;
TEST(EvolutionTest, SimpleTest) {
    // Create 10 random creatures
    shared_ptr<vector<EqSol>> pop = make_shared<vector<EqSol>>();
    *pop = vector<EqSol>();
    for (int i = 0; i < 10; i++) {
        pop->push_back(EqSol());
    }
    ASSERT_EQ(pop->size(), 10);

    // Instantiate Evolution
    Evolution<EqSol> ev = Evolution<EqSol>(pop);
    ASSERT_EQ(ev.popSize(), 10);

    // Sort
    ev.sortPopulation();
    for (int i = 0; i < int(pop->size()) - 1; i++) {
        ASSERT_TRUE((*pop)[i].GetCost() <= (*pop)[i + 1].GetCost());
    }
}

TEST(EvolutionTest, EndFittestTest) {
    // Create 10 random creatures
    shared_ptr<vector<EqSol>> pop = make_shared<vector<EqSol>>();
    *pop = vector<EqSol>();
    for (int i = 0; i < 10; i++) {
        pop->push_back(EqSol());
    }
    ASSERT_EQ(pop->size(), 10);

    // Instantiate Evolution
    Evolution<EqSol> ev = Evolution<EqSol>(pop);

    ASSERT_EQ(ev.endFittest(), 4);

    // Test other values of nKeep
    ev.nKeep = 0.33;
    ASSERT_EQ(ev.endFittest(), 2);
    ev.nKeep = 0.1;
    ASSERT_EQ(ev.endFittest(), 0);
}

TEST(EvolutionTest, fitnessTest) {
    // Create population with all 0.5 DNAs (which have zero cost)
    shared_ptr<vector<EqSol>> pop = make_shared<vector<EqSol>>();
    *pop = vector<EqSol>();
    Bounded x = Bounded();
    x.Set(0.5);
    // xNorm = 0.5 -> x = 0.0;
    for (int i = 0; i < 4; i++) {
        pop->push_back(EqSol(x, x));
    }
    ASSERT_EQ(pop->size(), 4);
    Evolution<EqSol> ev = Evolution<EqSol>(pop);

    // costs = [0,0,0,0,0]
    // fitness = [1/4,1/4,1/4,1/4]
    auto f = ev.fitness();
    ASSERT_EQ(f.size(), 4);
    ASSERT_DOUBLE_EQ(f[0], 1 / 4.0);
    ASSERT_DOUBLE_EQ(f[1], 1 / 4.0);
    ASSERT_DOUBLE_EQ(f[2], 1 / 4.0);
    ASSERT_DOUBLE_EQ(f[3], 1 / 4.0);

    // Test other values of dna
    // Set x0Norm = 0.25 -> x = -1
    (*pop)[0].dna[0].Set(0.25);
    // Set x1Norm = 0.75 -> x = 1
    // Set y1Norm = 0.75 -> y = 1
    (*pop)[1].dna[0].Set(0.75);
    (*pop)[1].dna[1].Set(0.75);
    // f(x,y) = x^2 + y^2 + 2x + y
    // cost = [1-2, 5, 0, 0]
    // cost = [-1, 5, 0, 0]
    ev.sortPopulation();
    // cost = [-1, 0, 0, 5]
    ASSERT_DOUBLE_EQ((*pop)[0].GetCost(), -1.0);
    ASSERT_DOUBLE_EQ((*pop)[1].GetCost(), 0.0);
    ASSERT_DOUBLE_EQ((*pop)[2].GetCost(), 0.0);
    ASSERT_DOUBLE_EQ((*pop)[3].GetCost(), 5.0);
    // fitness = [1,2,2,7] -> [1,1/2,1/2,1/7] -> [1/(1+1/2+1/2+1/7), ...]
    f = ev.fitness();
    ASSERT_EQ(f.size(), 4);
    ASSERT_DOUBLE_EQ(f[0], 1 / (1 + 1 / 2.0 + 1 / 2.0 + 1 / 7.0));
    ASSERT_DOUBLE_EQ(f[1], 1 / 2.0 / (1 + 1 / 2.0 + 1 / 2.0 + 1 / 7.0));
    ASSERT_DOUBLE_EQ(f[2], 1 / 2.0 / (1 + 1 / 2.0 + 1 / 2.0 + 1 / 7.0));
    ASSERT_DOUBLE_EQ(f[3], 1 / 7.0 / (1 + 1 / 2.0 + 1 / 2.0 + 1 / 7.0));
}

// TEST(EvolutionTest, getParentsTest) {
//     // Create 10 random creatures
//     shared_ptr<vector<EqSol>> pop = make_shared<vector<EqSol>>();
//     *pop = vector<EqSol>();
//     for (int i = 0; i < 10; i++) {
//         pop->push_back(EqSol());
//     }
//     ASSERT_EQ(pop->size(), 10);

//     // Instantiate Evolution
//     Evolution<EqSol> ev = Evolution<EqSol>(pop);

//     auto e = ev.getParents();

//     // ASSERT_EQ(ev.endFittest(), 4);

//     // // Test other values of nKeep
//     // ev.nKeep = 0.33;
//     // ASSERT_EQ(ev.endFittest(), 2);
//     // ev.nKeep = 0.1;
//     // ASSERT_EQ(ev.endFittest(), 0);
// }

TEST(EvolutionTest, MutateTest) {
    // Create 4 creatures with 0.0 values
    shared_ptr<vector<EqSol>> pop = make_shared<vector<EqSol>>();
    *pop = vector<EqSol>();
    Bounded x = Bounded();
    for (int i = 0; i < 4; i++) {
        pop->push_back(EqSol(x, x));
    }
    ASSERT_EQ(pop->size(), 4);

    // Instantiate Evolution
    Evolution<EqSol> ev = Evolution<EqSol>(pop);

    // Set invalid values for dna
    // This is done just to make sure which values were edited by the mutate
    // function, and should never be done outside tests.
    for (int i = 0; i < int(pop->size()); i++) {
        (*pop)[i].dna[0].val = -99.0;
        (*pop)[i].dna[1].val = -99.0;
    }

    ev.mutationRate = 0.5;
    ev.mutate();

    // Best is not changed
    ASSERT_DOUBLE_EQ((*pop)[0].dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ((*pop)[0].dna[1].val, -99.0);

    int changeCount = 0;
    for (int i = 0; i < int(pop->size()); i++) {
        if ((*pop)[i].dna[0].val != -99.0) {
            changeCount += 1;
        }
        if ((*pop)[i].dna[1].val != -99.0) {
            changeCount += 1;
        }
    }

    // 4 creatures with 2 dna slots each
    // 1 is not mutated (elitism)
    // 3 remain -> 6 slots
    // 0.5 mutation rate -> 3 changes
    ASSERT_EQ(changeCount, 3);

    // Repeat for different mutation rates: 2/3
    for (int i = 0; i < int(pop->size()); i++) {
        (*pop)[i].dna[0].val = -99.0;
        (*pop)[i].dna[1].val = -99.0;
    }
    ev.mutationRate = 2 / 3.0;
    ev.mutate();
    ASSERT_DOUBLE_EQ((*pop)[0].dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ((*pop)[0].dna[1].val, -99.0);
    changeCount = 0;
    for (int i = 0; i < int(pop->size()); i++) {
        if ((*pop)[i].dna[0].val != -99.0) {
            changeCount += 1;
        }
        if ((*pop)[i].dna[1].val != -99.0) {
            changeCount += 1;
        }
    }
    ASSERT_EQ(changeCount, 4);

    // Repeat for different mutation rates: 0%
    for (int i = 0; i < int(pop->size()); i++) {
        (*pop)[i].dna[0].val = -99.0;
        (*pop)[i].dna[1].val = -99.0;
    }
    ev.mutationRate = 0;
    ev.mutate();
    ASSERT_DOUBLE_EQ((*pop)[0].dna[0].val, -99.0);
    ASSERT_DOUBLE_EQ((*pop)[0].dna[1].val, -99.0);
    changeCount = 0;
    for (int i = 0; i < int(pop->size()); i++) {
        if ((*pop)[i].dna[0].val != -99.0) {
            changeCount += 1;
        }
        if ((*pop)[i].dna[1].val != -99.0) {
            changeCount += 1;
        }
    }
    ASSERT_EQ(changeCount, 0);
}
