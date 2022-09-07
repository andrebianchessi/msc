#include "creature_test.h"

#include <gtest/gtest.h>

#include <memory>

using namespace std;
TEST(CreatureTest, SimpleTest) {
    // Create creature with
    // (xNormalized,yNormalized) = (1,0)
    // (x,y) = (2,-2)
    Bounded x = Bounded();
    Bounded y = Bounded();
    x.Set(1.0);
    EXPECT_DOUBLE_EQ(x.Get(), 1.0);
    EXPECT_DOUBLE_EQ(y.Get(), 0.0);

    EqSol c0 = EqSol(x, y);

    EXPECT_DOUBLE_EQ(c0.dna[0].Get(), 1.0);
    EXPECT_DOUBLE_EQ(c0.dna[1].Get(), 0.0);

    // Verify creature is constructed by value
    x.Set(0.0);
    EXPECT_DOUBLE_EQ(c0.dna[0].Get(), 1.0);

    EXPECT_EQ(c0.DnaSize(), 2);

    // f(x,y) = x^2 + y^2 + 2x + y
    // f(2,-2) = 2^2 + (-2)^2 + 2*2 -2
    // f(2,-2) = 10
    EXPECT_DOUBLE_EQ(c0.GetCost(), 10.0);
}

TEST(CreatureTest, MateAndMutationTest) {
    Bounded x = Bounded();
    Bounded y = Bounded();
    x.Set(1.0);

    EqSol c0 = EqSol(x, y);
    x.Set(0.0);
    EqSol c1 = EqSol(x, y);
    EXPECT_DOUBLE_EQ(c0.dna[0].Get(), 1.0);
    EXPECT_DOUBLE_EQ(c0.dna[1].Get(), 0.0);
    EXPECT_DOUBLE_EQ(c1.dna[0].Get(), 0.0);
    EXPECT_DOUBLE_EQ(c1.dna[1].Get(), 0.0);

    // c0: [1.0, 0.0]
    // c1: [0.0, 0.0]
    //
    // Expected:
    // x in [0,1]
    // child0: [x, 0.0]
    // child1: [x, 0.0]
    EqSol child0 = EqSol(x, y);
    EqSol child1 = EqSol(x, y);

    c0.Mate(c1, &child0, &child1);

    EXPECT_TRUE(0.0 <= child0.dna[0].Get() && child0.dna[0].Get() <= 1.0);
    EXPECT_TRUE(0.0 <= child1.dna[0].Get() && child0.dna[0].Get() <= 1.0);
    EXPECT_TRUE(child1.dna[1].Get() == 0.0);
    EXPECT_TRUE(child0.dna[1].Get() == 0.0);

    c0.Mutate();
    // print("dna: ", c0.dna[0].Get(), c0.dna[1].Get());
    bool changed = c0.dna[0].Get() != 1.0 || c0.dna[1].Get() != 0.0;
    EXPECT_TRUE(changed);
}