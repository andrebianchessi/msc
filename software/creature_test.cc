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

    C c0 = C(x, y);

    EXPECT_DOUBLE_EQ(c0.dna[0].Get(), 1.0);
    EXPECT_DOUBLE_EQ(c0.dna[1].Get(), 0.0);

    // Verify creature is constructed by value
    x.Set(0.0);
    EXPECT_DOUBLE_EQ(c0.dna[0].Get(), 1.0);

    EXPECT_EQ(c0.DnaSize(), 2);

    // f(x,y) = x^2 + y^2 + 2x + y
    // f(2,-2) = 2^2 + (-2)^2 + 2*2 -2
    // f(2,-2) = 10
    ASSERT_FALSE(c0.hasGetCostCache);
    EXPECT_DOUBLE_EQ(c0.GetCost(), 10.0);
}

TEST(CreatureTest, MateAndMutationTest) {
    Bounded x = Bounded();
    Bounded y = Bounded();

    x.Set(1.0);
    y.Set(0.0);
    C c0 = C(x, y);

    x.Set(0.0);
    y.Set(0.0);
    C c1 = C(x, y);
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
    C child0 = C(x, y);
    C child1 = C(x, y);

    c0.Mate(c1, &child0, &child1);

    EXPECT_TRUE(0.0 <= child0.dna[0].Get() && child0.dna[0].Get() <= 1.0);
    EXPECT_TRUE(0.0 <= child1.dna[0].Get() && child0.dna[0].Get() <= 1.0);
    EXPECT_TRUE(child1.dna[1].Get() == 0.0);
    EXPECT_TRUE(child0.dna[1].Get() == 0.0);

    c0.Mutate(RandomInt(0, c0.DnaSize() - 1));
    // print("dna: ", c0.dna[0].Get(), c0.dna[1].Get());
    bool changed = c0.dna[0].Get() != 1.0 || c0.dna[1].Get() != 0.0;
    EXPECT_TRUE(changed);
}

TEST(CreatureTest, GetCostCacheTest) {
    Bounded x = Bounded();

    C c0 = C(x, x);
    ASSERT_FALSE(c0.hasGetCostCache);
    c0.GetCost();
    ASSERT_TRUE(c0.hasGetCostCache);

    C p0 = C(x, x);
    C p1 = C(x, x);
    p0.GetCost();
    p1.GetCost();

    C child0 = C(x, x);
    C child1 = C(x, x);

    child0.GetCost();
    child1.GetCost();

    ASSERT_TRUE(p0.hasGetCostCache);
    ASSERT_TRUE(p1.hasGetCostCache);
    ASSERT_TRUE(child0.hasGetCostCache);
    ASSERT_TRUE(child1.hasGetCostCache);
    p0.Mate(p1, &child0, &child1);
    ASSERT_TRUE(p0.hasGetCostCache);
    ASSERT_TRUE(p1.hasGetCostCache);
    ASSERT_FALSE(child0.hasGetCostCache);
    ASSERT_FALSE(child1.hasGetCostCache);

    p0.Mutate(RandomInt(0, c0.DnaSize()));
    ASSERT_FALSE(p0.hasGetCostCache);
}