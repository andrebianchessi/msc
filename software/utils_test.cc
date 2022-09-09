#include "utils.h"

#include <gtest/gtest.h>

#include "bounded.h"

TEST(UtilsTest, RandomTest) {
    for (int i = 0; i < 20; i++) {
        double r = Random();
        ASSERT_TRUE(r >= 0 && r <= 1);
    }
    for (int i = 0; i < 20; i++) {
        double r = Random(-2.0, -1.0);
        ASSERT_TRUE(r >= -2.0 && r <= 1.0);
    }
    for (int i = 0; i < 20; i++) {
        double r = Random(1.0, 10.0);
        ASSERT_TRUE(r >= 1 && r <= 10.0);
    }
    for (int i = 0; i < 20; i++) {
        double r = Random(-0.1, 0.1);
        ASSERT_TRUE(r >= -0.1 && r <= 0.1);
    }
    for (int i = 0; i < 20; i++) {
        int ri = RandomInt(-10, 10);
        ASSERT_TRUE(ri >= -10 && ri <= 10);
    }
}

TEST(UtilsTest, NormalizeTest) {
    EXPECT_TRUE(Normalize(0.0, 1.0, 2.0).isError);
    EXPECT_TRUE(Normalize(2.1, 1.0, 2.0).isError);

    auto e = Normalize(0.0, -1.0, 1.0);
    EXPECT_FALSE(e.isError);
    EXPECT_DOUBLE_EQ(e.val.Get(), 0.5);

    e = Normalize(1.0, 1.0, 11.0);
    EXPECT_FALSE(e.isError);
    EXPECT_DOUBLE_EQ(e.val.Get(), 0.0);

    e = Normalize(11.0, 1.0, 11.0);
    EXPECT_FALSE(e.isError);
    EXPECT_DOUBLE_EQ(e.val.Get(), 1.0);

    e = Normalize(0.0, -1.0, 2.0);
    EXPECT_FALSE(e.isError);
    EXPECT_DOUBLE_EQ(e.val.Get(), 1 / 3.0);
}

TEST(UtilsTest, UnnormalizeTest) {
    Bounded b = Bounded();

    b.Set(0.0);
    ASSERT_DOUBLE_EQ(Unnormalize(b, -1, 1), -1.0);

    b.Set(1.0);
    ASSERT_DOUBLE_EQ(Unnormalize(b, -1, 1), 1.0);

    b.Set(1 / 3.0);
    ASSERT_DOUBLE_EQ(Unnormalize(b, -2, 1), -1.0);
}