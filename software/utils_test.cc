#include "utils.h"

#include <gtest/gtest.h>

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
}