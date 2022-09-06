#include "bounded.h"

#include <gtest/gtest.h>

#include "maybe.h"

TEST(BoundedTest, SimpleTest) {
    Bounded b = Bounded();

    ASSERT_DOUBLE_EQ(b.Get(), 0.0);

    ASSERT_TRUE(b.Set(1.1).isError);
    ASSERT_TRUE(b.Set(-0.1).isError);

    b.Set(0.6);
    ASSERT_DOUBLE_EQ(b.Get(), 0.6);
}