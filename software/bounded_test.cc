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

    Maybe<Bounded> b1 = Bounded::CreateBounded(0.7);
    ASSERT_FALSE(b1.isError);
    ASSERT_DOUBLE_EQ(b1.val.Get(), 0.7);

    ASSERT_TRUE(Bounded::CreateBounded(1.01).isError);
    ASSERT_TRUE(Bounded::CreateBounded(-0.01).isError);
}