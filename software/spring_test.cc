#include <gtest/gtest.h>
#include "spring.h"
#include "problem.h"
#include "mass.h"

TEST(SpringTest, SimpleTest) {
    Problem* p = new Problem();
    p->AddMass(10.0,0.0,0.0);
    p->AddMass(20.0,1.0,1.0);

    p->AddSpring(0,1,30.0);
    
    auto e = p->GetSpring(0);
    ASSERT_FALSE(e.isError);

    auto s = e.val;
    
    auto M = s->GetM();
    ASSERT_DOUBLE_EQ(M(0,0),10.0);
    ASSERT_DOUBLE_EQ(M(0,1),0.0);
    ASSERT_DOUBLE_EQ(M(1,0),0.0);
    ASSERT_DOUBLE_EQ(M(1,1),20.0);

    auto K = s->GetK();
    ASSERT_DOUBLE_EQ(K(0,0),-30.0);
    ASSERT_DOUBLE_EQ(K(0,1),30.0);
    ASSERT_DOUBLE_EQ(K(1,0),30.0);
    ASSERT_DOUBLE_EQ(K(1,1),-30.0);
}

TEST(SpringTest, ErrorsTest) {
    Problem* p = new Problem();
    p->AddMass(10.0,0.0,0.0);
    p->AddMass(10.0,1.0,1.0);

    auto e = p->AddSpring(0,0,30.0);
    ASSERT_TRUE(e.isError);

    e = p->AddSpring(0,99,30.0);
    ASSERT_TRUE(e.isError);

    e = p->AddSpring(99,0,30.0);
    ASSERT_TRUE(e.isError);
}