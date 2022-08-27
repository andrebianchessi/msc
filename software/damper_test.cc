#include <gtest/gtest.h>
#include "damper.h"
#include "problem.h"
#include "mass.h"

TEST(DamperTest, SimpleTest) {
    Problem* p = new Problem();
    p->AddMass(10.0,0.0,0.0);
    p->AddMass(20.0,1.0,1.0);

    p->AddDamper(0,1,30.0);
    
    auto e = p->GetDamper(0);
    ASSERT_FALSE(e.isError);

    auto s = e.val;
    
    auto M = s->GetM();
    ASSERT_DOUBLE_EQ(M(0,0),10.0);
    ASSERT_DOUBLE_EQ(M(0,1),0.0);
    ASSERT_DOUBLE_EQ(M(1,0),0.0);
    ASSERT_DOUBLE_EQ(M(1,1),20.0);

    auto C = s->GetC();
    ASSERT_DOUBLE_EQ(C(0,0),-30.0);
    ASSERT_DOUBLE_EQ(C(0,1),30.0);
    ASSERT_DOUBLE_EQ(C(1,0),30.0);
    ASSERT_DOUBLE_EQ(C(1,1),-30.0);
}

TEST(DamperTest, ErrorsTest) {
    Problem* p = new Problem();
    p->AddMass(10.0,0.0,0.0);
    p->AddMass(10.0,1.0,1.0);

    auto e = p->AddDamper(0,0,30.0);
    ASSERT_TRUE(e.isError);

    e = p->AddDamper(0,99,30.0);
    ASSERT_TRUE(e.isError);

    e = p->AddDamper(99,0,30.0);
    ASSERT_TRUE(e.isError);
}