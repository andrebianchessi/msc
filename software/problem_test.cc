#include <gtest/gtest.h>
#include "problem.h"
#include <iostream>

TEST(ProblemTest, MassCreationTest) {
  Problem* p = new Problem();

  auto e = p->AddMass(0.0,0.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 2);
  EXPECT_TRUE(e.val == 0);

  e = p->AddMass(0.0,1.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 4);
  EXPECT_TRUE(e.val == 1);

  e = p->AddMass(0.0,1.0,1.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 6);
  EXPECT_TRUE(e.val == 2);

  e = p->AddMass(0.0,0.0,0.0);
  EXPECT_TRUE(e.isError);
  EXPECT_TRUE(p->GetDof() == 6);
}

TEST(ProblemTest, GetMassTest) {
  Problem* p = new Problem();

  p->AddMass(0.0,0.0,0.0);
  p->AddMass(0.0,1.0,1.0);
  p->AddMass(0.0,2.0,2.0);

  auto e = p->GetMass(0);
  EXPECT_FALSE(e.isError);
  auto m = e.val;
  EXPECT_TRUE(m->id == 0);
  EXPECT_TRUE(m->x == 0.0);

  e = p->GetMass(2);
  EXPECT_FALSE(e.isError);
  m = e.val;
  EXPECT_TRUE(m->id == 2);
  EXPECT_TRUE(m->x == 2.0);

  e = p->GetMass(3);
  EXPECT_TRUE(e.isError);
}