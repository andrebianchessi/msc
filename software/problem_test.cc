#include <gtest/gtest.h>
#include "problem.h"
#include <iostream>

TEST(ProblemTest, MassCreationTest) {
  Problem* p = new Problem();

  auto e = p->AddMass(0.0,0.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 1);
  EXPECT_TRUE(e.val == 0);

  e = p->AddMass(0.0,1.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 2);
  EXPECT_TRUE(e.val == 1);

  e = p->AddMass(0.0,1.0,1.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 3);
  EXPECT_TRUE(e.val == 2);

  e = p->AddMass(0.0,0.0,0.0);
  EXPECT_TRUE(e.isError);
  EXPECT_TRUE(p->GetDof() == 3);
}

TEST(ProblemTest, GetMassTest) {
  Problem* p = new Problem();

  p->AddMass(0.0,0.0,0.1);
  p->AddMass(0.0,1.0,1.1);
  p->AddMass(0.0,2.0,2.1);

  auto e = p->GetMass(0);
  EXPECT_FALSE(e.isError);
  auto m = e.val;
  EXPECT_TRUE(m->xIndex == 0);
  EXPECT_TRUE(m->x == 0.0);
  EXPECT_TRUE(m->y == 0.1);

  e = p->GetMass(2);
  EXPECT_FALSE(e.isError);
  m = e.val;
  EXPECT_TRUE(m->xIndex == 2);
  EXPECT_TRUE(m->x == 2.0);
  EXPECT_TRUE(m->y == 2.1);

  e = p->GetMass(3);
  EXPECT_TRUE(e.isError);
}

TEST(ProblemTest, SimpleBuildTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,9.0);
  p.Build();
  p.Build(); // Calling again to make sure matrices are reset

  auto MInv = p.MInv;
  EXPECT_EQ(MInv.size1(), 2);
  EXPECT_EQ(MInv.size2(), 2);
  EXPECT_DOUBLE_EQ(MInv(0,0), 1.0);
  EXPECT_DOUBLE_EQ(MInv(0,1), 0.0);
  EXPECT_DOUBLE_EQ(MInv(1,0), 0.0);
  EXPECT_DOUBLE_EQ(MInv(1,1), 0.5);

  auto K = p.K;
  EXPECT_EQ(K.size1(), 2);
  EXPECT_EQ(K.size2(), 2);
  EXPECT_DOUBLE_EQ(K(0,0), -9.0);
  EXPECT_DOUBLE_EQ(K(0,1), 9.0);
  EXPECT_DOUBLE_EQ(K(1,0), 9.0);
  EXPECT_DOUBLE_EQ(K(1,1), -9.0);

  auto X = p.X;
  EXPECT_EQ(X.size(), 2);
  EXPECT_DOUBLE_EQ(X(0), 0.0);
  EXPECT_DOUBLE_EQ(X(1), 0.0);

  auto XDot = p.XDot;
  EXPECT_EQ(XDot.size(), 2);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0);
  EXPECT_DOUBLE_EQ(XDot(1), 0.0);
}