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
  EXPECT_TRUE(m->px == 0.0);
  EXPECT_TRUE(m->py == 0.1);

  e = p->GetMass(2);
  EXPECT_FALSE(e.isError);
  m = e.val;
  EXPECT_TRUE(m->xIndex == 2);
  EXPECT_TRUE(m->px == 2.0);
  EXPECT_TRUE(m->py == 2.1);

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

  EXPECT_EQ(p.X.size(), 4);
  EXPECT_DOUBLE_EQ(p.X(0), 0.0);
  EXPECT_DOUBLE_EQ(p.X(1), 0.0);
  EXPECT_DOUBLE_EQ(p.X(2), 0.0);
  EXPECT_DOUBLE_EQ(p.X(3), 0.0);
}

TEST(ProblemTest, InitialConditionsTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,9.0);
  p.Build();

  p.SetInitialDisp(0,9.1);
  p.SetInitialDisp(1,9.2);
  p.SetInitialVel(0,99.1);
  p.SetInitialVel(1,99.2);
  EXPECT_EQ(p.X.size(), 4);
  EXPECT_DOUBLE_EQ(p.X(0), 9.1);
  EXPECT_DOUBLE_EQ(p.X(1), 9.2);
  EXPECT_DOUBLE_EQ(p.X(2), 99.1);
  EXPECT_DOUBLE_EQ(p.X(3), 99.2);

  p.SetInitialDisp(7.0);
  p.SetInitialVel(8.0);
  EXPECT_DOUBLE_EQ(p.X(0), 7.0);
  EXPECT_DOUBLE_EQ(p.X(1), 7.0);
  EXPECT_DOUBLE_EQ(p.X(2), 8.0);
  EXPECT_DOUBLE_EQ(p.X(3), 8.0);
}

TEST(ProblemTest, FixMassTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  
  Maybe<Void> e = p.FixMass(0);
  EXPECT_FALSE(e.isError);

  e = p.FixMass(1);
  EXPECT_FALSE(e.isError);

  EXPECT_EQ(p.fixedMasses.size(),2);
  EXPECT_TRUE(p.fixedMasses.find(p.GetMass(0).val) != p.fixedMasses.end());
  EXPECT_TRUE(p.fixedMasses.find(p.GetMass(1).val) != p.fixedMasses.end());

  e = p.FixMass(2);
  EXPECT_TRUE(e.isError);

  e = p.SetInitialVel(0, 10.0);
  EXPECT_TRUE(e.isError);
}

TEST(ProblemTest, GetDispAndVelTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,9.0);
  p.Build();

  p.SetInitialDisp(7.0);
  p.SetInitialVel(8.0);

  auto s = p.getDisp();
  EXPECT_EQ(s.size1(), 2);
  EXPECT_EQ(s.size2(), 1);
  EXPECT_DOUBLE_EQ(s(0,0), 7.0);
  EXPECT_DOUBLE_EQ(s(1,0), 7.0);

  auto v = p.getVel();
  EXPECT_EQ(v.size1(), 2);
  EXPECT_EQ(v.size2(), 1);
  EXPECT_DOUBLE_EQ(v(0,0), 8.0);
  EXPECT_DOUBLE_EQ(v(1,0), 8.0);
}

TEST(ProblemTest, XDotSimpleTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,1.0);
  p.Build();

  // Zero initial displacements and velocities:
  // Initial XDot is zero.
  auto XDot = p.XDot();
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0);
  EXPECT_DOUBLE_EQ(XDot(1), 0.0);
  EXPECT_DOUBLE_EQ(XDot(2), 0.0);
  EXPECT_DOUBLE_EQ(XDot(3), 0.0);
}

TEST(ProblemTest, XDotInitialVelocityTest) {
  Problem p = Problem();
  auto x0Dot = 13.0;
  auto x1Dot = 15.0;
  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,1.0);
  p.Build();

  p.SetInitialVel(0,x0Dot);
  p.SetInitialVel(1,x1Dot);

  // Zero initial displacements and non-zero Initial velocities:
  //    Initial accelerations are zero.
  auto XDot = p.XDot();
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), 0.0);
  EXPECT_DOUBLE_EQ(XDot(3), 0.0);
}

TEST(ProblemTest, XDotInitialDisplacementTest) {
  Problem p = Problem();
  auto m0 = 1.0;
  auto m1 = 2.0;
  auto k = 5.0;
  auto x0 = 9.0;
  auto x1 = 11.0;
  p.AddMass(m0,0.0,0.0);
  p.AddMass(m1,1.0,1.0);
  p.AddSpring(0,1,k);
  p.Build();

  p.SetInitialDisp(0,x0);
  p.SetInitialDisp(1,x1);

  // |x0DotDot| = |1/m0 0   | * |-k k| * |x0|
  // |x1DotDot|   |0    1/m1|   |k -k|   |x1|

  // Non-zero initial displacements and zero Initial velocities:
  //    Initial accelerations are non-zero.
  auto XDot = p.XDot();
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0);
  EXPECT_DOUBLE_EQ(XDot(1), 0.0);
  EXPECT_DOUBLE_EQ(XDot(2), 1/m0*(-k*x0+k*x1));
  EXPECT_DOUBLE_EQ(XDot(3), 1/m1*(k*x0-k*x1));
}

TEST(ProblemTest, XDotInitialDisplacementAndVelocityTest) {
  // XDotInitialVelocityTest and XDotInitialDisplacementTest combined
  Problem p = Problem();
  auto m0 = 1.0;
  auto m1 = 2.0;
  auto k = 5.0;
  auto x0 = 9.0;
  auto x1 = 11.0;
  auto x0Dot = 13.0;
  auto x1Dot = 15.0;
  p.AddMass(m0,0.0,0.0);
  p.AddMass(m1,1.0,1.0);
  p.AddSpring(0,1,k);
  p.Build();

  p.SetInitialDisp(0,x0);
  p.SetInitialDisp(1,x1);
  p.SetInitialVel(0,x0Dot);
  p.SetInitialVel(1,x1Dot);

  auto XDot = p.XDot();
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), 1/m0*(-k*x0+k*x1));
  EXPECT_DOUBLE_EQ(XDot(3), 1/m1*(k*x0-k*x1));
}

TEST(ProblemTest, XDotInitialDisplacementAndVelocityWithFixedMassTest) {
  Problem p = Problem();
  auto m0 = 1.0;
  auto m1 = 2.0;
  auto k = 5.0;
  auto x1 = 11.0;
  auto x1Dot = 15.0;
  p.AddMass(m0,0.0,0.0);
  p.AddMass(m1,1.0,1.0);
  p.AddSpring(0,1,k);
  p.Build();

  p.FixMass(0);
  p.SetInitialVel(0,99.0); // FixMas should take precedence over SetInitialVel

  p.SetInitialDisp(1,x1);
  p.SetInitialVel(1,x1Dot);

  auto XDot = p.XDot();
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), 0.0);
  EXPECT_DOUBLE_EQ(XDot(3), 1/m1*(k*0-k*x1));
}