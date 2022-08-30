#include <gtest/gtest.h>
#include "problem.h"
#include <iostream>
#include <math.h>

TEST(ProblemTest, MassCreationTest) {
  Problem* p = new Problem();

  auto e = p->AddMass(1.0,0.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 1);
  EXPECT_TRUE(e.val == 0);

  e = p->AddMass(1.0,1.0,0.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 2);
  EXPECT_TRUE(e.val == 1);

  e = p->AddMass(1.0,1.0,1.0);
  EXPECT_FALSE(e.isError);
  EXPECT_TRUE(p->GetDof() == 3);
  EXPECT_TRUE(e.val == 2);

  e = p->AddMass(1.0,0.0,0.0);
  EXPECT_TRUE(e.isError);
  EXPECT_TRUE(p->GetDof() == 3);

  EXPECT_TRUE(p->AddMass(0.0,9.0,9.0).isError); // 0.0 mass
}

TEST(ProblemTest, GetMassDispAndVelTest) {
  auto p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(1.0,1.0,1.0);
  p.AddMass(2.0,2.0,2.0);

  EXPECT_EQ(p.GetMassDispIndex(0),0);
  EXPECT_EQ(p.GetMassDispIndex(1),1);
  EXPECT_EQ(p.GetMassDispIndex(2),2);
  EXPECT_EQ(p.GetMassVelIndex(0),3);
  EXPECT_EQ(p.GetMassVelIndex(1),4);
  EXPECT_EQ(p.GetMassVelIndex(2),5);
}

TEST(ProblemTest, GetMassTest) {
  Problem* p = new Problem();

  p->AddMass(1.0,0.0,0.1);
  p->AddMass(1.0,1.0,1.1);
  p->AddMass(1.0,2.0,2.1);

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
  p.AddMass(1.0,1.0,1.0);
  p.AddMass(1.0,2.0,2.0);

  p.AddSpring(0,1,1.0);
  p.AddSpring(1,2,1.0);
  p.AddDamper(0,1,1.0);
  p.AddDamper(1,2,1.0);


  p.Build();
  p.Build(); // Calling again to make sure matrices are reset

  // |x0DotDot| = |1/m0 0       0| * (|-k1  k1       0| * |x0| + |-c1  c1        0| * |xDot0|)
  // |x1DotDot|   |0    1/m1    0|   (|k1   -k1-k2  k2|   |x1|   |c1   -c1-c2   c2| * |xDot1|)
  // |x2DotDot|   |0    0    1/m2|   (|0    k2     -k2|   |x2|   |0    c2      -c2| * |xDot1|)

  auto MInv = p.MInv;
  EXPECT_EQ(MInv.size1(), 3);
  EXPECT_EQ(MInv.size2(), 3);
  EXPECT_DOUBLE_EQ(MInv(0,0), 1.0);
  EXPECT_DOUBLE_EQ(MInv(0,1), 0.0);
  EXPECT_DOUBLE_EQ(MInv(0,2), 0.0);
  EXPECT_DOUBLE_EQ(MInv(1,0), 0.0);
  EXPECT_DOUBLE_EQ(MInv(1,1), 1.0);
  EXPECT_DOUBLE_EQ(MInv(1,2), 0.0);
  EXPECT_DOUBLE_EQ(MInv(2,0), 0.0);
  EXPECT_DOUBLE_EQ(MInv(2,1), 0.0);
  EXPECT_DOUBLE_EQ(MInv(2,2), 1.0);

  auto K = p.K;
  EXPECT_EQ(K.size1(), 3);
  EXPECT_EQ(K.size2(), 3);
  EXPECT_DOUBLE_EQ(K(0,0), -1.0);
  EXPECT_DOUBLE_EQ(K(0,1), 1.0);
  EXPECT_DOUBLE_EQ(K(0,2), 0.0);
  EXPECT_DOUBLE_EQ(K(1,0), 1.0);
  EXPECT_DOUBLE_EQ(K(1,1), -2.0);
  EXPECT_DOUBLE_EQ(K(1,2), 1.0);
  EXPECT_DOUBLE_EQ(K(2,0), 0.0);
  EXPECT_DOUBLE_EQ(K(2,1), 1.0);
  EXPECT_DOUBLE_EQ(K(2,2), -1.0);

  auto C = p.C;
  EXPECT_EQ(C.size1(), 3);
  EXPECT_EQ(C.size2(), 3);
  EXPECT_DOUBLE_EQ(C(0,0), -1.0);
  EXPECT_DOUBLE_EQ(C(0,1), 1.0);
  EXPECT_DOUBLE_EQ(C(0,2), 0.0);
  EXPECT_DOUBLE_EQ(C(1,0), 1.0);
  EXPECT_DOUBLE_EQ(C(1,1), -2.0);
  EXPECT_DOUBLE_EQ(C(1,2), 1.0);
  EXPECT_DOUBLE_EQ(C(2,0), 0.0);
  EXPECT_DOUBLE_EQ(C(2,1), 1.0);
  EXPECT_DOUBLE_EQ(C(2,2), -1.0);

  EXPECT_EQ(p.X.size(), 6);
  EXPECT_DOUBLE_EQ(p.X(0), 0.0);
  EXPECT_DOUBLE_EQ(p.X(1), 0.0);
  EXPECT_DOUBLE_EQ(p.X(2), 0.0);
  EXPECT_DOUBLE_EQ(p.X(3), 0.0);
  EXPECT_DOUBLE_EQ(p.X(4), 0.0);
  EXPECT_DOUBLE_EQ(p.X(5), 0.0);

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
  EXPECT_TRUE(e.isError); // Must first build

  p.Build();
  p.SetInitialVel(10.0);
  e = p.FixMass(0);
  EXPECT_FALSE(e.isError);

  EXPECT_DOUBLE_EQ(p.X[0], 0.0);
  EXPECT_DOUBLE_EQ(p.X[1], 0.0);
  EXPECT_DOUBLE_EQ(p.X[2], 0.0);
  EXPECT_DOUBLE_EQ(p.X[3], 10.0);

  EXPECT_EQ(p.fixedMasses.size(),1);
  EXPECT_TRUE(p.fixedMasses.find(p.GetMass(0).val) != p.fixedMasses.end());
  EXPECT_TRUE(p.fixedMasses.find(p.GetMass(1).val) == p.fixedMasses.end());

  e = p.FixMass(2);
  EXPECT_TRUE(e.isError);

  e = p.SetInitialVel(0, 10.0);
  EXPECT_TRUE(e.isError);
}

TEST(ProblemTest, GetDispAndVelTest) {

  vector<double> X = vector<double>(6);
  // System with 3 masses: 7.0 displacements and 8.0 velocities 
  X[0] = 7.0;
  X[1] = 7.0;
  X[2] = 7.0;
  X[3] = 8.0;
  X[4] = 8.0;
  X[5] = 8.0;

  auto s = Problem::getDisp(X, 3);
  EXPECT_EQ(s.size1(), 3);
  EXPECT_EQ(s.size2(), 1);
  EXPECT_DOUBLE_EQ(s(0,0), 7.0);
  EXPECT_DOUBLE_EQ(s(1,0), 7.0);
  EXPECT_DOUBLE_EQ(s(2,0), 7.0);

  auto v = Problem::getVel(X, 3);
  EXPECT_EQ(v.size1(), 3);
  EXPECT_EQ(v.size2(), 1);
  EXPECT_DOUBLE_EQ(v(0,0), 8.0);
  EXPECT_DOUBLE_EQ(v(1,0), 8.0);
  EXPECT_DOUBLE_EQ(v(2,0), 8.0);
}

TEST(ProblemTest, XDotSimpleTest) {
  Problem p = Problem();

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,1.0);
  p.Build();

  vector<double> XDot = vector<double>(p.X.size());
  p.SetXDot(p.X, XDot, 0.0);

  // Zero initial displacements and velocities:
  // Initial XDot is zero.
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

  vector<double> XDot = vector<double>(p.X.size());
  p.SetXDot(p.X, XDot, 0.0);
  
  // Zero initial displacements and non-zero Initial velocities:
  //    Initial accelerations are zero.
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
  auto c = 7.0;
  auto x0 = 9.0;
  auto x1 = 11.0;
  p.AddMass(m0,0.0,0.0);
  p.AddMass(m1,1.0,1.0);
  p.AddSpring(0,1,k);
  p.AddDamper(0,1,c);
  p.Build();

  p.SetInitialDisp(0,x0);
  p.SetInitialDisp(1,x1);

  vector<double> XDot = vector<double>(p.X.size());
  p.SetXDot(p.X, XDot, 0.0);

  // |x0DotDot| = |1/m0 0   | * (|-k k| * |x0| + |-c c| * |xDot0|)
  // |x1DotDot|   |0    1/m1|   (|k -k|   |x1|   |c -c| * |xDot1|)

  // Non-zero initial displacements and zero Initial velocities:
  //    Initial accelerations are non-zero.
  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0); // x0Dot
  EXPECT_DOUBLE_EQ(XDot(1), 0.0); // x1Dot
  EXPECT_DOUBLE_EQ(XDot(2), 1/m0*(-k*x0+k*x1)); // x0DotDot
  EXPECT_DOUBLE_EQ(XDot(3), 1/m1*(k*x0-k*x1));  // x1DotDot
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

  vector<double> XDot = vector<double>(p.X.size());
  p.SetXDot(p.X, XDot, 0.0);

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

  p.SetInitialVel(0,99.0);
  p.FixMass(0);

  p.SetInitialDisp(1,x1);
  p.SetInitialVel(1,x1Dot);

  vector<double> XDot = vector<double>(p.X.size());
  p.SetXDot(p.X, XDot, 0.0);

  EXPECT_EQ(XDot.size(), 4);
  EXPECT_DOUBLE_EQ(XDot(0), 0.0);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), 0.0);
  EXPECT_DOUBLE_EQ(XDot(3), 1/m1*(k*0-k*x1));
}

TEST(ProblemTest, GetXDotTest) {
  Problem p = Problem();
  auto m0 = 1.0;
  auto m1 = 2.0;
  auto k = 5.0;
  auto c = 7.0;
  p.AddMass(m0,0.0,-8.0); // random y values to test that they have no effect
  p.AddMass(m1,1.0,99.0);
  p.AddSpring(0,1,k);
  p.AddDamper(0,1,c);
  p.Build();

  ASSERT_EQ(p.X.size(), 4);
  vector<double> X = vector<double>(p.X.size());
  vector<double> XDot;

  // ZERO DISPLACEMENTS AND VELOCITIES TEST
  X[0] = 0.0;
  X[1] = 0.0;
  X[2] = 0.0;
  X[3] = 0.0;
  XDot = p.GetXDot(X, 0.0);
  EXPECT_EQ(XDot.size(), 4);

  // Expected values
  // |x0DotDot| = |1/m0 0   | * (|-k k| * |x0| + |-c c| * |xDot0|)
  // |x1DotDot|   |0    1/m1|   (|k -k|   |x1|   |c -c| * |xDot1|)
  double x0Dot = X[2];
  double x1Dot = X[3];
  double x0DotDot = 1/m0*(-k*X[0]+k*X[1]-c*X[2]+c*X[3]);
  double x1DotDot = 1/m1*(k*X[0]-k*X[1]+c*X[2]-c*X[3]);

  EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), x0DotDot);
  EXPECT_DOUBLE_EQ(XDot(3), x1DotDot);

  // INITIAL DISPLACEMENTS AND ZERO VELOCITIES TEST
  X[0] = 5.0;
  X[1] = 17.0;
  X[2] = 0.0;
  X[3] = 0.0;
  XDot = p.GetXDot(X, 0.0);
  EXPECT_EQ(XDot.size(), 4);
  x0Dot = X[2];
  x1Dot = X[3];
  x0DotDot = 1/m0*(-k*X[0]+k*X[1]-c*X[2]+c*X[3]);
  x1DotDot = 1/m1*(k*X[0]-k*X[1]+c*X[2]-c*X[3]);
  EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), x0DotDot);
  EXPECT_DOUBLE_EQ(XDot(3), x1DotDot);

  // INITIAL DISPLACEMENTS AND VELOCITIES TEST
  X[0] = 3.0;
  X[1] = 30.0;
  X[2] = 33.0;
  X[3] = 11.0;
  XDot = p.GetXDot(X, 0.0);
  EXPECT_EQ(XDot.size(), 4);
  x0Dot = X[2];
  x1Dot = X[3];
  x0DotDot = 1/m0*(-k*X[0]+k*X[1]-c*X[2]+c*X[3]);
  x1DotDot = 1/m1*(k*X[0]-k*X[1]+c*X[2]-c*X[3]);
  EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
  EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
  EXPECT_DOUBLE_EQ(XDot(2), x0DotDot);
  EXPECT_DOUBLE_EQ(XDot(3), x1DotDot);
}

TEST(ProblemTest,IntegrateStationaryTest) {
  Problem p = Problem();
  p.AddMass(1.0,0.0,0.0);
  p.AddMass(2.0,1.0,1.0);
  p.AddSpring(0,1,1.0);
  p.AddDamper(0,1,1.0);
  p.Build();

  p.FixMass(0);
  p.FixMass(1);

  p.Integrate(0.0, 1.0, 0.02);

  ASSERT_EQ(p.t.size(), p.XHistory.size());

  for (int i = 0; i < int(p.XHistory.size()); i++){

    ASSERT_EQ(p.XHistory[i].size(), 4);

    ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassDispIndex(0)],0.0);
    ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassDispIndex(1)],0.0);
    ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassVelIndex(0)],0.0);
    ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassVelIndex(1)],0.0);
  }
}

// Helper function that gets half the period of oscillation of the system
// Expects a setup in the following form:
//  p.AddMass(1.0,0.0,0.0);
//  p.AddMass(<M>,1.0,0.0);
//  p.AddSpring(0,1,<K>);
//  p.Build();
//  p.FixMass(0);
//  p.SetInitialDisp(1, <Number greater than 0>); --> IMPORTANT
//  p.Integrate(<t0>, <t1>, <tStep>);
double HalfPeriodOfOscillation(Problem p){

  // Find time in which speed becomes zero
  int tNeg; // t in which velocity is negative
  int tPos; // t in which velocity is positive
  for (int i = 0; i < int(p.XHistory.size()) -1 ; i++){
    if (p.XHistory[i][p.GetMassVelIndex(1)]<=0.0 && p.XHistory[i+1][p.GetMassVelIndex(1)]>0.0){
      tNeg = i;
      tPos = i+1;
      break;
    }
  }
  // Interpolate to find time in which velocity is zero
  double x0 = p.t[tNeg];
  double y0 = p.XHistory[tNeg][p.GetMassVelIndex(1)];
  double x1 = p.t[tPos];
  double y1 = p.XHistory[tPos][p.GetMassVelIndex(1)];
  // (y1-y0)/(x1-x0) = (0-y0)/(x-x0)
  // x = x0 -y0*(x1-x0)/(y1-y0)
  double halfT = x0 - y0*(x1-x0)/(y1-y0);

  return halfT;
}

// Helper function that tests if mass is stationary
// Should be used after p.Integrate is called
bool IsStationary(Problem p, int massId){
  // Verify first mass is stationary
  for (int i = 0; i < int(p.XHistory.size()); i++){
    if (p.XHistory[i][p.GetMassDispIndex(0)] != 0.0){
      return false;
    }
    if (p.XHistory[i][p.GetMassVelIndex(0)] != 0.0){
      return false;
    }
  }
  return true;
}

TEST(ProblemTest,HarmonicMotionTest) {
  // This test simulates a mass and spring system and automatically
  // validates the answer by comparing the numerical and theoretical
  // period of oscillation
  //
  // Expected natural frequency:
  // 1/(2pi)*sqrt(k/M) -> 1/2pi
  // Period of oscillation: 2pi
  // Half period: pi
  //
  // The mass's speed starts off negative (it has a positive initial disp).
  // Half period is how long it takes the mass to change direction, i.e.
  // For the speed to become zero.

  Problem p = Problem();
  p.AddMass(1.0,0.0,0.0);
  p.AddMass(5.0,1.0,0.0);
  p.AddSpring(0,1,5.0);
  p.Build();
  p.FixMass(0);
  p.SetInitialDisp(1, 0.1);

  p.Integrate(0.0, 4.0, 0.02);

  // Mass 0 should be stationary, and initial speed of mass 1 should be negative
  ASSERT_TRUE(IsStationary(p, 0));
  ASSERT_TRUE(p.XHistory[1][p.GetMassVelIndex(1)]<0);

  double halfT = HalfPeriodOfOscillation(p);

  double err = std::abs((halfT-M_PI)/M_PI);
  ASSERT_TRUE(err <= 0.002); // error smaller than 0.2%
}

TEST(ProblemTest,DampedOscillatorTest) {
  // This test simulates a mass, damper and spring system and automatically
  // validates the answer by comparing the numerical and theoretical
  // period of oscillation
  //
  // Expected natural frequency:
  // 1/(2pi)*sqrt(k/M - c^2/4m)
  // Period of oscillation: 2pi/sqrt(k/M - c^2/4m)
  // Half period: pi/sqrt(k/M - c^2/4m)
  //
  // for c=2, M=1 and k=5 -> Half period: pi/sqrt(5/1 - 2^2/4) = pi/2
  //
  // The mass's speed starts off negative (it has a positive initial disp).
  // Half period is how long it takes the mass to change direction, i.e.
  // For the speed to become zero.

  Problem p = Problem();
  p.AddMass(1.0,0.0,0.0);

  p.AddMass(1.0,1.0,0.0);
  p.AddSpring(0,1,5.0);
  p.AddDamper(0,1,2.0);
  p.Build();
  p.FixMass(0);
  p.SetInitialDisp(1, 10.0);

  p.Integrate(0.0, 5.0, 0.02);

  // Mass 0 should be stationary, and initial speed of mass 1 should be negative
  ASSERT_TRUE(IsStationary(p, 0));
  ASSERT_TRUE(p.XHistory[1][p.GetMassVelIndex(1)]<0);

  double halfT = HalfPeriodOfOscillation(p);

  double expectedHalfT  = M_PI/2.0;
  double err = std::abs((halfT-expectedHalfT)/expectedHalfT);
  ASSERT_TRUE(err <= 0.007);// error smaller than 0.7%
}

TEST(ProblemTest,DampedOscillatorPlotTest) {
  // This test simulates a mass, damper and spring system and prints a csv
  // to be visually tested against a plot available under:
  // http://spiff.rit.edu/classes/phys312/workshops/w5b/damped_theory.html#:~:text=A%20lightly%20damped%20harmonic%20oscillator,the%20decay%20happens%20more%20quickly.

  Problem p = Problem();
  p.AddMass(1.0,0.0,0.0);

  p.AddMass(20.0,1.0,0.0);
  p.AddSpring(0,1,30.0);
  p.AddDamper(0,1,2.9);
  p.Build();
  p.FixMass(0);
  p.SetInitialDisp(1, 1.0);

  p.Integrate(0.0, 40, 0.1);
  
  std::cout<<"DampedOscillatorPlotTest output:\n";
  p.PrintMassTimeHistory(1);
}

TEST(ProblemTest, MultiBodyBibliographyDataTest) {
  // This test simulates a system with 2 masses, similar to one found in the
  // bibliography and plots it's response so that we can compare it with the
  // results in the bibliography.
  // Source:
  // Mostafa, Marzbanrad Javad And. 2011. “A System Identification Algorithm for Vehicle Lumped.” International Journal of Modeling and Optimization 1 (January): 163–66.
  // Figure 5

  double m1 = 800;
  double m2 = 80;
  double c1 = 10000;
  double c2 = 1100;
  double k1 = 1000;
  double k2 = 160;
  double k3 = 2700;

  Problem p = Problem();
  EXPECT_FALSE(p.AddMass(1.0,0.0,0.0).isError);
  EXPECT_FALSE(p.AddMass(m1,1.0,0.0).isError); // m1
  EXPECT_FALSE(p.AddMass(m2,2.0,0.0).isError); // m2

  EXPECT_FALSE(p.AddSpring(0,1,k1).isError);
  EXPECT_FALSE(p.AddDamper(0,1,c1).isError);

  EXPECT_FALSE(p.AddSpring(1,2,k2).isError);
  EXPECT_FALSE(p.AddDamper(1,2,c2).isError);
  EXPECT_FALSE(p.AddSpring(1,2,k3).isError);

  p.Build();
  p.SetInitialVel(14.0);
  EXPECT_FALSE(p.FixMass(0).isError);

  EXPECT_EQ(p.X.size(), 6);
  EXPECT_DOUBLE_EQ(p.X[0], 0.0);
  EXPECT_DOUBLE_EQ(p.X[1], 0.0);
  EXPECT_DOUBLE_EQ(p.X[2], 0.0);
  EXPECT_DOUBLE_EQ(p.X[3], 0.0);
  EXPECT_DOUBLE_EQ(p.X[4], 14.0);
  EXPECT_DOUBLE_EQ(p.X[5], 14.0);

  auto XDot0 = p.GetXDot(p.X, 0.0);
  // Expected values
  // |x0DotDot| = |1/m0 0       0| * (|-k1  k1             0| * |x0| + |-c1  c1        0| * |xDot0|)
  // |x1DotDot|   |0    1/m1    0|   (|k1   -k1-k2-k3  k2+k3|   |x1|   |c1   -c1-c2   c2| * |xDot1|)
  // |x2DotDot|   |0    0    1/m2|   (|0    k2+k3     -k2-k3|   |x2|   |0    c2      -c2| * |xDot1|)
  EXPECT_EQ(XDot0.size(), 6);
  EXPECT_DOUBLE_EQ(XDot0[0], 0.0);  // xDot0
  EXPECT_DOUBLE_EQ(XDot0[1], 14.0); // xDot1
  EXPECT_DOUBLE_EQ(XDot0[2], 14.0); // xDot2
  EXPECT_DOUBLE_EQ(XDot0[3], 0.0);  // xDotDot0
  EXPECT_DOUBLE_EQ(XDot0[4], 1/m1*((-c1-c2)*14 + c2*14)); // xDotDot1
  EXPECT_DOUBLE_EQ(XDot0[5], 1/m2*(c2*14 -c2*14)); // xDotDot1

  EXPECT_FALSE(p.Integrate(0.0, 1.0, 0.001).isError);
  
  auto e = p.GetMassMinAccel(2);
  ASSERT_FALSE(e.isError);
  ASSERT_TRUE(-80<e.val && e.val<-60);
 
  e = p.GetMassMaxAccel(2);
  ASSERT_FALSE(e.isError);
  ASSERT_TRUE(0<e.val && e.val<10);

  ASSERT_DOUBLE_EQ(p.GetMassMaxAbsAccel(2).val, abs(p.GetMassMinAccel(2).val));
}

TEST(ProblemTest, MultiBodyBibliographyDataTest2) {
  // This test simulates a system with 5 masses, similar to one found in the
  // bibliography and plots it's response so that we can compare it with the
  // results in the bibliography.
  // Source:
  // Mostafa, Marzbanrad Javad And. 2011. “A System Identification Algorithm for Vehicle Lumped.” International Journal of Modeling and Optimization 1 (January): 163–66.
  // Figure 11


  Problem p = Problem();
  EXPECT_FALSE(p.AddMass(1.0,0.0,0.0).isError); // m0
  EXPECT_FALSE(p.AddMass(300,1.0,1.0).isError); // m1
  EXPECT_FALSE(p.AddMass(120,1.0,0.0).isError); // m2
  EXPECT_FALSE(p.AddMass(150,1.0,3.0).isError); // m3
  EXPECT_FALSE(p.AddMass(700,2.0,0.0).isError); // m4
  EXPECT_FALSE(p.AddMass(80 ,3.0,0.0).isError); // m5

  double k1 = 48.866;
  double k2 = 915522.58;
  double k3 = 1206875.43;
  double k4 = 1178694.84;
  double k5 = 36.7265;
  double k6 = 136.4661;
  double k7 = 38.6678;
  double k8 = 4761249.39;
  double k9 = 389232.42;
  double c4 = 0.8981;
  double c5 = 33114.21;
  double c6 = 1.7284;
  double c7 = 6764.6574;
  double c8 = 8648277.61;
  double c9 = 1595.52;

  EXPECT_FALSE(p.AddSpring(0,1,k1).isError);
  EXPECT_FALSE(p.AddSpring(1,2,k2).isError);
  EXPECT_FALSE(p.AddSpring(1,3,k3).isError);
  EXPECT_FALSE(p.AddSpring(1,4,k4).isError);
  EXPECT_FALSE(p.AddDamper(1,4,c4).isError);
  EXPECT_FALSE(p.AddSpring(0,2,k5).isError);
  EXPECT_FALSE(p.AddDamper(0,2,c5).isError);
  EXPECT_FALSE(p.AddSpring(2,4,k6).isError);
  EXPECT_FALSE(p.AddDamper(2,4,c6).isError);
  EXPECT_FALSE(p.AddSpring(0,3,k7).isError);
  EXPECT_FALSE(p.AddDamper(0,3,c7).isError);
  EXPECT_FALSE(p.AddSpring(3,4,k8).isError);
  EXPECT_FALSE(p.AddDamper(3,4,c8).isError);
  EXPECT_FALSE(p.AddSpring(4,5,k9).isError);
  EXPECT_FALSE(p.AddDamper(4,5,c9).isError);


  p.Build();
  p.SetInitialVel(14.0);
  EXPECT_FALSE(p.FixMass(0).isError);

  EXPECT_FALSE(p.Integrate(0.0, 0.15, 0.05).isError);

  // The velocities we get match with the paper. The acceleration follows the
  // same shape, but the paper doesn't show the units for the acceleration.
  // If our data is correct, it seems to be (10m/s^2). To validate the
  // accelerations we get is also correct, we do a single check using linear
  // interpolation at a random time instant.
  int i = p.XHistory.size()/2;
  double t0 = p.t[i];
  double v0 = p.XHistory[i][p.GetMassVelIndex(5)];
  double vDot0 = p.AccelHistory[i][5];
  double t1 = p.t[i+1];
  double v1 = p.XHistory[i+1][p.GetMassVelIndex(5)];
  double expectedV1 = v0 + vDot0*(t1-t0);
  double err = abs((expectedV1-v1)/v1);
  EXPECT_TRUE(err< 0.0001);

  p.PrintMassTimeHistory(5);

  auto e = p.GetMassMinAccel(5);
  ASSERT_FALSE(e.isError);
  ASSERT_TRUE(-400<e.val && e.val<-350);
 
  e = p.GetMassMaxAccel(5);
  ASSERT_FALSE(e.isError);
  ASSERT_TRUE(0<e.val && e.val<100);

  // p.PrintMassTimeHistory(5);
}

TEST(ProblemTest, MinMaxAccelTest){
  Problem p = Problem();

  auto e = p.GetMassMaxAbsAccel(0);
  ASSERT_TRUE(e.isError);
  e = p.GetMassMaxAccel(0);
  ASSERT_TRUE(e.isError);
  e = p.GetMassMinAccel(0);
  ASSERT_TRUE(e.isError);

  p.AddMass(1.0,0.0,0.0);
  p.AddMass(1.0,1.0,0.0);
  p.AddMass(1.0,2.0,0.0);
  
  // Called to prevent "problem must be integrated" error.
  // AccelHistory is then cleared and mocked for the test
  p.Build();
  p.Integrate(0.0,0.0,0.0);
  p.AccelHistory.clear();
  ASSERT_EQ(p.AccelHistory.size(), 0);

  vector<double> a = vector<double>(3);
  a[0] = -1.0;
  a[1] = 0.0;
  a[2] = -8.0;
  p.AccelHistory.push_back(a);
  a[0] = 0.0;
  a[1] = 1.0;
  a[2] = -90.0;
  p.AccelHistory.push_back(a);
  a[0] = 1.0;
  a[1] = 2.0;
  a[2] = -9.0;
  p.AccelHistory.push_back(a);

  ASSERT_DOUBLE_EQ(p.GetMassMaxAccel(0).val, 1.0);
  ASSERT_DOUBLE_EQ(p.GetMassMinAccel(0).val, -1.0);
  ASSERT_DOUBLE_EQ(p.GetMassMaxAbsAccel(0).val, 1.0);

  ASSERT_DOUBLE_EQ(p.GetMassMaxAccel(1).val, 2.0);
  ASSERT_DOUBLE_EQ(p.GetMassMinAccel(1).val, 0.0);
  ASSERT_DOUBLE_EQ(p.GetMassMaxAbsAccel(1).val, 2.0);

  ASSERT_DOUBLE_EQ(p.GetMassMaxAccel(2).val, -8.0);
  ASSERT_DOUBLE_EQ(p.GetMassMinAccel(2).val, -90.0);
  ASSERT_DOUBLE_EQ(p.GetMassMaxAbsAccel(2).val, 90.0);
}