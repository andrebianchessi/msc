#include "problem.h"

#include <gtest/gtest.h>
#include <math.h>

#include <iostream>

TEST(ProblemTest, MassCreationTest) {
    Problem p = Problem();

    auto e = p.AddMass(1.0, 0.0, 0.0);
    EXPECT_FALSE(e.isError);
    EXPECT_TRUE(p.GetDof() == 1);
    EXPECT_TRUE(e.val == 0);

    e = p.AddMass(1.0, 1.0, 0.0);
    EXPECT_FALSE(e.isError);
    EXPECT_TRUE(p.GetDof() == 2);
    EXPECT_TRUE(e.val == 1);

    e = p.AddMass(1.0, 1.0, 1.0);
    EXPECT_FALSE(e.isError);
    EXPECT_TRUE(p.GetDof() == 3);
    EXPECT_TRUE(e.val == 2);

    e = p.AddMass(1.0, 0.0, 0.0);
    EXPECT_TRUE(e.isError);
    EXPECT_TRUE(p.GetDof() == 3);

    EXPECT_TRUE(p.AddMass(0.0, 9.0, 9.0).isError);  // 0.0 mass
}

TEST(ProblemTest, GetMassDispAndVelTest) {
    auto p = Problem();

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(1.0, 1.0, 1.0);
    p.AddMass(2.0, 2.0, 2.0);

    EXPECT_EQ(p.GetMassDispIndex(0), 0);
    EXPECT_EQ(p.GetMassDispIndex(1), 1);
    EXPECT_EQ(p.GetMassDispIndex(2), 2);
    EXPECT_EQ(p.GetMassVelIndex(0), 3);
    EXPECT_EQ(p.GetMassVelIndex(1), 4);
    EXPECT_EQ(p.GetMassVelIndex(2), 5);
}

TEST(ProblemTest, GetMassTest) {
    Problem p = Problem();

    p.AddMass(1.0, 0.0, 0.1);
    p.AddMass(1.0, 1.0, 1.1);
    p.AddMass(1.0, 2.0, 2.1);

    auto e = p.GetMass(0);
    EXPECT_FALSE(e.isError);
    auto m = e.val;
    EXPECT_TRUE(m->xIndex == 0);
    EXPECT_TRUE(m->px == 0.0);
    EXPECT_TRUE(m->py == 0.1);

    e = p.GetMass(2);
    EXPECT_FALSE(e.isError);
    m = e.val;
    EXPECT_TRUE(m->xIndex == 2);
    EXPECT_TRUE(m->px == 2.0);
    EXPECT_TRUE(m->py == 2.1);

    e = p.GetMass(3);
    EXPECT_TRUE(e.isError);
}

TEST(ProblemTest, SimpleBuildTest) {
    Problem p = Problem();

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(1.0, 1.0, 1.0);
    p.AddMass(1.0, 2.0, 2.0);

    p.AddSpring(0, 1, 1.0);
    p.AddSpring(1, 2, 1.0);
    p.AddDamper(0, 1, 1.0);
    p.AddDamper(1, 2, 1.0);

    p.Build();
    p.Build();  // Calling again to make sure matrices are reset

    // |x0DotDot| = |1/m0 0       0| * (|-k1  k1       0| * |x0| + |-c1  c1 0| *
    // |xDot0|) |x1DotDot|   |0    1/m1    0|   (|k1   -k1-k2  k2|   |x1|   |c1
    // -c1-c2   c2| * |xDot1|) |x2DotDot|   |0    0    1/m2|   (|0    k2 -k2|
    // |x2|   |0    c2      -c2| * |xDot1|)

    auto MInv = p.MInv;
    EXPECT_EQ(MInv.size1(), 3);
    EXPECT_EQ(MInv.size2(), 3);
    EXPECT_DOUBLE_EQ(MInv(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(MInv(0, 1), 0.0);
    EXPECT_DOUBLE_EQ(MInv(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(MInv(1, 0), 0.0);
    EXPECT_DOUBLE_EQ(MInv(1, 1), 1.0);
    EXPECT_DOUBLE_EQ(MInv(1, 2), 0.0);
    EXPECT_DOUBLE_EQ(MInv(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(MInv(2, 1), 0.0);
    EXPECT_DOUBLE_EQ(MInv(2, 2), 1.0);

    auto K = p.K;
    EXPECT_EQ(K.size1(), 3);
    EXPECT_EQ(K.size2(), 3);
    EXPECT_DOUBLE_EQ(K(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(K(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(K(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(K(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(K(1, 1), -2.0);
    EXPECT_DOUBLE_EQ(K(1, 2), 1.0);
    EXPECT_DOUBLE_EQ(K(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(K(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(K(2, 2), -1.0);

    auto C = p.C;
    EXPECT_EQ(C.size1(), 3);
    EXPECT_EQ(C.size2(), 3);
    EXPECT_DOUBLE_EQ(C(0, 0), -1.0);
    EXPECT_DOUBLE_EQ(C(0, 1), 1.0);
    EXPECT_DOUBLE_EQ(C(0, 2), 0.0);
    EXPECT_DOUBLE_EQ(C(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(C(1, 1), -2.0);
    EXPECT_DOUBLE_EQ(C(1, 2), 1.0);
    EXPECT_DOUBLE_EQ(C(2, 0), 0.0);
    EXPECT_DOUBLE_EQ(C(2, 1), 1.0);
    EXPECT_DOUBLE_EQ(C(2, 2), -1.0);

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

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);
    p.AddSpring(0, 1, 9.0);
    p.Build();

    p.SetInitialDisp(0, 9.1);
    p.SetInitialDisp(1, 9.2);
    p.SetInitialVel(0, 99.1);
    p.SetInitialVel(1, 99.2);
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

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);

    Maybe<Void> e = p.FixMass(0);
    EXPECT_TRUE(e.isError);  // Must first build

    p.Build();
    p.SetInitialVel(10.0);
    e = p.FixMass(0);
    EXPECT_FALSE(e.isError);

    EXPECT_DOUBLE_EQ(p.X[0], 0.0);
    EXPECT_DOUBLE_EQ(p.X[1], 0.0);
    EXPECT_DOUBLE_EQ(p.X[2], 0.0);
    EXPECT_DOUBLE_EQ(p.X[3], 10.0);

    EXPECT_EQ(p.fixedMasses.size(), 1);
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
    EXPECT_DOUBLE_EQ(s(0, 0), 7.0);
    EXPECT_DOUBLE_EQ(s(1, 0), 7.0);
    EXPECT_DOUBLE_EQ(s(2, 0), 7.0);

    auto v = Problem::getVel(X, 3);
    EXPECT_EQ(v.size1(), 3);
    EXPECT_EQ(v.size2(), 1);
    EXPECT_DOUBLE_EQ(v(0, 0), 8.0);
    EXPECT_DOUBLE_EQ(v(1, 0), 8.0);
    EXPECT_DOUBLE_EQ(v(2, 0), 8.0);
}

TEST(ProblemTest, XDotSimpleTest) {
    Problem p = Problem();

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);
    p.AddSpring(0, 1, 1.0);
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
    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);
    p.AddSpring(0, 1, 1.0);
    p.Build();

    p.SetInitialVel(0, x0Dot);
    p.SetInitialVel(1, x1Dot);

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
    p.AddMass(m0, 0.0, 0.0);
    p.AddMass(m1, 1.0, 1.0);
    p.AddSpring(0, 1, k);
    p.AddDamper(0, 1, c);
    p.Build();

    p.SetInitialDisp(0, x0);
    p.SetInitialDisp(1, x1);

    vector<double> XDot = vector<double>(p.X.size());
    p.SetXDot(p.X, XDot, 0.0);

    // |x0DotDot| = |1/m0 0   | * (|-k k| * |x0| + |-c c| * |xDot0|)
    // |x1DotDot|   |0    1/m1|   (|k -k|   |x1|   |c -c| * |xDot1|)

    // Non-zero initial displacements and zero Initial velocities:
    //    Initial accelerations are non-zero.
    EXPECT_EQ(XDot.size(), 4);
    EXPECT_DOUBLE_EQ(XDot(0), 0.0);                          // x0Dot
    EXPECT_DOUBLE_EQ(XDot(1), 0.0);                          // x1Dot
    EXPECT_DOUBLE_EQ(XDot(2), 1 / m0 * (-k * x0 + k * x1));  // x0DotDot
    EXPECT_DOUBLE_EQ(XDot(3), 1 / m1 * (k * x0 - k * x1));   // x1DotDot
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
    p.AddMass(m0, 0.0, 0.0);
    p.AddMass(m1, 1.0, 1.0);
    p.AddSpring(0, 1, k);
    p.Build();

    p.SetInitialDisp(0, x0);
    p.SetInitialDisp(1, x1);
    p.SetInitialVel(0, x0Dot);
    p.SetInitialVel(1, x1Dot);

    vector<double> XDot = vector<double>(p.X.size());
    p.SetXDot(p.X, XDot, 0.0);

    EXPECT_EQ(XDot.size(), 4);
    EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
    EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
    EXPECT_DOUBLE_EQ(XDot(2), 1 / m0 * (-k * x0 + k * x1));
    EXPECT_DOUBLE_EQ(XDot(3), 1 / m1 * (k * x0 - k * x1));
}

TEST(ProblemTest, XDotInitialDisplacementAndVelocityWithFixedMassTest) {
    Problem p = Problem();
    auto m0 = 1.0;
    auto m1 = 2.0;
    auto k = 5.0;
    auto x1 = 11.0;
    auto x1Dot = 15.0;
    p.AddMass(m0, 0.0, 0.0);
    p.AddMass(m1, 1.0, 1.0);
    p.AddSpring(0, 1, k);
    p.Build();

    p.SetInitialVel(0, 99.0);
    p.FixMass(0);

    p.SetInitialDisp(1, x1);
    p.SetInitialVel(1, x1Dot);

    vector<double> XDot = vector<double>(p.X.size());
    p.SetXDot(p.X, XDot, 0.0);

    EXPECT_EQ(XDot.size(), 4);
    EXPECT_DOUBLE_EQ(XDot(0), 0.0);
    EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
    EXPECT_DOUBLE_EQ(XDot(2), 0.0);
    EXPECT_DOUBLE_EQ(XDot(3), 1 / m1 * (k * 0 - k * x1));
}

TEST(ProblemTest, GetXDotTest) {
    Problem p = Problem();
    auto m0 = 1.0;
    auto m1 = 2.0;
    auto k = 5.0;
    auto c = 7.0;
    p.AddMass(m0, 0.0,
              -8.0);  // random y values to test that they have no effect
    p.AddMass(m1, 1.0, 99.0);
    p.AddSpring(0, 1, k);
    p.AddDamper(0, 1, c);
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
    double x0DotDot = 1 / m0 * (-k * X[0] + k * X[1] - c * X[2] + c * X[3]);
    double x1DotDot = 1 / m1 * (k * X[0] - k * X[1] + c * X[2] - c * X[3]);

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
    x0DotDot = 1 / m0 * (-k * X[0] + k * X[1] - c * X[2] + c * X[3]);
    x1DotDot = 1 / m1 * (k * X[0] - k * X[1] + c * X[2] - c * X[3]);
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
    x0DotDot = 1 / m0 * (-k * X[0] + k * X[1] - c * X[2] + c * X[3]);
    x1DotDot = 1 / m1 * (k * X[0] - k * X[1] + c * X[2] - c * X[3]);
    EXPECT_DOUBLE_EQ(XDot(0), x0Dot);
    EXPECT_DOUBLE_EQ(XDot(1), x1Dot);
    EXPECT_DOUBLE_EQ(XDot(2), x0DotDot);
    EXPECT_DOUBLE_EQ(XDot(3), x1DotDot);
}

TEST(ProblemTest, GetAccelTest) {
    // Values taken from GetXDotTest
    Problem p = Problem();
    auto m0 = 1.0;
    auto m1 = 2.0;
    auto k = 5.0;
    auto c = 7.0;
    p.AddMass(m0, 0.0, -8.0);
    p.AddMass(m1, 1.0, 99.0);
    p.AddSpring(0, 1, k);
    p.AddDamper(0, 1, c);
    p.Build();

    ASSERT_EQ(p.X.size(), 4);
    vector<double> X = vector<double>(p.X.size());

    // INITIAL DISPLACEMENTS AND VELOCITIES TEST
    X[0] = 3.0;
    X[1] = 30.0;
    X[2] = 33.0;
    X[3] = 11.0;

    double x0DotDot = 1 / m0 * (-k * X[0] + k * X[1] - c * X[2] + c * X[3]);
    double x1DotDot = 1 / m1 * (k * X[0] - k * X[1] + c * X[2] - c * X[3]);

    auto accels = p.GetAccel(X, 0);
    ASSERT_FALSE(accels.isError);
    ASSERT_EQ(accels.val.size(), 2);
    ASSERT_DOUBLE_EQ(accels.val[0], x0DotDot);
    ASSERT_DOUBLE_EQ(accels.val[1], x1DotDot);

    // Error cases
    X = vector<double>(3);
    accels = p.GetAccel(X, 0);
    ASSERT_TRUE(accels.isError);
    X = vector<double>(5);
    accels = p.GetAccel(X, 0);
    ASSERT_TRUE(accels.isError);
}

TEST(ProblemTest, IntegrateStationaryTest) {
    Problem p = Problem();
    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);
    p.AddSpring(0, 1, 1.0);
    p.AddDamper(0, 1, 1.0);
    p.Build();

    p.FixMass(0);
    p.FixMass(1);

    p.Integrate(1.0);

    ASSERT_EQ(p.t.size(), p.XHistory.size());

    for (int i = 0; i < int(p.XHistory.size()); i++) {
        ASSERT_EQ(p.XHistory[i].size(), 4);

        ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassDispIndex(0)], 0.0);
        ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassDispIndex(1)], 0.0);
        ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassVelIndex(0)], 0.0);
        ASSERT_DOUBLE_EQ(p.XHistory[i][p.GetMassVelIndex(1)], 0.0);
    }
}

TEST(ProblemTest, IntegrateSignTest) {
    // m0 -- m1
    // m1 has initial positive displacement
    // Expected: m1 will have initial negative speed

    Problem p = Problem();
    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(2.0, 1.0, 1.0);
    p.AddSpring(0, 1, 1.0);
    p.Build();

    p.FixMass(0);
    p.SetInitialDisp(1, 10);

    p.Integrate(1.0);

    ASSERT_EQ(p.t.size(), p.XHistory.size());
    ASSERT_TRUE(p.XHistory[1][p.GetMassVelIndex(1)] < 0);
}

TEST(ProblemTest, AddSpringAndDamperOrderTest) {
    // Integrates problems that should have the same results. The only
    // difference is that the springs/damper/masses were added in different
    // order

    Problem p0 = Problem();
    int fixedMassId0 = p0.AddMass(1.0, 0.0, 0.0).val;        // 1
    int initialDispMassId0 = p0.AddMass(2.0, 1.0, 1.0).val;  // 2
    p0.AddSpring(0, 1, 1.0);                                 // 0,1
    p0.AddDamper(1, 0, 2.0);                                 // 1,0
    p0.Build();
    p0.FixMass(fixedMassId0);
    p0.SetInitialDisp(initialDispMassId0, 10);
    p0.Integrate(1.0);

    Problem p1 = Problem();
    int initialDispMassId1 = p1.AddMass(2.0, 1.0, 1.0).val;  // 2
    int fixedMassId1 = p1.AddMass(1.0, 0.0, 0.0).val;        // 1
    p1.AddSpring(0, 1, 1.0);                                 // 0,1
    p1.AddDamper(0, 1, 2.0);                                 // 0,1
    p1.Build();
    p1.FixMass(fixedMassId1);
    p1.SetInitialDisp(initialDispMassId1, 10);
    p1.Integrate(1.0);

    Problem p2 = Problem();
    int initialDispMassId2 = p2.AddMass(2.0, 1.0, 1.0).val;  // 2
    int fixedMassId2 = p2.AddMass(1.0, 0.0, 0.0).val;        // 1
    p2.AddSpring(1, 0, 1.0);                                 // 1,0
    p2.AddDamper(1, 0, 2.0);                                 // 1,0
    p2.Build();
    p2.FixMass(fixedMassId2);
    p2.SetInitialDisp(initialDispMassId2, 10);
    p2.Integrate(1.0);

    Problem p3 = Problem();
    int fixedMassId3 = p3.AddMass(1.0, 0.0, 0.0).val;        // 1
    int initialDispMassId3 = p3.AddMass(2.0, 1.0, 1.0).val;  // 2
    p3.AddSpring(1, 0, 1.0);                                 // 1,0
    p3.AddDamper(0, 1, 2.0);                                 // 0,1
    p3.Build();
    p3.FixMass(fixedMassId3);
    p3.SetInitialDisp(initialDispMassId3, 10);
    p3.Integrate(1.0);

    // Assert initial speed is negative
    ASSERT_TRUE(p0.XHistory[1][p0.GetMassVelIndex(1)] < 0);

    // Assert all have same history length
    ASSERT_TRUE(p0.XHistory.size() == p1.XHistory.size());
    ASSERT_TRUE(p0.XHistory.size() == p2.XHistory.size());
    ASSERT_TRUE(p0.XHistory.size() == p3.XHistory.size());

    for (int i = 0; i < int(p0.XHistory.size()) - 1; i++) {
        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassDispIndex(fixedMassId0)],
                         p1.XHistory[i][p1.GetMassDispIndex(fixedMassId1)]);
        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassVelIndex(fixedMassId0)],
                         p1.XHistory[i][p1.GetMassVelIndex(fixedMassId1)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassDispIndex(initialDispMassId0)],
            p1.XHistory[i][p1.GetMassDispIndex(initialDispMassId1)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassVelIndex(initialDispMassId0)],
            p1.XHistory[i][p1.GetMassVelIndex(initialDispMassId1)]);

        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassDispIndex(fixedMassId0)],
                         p2.XHistory[i][p2.GetMassDispIndex(fixedMassId2)]);
        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassVelIndex(fixedMassId0)],
                         p2.XHistory[i][p2.GetMassVelIndex(fixedMassId2)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassDispIndex(initialDispMassId0)],
            p2.XHistory[i][p2.GetMassDispIndex(initialDispMassId2)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassVelIndex(initialDispMassId0)],
            p2.XHistory[i][p2.GetMassVelIndex(initialDispMassId2)]);

        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassDispIndex(fixedMassId0)],
                         p3.XHistory[i][p3.GetMassDispIndex(fixedMassId3)]);
        ASSERT_DOUBLE_EQ(p0.XHistory[i][p0.GetMassVelIndex(fixedMassId0)],
                         p3.XHistory[i][p3.GetMassVelIndex(fixedMassId3)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassDispIndex(initialDispMassId0)],
            p3.XHistory[i][p3.GetMassDispIndex(initialDispMassId3)]);
        ASSERT_DOUBLE_EQ(
            p0.XHistory[i][p0.GetMassVelIndex(initialDispMassId0)],
            p3.XHistory[i][p3.GetMassVelIndex(initialDispMassId3)]);
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
double HalfPeriodOfOscillation(Problem p) {
    // Find time in which speed becomes zero
    int tNeg;  // t in which velocity is negative
    int tPos;  // t in which velocity is positive
    for (int i = 0; i < int(p.XHistory.size()) - 1; i++) {
        if (p.XHistory[i][p.GetMassVelIndex(1)] <= 0.0 &&
            p.XHistory[i + 1][p.GetMassVelIndex(1)] > 0.0) {
            tNeg = i;
            tPos = i + 1;
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
    double halfT = x0 - y0 * (x1 - x0) / (y1 - y0);

    return halfT;
}

// Helper function that tests if mass is stationary
// Should be used after p.Integrate is called
bool IsStationary(Problem p, int massId) {
    // Verify first mass is stationary
    for (int i = 0; i < int(p.XHistory.size()); i++) {
        if (p.XHistory[i][p.GetMassDispIndex(0)] != 0.0) {
            return false;
        }
        if (p.XHistory[i][p.GetMassVelIndex(0)] != 0.0) {
            return false;
        }
    }
    return true;
}

TEST(ProblemTest, HarmonicMotionTest) {
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
    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(5.0, 1.0, 0.0);
    p.AddSpring(0, 1, 5.0);
    p.Build();
    p.FixMass(0);
    p.SetInitialDisp(1, 0.1);

    p.Integrate(4.0);

    // Mass 0 should be stationary, and initial speed of mass 1 should be
    // negative
    ASSERT_TRUE(IsStationary(p, 0));
    ASSERT_TRUE(p.XHistory[1][p.GetMassVelIndex(1)] < 0);

    double halfT = HalfPeriodOfOscillation(p);

    double err = std::abs((halfT - M_PI) / M_PI);
    ASSERT_TRUE(err <= 0.002);  // error smaller than 0.2%
}

TEST(ProblemTest, DampedOscillatorTest) {
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
    p.AddMass(1.0, 0.0, 0.0);

    p.AddMass(1.0, 1.0, 0.0);
    p.AddSpring(0, 1, 5.0);
    p.AddDamper(0, 1, 2.0);
    p.Build();
    p.FixMass(0);
    p.SetInitialDisp(1, 10.0);

    p.Integrate(5.0);

    // Mass 0 should be stationary, and initial speed of mass 1 should be
    // negative
    ASSERT_TRUE(IsStationary(p, 0));
    ASSERT_TRUE(p.XHistory[1][p.GetMassVelIndex(1)] < 0);

    double halfT = HalfPeriodOfOscillation(p);

    double expectedHalfT = M_PI / 2.0;
    double err = std::abs((halfT - expectedHalfT) / expectedHalfT);
    ASSERT_TRUE(err <= 0.007);  // error smaller than 0.7%
}

TEST(ProblemTest, DampedOscillatorPlotTest) {
    // This test simulates a mass, damper and spring system and prints a csv
    // to be visually tested against a plot available under:
    // http://spiff.rit.edu/classes/phys312/workshops/w5b/damped_theory.html#:~:text=A%20lightly%20damped%20harmonic%20oscillator,the%20decay%20happens%20more%20quickly.

    Problem p = Problem();
    p.AddMass(1.0, 0.0, 0.0);

    p.AddMass(20.0, 1.0, 0.0);
    p.AddSpring(0, 1, 30.0);
    p.AddDamper(0, 1, 2.9);
    p.Build();
    p.FixMass(0);
    p.SetInitialDisp(1, 1.0);

    p.Integrate(40);

    std::cout << "DampedOscillatorPlotTest output:\n";
    p.PrintMassTimeHistory(1);
}

TEST(ProblemTest, MultiBodyBibliographyDataTest) {
    // This test simulates a system with 2 masses, similar to one found in the
    // bibliography and plots it's response so that we can compare it with the
    // results in the bibliography.
    // Source:
    // Mostafa, Marzbanrad Javad And. 2011. “A System Identification Algorithm
    // for Vehicle Lumped.” International Journal of Modeling and Optimization 1
    // (January): 163–66. Figure 5

    double m1 = 800;
    double m2 = 80;
    double c1 = 10000;
    double c2 = 1100;
    double k1 = 1000;
    double k2 = 160;
    double k3 = 2700;

    Problem p = Problem();
    EXPECT_FALSE(p.AddMass(1.0, 0.0, 0.0).isError);
    EXPECT_FALSE(p.AddMass(m1, 1.0, 0.0).isError);  // m1
    EXPECT_FALSE(p.AddMass(m2, 2.0, 0.0).isError);  // m2

    EXPECT_FALSE(p.AddSpring(0, 1, k1).isError);
    EXPECT_FALSE(p.AddDamper(0, 1, c1).isError);

    EXPECT_FALSE(p.AddSpring(1, 2, k2).isError);
    EXPECT_FALSE(p.AddDamper(1, 2, c2).isError);
    EXPECT_FALSE(p.AddSpring(1, 2, k3).isError);

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
    // |x0DotDot| = |1/m0 0       0| * (|-k1  k1             0| * |x0| + |-c1 c1
    // 0| * |xDot0|) |x1DotDot|   |0    1/m1    0|   (|k1   -k1-k2-k3  k2+k3|
    // |x1|   |c1   -c1-c2   c2| * |xDot1|) |x2DotDot|   |0    0    1/m2|   (|0
    // k2+k3     -k2-k3|   |x2|   |0    c2      -c2| * |xDot1|)
    EXPECT_EQ(XDot0.size(), 6);
    EXPECT_DOUBLE_EQ(XDot0[0], 0.0);   // xDot0
    EXPECT_DOUBLE_EQ(XDot0[1], 14.0);  // xDot1
    EXPECT_DOUBLE_EQ(XDot0[2], 14.0);  // xDot2
    EXPECT_DOUBLE_EQ(XDot0[3], 0.0);   // xDotDot0
    EXPECT_DOUBLE_EQ(XDot0[4],
                     1 / m1 * ((-c1 - c2) * 14 + c2 * 14));    // xDotDot1
    EXPECT_DOUBLE_EQ(XDot0[5], 1 / m2 * (c2 * 14 - c2 * 14));  // xDotDot1

    EXPECT_FALSE(p.Integrate(1.0).isError);

    auto e = p.GetMassMinAccel(2);
    ASSERT_FALSE(e.isError);
    ASSERT_TRUE(-80 < e.val && e.val < -60);

    e = p.GetMassMaxAccel(2);
    ASSERT_FALSE(e.isError);
    ASSERT_TRUE(0 < e.val && e.val < 10);

    ASSERT_DOUBLE_EQ(p.GetMassMaxAbsAccel(2).val,
                     abs(p.GetMassMinAccel(2).val));
}

TEST(ProblemTest, MultiBodyBibliographyDataTest2) {
    // This test simulates a system with 5 masses, similar to one found in the
    // bibliography and plots it's response so that we can compare it with the
    // results in the bibliography.
    // Source:
    // Mostafa, Marzbanrad Javad And. 2011. “A System Identification Algorithm
    // for Vehicle Lumped.” International Journal of Modeling and Optimization 1
    // (January): 163–66. Figure 11

    Problem p = Problem();
    EXPECT_FALSE(p.AddMass(1.0, 0.0, 0.0).isError);  // m0
    EXPECT_FALSE(p.AddMass(300, 1.0, 1.0).isError);  // m1
    EXPECT_FALSE(p.AddMass(120, 1.0, 0.0).isError);  // m2
    EXPECT_FALSE(p.AddMass(150, 1.0, 3.0).isError);  // m3
    EXPECT_FALSE(p.AddMass(700, 2.0, 0.0).isError);  // m4
    EXPECT_FALSE(p.AddMass(80, 3.0, 0.0).isError);   // m5

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

    EXPECT_FALSE(p.AddSpring(0, 1, k1).isError);
    EXPECT_FALSE(p.AddSpring(1, 2, k2).isError);
    EXPECT_FALSE(p.AddSpring(1, 3, k3).isError);
    EXPECT_FALSE(p.AddSpring(1, 4, k4).isError);
    EXPECT_FALSE(p.AddDamper(1, 4, c4).isError);
    EXPECT_FALSE(p.AddSpring(0, 2, k5).isError);
    EXPECT_FALSE(p.AddDamper(0, 2, c5).isError);
    EXPECT_FALSE(p.AddSpring(2, 4, k6).isError);
    EXPECT_FALSE(p.AddDamper(2, 4, c6).isError);
    EXPECT_FALSE(p.AddSpring(0, 3, k7).isError);
    EXPECT_FALSE(p.AddDamper(0, 3, c7).isError);
    EXPECT_FALSE(p.AddSpring(3, 4, k8).isError);
    EXPECT_FALSE(p.AddDamper(3, 4, c8).isError);
    EXPECT_FALSE(p.AddSpring(4, 5, k9).isError);
    EXPECT_FALSE(p.AddDamper(4, 5, c9).isError);

    p.Build();
    p.SetInitialVel(14.0);
    EXPECT_FALSE(p.FixMass(0).isError);

    EXPECT_FALSE(p.Integrate(0.15).isError);

    // The velocities we get match with the paper. The acceleration follows the
    // same shape, but the paper doesn't show the units for the acceleration.
    // If our data is correct, it seems to be (10m/s^2). To validate the
    // accelerations we get is also correct, we do a single check using linear
    // interpolation at a random time instant.
    int i = p.XHistory.size() / 2;
    double t0 = p.t[i];
    double v0 = p.XHistory[i][p.GetMassVelIndex(5)];
    double vDot0 = p.AccelHistory[i][5];
    double t1 = p.t[i + 1];
    double v1 = p.XHistory[i + 1][p.GetMassVelIndex(5)];
    double expectedV1 = v0 + vDot0 * (t1 - t0);
    double err = abs((expectedV1 - v1) / v1);
    EXPECT_TRUE(err < 0.0001);

    p.PrintMassTimeHistory(5);

    auto e = p.GetMassMinAccel(5);
    ASSERT_FALSE(e.isError);
    ASSERT_TRUE(-400 < e.val && e.val < -350);

    e = p.GetMassMaxAccel(5);
    ASSERT_FALSE(e.isError);
    ASSERT_TRUE(0 < e.val && e.val < 100);

    // p.PrintMassTimeHistory(5);
}

TEST(ProblemTest, MinMaxAccelTest) {
    Problem p = Problem();

    auto e = p.GetMassMaxAbsAccel(0);
    ASSERT_TRUE(e.isError);
    e = p.GetMassMaxAccel(0);
    ASSERT_TRUE(e.isError);
    e = p.GetMassMinAccel(0);
    ASSERT_TRUE(e.isError);

    p.AddMass(1.0, 0.0, 0.0);
    p.AddMass(1.0, 1.0, 0.0);
    p.AddMass(1.0, 2.0, 0.0);

    // Called to prevent "problem must be integrated" error.
    // AccelHistory is then cleared and mocked for the test
    p.Build();
    p.Integrate(0.0);
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

TEST(ProblemTest, ClearHistoryTest) {
    Problem p = Problem();

    // Create dummy problem
    ASSERT_FALSE(p.AddMass(1.0, 0.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 1.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 2.0, 0.0).isError);
    ASSERT_FALSE(p.AddSpring(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddSpring(1, 2, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(1, 2, 1.0).isError);

    p.Build();
    ASSERT_FALSE(p.SetInitialVel(10).isError);
    ASSERT_FALSE(p.FixMass(0).isError);

    ASSERT_FALSE(p.Integrate(1.0).isError);

    // Store some values of state vector and accelerations to compare latter
    ASSERT_TRUE(int(p.XHistory.size()) > 3);
    int XHistory0Size = int(p.XHistory.size());
    int AccelHistory0Size = int(p.AccelHistory.size());
    vector<double> X0 = vector<double>(p.XHistory[0].size());
    vector<double> XN = vector<double>(p.XHistory[0].size());
    vector<double> A0 = vector<double>(p.AccelHistory[0].size());
    vector<double> AN = vector<double>(p.AccelHistory[0].size());

    int n = int(p.XHistory.size()) / 2;
    std::copy(p.XHistory[0].begin(), p.XHistory[0].end(), X0.begin());
    std::copy(p.XHistory[n].begin(), p.XHistory[n].end(), XN.begin());
    std::copy(p.AccelHistory[0].begin(), p.AccelHistory[0].end(), A0.begin());
    std::copy(p.AccelHistory[n].begin(), p.AccelHistory[n].end(), AN.begin());

    p.ClearHistory();
    ASSERT_FALSE(p.isIntegrated);
    ASSERT_EQ(p.XHistory.size(), 0);
    ASSERT_EQ(p.AccelHistory.size(), 0);

    // Integrate again
    ASSERT_FALSE(p.Integrate(1.0).isError);
    ASSERT_TRUE(p.isIntegrated);

    // Compare values
    ASSERT_EQ(p.XHistory.size(), XHistory0Size);
    ASSERT_EQ(p.AccelHistory.size(), AccelHistory0Size);
    for (int i = 0; i < int(p.XHistory[0].size()); i++) {
        ASSERT_EQ(p.XHistory[0][i], X0[i]);
        ASSERT_EQ(p.XHistory[n][i], XN[i]);
    }
    for (int i = 0; i < int(p.AccelHistory[0].size()); i++) {
        ASSERT_EQ(p.AccelHistory[0][i], A0[i]);
        ASSERT_EQ(p.AccelHistory[n][i], AN[i]);
    }
}

TEST(ProblemTest, AssignmentTest) {
    Problem p = Problem();
    ASSERT_FALSE(p.AddMass(1.0, 0.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 1.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 2.0, 0.0).isError);
    ASSERT_FALSE(p.AddSpring(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddSpring(1, 2, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(1, 2, 1.0).isError);
    p.Build();
    ASSERT_FALSE(p.SetInitialVel(10).isError);
    ASSERT_FALSE(p.FixMass(0).isError);
    ASSERT_FALSE(p.Integrate(1.0).isError);

    Problem pA = Problem();
    pA = p;

    ASSERT_TRUE(p.masses.size() == pA.masses.size());
    for (int i = 0; i < int(p.masses.size()); i++) {
        ASSERT_TRUE(p.masses[i].m == pA.masses[i].m);
        ASSERT_TRUE(p.masses[i].px == pA.masses[i].px);
        ASSERT_TRUE(p.masses[i].py == pA.masses[i].py);
        ASSERT_TRUE(p.masses[i].xIndex == pA.masses[i].xIndex);
    }

    ASSERT_TRUE(p.springs.size() == pA.springs.size());
    for (int i = 0; i < int(p.springs.size()); i++) {
        ASSERT_TRUE(p.springs[i].k == pA.springs[i].k);
        ASSERT_TRUE(p.springs[i].m0->xIndex == pA.springs[i].m0->xIndex);
        ASSERT_TRUE(p.springs[i].m1->xIndex == pA.springs[i].m1->xIndex);
    }

    ASSERT_TRUE(p.dampers.size() == pA.dampers.size());
    for (int i = 0; i < int(p.dampers.size()); i++) {
        ASSERT_TRUE(p.dampers[i].c == pA.dampers[i].c);
        ASSERT_TRUE(p.dampers[i].m0->xIndex == pA.dampers[i].m0->xIndex);
        ASSERT_TRUE(p.dampers[i].m1->xIndex == pA.dampers[i].m1->xIndex);
    }

    ASSERT_TRUE(p.MInv.size1() == pA.MInv.size1());
    ASSERT_TRUE(p.MInv.size2() == pA.MInv.size2());
    for (int i = 0; i < int(p.MInv.size1()); i++) {
        for (int j = 0; j < int(p.MInv.size2()); j++) {
            ASSERT_TRUE(p.MInv(i, j) == pA.MInv(i, j));
        }
    }

    ASSERT_TRUE(p.K.size1() == pA.K.size1());
    ASSERT_TRUE(p.K.size2() == pA.K.size2());
    for (int i = 0; i < int(p.K.size1()); i++) {
        for (int j = 0; j < int(p.K.size2()); j++) {
            ASSERT_TRUE(p.K(i, j) == pA.K(i, j));
        }
    }

    ASSERT_TRUE(p.C.size1() == pA.C.size1());
    ASSERT_TRUE(p.C.size2() == pA.C.size2());
    for (int i = 0; i < int(p.C.size1()); i++) {
        for (int j = 0; j < int(p.C.size2()); j++) {
            ASSERT_TRUE(p.C(i, j) == pA.C(i, j));
        }
    }

    ASSERT_TRUE(p.X.size() == pA.X.size());
    for (int i = 0; i < int(p.X.size()); i++) {
        ASSERT_TRUE(p.X[i] == pA.X[i]);
    }

    ASSERT_TRUE(p.t.size() == pA.t.size());
    for (int i = 0; i < int(p.t.size()); i++) {
        ASSERT_TRUE(p.t[i] == pA.t[i]);
    }

    ASSERT_TRUE(p.XHistory.size() == pA.XHistory.size());
    for (int i = 0; i < int(p.XHistory.size()); i++) {
        ASSERT_TRUE(p.XHistory[i].size() == pA.XHistory[i].size());
        for (int j = 0; j < int(p.XHistory[i].size()); j++) {
            ASSERT_TRUE(p.XHistory[i][j] == pA.XHistory[i][j]);
        }
    }

    ASSERT_TRUE(p.AccelHistory.size() == pA.AccelHistory.size());
    for (int i = 0; i < int(p.AccelHistory.size()); i++) {
        ASSERT_TRUE(p.AccelHistory[i].size() == pA.AccelHistory[i].size());
        for (int j = 0; j < int(p.AccelHistory[i].size()); j++) {
            ASSERT_TRUE(p.AccelHistory[i][j] == pA.AccelHistory[i][j]);
        }
    }

    for (auto fm : p.fixedMasses) {
        Maybe<Mass*> mA = pA.GetMass(fm->xIndex);
        ASSERT_FALSE(mA.isError);
        ASSERT_TRUE(pA.fixedMasses.find(mA.val) != pA.fixedMasses.end());
    }
    for (auto fm : pA.fixedMasses) {
        Maybe<Mass*> m = p.GetMass(fm->xIndex);
        ASSERT_FALSE(m.isError);
        ASSERT_TRUE(p.fixedMasses.find(m.val) != p.fixedMasses.end());
    }

    ASSERT_TRUE(p.isBuilt == pA.isBuilt);
    ASSERT_TRUE(p.isIntegrated == pA.isIntegrated);
}

TEST(ProblemTest, CopyConstructorTest) {
    // Same as AssignmentTest but using copy constructor
    Problem p = Problem();
    ASSERT_FALSE(p.AddMass(1.0, 0.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 1.0, 0.0).isError);
    ASSERT_FALSE(p.AddMass(1.0, 2.0, 0.0).isError);
    ASSERT_FALSE(p.AddSpring(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddSpring(1, 2, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(0, 1, 1.0).isError);
    ASSERT_FALSE(p.AddDamper(1, 2, 1.0).isError);
    p.Build();
    ASSERT_FALSE(p.SetInitialVel(10).isError);
    ASSERT_FALSE(p.FixMass(0).isError);
    ASSERT_FALSE(p.Integrate(1.0).isError);

    Problem pC = p;

    ASSERT_TRUE(p.masses.size() == pC.masses.size());
    for (int i = 0; i < int(p.masses.size()); i++) {
        ASSERT_TRUE(p.masses[i].m == pC.masses[i].m);
        ASSERT_TRUE(p.masses[i].px == pC.masses[i].px);
        ASSERT_TRUE(p.masses[i].py == pC.masses[i].py);
        ASSERT_TRUE(p.masses[i].xIndex == pC.masses[i].xIndex);
    }

    ASSERT_TRUE(p.springs.size() == pC.springs.size());
    for (int i = 0; i < int(p.springs.size()); i++) {
        ASSERT_TRUE(p.springs[i].k == pC.springs[i].k);
        ASSERT_TRUE(p.springs[i].m0->xIndex == pC.springs[i].m0->xIndex);
        ASSERT_TRUE(p.springs[i].m1->xIndex == pC.springs[i].m1->xIndex);
    }

    ASSERT_TRUE(p.dampers.size() == pC.dampers.size());
    for (int i = 0; i < int(p.dampers.size()); i++) {
        ASSERT_TRUE(p.dampers[i].c == pC.dampers[i].c);
        ASSERT_TRUE(p.dampers[i].m0->xIndex == pC.dampers[i].m0->xIndex);
        ASSERT_TRUE(p.dampers[i].m1->xIndex == pC.dampers[i].m1->xIndex);
    }

    ASSERT_TRUE(p.MInv.size1() == pC.MInv.size1());
    ASSERT_TRUE(p.MInv.size2() == pC.MInv.size2());
    for (int i = 0; i < int(p.MInv.size1()); i++) {
        for (int j = 0; j < int(p.MInv.size2()); j++) {
            ASSERT_TRUE(p.MInv(i, j) == pC.MInv(i, j));
        }
    }

    ASSERT_TRUE(p.K.size1() == pC.K.size1());
    ASSERT_TRUE(p.K.size2() == pC.K.size2());
    for (int i = 0; i < int(p.K.size1()); i++) {
        for (int j = 0; j < int(p.K.size2()); j++) {
            ASSERT_TRUE(p.K(i, j) == pC.K(i, j));
        }
    }

    ASSERT_TRUE(p.C.size1() == pC.C.size1());
    ASSERT_TRUE(p.C.size2() == pC.C.size2());
    for (int i = 0; i < int(p.C.size1()); i++) {
        for (int j = 0; j < int(p.C.size2()); j++) {
            ASSERT_TRUE(p.C(i, j) == pC.C(i, j));
        }
    }

    ASSERT_TRUE(p.X.size() == pC.X.size());
    for (int i = 0; i < int(p.X.size()); i++) {
        ASSERT_TRUE(p.X[i] == pC.X[i]);
    }

    ASSERT_TRUE(p.t.size() == pC.t.size());
    for (int i = 0; i < int(p.t.size()); i++) {
        ASSERT_TRUE(p.t[i] == pC.t[i]);
    }

    ASSERT_TRUE(p.XHistory.size() == pC.XHistory.size());
    for (int i = 0; i < int(p.XHistory.size()); i++) {
        ASSERT_TRUE(p.XHistory[i].size() == pC.XHistory[i].size());
        for (int j = 0; j < int(p.XHistory[i].size()); j++) {
            ASSERT_TRUE(p.XHistory[i][j] == pC.XHistory[i][j]);
        }
    }

    ASSERT_TRUE(p.AccelHistory.size() == pC.AccelHistory.size());
    for (int i = 0; i < int(p.AccelHistory.size()); i++) {
        ASSERT_TRUE(p.AccelHistory[i].size() == pC.AccelHistory[i].size());
        for (int j = 0; j < int(p.AccelHistory[i].size()); j++) {
            ASSERT_TRUE(p.AccelHistory[i][j] == pC.AccelHistory[i][j]);
        }
    }

    for (auto fm : p.fixedMasses) {
        Maybe<Mass*> mA = pC.GetMass(fm->xIndex);
        ASSERT_FALSE(mA.isError);
        ASSERT_TRUE(pC.fixedMasses.find(mA.val) != pC.fixedMasses.end());
    }
    for (auto fm : pC.fixedMasses) {
        Maybe<Mass*> m = p.GetMass(fm->xIndex);
        ASSERT_FALSE(m.isError);
        ASSERT_TRUE(p.fixedMasses.find(m.val) != p.fixedMasses.end());
    }

    ASSERT_TRUE(p.isBuilt == pC.isBuilt);
    ASSERT_TRUE(p.isIntegrated == pC.isIntegrated);
}