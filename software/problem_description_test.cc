#include "problem_description.h"

#include <gtest/gtest.h>

#include <memory>

#include "utils.h"

TEST(ProblemDescriptionTest, MassSpringDamperCreationTest) {
    // Instance to check errors returned by BuildRandom
    ProblemDescription pdErr = ProblemDescription();

    // Invalid damper
    pdErr.AddDamper(0, 0, 0.0, 0.0);
    ASSERT_TRUE(pdErr.BuildRandom().isError);

    // Invalid string
    pdErr = ProblemDescription();
    pdErr.AddSpring(0, 0, 0.0, 0.0);
    ASSERT_TRUE(pdErr.BuildRandom().isError);

    // Invalid mass (0.0 mass)
    pdErr = ProblemDescription();
    pdErr.AddMass(0.0, 1, 1);
    ASSERT_TRUE(pdErr.BuildRandom().isError);
    // Invalid mass (duplicate location)
    pdErr = ProblemDescription();
    pdErr.AddMass(1.0, 1.0, 1.0);
    pdErr.AddMass(1.0, 1.0, 1.0);
    ASSERT_TRUE(pdErr.BuildRandom().isError);

    // Instance stationary problem to check success
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(1.0, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    auto e0 = pd.BuildRandom();
    ASSERT_FALSE(e0.isError);

    Problem p = (*e0.val);
    ASSERT_FALSE(p.Integrate(0.0, 0.1, 0.01).isError);
    auto maxAccel0 = p.GetMassMaxAbsAccel(0);
    ASSERT_FALSE(maxAccel0.isError);
    ASSERT_DOUBLE_EQ(maxAccel0.val, 0.0);
    auto maxAccel1 = p.GetMassMaxAbsAccel(1);
    ASSERT_FALSE(maxAccel1.isError);
    ASSERT_DOUBLE_EQ(maxAccel1.val, 0.0);
}

TEST(ProblemDescriptionTest, DynamicTestWithInitialVel) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(1.0, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    pd.SetFixedMass(0);
    pd.AddInitialVel(1, 10.0);
    auto e0 = pd.BuildRandom();
    ASSERT_FALSE(e0.isError);

    Problem p = (*e0.val);
    ASSERT_FALSE(p.Integrate(0.0, 0.1, 0.01).isError);
    auto maxAccel0 = p.GetMassMaxAbsAccel(0);
    ASSERT_FALSE(maxAccel0.isError);
    ASSERT_DOUBLE_EQ(maxAccel0.val, 0.0);
    auto maxAccel1 = p.GetMassMinAccel(1);
    ASSERT_FALSE(maxAccel1.isError);
    ASSERT_TRUE(maxAccel1.val < 0.0);
}

TEST(ProblemDescriptionTest, DynamicTestWithInitialDisp) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(1.0, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.5, 1.0);
    pd.AddDamper(0, 1, 0.5, 1.0);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, 10.0);
    auto e0 = pd.BuildRandom();
    ASSERT_FALSE(e0.isError);

    Problem p = (*e0.val);
    ASSERT_FALSE(p.Integrate(0.0, 0.1, 0.01).isError);
    auto maxAccel0 = p.GetMassMaxAbsAccel(0);
    ASSERT_FALSE(maxAccel0.isError);
    ASSERT_DOUBLE_EQ(maxAccel0.val, 0.0);
    auto maxAccel1 = p.GetMassMinAccel(1);
    ASSERT_FALSE(maxAccel1.isError);
    ASSERT_TRUE(maxAccel1.val < 0.0);
}

TEST(ProblemDescriptionTest, MultiBodyBibliographyDataTest2) {
    // Equivalent of problem_test.cc TEST(ProblemTest,
    // MultiBodyBibliographyDataTest2)

    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(300, 1.0, 1.0);
    pd.AddMass(120, 1.0, 0.0);
    pd.AddMass(150, 1.0, 3.0);
    pd.AddMass(700, 2.0, 0.0);
    pd.AddMass(80, 3.0, 0.0);
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
    pd.AddSpring(0, 1, k1, k1);
    pd.AddSpring(1, 2, k2, k2);
    pd.AddSpring(1, 3, k3, k3);
    pd.AddSpring(1, 4, k4, k4);
    pd.AddDamper(1, 4, c4, c4);
    pd.AddSpring(0, 2, k5, k5);
    pd.AddDamper(0, 2, c5, c5);
    pd.AddSpring(2, 4, k6, k6);
    pd.AddDamper(2, 4, c6, c6);
    pd.AddSpring(0, 3, k7, k7);
    pd.AddDamper(0, 3, c7, c7);
    pd.AddSpring(3, 4, k8, k8);
    pd.AddDamper(3, 4, c8, c8);
    pd.AddSpring(4, 5, k9, k9);
    pd.AddDamper(4, 5, c9, c9);
    pd.AddInitialVel(14.0);
    pd.SetFixedMass(0);
    auto e = pd.BuildRandom();
    EXPECT_FALSE(e.isError);

    auto p = e.val;
    EXPECT_FALSE(p->Integrate(0.0, 0.15, 0.05).isError);

    auto e2 = p->GetMassMinAccel(5);
    ASSERT_FALSE(e2.isError);
    ASSERT_TRUE(-400 < e2.val && e2.val < -350);

    e2 = p->GetMassMaxAccel(5);
    ASSERT_FALSE(e2.isError);
    ASSERT_TRUE(0 < e2.val && e2.val < 100);
}

TEST(ProblemDescriptionTest, CreateMultipleTest) {
    // Tests creation of multiple Problem instances

    // Create arbitrary pd with a large range for spring and damper constants
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);
    pd.AddMass(300, 1.0, 1.0);
    pd.AddMass(120, 1.0, 0.0);
    pd.AddMass(150, 1.0, 3.0);
    pd.AddMass(700, 2.0, 0.0);
    pd.AddMass(80, 3.0, 0.0);
    pd.AddSpring(0, 1, 0.0, 1000000);
    pd.AddSpring(1, 2, 0.0, 1000000);
    pd.AddSpring(1, 3, 0.0, 1000000);
    pd.AddSpring(1, 4, 0.0, 1000000);
    pd.AddDamper(1, 4, 0.0, 1000000);
    pd.AddSpring(0, 2, 0.0, 1000000);
    pd.AddDamper(0, 2, 0.0, 1000000);
    pd.AddSpring(2, 4, 0.0, 1000000);
    pd.AddDamper(2, 4, 0.0, 1000000);
    pd.AddSpring(0, 3, 0.0, 1000000);
    pd.AddDamper(0, 3, 0.0, 1000000);
    pd.AddSpring(3, 4, 0.0, 1000000);
    pd.AddDamper(3, 4, 0.0, 1000000);
    pd.AddSpring(4, 5, 0.0, 1000000);
    pd.AddDamper(4, 5, 0.0, 1000000);
    pd.AddInitialVel(14.0);
    pd.SetFixedMass(0);

    // Get multiple problem instances and check they have different values
    auto e = pd.BuildRandom();
    EXPECT_FALSE(e.isError);
    auto p0 = e.val;

    e = pd.BuildRandom();
    EXPECT_FALSE(e.isError);
    auto p1 = e.val;

    e = pd.BuildRandom();
    EXPECT_FALSE(e.isError);
    auto p2 = e.val;

    EXPECT_FALSE(p0->Integrate(0.0, 0.15, 0.05).isError);
    EXPECT_FALSE(p1->Integrate(0.0, 0.15, 0.05).isError);
    EXPECT_FALSE(p2->Integrate(0.0, 0.15, 0.05).isError);

    auto a0 = p0->GetMassMaxAbsAccel(5);
    auto a1 = p1->GetMassMaxAbsAccel(5);
    auto a2 = p2->GetMassMaxAbsAccel(5);
    EXPECT_FALSE(a0.isError);
    EXPECT_FALSE(a1.isError);
    EXPECT_FALSE(a2.isError);

    EXPECT_TRUE(a0.val != a1.val);
    EXPECT_TRUE(a1.val != a2.val);
    EXPECT_TRUE(a0.val != a2.val);

    print("a0.val", a0.val);
    print("a1.val", a1.val);
    print("a2.val", a2.val);
}