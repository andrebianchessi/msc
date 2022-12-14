#include "pimodel.h"

#include <gtest/gtest.h>

#include <vector>

#include "bounded.h"
#include "problem.h"
#include "utils.h"

const double m = 3.0;
const double kMin = 1.0;
const double kMax = 2.0;
const double cMin = 10.0;
const double cMax = 20.0;
const double initialDisplacement = 10.0;

class PimodelTest : public testing::Test {
   public:
    ProblemDescription pd;
    double x;
    double y;
    std::vector<double> xy;

    void SetUp() {
        // Called before every TEST_F
        // Create problem description of two masses,
        // connected with a spring and damper. One is fixed
        // and the other has initial displacement.
        this->pd = ProblemDescription();
        this->pd.AddMass(m, 0.0, 0.0);
        this->pd.AddMass(m, 1.0, 0.0);
        this->pd.AddSpring(0, 1, kMin, kMax);
        this->pd.AddDamper(0, 1, cMin, cMax);
        this->pd.SetFixedMass(0);
        this->pd.AddInitialDisp(1, initialDisplacement);
        std::vector<Bounded> dna = std::vector<Bounded>(2);
        auto e0 = this->pd.BuildFromDNA(dna);
        ASSERT_FALSE(e0.isError);

        x = 19354;
        y = 1235;
        xy = std::vector<double>{x, y};
    }
};

TEST_F(PimodelTest, ConstructorTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 10, 2);

    // There should be one polynomial for each mass
    ASSERT_EQ(model.polys.size1(), 2);
    ASSERT_EQ(model.polys.size2(), 1);
}

TEST_F(PimodelTest, OperatorTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 10, 2);

    // Test the model's prediction for values of time, spring and damper
    std::vector<double> tkc = {1.0, 2.0, 3.0};
    auto eval = model(&tkc);
    ASSERT_FALSE(eval.isError);

    // There are two masses; the returned vector should contain
    // the position of both of them
    ASSERT_EQ(eval.val.size(), 2);

    // Given that the polynomials all start with parameters = 0,
    // the initial prediction (before training) should be position = 0
    // for all masses
    ASSERT_DOUBLE_EQ(eval.val[0], 0);
    ASSERT_DOUBLE_EQ(eval.val[1], 0);

    // Test error cases
    // tkc too large
    tkc = {1.1, 1.0, 2.0, 3.0};
    eval = model(&tkc);
    ASSERT_TRUE(eval.isError);

    // tkc too small
    tkc = {1.1, 1.0};
    eval = model(&tkc);
    ASSERT_TRUE(eval.isError);
}

TEST_F(PimodelTest, SetParametersTest) {
    // The polynomial of each mass will be of order 1:
    // x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    // x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 10, 1);

    std::vector<double> params = {1, 2, 3, 4, 5, 6, 7, 8};
    auto r = model.SetParameters(&params);
    ASSERT_FALSE(r.isError);

    double t = 99;
    double k = 98;
    double c = 97;
    std::vector<double> tkc = {t, k, c};
    auto eval = model(&tkc);
    ASSERT_FALSE(eval.isError);

    ASSERT_EQ(eval.val.size(), 2);

    ASSERT_DOUBLE_EQ(eval.val[0], 1 * t + 2 * k + 3 * c + 4);
    ASSERT_DOUBLE_EQ(eval.val[1], 5 * t + 6 * k + 7 * c + 8);

    // Test error cases
    params = {1, 2, 3, 4, 5, 6, 7};
    r = model.SetParameters(&params);
    ASSERT_TRUE(r.isError);
    params = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    r = model.SetParameters(&params);
    ASSERT_TRUE(r.isError);
}

TEST_F(PimodelTest, GetParametersTest) {
    // The polynomial of each mass will be of order 1:
    // x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    // x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 10, 1);

    std::vector<double> params = std::vector<double>(8);
    auto r = model.GetParameters(&params);
    ASSERT_FALSE(r.isError);

    ASSERT_DOUBLE_EQ(params[0], 0.0);
    ASSERT_DOUBLE_EQ(params[1], 0.0);
    ASSERT_DOUBLE_EQ(params[2], 0.0);
    ASSERT_DOUBLE_EQ(params[3], 0.0);
    ASSERT_DOUBLE_EQ(params[4], 0.0);
    ASSERT_DOUBLE_EQ(params[5], 0.0);
    ASSERT_DOUBLE_EQ(params[6], 0.0);
    ASSERT_DOUBLE_EQ(params[7], 0.0);

    params = {1, 2, 3, 4, 5, 6, 7, 8};
    r = model.SetParameters(&params);
    ASSERT_FALSE(r.isError);
    r = model.GetParameters(&params);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(params[0], 1.0);
    ASSERT_DOUBLE_EQ(params[1], 2.0);
    ASSERT_DOUBLE_EQ(params[2], 3.0);
    ASSERT_DOUBLE_EQ(params[3], 4.0);
    ASSERT_DOUBLE_EQ(params[4], 5.0);
    ASSERT_DOUBLE_EQ(params[5], 6.0);
    ASSERT_DOUBLE_EQ(params[6], 7.0);
    ASSERT_DOUBLE_EQ(params[7], 8.0);

    // Test error cases
    params = std::vector<double>(7);
    r = model.GetParameters(&params);
    ASSERT_TRUE(r.isError);
    params = std::vector<double>(9);
    r = model.GetParameters(&params);
    ASSERT_TRUE(r.isError);
}
TEST_F(PimodelTest, DaTest) {
    Poly p0;
    p0.Build(2, 1, 1.0);
    // p0 = 1.0*x + 1.0*y + 1.0*1

    Poly p1;
    p1.Build(2, 2, 1.0);
    // p1 = 1.0*x^2 + 1.0*xy + 1.0*x + 1.0*y^2 + 1.0*y + 1.0*1

    std::vector<double> grad = std::vector<double>(9);

    std::vector<Poly> v = std::vector<Poly>{p0, p1};
    ASSERT_FALSE(Da(&v, &xy, &grad).isError);
    ASSERT_DOUBLE_EQ(grad[0], x);
    ASSERT_DOUBLE_EQ(grad[1], y);
    ASSERT_DOUBLE_EQ(grad[2], 1);
    ASSERT_DOUBLE_EQ(grad[3], x * x);
    ASSERT_DOUBLE_EQ(grad[4], x * y);
    ASSERT_DOUBLE_EQ(grad[5], x);
    ASSERT_DOUBLE_EQ(grad[6], y * y);
    ASSERT_DOUBLE_EQ(grad[7], y);
    ASSERT_DOUBLE_EQ(grad[8], 1);

    boost::numeric::ublas::matrix<Poly> m =
        boost::numeric::ublas::matrix<Poly>(2, 1);
    m(0, 0) = p0;
    m(1, 0) = p1;
    grad = std::vector<double>(9);
    ASSERT_FALSE(Da(&m, &xy, &grad).isError);
    ASSERT_DOUBLE_EQ(grad[0], x);
    ASSERT_DOUBLE_EQ(grad[1], y);
    ASSERT_DOUBLE_EQ(grad[2], 1);
    ASSERT_DOUBLE_EQ(grad[3], x * x);
    ASSERT_DOUBLE_EQ(grad[4], x * y);
    ASSERT_DOUBLE_EQ(grad[5], x);
    ASSERT_DOUBLE_EQ(grad[6], y * y);
    ASSERT_DOUBLE_EQ(grad[7], y);
    ASSERT_DOUBLE_EQ(grad[8], 1);
}

TEST_F(PimodelTest, InitialConditionsLossTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 1, 1, 1);
    std::vector<double> params = std::vector<double>(8);
    params = {1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_FALSE(model.SetParameters(&params).isError);
    // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    // modelXDot0(t,k,c) = 1
    // modelXDot1(t,k,c) = 5

    double expectedLoss = 0;

    // Initial values for the problem created in SetUp()
    double initialX0 = 0;
    double initialX1 = initialDisplacement;
    double initialXDot0 = 0;
    double initialXDot1 = 0;
    double t = 0;
    // Initial conditions loss:
    for (double k : std::vector<double>{kMin, kMax}) {
        for (double c : std::vector<double>{cMin, cMax}) {
            double modelX0 = 1 * t + 2 * k + 3 * c + 4 * 1;
            double modelX1 = 5 * t + 6 * k + 7 * c + 8 * 1;
            double modelX0Dot = 1;
            double modelX1Dot = 5;
            expectedLoss += pow(modelX0 - initialX0, 2);
            expectedLoss += pow(modelX1 - initialX1, 2);
            expectedLoss += pow(modelX0Dot - initialXDot0, 2);
            expectedLoss += pow(modelX1Dot - initialXDot1, 2);
        }
    }

    std::vector<double> tkc = std::vector<double>(model.inputSize());
    tkc[0] = 0.0;
    double loss = 0;
    model.InitialConditionsLossDfs(&tkc, 1, &loss);

    ASSERT_DOUBLE_EQ(expectedLoss, loss);
}

TEST_F(PimodelTest, getAccelsFromDiffEqTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 1, 1, 1);
    std::vector<double> params = std::vector<double>(8);
    params = {1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_FALSE(model.SetParameters(&params).isError);
    // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    // modelXDot0(t,k,c) = 1*t^0 + 0*k + 0*c + 0*1
    // modelXDot1(t,k,c) = 5*t^0 + 0*k + 0*c + 0*1

    Maybe<Problem> mProblem =
        this->pd.BuildFromVector(std::vector<double>{kMin, cMin});
    ASSERT_FALSE(mProblem.isError);
    Problem problem = mProblem.val;

    auto A = model.getAccelsFromDiffEq(&problem);
    ASSERT_EQ(A.size1(), 2);
    ASSERT_EQ(A.size2(), 1);

    // expectedAccel0(t,k,c) = 1 / m * (-kMin * x0 + kMin * x1 - cMin * dx0dt +
    // cMin * dx1dt);
    // expectedAccel0 =
    //                 1 / m * (-kMin * (1*t + 2*k + 3*c + 4*1) + kMin * (5*t +
    //                 6*k + 7*c + 8*1) - cMin * 1 + cMin * 5);
    // expectedAccel0 = (4*kMin/m)*t + (4*kMin)/m*k + (4 *kMin)/m*c  +
    // (4*kMin+4*cMin)/m*1
    std::vector<double> coefs = std::vector<double>(4);
    ASSERT_FALSE(A(0, 0).GetCoefficients(&coefs).isError);
    ASSERT_DOUBLE_EQ(coefs[0], 4 * kMin / m);
    ASSERT_DOUBLE_EQ(coefs[1], 4 * kMin / m);
    ASSERT_DOUBLE_EQ(coefs[2], 4 * kMin / m);
    ASSERT_DOUBLE_EQ(coefs[3], (4 * kMin + 4 * cMin) / m);
    double t = 0;
    double k = 923847;
    double c = 74437;
    double expectedEval =
        1 / m *
        (-kMin * (1 * t + 2 * k + 3 * c + 4 * 1) +
         kMin * (5 * t + 6 * k + 7 * c + 8 * 1) - cMin * 1 + cMin * 5);
    std::vector<double> tkc = {t, k, c};
    ASSERT_DOUBLE_EQ(A(0, 0)(&tkc).val, expectedEval);

    // expectedAccel1 = -1*expectedAccel0
    coefs = std::vector<double>(4);
    ASSERT_FALSE(A(1, 0).GetCoefficients(&coefs).isError);
    EXPECT_DOUBLE_EQ(coefs[0], -4 * kMin / m);
    EXPECT_DOUBLE_EQ(coefs[1], -4 * kMin / m);
    EXPECT_DOUBLE_EQ(coefs[2], -4 * kMin / m);
    EXPECT_DOUBLE_EQ(coefs[3], -(4 * kMin + 4 * cMin) / m);
}

TEST_F(PimodelTest, PhysicsLossTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 1, 1, 1);
    std::vector<double> params = std::vector<double>(8);
    params = {1, 2, 3, 4, 5, 6, 7, 8};
    ASSERT_FALSE(model.SetParameters(&params).isError);
    // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    // modelXDot0(t,k,c) = 1
    // modelXDot1(t,k,c) = 5
    // modelXDotDot0(t,k,c) = 0
    // modelXDotDot1(t,k,c) = 0

    double expectedLoss = 0;

    // Physics loss:
    for (double t : std::vector<double>{0, 1.0}) {
        for (double k : std::vector<double>{kMin, kMax}) {
            for (double c : std::vector<double>{cMin, cMax}) {
                double modelX0 = 1 * t + 2 * k + 3 * c + 4 * 1;
                double modelX1 = 5 * t + 6 * k + 7 * c + 8 * 1;
                double modelX0Dot = 1;
                double modelX1Dot = 5;
                double modelX0DotDot = 0;
                double modelX1DotDot = 0;

                double x0DotDot = 1 / m *
                                  (-k * modelX0 + k * modelX1 - c * modelX0Dot +
                                   c * modelX1Dot);
                double x1DotDot = 1 / m *
                                  (k * modelX0 - k * modelX1 + c * modelX0Dot -
                                   c * modelX1Dot);
                expectedLoss += pow(modelX0DotDot - x0DotDot, 2);
                expectedLoss += pow(modelX1DotDot - x1DotDot, 2);
            }
        }
    }

    std::vector<double> tkc = std::vector<double>(model.inputSize());
    double loss = 0;
    model.PhysicsLossDfs(&tkc, 0, &loss);

    ASSERT_DOUBLE_EQ(expectedLoss, loss);
}

TEST_F(PimodelTest, LossTest2) {
    // Similar to LossTest, but with second order polynomials
    // and more refined discretization
    Pimodel model = Pimodel(&this->pd, 1.0, 2, 2, 2);
    std::vector<double> params = std::vector<double>(8);
    params = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
              11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    ASSERT_FALSE(model.SetParameters(&params).isError);
    // modelX0(t,k,c) = 1*t^2 + 2*tk + 3*tc + 4*t + 5*k^2 + 6*kc + 7*k
    // + 8*c^2 + 9*c + 10*1
    // modelX1(t,k,c) = 11*t^2 + 12*tk + 13*tc + 14*t + 15*k^2 + 16*kc + 17*k
    // + 18*c^2 + 19*c + 20*1
    // modelX0Dot(t,k,c) = 2t + 2k + 3c + 4
    // modelX1Dot(t,k,c) = 22t + 12k + 13c + 14
    // modelX0DotDot(t,k,c) = 2
    // modelX1DotDot(t,k,c) = 22

    double expectedLoss = 0;

    // Initial values for the problem created in SetUp()
    double initialX0 = 0;
    double initialX1 = initialDisplacement;
    double initialXDot0 = 0;
    double initialXDot1 = 0;
    double t = 0;
    // Initial conditions loss:
    for (double k : std::vector<double>{kMin, (kMin + kMax) / 2, kMax}) {
        for (double c : std::vector<double>{cMin, (cMin + cMax) / 2, cMax}) {
            double modelX0 = 1 * t * t + 2 * t * k + 3 * t * c + 4 * t +
                             5 * k * k + 6 * k * c + 7 * k + 8 * c * c + 9 * c +
                             10 * 1;
            double modelX1 = 11 * t * t + 12 * t * k + 13 * t * c + 14 * t +
                             15 * k * k + 16 * k * c + 17 * k + 18 * c * c +
                             19 * c + 20 * 1;
            double modelX0Dot = 2 * t + 2 * k + 3 * c + 4;
            double modelX1Dot = 22 * t + 12 * k + 13 * c + 14;
            expectedLoss += pow(modelX0 - initialX0, 2);
            expectedLoss += pow(modelX1 - initialX1, 2);
            expectedLoss += pow(modelX0Dot - initialXDot0, 2);
            expectedLoss += pow(modelX1Dot - initialXDot1, 2);
        }
    }

    // Physics loss:
    for (double t : std::vector<double>{0, 0.5, 1.0}) {
        for (double k : std::vector<double>{kMin, (kMin + kMax) / 2, kMax}) {
            for (double c :
                 std::vector<double>{cMin, (cMin + cMax) / 2, cMax}) {
                double modelX0 = 1 * t * t + 2 * t * k + 3 * t * c + 4 * t +
                                 5 * k * k + 6 * k * c + 7 * k + 8 * c * c +
                                 9 * c + 10 * 1;
                double modelX1 = 11 * t * t + 12 * t * k + 13 * t * c + 14 * t +
                                 15 * k * k + 16 * k * c + 17 * k + 18 * c * c +
                                 19 * c + 20 * 1;
                double modelX0Dot = 2 * t + 2 * k + 3 * c + 4;
                double modelX1Dot = 22 * t + 12 * k + 13 * c + 14;
                double modelX0DotDot = 2;
                double modelX1DotDot = 22;

                double x0DotDot = 1 / m *
                                  (-k * modelX0 + k * modelX1 - c * modelX0Dot +
                                   c * modelX1Dot);
                double x1DotDot = 1 / m *
                                  (k * modelX0 - k * modelX1 + c * modelX0Dot -
                                   c * modelX1Dot);
                expectedLoss += pow(modelX0DotDot - x0DotDot, 2);
                expectedLoss += pow(modelX1DotDot - x1DotDot, 2);
            }
        }
    }

    double loss = model.Loss();

    ASSERT_DOUBLE_EQ(expectedLoss, loss);
}

// TEST_F(PimodelTest, LossGradientTest) {
//     Pimodel model = Pimodel(&this->pd, 1.0, 1, 1, 1);
//     std::vector<double> params = std::vector<double>(8);
//     params = {1, 2, 3, 4, 5, 6, 7, 8};
//     ASSERT_FALSE(model.SetParameters(&params).isError);
//     // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
//     // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
//     // modelXDot0(t,k,c) = 1
//     // modelXDot1(t,k,c) = 5
//     // modelXDotDot0(t,k,c) = 0
//     // modelXDotDot1(t,k,c) = 0

//     std::vector<double> expectedGrad = std::vector<double>(8);

//     // Initial values for the problem created in SetUp()
//     double initialX0 = 0;
//     double initialX1 = initialDisplacement;
//     double initialX0Dot = 0;
//     double initialX1Dot = 0;

//     double t = 0;
//     // Initial conditions loss:
//     for (double k : std::vector<double>{kMin, kMax}) {
//         for (double c : std::vector<double>{cMin, cMax}) {
//             double modelX0 = 1 * t + 2 * k + 3 * c + 4 * 1;
//             double modelX1 = 5 * t + 6 * k + 7 * c + 8 * 1;
//             double modelX0Dot = 1;
//             double modelX1Dot = 5;

//             double d_da0_modelX0 = t;
//             double d_da1_modelX0 = k;
//             double d_da2_modelX0 = c;
//             double d_da3_modelX0 = 1;
//             expectedGrad[0] += 2 * (modelX0 - initialX0) * d_da0_modelX0;
//             expectedGrad[0] += 2 * (modelX0Dot - initialX0Dot) *
//             d_da0_modelX0; expectedGrad[1] += 2 * (modelX0 - initialX0) *
//             d_da1_modelX0; expectedGrad[1] += 2 * (modelX0Dot - initialX0Dot)
//             * d_da1_modelX0; expectedGrad[2] += 2 * (modelX0 - initialX0) *
//             d_da2_modelX0; expectedGrad[2] += 2 * (modelX0Dot - initialX0Dot)
//             * d_da2_modelX0; expectedGrad[3] += 2 * (modelX0 - initialX0) *
//             d_da3_modelX0; expectedGrad[3] += 2 * (modelX0Dot - initialX0Dot)
//             * d_da3_modelX0;

//             double d_da0_modelX1 = t;
//             double d_da1_modelX1 = k;
//             double d_da2_modelX1 = c;
//             double d_da3_modelX1 = 1;
//             expectedGrad[4] += 2 * (modelX1 - initialX1) * d_da0_modelX1;
//             expectedGrad[4] += 2 * (modelX1Dot - initialX1Dot) *
//             d_da0_modelX1; expectedGrad[5] += 2 * (modelX1 - initialX1) *
//             d_da1_modelX1; expectedGrad[5] += 2 * (modelX1Dot - initialX1Dot)
//             * d_da1_modelX1; expectedGrad[6] += 2 * (modelX1 - initialX1) *
//             d_da2_modelX1; expectedGrad[6] += 2 * (modelX1Dot - initialX1Dot)
//             * d_da2_modelX1; expectedGrad[7] += 2 * (modelX1 - initialX1) *
//             d_da3_modelX1; expectedGrad[7] += 2 * (modelX1Dot - initialX1Dot)
//             * d_da3_modelX1;
//         }
//     }

//     // Physics loss:
//     for (double t : std::vector<double>{0, 1.0}) {
//         for (double k : std::vector<double>{kMin, kMax}) {
//             for (double c : std::vector<double>{cMin, cMax}) {
//                 double modelX0 = 1 * t + 2 * k + 3 * c + 4 * 1;
//                 double modelX1 = 5 * t + 6 * k + 7 * c + 8 * 1;
//                 double modelX0Dot = 1;
//                 double modelX1Dot = 5;
//                 double modelX0DotDot = 0;
//                 double modelX1DotDot = 0;

//                 double x0DotDot = 0;  // mass 0 is fixed
//                 double x1DotDot = 1 / m *
//                                   (k * modelX0 - k * modelX1 + c * modelX0Dot
//                                   -
//                                    c * modelX1Dot);

//                 double d_da0_modelX0 = t;
//                 double d_da1_modelX0 = k;
//                 double d_da2_modelX0 = c;
//                 double d_da3_modelX0 = 1;
//                 expectedGrad[0] +=
//                     2 * (modelX0DotDot - x0DotDot) * d_da0_modelX0;
//                 expectedGrad[1] +=
//                     2 * (modelX0DotDot - x0DotDot) * d_da1_modelX0;
//                 expectedGrad[2] +=
//                     2 * (modelX0DotDot - x0DotDot) * d_da2_modelX0;
//                 expectedGrad[3] +=
//                     2 * (modelX0DotDot - x0DotDot) * d_da3_modelX0;

//                 double d_da0_modelX1 = t;
//                 double d_da1_modelX1 = k;
//                 double d_da2_modelX1 = c;
//                 double d_da3_modelX1 = 1;
//                 expectedGrad[4] +=
//                     2 * (modelX1DotDot - x1DotDot) * d_da0_modelX1;
//                 expectedGrad[5] +=
//                     2 * (modelX1DotDot - x1DotDot) * d_da1_modelX1;
//                 expectedGrad[6] +=
//                     2 * (modelX1DotDot - x1DotDot) * d_da2_modelX1;
//                 expectedGrad[7] +=
//                     2 * (modelX1DotDot - x1DotDot) * d_da3_modelX1;
//             }
//         }
//     }

//     std::vector<double> grad = model.LossGradient();

//     ASSERT_EQ(grad.size(), 8);
//     ASSERT_DOUBLE_EQ(grad[0], expectedGrad[0]);
//     ASSERT_DOUBLE_EQ(grad[1], expectedGrad[1]);
//     ASSERT_DOUBLE_EQ(grad[2], expectedGrad[2]);
//     ASSERT_DOUBLE_EQ(grad[3], expectedGrad[3]);
//     ASSERT_DOUBLE_EQ(grad[4], expectedGrad[4]);
//     ASSERT_DOUBLE_EQ(grad[5], expectedGrad[5]);
//     ASSERT_DOUBLE_EQ(grad[6], expectedGrad[6]);
//     ASSERT_DOUBLE_EQ(grad[7], expectedGrad[7]);
// }

// // Training is not showing good results. Better implement more unit-tests
// // to check if everything is working as intended.
// TEST_F(PimodelTest, TrainTest) {
//     double T = 15.0;
//     int timeDiscretization = 10;
//     int kcDiscretization = 1;
//     int order = 10;

//     // Train model
//     Pimodel model =
//         Pimodel(&this->pd, T, timeDiscretization, kcDiscretization, order);
//     double lr = 0.001;
//     model.Train(lr, lr / 4, true);

//     // Get problem using intermediate value for k and c, and integrate it
//     double k = (kMin + kMax) / 2;
//     double c = (cMin + cMax) / 2;
//     Maybe<Problem> mP = this->pd.BuildFromVector(std::vector<double>{k, c});
//     ASSERT_FALSE(mP.isError);
//     Problem p = mP.val;
//     p.Integrate(T);

//     // Compare the model's prediction with the problem's result
//     std::vector<double> tkc = std::vector<double>{0.0, k, c};
//     Maybe<std::vector<double>> X;
//     std::cout << "t,x0,modelX0,x1,modelX1" << std::endl;
//     for (int i = 0; i < int(p.t.size()); i++) {
//         tkc[0] = p.t[i];
//         X = model(&tkc);
//         ASSERT_FALSE(X.isError);

//         std::cout << p.t[i] << ",";
//         std::cout << p.XHistory[i][0] << ",";
//         std::cout << X.val[0] << ",";
//         std::cout << p.XHistory[i][1] << ",";
//         std::cout << X.val[1] << std::endl;
//     }
// }
