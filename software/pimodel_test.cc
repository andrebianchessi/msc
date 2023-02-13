#include "pimodel.h"

#include <gtest/gtest.h>

#include <vector>

#include "bounded.h"
#include "problem.h"
#include "utils.h"

const double m = 3.0;
const double tMax = 0.5;
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

    // First order model (polynomials of order 1) with time and kc
    // discretization 1.
    // i.e. in loss function we considered:
    // t = 0 and t = tMax
    // k = kMin and kMax
    // c = cMin and cMax
    //
    // The polynomials that describe each mass are:
    // x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    // x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    Pimodel* simpleModel;

    // Second order model (polynomials of order 2) with time and kc
    // discretization 2.
    // i.e. in loss function we considered:
    // t = 0, t=tMax/2 and t = tMax
    // k = kMin, k = (kMin+kMax)/2 and k = kMax
    // c = cMin, c = (cMin+cMax)/2 and c = cMax
    //
    // The polynomials (with coefficients a0, a1... omitted) that describe each
    // mass are in the form:
    // x0(t,k,c) = t^2 + tk + tc + t + k^2 + kc + k + c^2 + c + 1
    // x1(t,k,c) = t^2 + tk + tc + t + k^2 + kc + k + c^2 + c + 1
    Pimodel* secondOrderModel;

    // Same as simpleModel, but separates time into 2 buckets
    Pimodel* twoBucketsModel;

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

        this->simpleModel = new Pimodel(&this->pd, tMax, 1, 1, 1, 1);
        this->secondOrderModel = new Pimodel(&this->pd, tMax, 1, 2, 2, 2);

        this->twoBucketsModel = new Pimodel(&this->pd, tMax, 2, 1, 1, 1);

        x = 19354;
        y = 1235;
        xy = std::vector<double>{x, y};
    }

    void SetParameters() {
        // Convenient helper to set parameters of the model to increasing values
        std::vector<double> params = {1, 2, 3, 4, 5, 6, 7, 8};
        ASSERT_FALSE(this->simpleModel->SetParameters(&params).isError);

        params = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
        ASSERT_FALSE(twoBucketsModel->SetParameters(&params).isError);

        params = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
        ASSERT_FALSE(secondOrderModel->SetParameters(&params).isError);
    }
};

TEST_F(PimodelTest, ConstructorTest) {
    auto t = std::vector<double>{0, tMax};
    ASSERT_EQ(simpleModel->timeBuckets, t);
    ASSERT_EQ(secondOrderModel->timeBuckets, t);

    t = std::vector<double>{0, tMax / 2, tMax};
    ASSERT_EQ(twoBucketsModel->timeBuckets, t);

    ASSERT_EQ(simpleModel->models.size(), 1);
    ASSERT_EQ(simpleModel->models[0].size1(), 2);
    ASSERT_EQ(simpleModel->models[0].size2(), 1);
    ASSERT_EQ(simpleModel->models[0](0, 0).id, 0);
    ASSERT_EQ(simpleModel->models[0](1, 0).id, 1);
}

TEST_F(PimodelTest, TimeBucketTest) {
    ASSERT_DEATH({ simpleModel->timeBucket(-0.1); }, "");
    ASSERT_EQ(simpleModel->timeBucket(0.0), 0);
    ASSERT_EQ(simpleModel->timeBucket(tMax / 4), 0);
    ASSERT_EQ(simpleModel->timeBucket(tMax / 2), 0);
    ASSERT_EQ(simpleModel->timeBucket(0.99 * tMax), 0);
    ASSERT_EQ(simpleModel->timeBucket(tMax), 0);

    ASSERT_EQ(twoBucketsModel->timeBucket(0.0), 0);
    ASSERT_EQ(twoBucketsModel->timeBucket(tMax / 8), 0);
    ASSERT_EQ(twoBucketsModel->timeBucket(tMax / 4), 0);
    ASSERT_EQ(twoBucketsModel->timeBucket(tMax / 2), 1);
    ASSERT_EQ(twoBucketsModel->timeBucket(3 * tMax / 4), 1);
    ASSERT_EQ(twoBucketsModel->timeBucket(0.99 * tMax), 1);
    ASSERT_EQ(twoBucketsModel->timeBucket(tMax), 1);

    auto pd = ProblemDescription();
    // tFinal = 100, 10 time buckets
    // buckets = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    Pimodel m = Pimodel(&pd, 100, 10, 1, 1, 1);
    ASSERT_DEATH({ m.timeBucket(-0.1); }, "");
    ASSERT_EQ(m.timeBucket(0.0), 0);
    ASSERT_EQ(m.timeBucket(1), 0);
    ASSERT_EQ(m.timeBucket(10), 1);
    ASSERT_EQ(m.timeBucket(10.1), 1);
    ASSERT_EQ(m.timeBucket(15), 1);
    ASSERT_EQ(m.timeBucket(19.9), 1);
    ASSERT_EQ(m.timeBucket(20), 2);
    ASSERT_EQ(m.timeBucket(50), 5);
    ASSERT_EQ(m.timeBucket(55), 5);
    ASSERT_EQ(m.timeBucket(89.9), 8);
    ASSERT_EQ(m.timeBucket(90), 9);
    ASSERT_EQ(m.timeBucket(99.9), 9);
    ASSERT_EQ(m.timeBucket(100), 9);
    ASSERT_DEATH({ m.timeBucket(100.1); }, "");
}

TEST_F(PimodelTest, OperatorTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 2, 10, 10, 2);

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

    // invalid t
    tkc = {1.1, kMin, cMin};
    eval = model(&tkc);
    ASSERT_TRUE(eval.isError);
}

TEST_F(PimodelTest, nParametersTest) {
    // See PimodelTest class
    ASSERT_EQ(simpleModel->nParameters(), 8);

    // See PimodelTest class
    ASSERT_EQ(secondOrderModel->nParameters(), 20);

    // See PimodelTest class
    ASSERT_EQ(twoBucketsModel->nParameters(), 16);
}

TEST_F(PimodelTest, SetParametersTest) {
    std::vector<double> params = {1, 2, 3, 4, 5, 6, 7, 8};
    auto r = simpleModel->SetParameters(&params);
    ASSERT_FALSE(r.isError);

    // simpleModel->modelsCoefficients[timeBucket][mass][coeff]
    ASSERT_EQ(simpleModel->modelsCoefficients[0][0][0], 1);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][0][1], 2);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][0][2], 3);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][0][3], 4);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][1][0], 5);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][1][1], 6);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][1][2], 7);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][1][3], 8);

    params = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16};
    r = twoBucketsModel->SetParameters(&params);
    ASSERT_FALSE(r.isError);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][0][0], 1);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][0][1], 2);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][0][2], 3);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][0][3], 4);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][1][0], 5);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][1][1], 6);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][1][2], 7);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[0][1][3], 8);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][0][0], 9);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][0][1], 10);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][0][2], 11);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][0][3], 12);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][1][0], 13);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][1][1], 14);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][1][2], 15);
    ASSERT_EQ(twoBucketsModel->modelsCoefficients[1][1][3], 16);

    // Test error cases
    params = {1, 2, 3, 4, 5, 6, 7};
    r = simpleModel->SetParameters(&params);
    ASSERT_TRUE(r.isError);
    params = {1, 2, 3, 4, 5, 6, 7, 8, 9};
    r = simpleModel->SetParameters(&params);
    ASSERT_TRUE(r.isError);
}

TEST_F(PimodelTest, GetParametersTest) {
    std::vector<double> params = std::vector<double>(8);
    ASSERT_FALSE(simpleModel->GetParameters(&params).isError);

    ASSERT_DOUBLE_EQ(params[0], 0.0);
    ASSERT_DOUBLE_EQ(params[1], 0.0);
    ASSERT_DOUBLE_EQ(params[2], 0.0);
    ASSERT_DOUBLE_EQ(params[3], 0.0);
    ASSERT_DOUBLE_EQ(params[4], 0.0);
    ASSERT_DOUBLE_EQ(params[5], 0.0);
    ASSERT_DOUBLE_EQ(params[6], 0.0);
    ASSERT_DOUBLE_EQ(params[7], 0.0);

    this->SetParameters();

    ASSERT_FALSE(simpleModel->GetParameters(&params).isError);
    ASSERT_DOUBLE_EQ(params[0], 1.0);
    ASSERT_DOUBLE_EQ(params[1], 2.0);
    ASSERT_DOUBLE_EQ(params[2], 3.0);
    ASSERT_DOUBLE_EQ(params[3], 4.0);
    ASSERT_DOUBLE_EQ(params[4], 5.0);
    ASSERT_DOUBLE_EQ(params[5], 6.0);
    ASSERT_DOUBLE_EQ(params[6], 7.0);
    ASSERT_DOUBLE_EQ(params[7], 8.0);

    params = std::vector<double>(16);
    ASSERT_FALSE(twoBucketsModel->GetParameters(&params).isError);
    ASSERT_DOUBLE_EQ(params[0], 1.0);
    ASSERT_DOUBLE_EQ(params[1], 2.0);
    ASSERT_DOUBLE_EQ(params[2], 3.0);
    ASSERT_DOUBLE_EQ(params[3], 4.0);
    ASSERT_DOUBLE_EQ(params[4], 5.0);
    ASSERT_DOUBLE_EQ(params[5], 6.0);
    ASSERT_DOUBLE_EQ(params[6], 7.0);
    ASSERT_DOUBLE_EQ(params[7], 8.0);
    ASSERT_DOUBLE_EQ(params[8], 9.0);
    ASSERT_DOUBLE_EQ(params[9], 10.0);
    ASSERT_DOUBLE_EQ(params[10], 11.0);
    ASSERT_DOUBLE_EQ(params[11], 12.0);
    ASSERT_DOUBLE_EQ(params[12], 13.0);
    ASSERT_DOUBLE_EQ(params[13], 14.0);
    ASSERT_DOUBLE_EQ(params[14], 15.0);
    ASSERT_DOUBLE_EQ(params[15], 16.0);

    // Test error cases
    params = std::vector<double>(7);
    ASSERT_TRUE(simpleModel->GetParameters(&params).isError);
    params = std::vector<double>(9);
    ASSERT_TRUE(simpleModel->GetParameters(&params).isError);
}

TEST_F(PimodelTest, ProblemFromTkcTest) {
    double k = kMin;
    double c = cMin;
    auto tkc = std::vector<double>{1.0, k, c};
    Problem p = simpleModel->problemFromTkc(&tkc);

    ASSERT_EQ(p.masses.size(), 2);
    ASSERT_EQ(p.masses[0].m, m);
    ASSERT_EQ(p.masses[1].m, m);

    ASSERT_EQ(p.springs.size(), 1);
    ASSERT_EQ(p.springs[0].Get_k(), k);
    ASSERT_EQ(p.dampers.size(), 1);
    ASSERT_EQ(p.dampers[0].Get_c(), c);
}

TEST_F(PimodelTest, getXModelTest) {
    // x0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // x1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    this->SetParameters();

    double k = kMin;
    double c = cMin;
    double t = tMax;
    auto tkc = std::vector<double>{t, k, c};

    std::vector<double> X = simpleModel->getXModel(&tkc);

    ASSERT_EQ(X.size(), 4);
    ASSERT_EQ(X[0], 1 * t + 2 * k + 3 * c + 4 * 1);
    ASSERT_EQ(X[1], 5 * t + 6 * k + 7 * c + 8 * 1);
    ASSERT_EQ(X[2], 1);  // dx0/dt = 1
    ASSERT_EQ(X[3], 5);  // dx1/dt = 5

    t = 0;
    k = kMin;
    c = cMin;
    tkc = std::vector<double>{t, k, c};
    X = twoBucketsModel->getXModel(&tkc);
    ASSERT_EQ(X.size(), 4);
    ASSERT_EQ(X[0], 1 * t + 2 * k + 3 * c + 4 * 1);
    ASSERT_EQ(X[1], 5 * t + 6 * k + 7 * c + 8 * 1);
    ASSERT_EQ(X[2], 1);  // dx0/dt = 1
    ASSERT_EQ(X[3], 5);  // dx1/dt = 5

    t = tMax / 2;
    k = kMin;
    c = cMin;
    tkc = std::vector<double>{t, k, c};
    X = twoBucketsModel->getXModel(&tkc);
    ASSERT_EQ(X.size(), 4);
    ASSERT_EQ(X[0], 9 * t + 10 * k + 11 * c + 12 * 1);
    ASSERT_EQ(X[1], 13 * t + 14 * k + 15 * c + 16 * 1);
    ASSERT_EQ(X[2], 9);   // dx0/dt = 9
    ASSERT_EQ(X[3], 13);  // dx1/dt = 13
}

TEST_F(PimodelTest, getAccelsFromDiffEqTest) {
    this->SetParameters();
    // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    // modelXDot0(t,k,c) = 1
    // modelXDot1(t,k,c) = 5

    double t = tMax;
    double k = 923847;
    double c = 74437;
    std::vector<double> tkc = {t, k, c};

    Maybe<Problem> mProblem =
        this->pd.BuildFromVector(std::vector<double>{kMin, cMin});
    ASSERT_FALSE(mProblem.isError);
    Problem problem = mProblem.val;

    auto A = simpleModel->getAccelsFromDiffEq(&problem, tkc);
    ASSERT_EQ(A.size1(), 2);
    ASSERT_EQ(A.size2(), 1);

    // expectedAccel0(t,k,c) = 1 / m * (-kMin * x0 + kMin * x1 - cMin * dx0dt
    // + cMin * dx1dt);
    // expectedAccel0 = 1 / m * (-kMin * (1*t + 2*k + 3*c + 4*1) + kMin *(5*t
    // + 6*k + 7*c + 8*1) - cMin * 1 + cMin * 5);
    // expectedAccel0 = (4*kMin/m)*t + (4*kMin)/m*k + (4 *kMin)/m*c  +
    // (4*kMin+4*cMin)/m*1
    double expectedEval =
        1 / m *
        (-kMin * (1 * t + 2 * k + 3 * c + 4 * 1) +
         kMin * (5 * t + 6 * k + 7 * c + 8 * 1) - cMin * 1 + cMin * 5);

    std::vector<std::vector<double>> coefs = std::vector<std::vector<double>>{
        std::vector<double>{1, 2, 3, 4}, std::vector<double>{5, 6, 7, 8}};
    ASSERT_DOUBLE_EQ(A(0, 0)(coefs).val, expectedEval);
    ASSERT_DOUBLE_EQ(A(1, 0)(coefs).val, -expectedEval);

    // Test for two time buckets
    // For t=tMax, the model from twoBucketsModel is the second one
    A = twoBucketsModel->getAccelsFromDiffEq(&problem, tkc);
    ASSERT_EQ(A.size1(), 2);
    ASSERT_EQ(A.size2(), 1);
    expectedEval =
        1 / m *
        (-kMin * (9 * t + 10 * k + 11 * c + 12 * 1) +
         kMin * (13 * t + 14 * k + 15 * c + 16 * 1) - cMin * 1 + cMin * 5);
    coefs =
        std::vector<std::vector<double>>{std::vector<double>{9, 10, 11, 12},
                                         std::vector<double>{13, 14, 15, 16}};
    ASSERT_DOUBLE_EQ(A(0, 0)(coefs).val, expectedEval);
    ASSERT_DOUBLE_EQ(A(1, 0)(coefs).val, -expectedEval);
}

TEST_F(PimodelTest, getInitialXTest) {
    auto pd = ProblemDescription();
    pd.AddMass(m, 0.0, 0.0);
    pd.AddMass(m, 1.0, 0.0);
    pd.AddMass(m, 2.0, 0.0);
    pd.AddSpring(0, 1, kMin, kMax);
    pd.AddDamper(0, 1, cMin, cMax);
    pd.AddDamper(1, 2, cMin, cMax);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, 0.1);
    pd.AddInitialVel(1, 0.2);
    pd.AddInitialDisp(2, 0.3);
    pd.AddInitialVel(2, 0.4);

    auto model = Pimodel(&pd, tMax, 1, 1, 1, 1);

    std::vector<double> X = model.getInitialX();
    ASSERT_EQ(X.size(), 6);
    ASSERT_EQ(X[0], 0);
    ASSERT_EQ(X[1], 0.1);
    ASSERT_EQ(X[2], 0.3);
    ASSERT_EQ(X[3], 0);
    ASSERT_EQ(X[4], 0.2);
    ASSERT_EQ(X[5], 0.4);
}

TEST_F(PimodelTest, InitialConditionsResiduesTkcTest) {
    // Note: AddResidues is already called at the model constructor

    ASSERT_EQ(simpleModel->initialConditionsResiduesTkc.size(), 4);

    // kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][2], cMin);

    // kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][2], cMax);

    // kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][2], cMin);

    // kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][2], cMax);
}

TEST_F(PimodelTest, PhysicsResiduesTkcTest) {
    // Note: AddResidues is already called at the model constructor

    ASSERT_EQ(simpleModel->physicsResiduesTkc.size(), 8);

    // t0, kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][2], cMin);

    // t0, kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][2], cMax);

    // t0, kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][2], cMin);

    // t0, kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][0], 0);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][2], cMax);

    // tMax, kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][0], tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][2], cMin);

    // tMax, kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][0], tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][1], kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][2], cMax);

    // tMax, kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][0], tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][2], cMin);

    // tMax, kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][0], tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][1], kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][2], cMax);
}

TEST_F(PimodelTest, InitialConditionsResiduesTest) {
    // Note: AddResidues is already called at the model constructor

    // 4 total values of (k,c) per mass (see InitialConditionsResiduesTkcTest)
    // 2 masses
    // -> 8 residues in total
    ASSERT_EQ(simpleModel->initialDispResidues.size(), 8);
    ASSERT_EQ(simpleModel->initialVelResidues.size(), 8);

    std::vector<double> tkc;
    double initialX0 = 0;
    double initialX1 = initialDisplacement;

#define TEST_DISP_RESIDUES(residue, mass, initialX, expectedT, expectedK, \
                           expectedC)                                     \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].polys.size(), 1); \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].k.size(), 1);     \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].k[0], 1.0);       \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].polys[0],         \
              this->simpleModel->models[0](mass, 0));                     \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].plus, -initialX); \
    tkc = simpleModel->initialDispResidues[residue].polys[0].GetX();      \
    ASSERT_EQ(tkc.size(), 3);                                             \
    ASSERT_DOUBLE_EQ(tkc[0], expectedT);                                  \
    ASSERT_DOUBLE_EQ(tkc[1], expectedK);                                  \
    ASSERT_DOUBLE_EQ(tkc[2], expectedC);

    TEST_DISP_RESIDUES(0, 0, initialX0, 0, kMin, cMin);
    TEST_DISP_RESIDUES(1, 0, initialX0, 0, kMin, cMax);
    TEST_DISP_RESIDUES(2, 0, initialX0, 0, kMax, cMin);
    TEST_DISP_RESIDUES(3, 0, initialX0, 0, kMax, cMax);
    TEST_DISP_RESIDUES(4, 1, initialX1, 0, kMin, cMin);
    TEST_DISP_RESIDUES(5, 1, initialX1, 0, kMin, cMax);
    TEST_DISP_RESIDUES(6, 1, initialX1, 0, kMax, cMin);
    TEST_DISP_RESIDUES(7, 1, initialX1, 0, kMax, cMax);

    double initialX0Dot = 0;
    double initialX1Dot = 0;
    Poly model;
#define TEST_VEL_RESIDUES(residue, mass, initialXDot, expectedT, expectedK, \
                          expectedC)                                        \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].polys.size(), 1);    \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].k.size(), 1);        \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].k[0], 1.0);          \
    model = this->simpleModel->models[0](mass, 0);                          \
    model.Dxi(0);                                                           \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].polys[0], model);    \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].plus, -initialXDot); \
    tkc = simpleModel->initialVelResidues[residue].polys[0].GetX();         \
    ASSERT_EQ(tkc.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkc[0], expectedT);                                    \
    ASSERT_DOUBLE_EQ(tkc[1], expectedK);                                    \
    ASSERT_DOUBLE_EQ(tkc[2], expectedC);

    TEST_VEL_RESIDUES(0, 0, initialX0Dot, 0, kMin, cMin);
    TEST_VEL_RESIDUES(1, 0, initialX0Dot, 0, kMin, cMax);
    TEST_VEL_RESIDUES(2, 0, initialX0Dot, 0, kMax, cMin);
    TEST_VEL_RESIDUES(3, 0, initialX0Dot, 0, kMax, cMax);
    TEST_VEL_RESIDUES(4, 1, initialX1Dot, 0, kMin, cMin);
    TEST_VEL_RESIDUES(5, 1, initialX1Dot, 0, kMin, cMax);
    TEST_VEL_RESIDUES(6, 1, initialX1Dot, 0, kMax, cMin);
    TEST_VEL_RESIDUES(7, 1, initialX1Dot, 0, kMax, cMax);
}

TEST_F(PimodelTest, PhysicsResiduesTest) {
    // Note: AddResidues is already called at the model constructor

    // 8 total values of (t, k,c) per mass (see
    // InitialConditionsResiduesTkcTest) 2 masses
    // -> 16 residues in total
    ASSERT_EQ(simpleModel->physicsResidues.size(), 16);

    std::vector<double> tkc;

    Poly x0, x1, x0Dot, x1Dot, x0DotDot, x1DotDot;

#define TEST_MASS_0_PHYSICS_RESIDUE(residue, kValue, cValue, tValue)     \
    ASSERT_DOUBLE_EQ(simpleModel->physicsResidues[residue].plus, -0);    \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys.size(), 1);    \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k.size(), 1);        \
                                                                         \
    /* x0DotDot */                                                       \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[0], 1.0);          \
    x0DotDot = this->simpleModel->models[0](0, 0);                       \
    x0DotDot.Dxi(0);                                                     \
    x0DotDot.Dxi(0);                                                     \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[0], x0DotDot); \
                                                                         \
    tkc = simpleModel->physicsResidues[residue].polys[0].GetX();         \
    ASSERT_EQ(tkc.size(), 3);                                            \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                    \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                    \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);

    TEST_MASS_0_PHYSICS_RESIDUE(0, kMin, cMin, 0);
    TEST_MASS_0_PHYSICS_RESIDUE(1, kMin, cMax, 0);
    TEST_MASS_0_PHYSICS_RESIDUE(2, kMax, cMin, 0);
    TEST_MASS_0_PHYSICS_RESIDUE(3, kMax, cMax, 0);
    TEST_MASS_0_PHYSICS_RESIDUE(4, kMin, cMin, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(5, kMin, cMax, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(6, kMax, cMin, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(7, kMax, cMax, tMax);

#define TEST_MASS_1_PHYSICS_RESIDUE(residue, kValue, cValue, tValue)           \
    ASSERT_DOUBLE_EQ(simpleModel->physicsResidues[residue].plus, 0);           \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys.size(), 5);          \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k.size(), 5);              \
                                                                               \
    /* x1DotDot */                                                             \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[0], 1.0);                \
    x1DotDot = this->simpleModel->models[0](1, 0);                             \
    x1DotDot.Dxi(0);                                                           \
    x1DotDot.Dxi(0);                                                           \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[0], x1DotDot);       \
                                                                               \
    /* -1/m*(k*x0) */                                                          \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[1], -1 / m * (kValue));  \
    x0 = this->simpleModel->models[0](0, 0);                                   \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[1], x0);             \
                                                                               \
    /* -1/m*(-k*x1) */                                                         \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[2], -1 / m * (-kValue)); \
    x1 = this->simpleModel->models[0](1, 0);                                   \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[2], x1);             \
                                                                               \
    /* -1/m*(c*x0Dot) */                                                       \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[3], -1 / m * (cValue));  \
    x0Dot = this->simpleModel->models[0](0, 0);                                \
    x0Dot.Dxi(0);                                                              \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[3], x0Dot);          \
                                                                               \
    /* -1/m*(-c*x1Dot) */                                                      \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[4], -1 / m * (-cValue)); \
    x1Dot = this->simpleModel->models[0](1, 0);                                \
    x1Dot.Dxi(0);                                                              \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[4], x1Dot);          \
                                                                               \
    tkc = simpleModel->physicsResidues[residue].polys[0].GetX();               \
    ASSERT_EQ(tkc.size(), 3);                                                  \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);                                          \
    tkc = simpleModel->physicsResidues[residue].polys[1].GetX();               \
    ASSERT_EQ(tkc.size(), 3);                                                  \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);                                          \
    tkc = simpleModel->physicsResidues[residue].polys[2].GetX();               \
    ASSERT_EQ(tkc.size(), 3);                                                  \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);                                          \
    tkc = simpleModel->physicsResidues[residue].polys[3].GetX();               \
    ASSERT_EQ(tkc.size(), 3);                                                  \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);                                          \
    tkc = simpleModel->physicsResidues[residue].polys[4].GetX();               \
    ASSERT_EQ(tkc.size(), 3);                                                  \
    ASSERT_DOUBLE_EQ(tkc[0], tValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[1], kValue);                                          \
    ASSERT_DOUBLE_EQ(tkc[2], cValue);

    TEST_MASS_1_PHYSICS_RESIDUE(8, kMin, cMin, 0);
    TEST_MASS_1_PHYSICS_RESIDUE(9, kMin, cMax, 0);
    TEST_MASS_1_PHYSICS_RESIDUE(10, kMax, cMin, 0);
    TEST_MASS_1_PHYSICS_RESIDUE(11, kMax, cMax, 0);
    TEST_MASS_1_PHYSICS_RESIDUE(12, kMin, cMin, tMax);
    TEST_MASS_1_PHYSICS_RESIDUE(13, kMin, cMax, tMax);
    TEST_MASS_1_PHYSICS_RESIDUE(14, kMax, cMin, tMax);
    TEST_MASS_1_PHYSICS_RESIDUE(15, kMax, cMax, tMax);
}

TEST_F(PimodelTest, LossTest) {
    // Full loss test
    std::vector<double> params = std::vector<double>(8);
    params = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
              11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
    ASSERT_FALSE(secondOrderModel->SetParameters(&params).isError);
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
    for (double t : std::vector<double>{0, tMax / 2, tMax}) {
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

                // Mass 0 is fixed
                double x0DotDot = 0;
                double x1DotDot = 1 / m *
                                  (k * modelX0 - k * modelX1 + c * modelX0Dot -
                                   c * modelX1Dot);
                expectedLoss += pow(modelX0DotDot - x0DotDot, 2);
                expectedLoss += pow(modelX1DotDot - x1DotDot, 2);
            }
        }
    }

    double loss = secondOrderModel->Loss();

    ASSERT_DOUBLE_EQ(expectedLoss, loss);
}

TEST_F(PimodelTest, LossGradientTest) {
    Pimodel model = Pimodel(&this->pd, tMax, 1, 1, 1, 1);
    std::vector<double> params = std::vector<double>(8);
    double a0 = 13.0;
    double a1 = 26.0;
    double a2 = 332.0;
    double a3 = 44.0;
    double a4 = 56.0;
    double a5 = 61.0;
    double a6 = 722.0;
    double a7 = 58.0;
    params = {a0, a1, a2, a3, a4, a5, a6, a7};
    ASSERT_FALSE(model.SetParameters(&params).isError);
    // modelX0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    // modelX1(t,k,c) = a4*t + a5*k + a6*c + a7*1
    // modelXDot0(t,k,c) = a0
    // modelXDot1(t,k,c) = a4
    // modelXDotDot0(t,k,c) = 0
    // modelXDotDot1(t,k,c) = 0

    std::vector<double> expectedGrad = std::vector<double>(8);

    // Initial values for the problem created in SetUp()
    double initialX0 = 0;
    double initialX1 = initialDisplacement;
    double initialX0Dot = 0;
    double initialX1Dot = 0;

    double t = 0;
    // Initial conditions loss:
    for (double k : std::vector<double>{kMin, kMax}) {
        for (double c : std::vector<double>{cMin, cMax}) {
            double modelX0 = a0 * t + a1 * k + a2 * c + a3 * 1;
            double modelX1 = a4 * t + a5 * k + a6 * c + a7 * 1;
            double modelX0Dot = a0 * 1;
            double modelX1Dot = a4 * 1;

            double d_da0_modelX0 = t;
            double d_da1_modelX0 = k;
            double d_da2_modelX0 = c;
            double d_da3_modelX0 = 1;
            double d_da0_modelX0Dot = 1;
            double d_da1_modelX0Dot = 0;
            double d_da2_modelX0Dot = 0;
            double d_da3_modelX0Dot = 0;

            double d_da4_modelX1 = t;
            double d_da5_modelX1 = k;
            double d_da6_modelX1 = c;
            double d_da7_modelX1 = 1;
            double d_da4_modelX1Dot = 1;
            double d_da5_modelX1Dot = 0;
            double d_da6_modelX1Dot = 0;
            double d_da7_modelX1Dot = 0;

            expectedGrad[0] += (modelX0 - initialX0) * d_da0_modelX0;
            expectedGrad[0] += (modelX0Dot - initialX0Dot) * d_da0_modelX0Dot;

            expectedGrad[1] += (modelX0 - initialX0) * d_da1_modelX0;
            expectedGrad[1] += (modelX0Dot - initialX0Dot) * d_da1_modelX0Dot;

            expectedGrad[2] += (modelX0 - initialX0) * d_da2_modelX0;
            expectedGrad[2] += (modelX0Dot - initialX0Dot) * d_da2_modelX0Dot;

            expectedGrad[3] += (modelX0 - initialX0) * d_da3_modelX0;
            expectedGrad[3] += (modelX0Dot - initialX0Dot) * d_da3_modelX0Dot;

            expectedGrad[4] += (modelX1 - initialX1) * d_da4_modelX1;
            expectedGrad[4] += (modelX1Dot - initialX1Dot) * d_da4_modelX1Dot;

            expectedGrad[5] += (modelX1 - initialX1) * d_da5_modelX1;
            expectedGrad[5] += (modelX1Dot - initialX1Dot) * d_da5_modelX1Dot;

            expectedGrad[6] += (modelX1 - initialX1) * d_da6_modelX1;
            expectedGrad[6] += (modelX1Dot - initialX1Dot) * d_da6_modelX1Dot;

            expectedGrad[7] += (modelX1 - initialX1) * d_da7_modelX1;
            expectedGrad[7] += (modelX1Dot - initialX1Dot) * d_da7_modelX1Dot;
        }
    }

    // Physics loss:
    for (double t : std::vector<double>{0, tMax}) {
        for (double k : std::vector<double>{kMin, kMax}) {
            for (double c : std::vector<double>{cMin, cMax}) {
                double modelX0 = a0 * t + a1 * k + a2 * c + a3 * 1;
                double modelX1 = a4 * t + a5 * k + a6 * c + a7 * 1;
                double modelX0Dot = a0;
                double modelX1Dot = a4;
                double modelX0DotDot = 0;
                double modelX1DotDot = 0;

                // Mass 0 is fixed
                double X0DotDot = 0;
                double X1DotDot = 1 / m *
                                  (k * modelX0 - k * modelX1 + c * modelX0Dot -
                                   c * modelX1Dot);

                double d_da0_modelX0 = t;
                double d_da1_modelX0 = k;
                double d_da2_modelX0 = c;
                double d_da3_modelX0 = 1;
                double d_da4_modelX0 = 0;
                double d_da5_modelX0 = 0;
                double d_da6_modelX0 = 0;
                double d_da7_modelX0 = 0;
                double d_da0_modelX1 = 0;
                double d_da1_modelX1 = 0;
                double d_da2_modelX1 = 0;
                double d_da3_modelX1 = 0;
                double d_da4_modelX1 = t;
                double d_da5_modelX1 = k;
                double d_da6_modelX1 = c;
                double d_da7_modelX1 = 1;

                double d_da0_modelX0Dot = 1;
                double d_da1_modelX0Dot = 0;
                double d_da2_modelX0Dot = 0;
                double d_da3_modelX0Dot = 0;
                double d_da4_modelX0Dot = 0;
                double d_da5_modelX0Dot = 0;
                double d_da6_modelX0Dot = 0;
                double d_da7_modelX0Dot = 0;
                double d_da0_modelX1Dot = 0;
                double d_da1_modelX1Dot = 0;
                double d_da2_modelX1Dot = 0;
                double d_da3_modelX1Dot = 0;
                double d_da4_modelX1Dot = 1;
                double d_da5_modelX1Dot = 0;
                double d_da6_modelX1Dot = 0;
                double d_da7_modelX1Dot = 0;

                // modelX0DotDot = modelX1DotDot = 0
                // -> All derivatives are 0
                double d_dai_modelX0DotDot = 0;
                double d_dai_modelX1DotDot = 0;

                double d_da0_X0DotDot = 0;
                double d_da0_X1DotDot =
                    1 / m *
                    (k * d_da0_modelX0 - k * d_da0_modelX1 +
                     c * d_da0_modelX0Dot - c * d_da0_modelX1Dot);
                double d_da1_X0DotDot = 0;
                double d_da1_X1DotDot =
                    1 / m *
                    (k * d_da1_modelX0 - k * d_da1_modelX1 +
                     c * d_da1_modelX0Dot - c * d_da1_modelX1Dot);
                double d_da2_X0DotDot = 0;
                double d_da2_X1DotDot =
                    1 / m *
                    (k * d_da2_modelX0 - k * d_da2_modelX1 +
                     c * d_da2_modelX0Dot - c * d_da2_modelX1Dot);
                double d_da3_X0DotDot = 0;
                double d_da3_X1DotDot =
                    1 / m *
                    (k * d_da3_modelX0 - k * d_da3_modelX1 +
                     c * d_da3_modelX0Dot - c * d_da3_modelX1Dot);
                double d_da4_X0DotDot = 0;
                double d_da4_X1DotDot =
                    1 / m *
                    (k * d_da4_modelX0 - k * d_da4_modelX1 +
                     c * d_da4_modelX0Dot - c * d_da4_modelX1Dot);
                double d_da5_X0DotDot = 0;
                double d_da5_X1DotDot =
                    1 / m *
                    (k * d_da5_modelX0 - k * d_da5_modelX1 +
                     c * d_da5_modelX0Dot - c * d_da5_modelX1Dot);
                double d_da6_X0DotDot = 0;
                double d_da6_X1DotDot =
                    1 / m *
                    (k * d_da6_modelX0 - k * d_da6_modelX1 +
                     c * d_da6_modelX0Dot - c * d_da6_modelX1Dot);
                double d_da7_X0DotDot = 0;
                double d_da7_X1DotDot =
                    1 / m *
                    (k * d_da7_modelX0 - k * d_da7_modelX1 +
                     c * d_da7_modelX0Dot - c * d_da7_modelX1Dot);

                expectedGrad[0] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da0_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da0_X1DotDot);
                expectedGrad[1] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da1_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da1_X1DotDot);
                expectedGrad[2] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da2_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da2_X1DotDot);
                expectedGrad[3] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da3_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da3_X1DotDot);
                expectedGrad[4] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da4_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da4_X1DotDot);
                expectedGrad[5] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da5_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da5_X1DotDot);
                expectedGrad[6] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da6_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da6_X1DotDot);
                expectedGrad[7] += (modelX0DotDot - X0DotDot) *
                                       (d_dai_modelX0DotDot - d_da7_X0DotDot) +
                                   (modelX1DotDot - X1DotDot) *
                                       (d_dai_modelX1DotDot - d_da7_X1DotDot);
            }
        }
    }

    std::vector<double> grad = model.LossGradient();

    ASSERT_EQ(grad.size(), 8);
    EXPECT_DOUBLE_EQ(grad[0], expectedGrad[0]);
    EXPECT_DOUBLE_EQ(grad[1], expectedGrad[1]);
    EXPECT_DOUBLE_EQ(grad[2], expectedGrad[2]);
    EXPECT_DOUBLE_EQ(grad[3], expectedGrad[3]);
    EXPECT_DOUBLE_EQ(grad[4], expectedGrad[4]);
    EXPECT_DOUBLE_EQ(grad[5], expectedGrad[5]);
    EXPECT_DOUBLE_EQ(grad[6], expectedGrad[6]);
    EXPECT_DOUBLE_EQ(grad[7], expectedGrad[7]);
}

TEST(PimodelTrainingTest, TrainTest) {
    double mass = 1;
    double kMin = 0.8;
    double kMax = 1.0;
    double cMin = 0.0;
    double cMax = 0.2;
    double tMax = 2.0;
    double initialDisp = 1.0;

    auto pd = ProblemDescription();
    pd.AddMass(mass, 0.0, 0.0);
    pd.AddMass(mass, 1.0, 0.0);
    pd.AddSpring(0, 1, kMin, kMax);
    pd.AddDamper(0, 1, cMin, cMax);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, initialDisp);

    int timeBuckets = 3;
    int timeDiscretization = 3;
    int kcDiscretization = 1;
    int order = 2;
    double learningRate = 0.0001;
    int maxSteps = 10;
    bool log = true;

    // Train model
    Pimodel model = Pimodel(&pd, tMax, timeBuckets, timeDiscretization,
                            kcDiscretization, order);

    double initialLoss = model.Loss();
    model.Train(learningRate, maxSteps, log);
    ASSERT_TRUE(model.Loss() < initialLoss);

    // Get problem using intermediate value for k and c, and integrate it.
    // Then, compare the model's prediction with the problem's result.
    double k = (kMin + kMax) / 2;
    double c = (cMin + cMax) / 2;
    Maybe<Problem> mP = pd.BuildFromVector(std::vector<double>{k, c});
    ASSERT_FALSE(mP.isError);
    Problem p = mP.val;
    p.Integrate(tMax);

    std::vector<double> tkc = std::vector<double>{0.0, k, c};
    Maybe<std::vector<double>> X;
    std::cout << "t,x0,modelX0,x1,modelX1" << std::endl;
    for (int i = 0; i < int(p.t.size()); i += 1) {
        tkc[0] = p.t[i];
        X = model(&tkc);
        ASSERT_FALSE(X.isError);

        std::cout << p.t[i] << ",";
        std::cout << p.XHistory[i][0] << ",";
        std::cout << X.val[0] << ",";
        std::cout << p.XHistory[i][1] << ",";
        std::cout << X.val[1] << std::endl;
    }
}
