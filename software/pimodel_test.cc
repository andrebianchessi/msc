#include "pimodel.h"

#include <gtest/gtest.h>

#include <vector>

#include "bounded.h"
#include "problem.h"
#include "utils.h"

const double m = 3.0;
const double TMin = 0.0;
const double TMax = 0.5;
const double KMin = 1.0;
const double KMax = 2.0;
const double CMin = 10.0;
const double CMax = 20.0;
const double initialDisplacement = 10.0;

class PimodelTest : public testing::Test {
   public:
    ProblemDescription pd;

    // Normalized variables that are used to evaluate the models
    Bounded tB;
    Bounded kB;
    Bounded cB;
    double t;
    double k;
    double c;
    std::vector<Bounded> tkc;
    std::vector<double> tkc_;

    // Auxiliary variables
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

    void SetUp() {
        // Called before every TEST_F
        // Create problem description of two masses,
        // connected with a spring and damper. One is fixed
        // and the other has initial displacement.
        this->pd = ProblemDescription();
        this->pd.AddMass(m, 0.0, 0.0);
        this->pd.AddMass(m, 1.0, 0.0);
        this->pd.AddSpring(0, 1, KMin, KMax);
        this->pd.AddDamper(0, 1, CMin, CMax);
        this->pd.SetFixedMass(0);
        this->pd.AddInitialDisp(1, initialDisplacement);
        std::vector<Bounded> dna = std::vector<Bounded>(2);
        auto e0 = this->pd.BuildFromDNA(dna);
        ASSERT_FALSE(e0.isError);

        this->simpleModel = new Pimodel(this->pd, TMin, TMax, 1, 1, 1);
        this->simpleModel->AddResidues();

        this->secondOrderModel = new Pimodel(this->pd, TMin, TMax, 2, 2, 2);
        this->secondOrderModel->AddResidues();

        tB = RandomB();
        kB = RandomB();
        cB = RandomB();
        t = tB.Get();
        k = kB.Get();
        c = cB.Get();
        tkc = std::vector<Bounded>{tB, kB, cB};
        tkc_ = std::vector<double>(3);
        tkc_[0] = Unnormalize(tB, TMin, TMax);
        tkc_[1] = Unnormalize(kB, KMin, KMax);
        tkc_[2] = Unnormalize(cB, CMin, CMax);

        x = Random(10, 1000);
        y = Random(10, 1000);
        xy = std::vector<double>{x, y};
    }

    void SetParameters() {
        // Convenient helper to set parameters of the model to increasing values
        std::vector<double> params = {1, 2, 3, 4, 5, 6, 7, 8};
        ASSERT_FALSE(this->simpleModel->SetParameters(&params).isError);

        params = {1,  2,  3,  4,  5,  6,  7,  8,  9,  10,
                  11, 12, 13, 14, 15, 16, 17, 18, 19, 20};
        ASSERT_FALSE(secondOrderModel->SetParameters(&params).isError);
    }
};

TEST_F(PimodelTest, ConstructorTest) {
    ASSERT_EQ(simpleModel->models.size1(), 2);
    ASSERT_EQ(simpleModel->models.size2(), 1);
    ASSERT_EQ(simpleModel->models(0, 0).id, 0);
    ASSERT_EQ(simpleModel->models(1, 0).id, 1);

    ASSERT_EQ(simpleModel->modelsD.size(), 2);
    ASSERT_EQ(simpleModel->modelsD[0].id, 0);
    ASSERT_EQ(simpleModel->modelsD[1].id, 1);

    ASSERT_EQ(simpleModel->modelsDD.size(), 2);
    ASSERT_EQ(simpleModel->modelsDD[0].id, 0);
    ASSERT_EQ(simpleModel->modelsDD[1].id, 1);
}

TEST_F(PimodelTest, OperatorTest) {
    Pimodel& model = (*this->secondOrderModel);

    // x0(t,k,c) = a0*t^2 + a1*tk + a2*tc + a3*t + a4*k^2 + a5*kc + a6*k +
    //             a7*c^2 + a8*c + a9*1
    // dx0Dt(t,k,c) = 2*a0*t + a1*k + a2*c + a3
    // x1(t,k,c) = b0*t^2 + b1*tk + b2*tc + b3*t + b4*k^2 +
    //             b5*kc + b6*k + b7*c^2 + b8*c + b9*1
    // dx1Dt(t,k,c) = 2*b0*t + b1*k + b2*c + b3
    std::vector<double> a =
        std::vector<double>{Random(), Random(), Random(), Random(), Random(),
                            Random(), Random(), Random(), Random(), Random()};
    std::vector<double> b =
        std::vector<double>{Random(), Random(), Random(), Random(), Random(),
                            Random(), Random(), Random(), Random(), Random()};
    std::vector<double> parameters = std::vector<double>{
        a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8], a[9],
        b[0], b[1], b[2], b[3], b[4], b[5], b[6], b[7], b[8], b[9]};
    model.SetParameters(&parameters);

    auto eval = model(&tkc_);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val.size(), 2);
    ASSERT_DOUBLE_EQ(eval.val[0], a[0] * t * t + a[1] * t * k + a[2] * t * c +
                                      a[3] * t + a[4] * k * k + a[5] * k * c +
                                      a[6] * k + a[7] * c * c + a[8] * c +
                                      a[9] * 1);
    ASSERT_DOUBLE_EQ(eval.val[1], b[0] * t * t + b[1] * t * k + b[2] * t * c +
                                      b[3] * t + b[4] * k * k + b[5] * k * c +
                                      b[6] * k + b[7] * c * c + b[8] * c +
                                      b[9] * 1);
    eval = model.GetVelocities(&tkc_);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val.size(), 2);
    ASSERT_DOUBLE_EQ(eval.val[0], 2 * a[0] * t + a[1] * k + a[2] * c + a[3]);
    ASSERT_DOUBLE_EQ(eval.val[1], 2 * b[0] * t + b[1] * k + b[2] * c + b[3]);

    // Test error cases
    // tkc too large
    tkc_ = std::vector<double>{1.1, 1.0, 2.0, 3.0};
    eval = model(&tkc_);
    ASSERT_TRUE(eval.isError);

    // tkc too small
    tkc_ = std::vector<double>{1.1, 1.0};
    eval = model(&tkc_);
    ASSERT_TRUE(eval.isError);
}

TEST_F(PimodelTest, nParametersTest) {
    // See PimodelTest class
    ASSERT_EQ(simpleModel->nParameters(), 8);

    // See PimodelTest class
    ASSERT_EQ(secondOrderModel->nParameters(), 20);
}

TEST_F(PimodelTest, SetParametersTest) {
    std::vector<double> params = {1, 2, 3, 4, 5, 6, 7, 8};
    auto r = simpleModel->SetParameters(&params);
    ASSERT_FALSE(r.isError);

    // simpleModel->modelsCoefficients[timeBucket][mass][coeff]
    ASSERT_EQ(simpleModel->modelsCoefficients[0][0], 1);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][1], 2);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][2], 3);
    ASSERT_EQ(simpleModel->modelsCoefficients[0][3], 4);
    ASSERT_EQ(simpleModel->modelsCoefficients[1][0], 5);
    ASSERT_EQ(simpleModel->modelsCoefficients[1][1], 6);
    ASSERT_EQ(simpleModel->modelsCoefficients[1][2], 7);
    ASSERT_EQ(simpleModel->modelsCoefficients[1][3], 8);

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

    // Test error cases
    params = std::vector<double>(7);
    ASSERT_TRUE(simpleModel->GetParameters(&params).isError);
    params = std::vector<double>(9);
    ASSERT_TRUE(simpleModel->GetParameters(&params).isError);
}

TEST_F(PimodelTest, ProblemFromTkcTest) {
    Problem p = simpleModel->problemFromTkc(&tkc);

    ASSERT_EQ(p.masses.size(), 2);
    ASSERT_EQ(p.masses[0].m, m);
    ASSERT_EQ(p.masses[1].m, m);

    ASSERT_EQ(p.springs.size(), 1);
    ASSERT_EQ(p.springs[0].Get_k(), Unnormalize(tkc[1], KMin, KMax));
    ASSERT_EQ(p.dampers.size(), 1);
    ASSERT_EQ(p.dampers[0].Get_c(), Unnormalize(tkc[2], CMin, CMax));
}

TEST_F(PimodelTest, getXModelTest) {
    // x0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // x1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    this->SetParameters();

    std::vector<double> X = simpleModel->getXModel(&tkc);

    ASSERT_EQ(X.size(), 4);
    ASSERT_EQ(X[0], 1 * t + 2 * k + 3 * c + 4 * 1);
    ASSERT_EQ(X[1], 5 * t + 6 * k + 7 * c + 8 * 1);
    ASSERT_EQ(X[2], 1);  // dx0/dt = 1
    ASSERT_EQ(X[3], 5);  // dx1/dt = 5
}

TEST_F(PimodelTest, getAccelsFromDiffEqTest) {
    this->SetParameters();
    // modelX0(t,k,c) = 1*t + 2*k + 3*c + 4*1
    // modelX1(t,k,c) = 5*t + 6*k + 7*c + 8*1
    // modelXDot0(t,k,c) = 1
    // modelXDot1(t,k,c) = 5

    Maybe<Problem> mProblem =
        this->pd.BuildFromVector(std::vector<double>{KMin, CMin});
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
        (-KMin * (1 * t + 2 * k + 3 * c + 4 * 1) +
         KMin * (5 * t + 6 * k + 7 * c + 8 * 1) - CMin * 1 + CMin * 5);

    std::vector<std::vector<double>> coefs = std::vector<std::vector<double>>{
        std::vector<double>{1, 2, 3, 4}, std::vector<double>{5, 6, 7, 8}};
    ASSERT_DOUBLE_EQ(A(0, 0)(coefs).val, expectedEval);
    ASSERT_DOUBLE_EQ(A(1, 0)(coefs).val, -expectedEval);
}

TEST_F(PimodelTest, getInitialXTest) {
    auto pd = ProblemDescription();
    pd.AddMass(m, 0.0, 0.0);
    pd.AddMass(m, 1.0, 0.0);
    pd.AddMass(m, 2.0, 0.0);
    pd.AddSpring(0, 1, KMin, KMax);
    pd.AddDamper(0, 1, CMin, CMax);
    pd.AddDamper(1, 2, CMin, CMax);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, 0.1);
    pd.AddInitialVel(1, 0.2);
    pd.AddInitialDisp(2, 0.3);
    pd.AddInitialVel(2, 0.4);

    auto model = Pimodel(pd, TMin, TMax, 1, 1, 1);

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

    double tMin = 0, kMin = 0, cMin = 0;
    double kMax = 1, cMax = 1;

    // kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][0].Get(),
                     tMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][1].Get(),
                     kMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[0][2].Get(),
                     cMin);

    // kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][0].Get(),
                     tMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][1].Get(),
                     kMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[1][2].Get(),
                     cMax);

    // kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][0].Get(),
                     tMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][1].Get(),
                     kMax);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[2][2].Get(),
                     cMin);

    // kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][0].Get(),
                     tMin);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][1].Get(),
                     kMax);
    ASSERT_DOUBLE_EQ(simpleModel->initialConditionsResiduesTkc[3][2].Get(),
                     cMax);
}

TEST_F(PimodelTest, PhysicsResiduesTkcTest) {
    // Note: AddResidues is already called at the model constructor

    ASSERT_EQ(simpleModel->physicsResiduesTkc.size(), 8);

    double tMin = 0, kMin = 0, cMin = 0;
    double tMax = 1, kMax = 1, cMax = 1;

    // t0, kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][0].Get(), tMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][1].Get(), kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[0][2].Get(), cMin);

    // t0, kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][0].Get(), tMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][1].Get(), kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[1][2].Get(), cMax);

    // t0, kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][0].Get(), tMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][1].Get(), kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[2][2].Get(), cMin);

    // t0, kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][0].Get(), tMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][1].Get(), kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[3][2].Get(), cMax);

    // tMax, kMin, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][0].Get(), tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][1].Get(), kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[4][2].Get(), cMin);

    // tMax, kMin, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][0].Get(), tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][1].Get(), kMin);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[5][2].Get(), cMax);

    // tMax, kMax, cMin
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][0].Get(), tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][1].Get(), kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[6][2].Get(), cMin);

    // tMax, kMax, cMax
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7].size(), 3);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][0].Get(), tMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][1].Get(), kMax);
    ASSERT_DOUBLE_EQ(simpleModel->physicsResiduesTkc[7][2].Get(), cMax);
}

TEST_F(PimodelTest, InitialConditionsResiduesTest) {
    // Note: AddResidues is already called at the model constructor

    // 4 total values of (k,c) per mass (see InitialConditionsResiduesTkcTest)
    // 2 masses
    // -> 8 residues in total
    ASSERT_EQ(simpleModel->initialDispResidues.size(), 8);
    ASSERT_EQ(simpleModel->initialVelResidues.size(), 8);

    double tMin = 0, kMin = 0, cMin = 0;
    double kMax = 1, cMax = 1;

    std::vector<double> tkcSet;
    double initialX0 = 0;
    double initialX1 = initialDisplacement;

#define TEST_DISP_RESIDUES(residue, mass, initialX, expectedT, expectedK, \
                           expectedC)                                     \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].polys.size(), 1); \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].k.size(), 1);     \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].k[0], 1.0);       \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].polys[0],         \
              this->simpleModel->models(mass, 0));                        \
    ASSERT_EQ(simpleModel->initialDispResidues[residue].plus, -initialX); \
    tkcSet = simpleModel->initialDispResidues[residue].polys[0].GetX();   \
    ASSERT_EQ(tkcSet.size(), 3);                                          \
    ASSERT_DOUBLE_EQ(tkcSet[0], expectedT);                               \
    ASSERT_DOUBLE_EQ(tkcSet[1], expectedK);                               \
    ASSERT_DOUBLE_EQ(tkcSet[2], expectedC);

    TEST_DISP_RESIDUES(0, 0, initialX0, tMin, kMin, cMin);
    TEST_DISP_RESIDUES(1, 0, initialX0, tMin, kMin, cMax);
    TEST_DISP_RESIDUES(2, 0, initialX0, tMin, kMax, cMin);
    TEST_DISP_RESIDUES(3, 0, initialX0, tMin, kMax, cMax);
    TEST_DISP_RESIDUES(4, 1, initialX1, tMin, kMin, cMin);
    TEST_DISP_RESIDUES(5, 1, initialX1, tMin, kMin, cMax);
    TEST_DISP_RESIDUES(6, 1, initialX1, tMin, kMax, cMin);
    TEST_DISP_RESIDUES(7, 1, initialX1, tMin, kMax, cMax);

    double initialX0Dot = 0;
    double initialX1Dot = 0;
    Poly model;
#define TEST_VEL_RESIDUES(residue, mass, initialXDot, expectedT, expectedK, \
                          expectedC)                                        \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].polys.size(), 1);    \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].k.size(), 1);        \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].k[0], 1.0);          \
    model = this->simpleModel->models(mass, 0);                             \
    model.Dxi(0);                                                           \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].polys[0], model);    \
    ASSERT_EQ(simpleModel->initialVelResidues[residue].plus, -initialXDot); \
    tkcSet = simpleModel->initialVelResidues[residue].polys[0].GetX();      \
    ASSERT_EQ(tkcSet.size(), 3);                                            \
    ASSERT_DOUBLE_EQ(tkcSet[0], expectedT);                                 \
    ASSERT_DOUBLE_EQ(tkcSet[1], expectedK);                                 \
    ASSERT_DOUBLE_EQ(tkcSet[2], expectedC);

    TEST_VEL_RESIDUES(0, 0, initialX0Dot, tMin, kMin, cMin);
    TEST_VEL_RESIDUES(1, 0, initialX0Dot, tMin, kMin, cMax);
    TEST_VEL_RESIDUES(2, 0, initialX0Dot, tMin, kMax, cMin);
    TEST_VEL_RESIDUES(3, 0, initialX0Dot, tMin, kMax, cMax);
    TEST_VEL_RESIDUES(4, 1, initialX1Dot, tMin, kMin, cMin);
    TEST_VEL_RESIDUES(5, 1, initialX1Dot, tMin, kMin, cMax);
    TEST_VEL_RESIDUES(6, 1, initialX1Dot, tMin, kMax, cMin);
    TEST_VEL_RESIDUES(7, 1, initialX1Dot, tMin, kMax, cMax);
}

TEST_F(PimodelTest, PhysicsResiduesTest) {
    // Note: AddResidues is already called at the model constructor

    // 8 total values of (t, k,c) per mass (see
    // InitialConditionsResiduesTkcTest) 2 masses
    // -> 16 residues in total
    ASSERT_EQ(simpleModel->physicsResidues.size(), 16);

    std::vector<double> tkcSet;
    double tMin = 0, kMin = 0, cMin = 0;
    double tMax = 1, kMax = 1, cMax = 1;

    Poly x0, x1, x0Dot, x1Dot, x0DotDot, x1DotDot;

#define TEST_MASS_0_PHYSICS_RESIDUE(residue, kValue, cValue, tValue)     \
    ASSERT_DOUBLE_EQ(simpleModel->physicsResidues[residue].plus, -0);    \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys.size(), 1);    \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k.size(), 1);        \
                                                                         \
    /* x0DotDot */                                                       \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[0], 1.0);          \
    x0DotDot = this->simpleModel->models(0, 0);                          \
    x0DotDot.Dxi(0);                                                     \
    x0DotDot.Dxi(0);                                                     \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[0], x0DotDot); \
                                                                         \
    tkcSet = simpleModel->physicsResidues[residue].polys[0].GetX();      \
    ASSERT_EQ(tkcSet.size(), 3);                                         \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                 \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                 \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);

    TEST_MASS_0_PHYSICS_RESIDUE(0, kMin, cMin, tMin);
    TEST_MASS_0_PHYSICS_RESIDUE(1, kMin, cMax, tMin);
    TEST_MASS_0_PHYSICS_RESIDUE(2, kMax, cMin, tMin);
    TEST_MASS_0_PHYSICS_RESIDUE(3, kMax, cMax, tMin);
    TEST_MASS_0_PHYSICS_RESIDUE(4, kMin, cMin, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(5, kMin, cMax, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(6, kMax, cMin, tMax);
    TEST_MASS_0_PHYSICS_RESIDUE(7, kMax, cMax, tMax);

#define TEST_MASS_1_PHYSICS_RESIDUE(residue, kValue, cValue, tValue, KValue,   \
                                    CValue, TValue)                            \
    ASSERT_DOUBLE_EQ(simpleModel->physicsResidues[residue].plus, 0);           \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys.size(), 5);          \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k.size(), 5);              \
                                                                               \
    /* x1DotDot */                                                             \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[0], 1.0);                \
    x1DotDot = this->simpleModel->models(1, 0);                                \
    x1DotDot.Dxi(0);                                                           \
    x1DotDot.Dxi(0);                                                           \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[0], x1DotDot);       \
                                                                               \
    /* -1/m*(k*x0) */                                                          \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[1], -1 / m * (KValue));  \
    x0 = this->simpleModel->models(0, 0);                                      \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[1], x0);             \
                                                                               \
    /* -1/m*(-k*x1) */                                                         \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[2], -1 / m * (-KValue)); \
    x1 = this->simpleModel->models(1, 0);                                      \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[2], x1);             \
                                                                               \
    /* -1/m*(c*x0Dot) */                                                       \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[3], -1 / m * (CValue));  \
    x0Dot = this->simpleModel->models(0, 0);                                   \
    x0Dot.Dxi(0);                                                              \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[3], x0Dot);          \
                                                                               \
    /* -1/m*(-c*x1Dot) */                                                      \
    ASSERT_EQ(simpleModel->physicsResidues[residue].k[4], -1 / m * (-CValue)); \
    x1Dot = this->simpleModel->models(1, 0);                                   \
    x1Dot.Dxi(0);                                                              \
    ASSERT_EQ(simpleModel->physicsResidues[residue].polys[4], x1Dot);          \
                                                                               \
    tkcSet = simpleModel->physicsResidues[residue].polys[0].GetX();            \
    ASSERT_EQ(tkcSet.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);                                       \
    tkcSet = simpleModel->physicsResidues[residue].polys[1].GetX();            \
    ASSERT_EQ(tkcSet.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);                                       \
    tkcSet = simpleModel->physicsResidues[residue].polys[2].GetX();            \
    ASSERT_EQ(tkcSet.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);                                       \
    tkcSet = simpleModel->physicsResidues[residue].polys[3].GetX();            \
    ASSERT_EQ(tkcSet.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);                                       \
    tkcSet = simpleModel->physicsResidues[residue].polys[4].GetX();            \
    ASSERT_EQ(tkcSet.size(), 3);                                               \
    ASSERT_DOUBLE_EQ(tkcSet[0], tValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[1], kValue);                                       \
    ASSERT_DOUBLE_EQ(tkcSet[2], cValue);

    TEST_MASS_1_PHYSICS_RESIDUE(8, kMin, cMin, tMin, KMin, CMin, TMin);
    TEST_MASS_1_PHYSICS_RESIDUE(9, kMin, cMax, tMin, KMin, CMax, TMin);
    TEST_MASS_1_PHYSICS_RESIDUE(10, kMax, cMin, tMin, KMax, CMin, TMin);
    TEST_MASS_1_PHYSICS_RESIDUE(11, kMax, cMax, tMin, KMax, CMax, TMin);
    TEST_MASS_1_PHYSICS_RESIDUE(12, kMin, cMin, tMax, KMin, CMin, TMax);
    TEST_MASS_1_PHYSICS_RESIDUE(13, kMin, cMax, tMax, KMin, CMax, TMax);
    TEST_MASS_1_PHYSICS_RESIDUE(14, kMax, cMin, tMax, KMax, CMin, TMax);
    TEST_MASS_1_PHYSICS_RESIDUE(15, kMax, cMax, tMax, KMax, CMax, TMax);
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

    int nInitialDispLoss = 3 * 3;  // 3k * 3c
    int nInitialVelLoss = 3 * 3;   // 3k * 3c
    int nInitialCondLoss = nInitialDispLoss + nInitialVelLoss;
    int nPhysicsLoss = 3 * 3 * 3;  // 3t * 3k * 3c
    double totalLossTerms = nInitialDispLoss + nInitialDispLoss + nPhysicsLoss;

    double initialConditionsLossWeight = nPhysicsLoss / totalLossTerms;
    double physicsLossWeight = nInitialCondLoss / totalLossTerms;

    double expectedLoss = 0;

    // Initial values for the problem created in SetUp()
    double initialX0 = 0;
    double initialX1 = initialDisplacement;
    double initialXDot0 = 0;
    double initialXDot1 = 0;

    double tMin = 0, kMin = 0, cMin = 0;
    double tMax = 1, kMax = 1, cMax = 1;
    // Initial conditions loss:
    for (double k : std::vector<double>{kMin, (kMin + kMax) / 2, kMax}) {
        for (double c : std::vector<double>{cMin, (cMin + cMax) / 2, cMax}) {
            double t = 0;
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
    expectedLoss *= initialConditionsLossWeight;

    // Physics loss:
    for (double t : std::vector<double>{tMin, (tMax + tMin) / 2, tMax}) {
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
                double K = KMin + k * (KMax - KMin);
                double C = CMin + c * (CMax - CMin);
                double x1DotDot = 1 / m *
                                  (K * modelX0 - K * modelX1 + C * modelX0Dot -
                                   C * modelX1Dot);
                expectedLoss +=
                    physicsLossWeight * pow(modelX0DotDot - x0DotDot, 2);
                expectedLoss +=
                    physicsLossWeight * pow(modelX1DotDot - x1DotDot, 2);
            }
        }
    }

    double loss = secondOrderModel->Loss();

    ASSERT_DOUBLE_EQ(expectedLoss, loss);
}

TEST_F(PimodelTest, LossGradientTest) {
    Pimodel model = Pimodel(this->pd, TMin, TMax, 1, 1, 1);
    model.AddResidues();
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

    int nInitialDispLoss = 2 * 2;  // 2k * 2c
    int nInitialVelLoss = 2 * 2;   // 2k * 2c
    int nInitialCondLoss = nInitialDispLoss + nInitialVelLoss;
    int nPhysicsLoss = 2 * 2 * 2;  // 2t * 2k * 2c
    double totalLossTerms = nInitialDispLoss + nInitialDispLoss + nPhysicsLoss;

    double initialConditionsLossWeight = nPhysicsLoss / totalLossTerms;
    double physicsLossWeight = nInitialCondLoss / totalLossTerms;

    double tMin = 0, kMin = 0, cMin = 0;
    double tMax = 1, kMax = 1, cMax = 1;

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

            expectedGrad[0] += initialConditionsLossWeight *
                               (modelX0 - initialX0) * d_da0_modelX0;
            expectedGrad[0] += initialConditionsLossWeight *
                               (modelX0Dot - initialX0Dot) * d_da0_modelX0Dot;

            expectedGrad[1] += initialConditionsLossWeight *
                               (modelX0 - initialX0) * d_da1_modelX0;
            expectedGrad[1] += initialConditionsLossWeight *
                               (modelX0Dot - initialX0Dot) * d_da1_modelX0Dot;

            expectedGrad[2] += initialConditionsLossWeight *
                               (modelX0 - initialX0) * d_da2_modelX0;
            expectedGrad[2] += initialConditionsLossWeight *
                               (modelX0Dot - initialX0Dot) * d_da2_modelX0Dot;

            expectedGrad[3] += initialConditionsLossWeight *
                               (modelX0 - initialX0) * d_da3_modelX0;
            expectedGrad[3] += initialConditionsLossWeight *
                               (modelX0Dot - initialX0Dot) * d_da3_modelX0Dot;

            expectedGrad[4] += initialConditionsLossWeight *
                               (modelX1 - initialX1) * d_da4_modelX1;
            expectedGrad[4] += initialConditionsLossWeight *
                               (modelX1Dot - initialX1Dot) * d_da4_modelX1Dot;

            expectedGrad[5] += initialConditionsLossWeight *
                               (modelX1 - initialX1) * d_da5_modelX1;
            expectedGrad[5] += initialConditionsLossWeight *
                               (modelX1Dot - initialX1Dot) * d_da5_modelX1Dot;

            expectedGrad[6] += initialConditionsLossWeight *
                               (modelX1 - initialX1) * d_da6_modelX1;
            expectedGrad[6] += initialConditionsLossWeight *
                               (modelX1Dot - initialX1Dot) * d_da6_modelX1Dot;

            expectedGrad[7] += initialConditionsLossWeight *
                               (modelX1 - initialX1) * d_da7_modelX1;
            expectedGrad[7] += initialConditionsLossWeight *
                               (modelX1Dot - initialX1Dot) * d_da7_modelX1Dot;
        }
    }

    // Physics loss:
    for (double t : std::vector<double>{tMin, tMax}) {
        for (double k : std::vector<double>{kMin, kMax}) {
            for (double c : std::vector<double>{cMin, cMax}) {
                double modelX0 = a0 * t + a1 * k + a2 * c + a3 * 1;
                double modelX1 = a4 * t + a5 * k + a6 * c + a7 * 1;
                double modelX0Dot = a0;
                double modelX1Dot = a4;
                double modelX0DotDot = 0;
                double modelX1DotDot = 0;

                double K = KMin + k * (KMax - KMin);
                double C = CMin + c * (CMax - CMin);
                // Mass 0 is fixed
                double X0DotDot = 0;
                double X1DotDot = 1 / m *
                                  (K * modelX0 - K * modelX1 + C * modelX0Dot -
                                   C * modelX1Dot);

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
                    (K * d_da0_modelX0 - K * d_da0_modelX1 +
                     C * d_da0_modelX0Dot - C * d_da0_modelX1Dot);
                double d_da1_X0DotDot = 0;
                double d_da1_X1DotDot =
                    1 / m *
                    (K * d_da1_modelX0 - K * d_da1_modelX1 +
                     C * d_da1_modelX0Dot - C * d_da1_modelX1Dot);
                double d_da2_X0DotDot = 0;
                double d_da2_X1DotDot =
                    1 / m *
                    (K * d_da2_modelX0 - K * d_da2_modelX1 +
                     C * d_da2_modelX0Dot - C * d_da2_modelX1Dot);
                double d_da3_X0DotDot = 0;
                double d_da3_X1DotDot =
                    1 / m *
                    (K * d_da3_modelX0 - K * d_da3_modelX1 +
                     C * d_da3_modelX0Dot - C * d_da3_modelX1Dot);
                double d_da4_X0DotDot = 0;
                double d_da4_X1DotDot =
                    1 / m *
                    (K * d_da4_modelX0 - K * d_da4_modelX1 +
                     C * d_da4_modelX0Dot - C * d_da4_modelX1Dot);
                double d_da5_X0DotDot = 0;
                double d_da5_X1DotDot =
                    1 / m *
                    (K * d_da5_modelX0 - K * d_da5_modelX1 +
                     C * d_da5_modelX0Dot - C * d_da5_modelX1Dot);
                double d_da6_X0DotDot = 0;
                double d_da6_X1DotDot =
                    1 / m *
                    (K * d_da6_modelX0 - K * d_da6_modelX1 +
                     C * d_da6_modelX0Dot - C * d_da6_modelX1Dot);
                double d_da7_X0DotDot = 0;
                double d_da7_X1DotDot =
                    1 / m *
                    (K * d_da7_modelX0 - K * d_da7_modelX1 +
                     C * d_da7_modelX0Dot - C * d_da7_modelX1Dot);

                expectedGrad[0] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da0_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da0_X1DotDot));
                expectedGrad[1] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da1_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da1_X1DotDot));
                expectedGrad[2] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da2_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da2_X1DotDot));
                expectedGrad[3] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da3_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da3_X1DotDot));
                expectedGrad[4] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da4_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da4_X1DotDot));
                expectedGrad[5] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da5_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da5_X1DotDot));
                expectedGrad[6] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da6_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da6_X1DotDot));
                expectedGrad[7] += physicsLossWeight *
                                   ((modelX0DotDot - X0DotDot) *
                                        (d_dai_modelX0DotDot - d_da7_X0DotDot) +
                                    (modelX1DotDot - X1DotDot) *
                                        (d_dai_modelX1DotDot - d_da7_X1DotDot));
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

class PimodelTrainTest : public testing::Test {
   public:
    double mass;
    double kMin;
    double kMax;
    double cMin;
    double cMax;
    double tMin;
    double tMax;
    double initialDisp;
    ProblemDescription pd;

    Pimodel model;

    int timeDiscretization;
    int kcDiscretization;
    int order;
    double learningRate;
    int maxSteps;
    bool log;

    void SetUp() {
        this->mass = 1;
        this->kMin = 0.8;
        this->kMax = 1.0;
        this->cMin = 0.0;
        this->cMax = 0.2;
        this->tMin = 0.0;
        this->tMax = 1.0;
        this->timeDiscretization = 2;
        this->kcDiscretization = 1;
        this->order = 3;
        this->learningRate = 0.01;
        this->maxSteps = 500;
        this->log = true;
        this->pd.AddMass(mass, 0.0, 0.0);
        this->pd.AddMass(mass, 1.0, 0.0);
        this->pd.AddSpring(0, 1, kMin, kMax);
        this->pd.AddDamper(0, 1, cMin, cMax);
        this->pd.SetFixedMass(0);
    }

    void Train() {
        this->model =
            Pimodel(this->pd, this->tMin, this->tMax, this->timeDiscretization,
                    this->kcDiscretization, this->order);
        this->model.AddResidues();

        double initialLoss = this->model.Loss();
        this->model.Train(this->learningRate, this->maxSteps, this->log);
        ASSERT_TRUE(this->model.Loss() < initialLoss);
    }

    void Plot() {
        // Get problem using intermediate value for k and c, and integrate it.
        // Then, compare the model's prediction with the problem's result.
        double k = (this->kMin + this->kMax) / 2;
        double c = (this->cMin + this->cMax) / 2;
        Maybe<Problem> mP = pd.BuildFromVector(std::vector<double>{k, c});
        ASSERT_FALSE(mP.isError);
        Problem p = mP.val;
        p.Integrate(this->tMax);

        std::vector<double> tkc = std::vector<double>{0.0, k, c};
        Maybe<std::vector<double>> X;
        Maybe<std::vector<double>> XDot;
        std::cout << "t,x0,x0Dot,modelX0,modelX0Dot,x1,x1Dot,modelX1,modelX1Dot"
                  << std::endl;

        // Plot the analytical values and the model predictions
        for (int i = 0; i < int(p.t.size()); i += 1) {
            tkc[0] = p.t[i];
            X = model(&tkc);
            ASSERT_FALSE(X.isError);
            XDot = model.GetVelocities(&tkc);
            ASSERT_FALSE(XDot.isError);

            std::cout << p.t[i] << ",";             // t,
            std::cout << p.XHistory[i][0] << ",";   // x0,
            std::cout << p.XHistory[i][2] << ",";   // x0Dot,
            std::cout << X.val[0] << ",";           // modelX0,
            std::cout << XDot.val[0] << ",";        // modelX0Dot,
            std::cout << p.XHistory[i][1] << ",";   // x1,
            std::cout << p.XHistory[i][3] << ",";   // x1Dot,
            std::cout << X.val[1] << ",";           // modelX1,
            std::cout << XDot.val[1] << std::endl;  // modelX1Dot
        }
    }
};

TEST_F(PimodelTrainTest, InitialDispTrainTest) {
    this->log = false;
    this->pd.AddInitialDisp(1, 1.0);
    this->Train();
    this->Plot();
}

TEST_F(PimodelTrainTest, InitialVelTrainTest) {
    this->log = false;
    this->pd.AddInitialVel(1, 1.0);
    this->Train();
    this->Plot();
}

TEST_F(PimodelTrainTest, InitialDispAndVelTrainTest) {
    this->log = false;
    this->pd.AddInitialDisp(1, 1.0);
    this->pd.AddInitialVel(1, 1.0);
    this->Train();
    this->Plot();
}

class PimodelsTest : public testing::Test {
   public:
    ProblemDescription pd;

    void SetUp() {
        // Called before every TEST_F
        // Create problem description of two masses,
        // connected with a spring and damper. No fixed or initial
        // conditions are set
        this->pd = ProblemDescription();
        this->pd.AddMass(m, 0.0, 0.0);
        this->pd.AddMass(m, 1.0, 0.0);
        this->pd.AddSpring(0, 1, KMin, KMax);
        this->pd.AddDamper(0, 1, CMin, CMax);
    }
};

TEST_F(PimodelsTest, ConstructorTest) {
    double t = 100;
    // Only one model
    Pimodels pimodels = Pimodels(this->pd, t, 1, 2, 2, 2);
    ASSERT_EQ(pimodels.timeBuckets.size(), 2);
    ASSERT_EQ(pimodels.timeBuckets[0], 0.0);
    ASSERT_EQ(pimodels.timeBuckets[1], t);
    ASSERT_EQ(pimodels.pimodels.size(), 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t0, 0.0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t1, t);

    // Two models
    pimodels = Pimodels(this->pd, t, 2, 1, 1, 1);
    ASSERT_EQ(pimodels.timeBuckets.size(), 3);
    ASSERT_EQ(pimodels.timeBuckets[0], 0.0);
    ASSERT_DOUBLE_EQ(pimodels.timeBuckets[1], t / 2);
    ASSERT_DOUBLE_EQ(pimodels.timeBuckets[2], t);
    ASSERT_EQ(pimodels.pimodels.size(), 2);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t0, 0.0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t1, t / 2);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].t0, t / 2);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].t1, t);

    // Three models
    pimodels = Pimodels(this->pd, t, 3, 1, 1, 1);
    ASSERT_EQ(pimodels.timeBuckets.size(), 4);
    ASSERT_EQ(pimodels.timeBuckets[0], 0.0);
    ASSERT_DOUBLE_EQ(pimodels.timeBuckets[1], t / 3);
    ASSERT_DOUBLE_EQ(pimodels.timeBuckets[2], 2 * t / 3);
    ASSERT_DOUBLE_EQ(pimodels.timeBuckets[3], t);

    ASSERT_EQ(pimodels.pimodels.size(), 3);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t0, 0.0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[0].t1, t / 3);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].t0, t / 3);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].t1, 2 * t / 3);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].t0, 2 * t / 3);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].t1, t);

    // Errors
    ASSERT_DEATH({ Pimodels(this->pd, -1.0, 1, 1, 1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, -1, 1, 1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 0, 1, 1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, -1, 1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, 0, 1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, 1, -1, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, 1, 0, 1); }, "");
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, 1, 1, -1); }, "");
}

TEST_F(PimodelsTest, GetTimeBucketTest) {
    double t = 100;

    // Only one model
    Pimodels pimodels = Pimodels(this->pd, t, 1, 2, 2, 2);
    ASSERT_EQ(pimodels.getTimeBucket(0.0), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t / 2), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t), 0);

    // Two models
    // time_buckets = [0, t/2, t]
    // b0 = [0,t/2)
    // b1 = [t/2,t]
    pimodels = Pimodels(this->pd, t, 2, 2, 2, 2);
    ASSERT_EQ(pimodels.getTimeBucket(0.0), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t / 4), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t / 2 * 0.99), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t / 2), 1);
    ASSERT_EQ(pimodels.getTimeBucket((t / 2 + t) / 2), 1);
    ASSERT_EQ(pimodels.getTimeBucket(t), 1);

    // Three models
    // time_buckets = [0, t/3, 2t/3, t]
    // b0 = [0,t/3)
    // b1 = [t/3,2t/3)
    // b2 = [2t/3, t]
    pimodels = Pimodels(this->pd, t, 3, 2, 2, 2);
    double t0 = 0;
    double t1 = t / 3;
    double t2 = 2 * t / 3;
    double t3 = t;

    ASSERT_EQ(pimodels.getTimeBucket(t0), 0);
    ASSERT_EQ(pimodels.getTimeBucket((t0 + t1) / 2), 0);
    ASSERT_EQ(pimodels.getTimeBucket(t1 * 0.99), 0);

    ASSERT_EQ(pimodels.getTimeBucket(t1), 1);
    ASSERT_EQ(pimodels.getTimeBucket((t1 + t2) / 2), 1);
    ASSERT_EQ(pimodels.getTimeBucket(t2 * 0.99), 1);

    ASSERT_EQ(pimodels.getTimeBucket(t2), 2);
    ASSERT_EQ(pimodels.getTimeBucket((t2 + t3) / 2), 2);
    ASSERT_EQ(pimodels.getTimeBucket(t3 * 0.99), 2);
    ASSERT_EQ(pimodels.getTimeBucket(t3), 2);

    // Errors
    ASSERT_DEATH({ pimodels.getTimeBucket(-0.1); }, "");
    ASSERT_DEATH({ pimodels.getTimeBucket(t * 1.1); }, "");
}

TEST_F(PimodelsTest, getContinuityTkcTest) {
    double k0Min = Random(1, 100);
    double k0Max = k0Min * (1 + Random());
    double k1Min = Random(1, 100);
    double k1Max = k1Min * (1 + Random());
    double k2Min = Random(1, 100);
    double k2Max = k2Min * (1 + Random());
    double c0Min = Random(1, 100);
    double c0Max = c0Min * (1 + Random());
    double c1Min = Random(1, 100);
    double c1Max = c1Min * (1 + Random());
    double c2Min = Random(1, 100);
    double c2Max = c2Min * (1 + Random());

    auto pd = ProblemDescription();
    pd.AddMass(Random(), 0.0, 0.0);
    pd.AddMass(Random(), 1.0, 0.0);
    pd.AddMass(Random(), 2.0, 0.0);
    pd.AddSpring(0, 1, k0Min, k0Max);
    pd.AddDamper(0, 1, c0Min, c0Max);
    pd.AddSpring(1, 2, k1Min, k1Max);
    pd.AddDamper(1, 2, c1Min, c1Max);
    pd.AddSpring(0, 2, k2Min, k2Max);
    pd.AddDamper(0, 2, c2Min, c2Max);

    auto models = Pimodels(pd, 1.0, 2, 2, 2, 2);
    auto tkc = models.continuityTkc();

    ASSERT_EQ(tkc.size(), 7);
    ASSERT_DOUBLE_EQ(tkc[0], 0);
    ASSERT_DOUBLE_EQ(tkc[1], (k0Min + k0Max) / 2);
    ASSERT_DOUBLE_EQ(tkc[2], (k1Min + k1Max) / 2);
    ASSERT_DOUBLE_EQ(tkc[3], (k2Min + k2Max) / 2);
    ASSERT_DOUBLE_EQ(tkc[4], (c0Min + c0Max) / 2);
    ASSERT_DOUBLE_EQ(tkc[5], (c1Min + c1Max) / 2);
    ASSERT_DOUBLE_EQ(tkc[6], (c2Min + c2Max) / 2);
}

TEST_F(PimodelsTest, setContinuityTest) {
    auto a = std::vector<double>{Random(), Random(), Random(), Random()};
    auto b = std::vector<double>{Random(), Random(), Random(), Random()};
    auto ab =
        std::vector<double>{a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3]};
    auto c = std::vector<double>{Random(), Random(), Random(), Random()};
    auto d = std::vector<double>{Random(), Random(), Random(), Random()};
    auto cd =
        std::vector<double>{c[0], c[1], c[2], c[3], d[0], d[1], d[2], d[3]};
    auto e = std::vector<double>{Random(), Random(), Random(), Random()};
    auto f = std::vector<double>{Random(), Random(), Random(), Random()};
    auto ef =
        std::vector<double>{e[0], e[1], e[2], e[3], f[0], f[1], f[2], f[3]};

    double totalT = 100 + Random();
    double K = Random(KMin, KMax);
    double C = Random(CMin, CMax);
    std::vector<double> TKC = std::vector<double>{0, K, C};
    double k_ = Normalize(K, KMin, KMax).val.Get();
    double c_ = Normalize(C, CMin, CMax).val.Get();
    double t_ = 1;  // In the point where we set continuity, t(local) = 1

    // Single Pimodel
    Pimodels pimodels = Pimodels(this->pd, totalT, 1, 1, 1, 1);
    ASSERT_DEATH({ pimodels.setContinuity(0, TKC); }, "");

    // Two Pimodels
    // t in [0, tMax/2):
    //  x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    //  x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    // t in [tMax/2, tMax]:
    //  x0(t,k,c) = c0*t + c1*k + c2*c + c3*1
    //  x1(t,k,c) = d0*t + d1*k + d2*c + d3*1
    pimodels = Pimodels(this->pd, totalT, 2, 1, 1, 1);
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&ab).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&cd).isError);
    ASSERT_DEATH({ pimodels.setContinuity(0, TKC); }, "");
    pimodels.setContinuity(1, TKC);
    // No initial conditions are added to the first pimodel
    ASSERT_EQ(pimodels.pimodels[0].p.initialDisps.size(), 0);
    ASSERT_EQ(pimodels.pimodels[0].p.initialVels.size(), 0);
    // One initial condition (disp and vel) is added for each mass
    // for the second pimodel
    ASSERT_EQ(pimodels.pimodels[1].p.initialDisps.size(), 2);
    ASSERT_EQ(pimodels.pimodels[1].p.initialVels.size(), 2);

    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].val,
                     a[0] * t_ + a[1] * k_ + a[2] * c_ + a[3]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].val,
                     b[0] * t_ + b[1] * k_ + b[2] * c_ + b[3]);

    // C1 continuity
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].val, a[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].val, b[0]);

    // Three Pimodels
    // t in [0, tMax/3):
    //  x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    //  x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    // t in [tMax/3, 2*tMax/3):
    //  x0(t,k,c) = c0*t + c1*k + c2*c + c3*1
    //  x1(t,k,c) = d0*t + d1*k + d2*c + d3*1
    // t in [2*tMax/3, tMax]:
    //  x0(t,k,c) = e0*t + e1*k + e2*e + e3*1
    //  x1(t,k,c) = f0*t + f1*k + f2*c + f3*1
    pimodels = Pimodels(this->pd, totalT, 3, 1, 1, 1);
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&ab).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&cd).isError);
    ASSERT_FALSE(pimodels.pimodels[2].SetParameters(&ef).isError);
    pimodels.setContinuity(1, TKC);
    pimodels.setContinuity(2, TKC);
    ASSERT_DEATH({ pimodels.setContinuity(3, TKC); }, "");

    ASSERT_EQ(pimodels.pimodels[0].p.initialDisps.size(), 0);
    ASSERT_EQ(pimodels.pimodels[0].p.initialVels.size(), 0);
    ASSERT_EQ(pimodels.pimodels[1].p.initialDisps.size(), 2);
    ASSERT_EQ(pimodels.pimodels[1].p.initialVels.size(), 2);
    ASSERT_EQ(pimodels.pimodels[2].p.initialDisps.size(), 2);
    ASSERT_EQ(pimodels.pimodels[2].p.initialVels.size(), 2);

    // Pimodel 1
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].val,
                     a[0] * t_ + a[1] * k_ + a[2] * c_ + a[3]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].val,
                     b[0] * t_ + b[1] * k_ + b[2] * c_ + b[3]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].val, a[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].val, b[0]);

    // Pimodel 2
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[0].val,
                     c[0] * t_ + c[1] * k_ + c[2] * c_ + c[3]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[1].val,
                     d[0] * t_ + d[1] * k_ + d[2] * c_ + d[3]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[0].val, c[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[1].val, d[0]);
}

TEST_F(PimodelsTest, OperatorTest) {
    auto a = std::vector<double>{Random(), Random(), Random(), Random()};
    auto b = std::vector<double>{Random(), Random(), Random(), Random()};
    auto ab =
        std::vector<double>{a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3]};
    auto c = std::vector<double>{Random(), Random(), Random(), Random()};
    auto d = std::vector<double>{Random(), Random(), Random(), Random()};
    auto cd =
        std::vector<double>{c[0], c[1], c[2], c[3], d[0], d[1], d[2], d[3]};
    auto e = std::vector<double>{Random(), Random(), Random(), Random()};
    auto f = std::vector<double>{Random(), Random(), Random(), Random()};
    auto ef =
        std::vector<double>{e[0], e[1], e[2], e[3], f[0], f[1], f[2], f[3]};

    double totalT = 100 + Random();
    double K = Random(KMin, KMax);
    double C = Random(CMin, CMax);
    std::vector<double> TKC = std::vector<double>{0, K, C};
    double k_ = Normalize(K, KMin, KMax).val.Get();
    double c_ = Normalize(C, CMin, CMax).val.Get();

    // Three Pimodels
    // t in [0, tMax/3):
    //  x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    //  x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    // t in [tMax/3, 2*tMax/3):
    //  x0(t,k,c) = c0*t + c1*k + c2*c + c3*1
    //  x1(t,k,c) = d0*t + d1*k + d2*c + d3*1
    // t in [2*tMax/3, tMax]:
    //  x0(t,k,c) = e0*t + e1*k + e2*e + e3*1
    //  x1(t,k,c) = f0*t + f1*k + f2*c + f3*1
    Pimodels pimodels = Pimodels(this->pd, totalT, 3, 1, 1, 1);
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&ab).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&cd).isError);
    ASSERT_FALSE(pimodels.pimodels[2].SetParameters(&ef).isError);

    double T;
    double t_;
    // Pimodel 0
    T = 0;
    TKC[0] = T;
    t_ = Normalize(T, 0, totalT / 3).val.Get();
    auto eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val[0], a[0] * t_ + a[1] * k_ + a[2] * c_ + a[3]);
    ASSERT_DOUBLE_EQ(eval.val[1], b[0] * t_ + b[1] * k_ + b[2] * c_ + b[3]);

    // Pimodel 1
    T = totalT / 3;
    TKC[0] = T;
    t_ = Normalize(T, totalT / 3, 2 * totalT / 3).val.Get();
    eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val[0], c[0] * t_ + c[1] * k_ + c[2] * c_ + c[3]);
    ASSERT_DOUBLE_EQ(eval.val[1], d[0] * t_ + d[1] * k_ + d[2] * c_ + d[3]);

    // Pimodel 2
    T = 2 * totalT / 3;
    TKC[0] = T;
    t_ = Normalize(T, 2 * totalT / 3, totalT).val.Get();
    eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val[0], e[0] * t_ + e[1] * k_ + e[2] * c_ + e[3]);
    ASSERT_DOUBLE_EQ(eval.val[1], f[0] * t_ + f[1] * k_ + f[2] * c_ + f[3]);
}

TEST(PimodelsTrainingTest, TrainTest) {
    double mass = 1;
    double kMin = 0.8;
    double kMax = 1.0;
    double cMin = 0.0;
    double cMax = 0.05;
    double tMax = 8.0;
    double initialDisp = 1.0;

    auto pd = ProblemDescription();
    pd.AddMass(mass, 0.0, 0.0);
    pd.AddMass(mass, 1.0, 0.0);
    pd.AddSpring(0, 1, kMin, kMax);
    pd.AddDamper(0, 1, cMin, cMax);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, initialDisp);

    int nModels = 8;
    int timeDiscretization = 2;
    int kcDiscretization = 1;
    int order = 3;
    double learningRate = 0.01;
    int maxSteps = 500;
    bool log = false;

    // Train all models
    Pimodels models = Pimodels(pd, tMax, nModels, timeDiscretization,
                               kcDiscretization, order);
    models.Train(learningRate, maxSteps, log);

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
    Maybe<std::vector<double>> XDot;
    std::cout << "t,x0,x0Dot,modelX0,modelX0Dot,x1,x1Dot,modelX1,modelX1Dot"
              << std::endl;

    // Plot the analytical values and the model predictions
    for (int i = 0; i < int(p.t.size()); i += 1) {
        tkc[0] = p.t[i];
        X = models(&tkc);
        ASSERT_FALSE(X.isError);
        XDot = models.GetVelocities(&tkc);
        ASSERT_FALSE(XDot.isError);

        std::cout << p.t[i] << ",";             // t,
        std::cout << p.XHistory[i][0] << ",";   // x0,
        std::cout << p.XHistory[i][2] << ",";   // x0Dot,
        std::cout << X.val[0] << ",";           // modelX0,
        std::cout << XDot.val[0] << ",";        // modelX0Dot,
        std::cout << p.XHistory[i][1] << ",";   // x1,
        std::cout << p.XHistory[i][3] << ",";   // x1Dot,
        std::cout << X.val[1] << ",";           // modelX1,
        std::cout << XDot.val[1] << std::endl;  // modelX1Dot
    }
}
