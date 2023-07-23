#include "pimodel.h"

#include <gtest/gtest.h>

#include <vector>

#include "bounded.h"
#include "polynomial.h"
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
const double initialVelocity = 7.0;

class PimodelTest : public testing::Test {
   public:
    ProblemDescription pd;

    // Variables that are used to evaluate the models
    Bounded tB;
    Bounded kB;
    Bounded cB;
    std::vector<Bounded> tkc;  // Normalized
    std::vector<double> TKC;   // Un-normalized

    Pimodel piModel;
    Poly p0;  // model that describes displacement of mass 0
    Poly p1;  // model that describes displacement of mass 1
    const int timeDiscretization = 2;
    const int kcDiscretization = 2;
    const int order = 2;

    int parametersPerModel;
    std::vector<double> p0Parameters;
    std::vector<double> p1Parameters;
    std::vector<double> piModelParameters;

    double x;
    double y;
    std::vector<double> X;

    void SetUp() {
        // Called before every TEST_F
        // Create problem description of two masses,
        // connected with a spring and damper. One is fixed
        // and the other has initial displacement.
        pd = ProblemDescription();
        pd.AddMass(m, 0.0, 0.0);
        pd.AddMass(m, 1.0, 0.0);
        pd.AddSpring(0, 1, KMin, KMax);
        pd.AddDamper(0, 1, CMin, CMax);
        pd.SetFixedMass(0);
        pd.AddInitialDisp(1, initialDisplacement);
        pd.AddInitialVel(1, initialVelocity);
        std::vector<Bounded> dna = std::vector<Bounded>(2);
        auto e0 = pd.BuildFromDNA(dna);
        ASSERT_FALSE(e0.isError);

        piModel = Pimodel(pd, TMin, TMax, timeDiscretization, kcDiscretization,
                          order);
        piModel.SetResidues(true, true);
        parametersPerModel = piModel.nParameters() / 2;

        // 3 input variables:
        // t, k, c
        p0.Build(3, order, true, 0);
        p1.Build(3, order, true, 1);

        tB = RandomB();
        kB = RandomB();
        cB = RandomB();
        tkc = std::vector<Bounded>{tB, kB, cB};
        TKC = std::vector<double>(3);
        TKC[0] = Unnormalize(tB, TMin, TMax);
        TKC[1] = Unnormalize(kB, KMin, KMax);
        TKC[2] = Unnormalize(cB, CMin, CMax);

        p0.SetX(Bounded::Get(tkc));
        p1.SetX(Bounded::Get(tkc));

        p0Parameters = std::vector<double>(parametersPerModel);
        p1Parameters = std::vector<double>(parametersPerModel);
        for (int i = 0; i < parametersPerModel; i++) {
            p0Parameters[i] = Random();
            p1Parameters[i] = Random();
        }
        piModelParameters = std::vector<double>(2 * parametersPerModel);
        int param = 0;
        for (int i = 0; i < parametersPerModel; i++) {
            piModelParameters[param] = p0Parameters[i];
            param++;
        }
        for (int i = 0; i < parametersPerModel; i++) {
            piModelParameters[param] = p1Parameters[i];
            param++;
        }
        piModel.SetParameters(&piModelParameters);

        x = Random(10, 1000);
        y = Random(10, 1000);
        X = std::vector<double>{x, y};
    }
};

TEST_F(PimodelTest, ConstructorTest) {
    ASSERT_EQ(piModel.models.size1(), 2);
    ASSERT_EQ(piModel.models.size2(), 1);
    ASSERT_EQ(piModel.models(0, 0).id, 0);
    ASSERT_EQ(piModel.models(1, 0).id, 1);

    ASSERT_EQ(piModel.modelsD.size(), 2);
    ASSERT_EQ(piModel.modelsD[0].id, 0);
    ASSERT_EQ(piModel.modelsD[1].id, 1);

    ASSERT_EQ(piModel.modelsDD.size(), 2);
    ASSERT_EQ(piModel.modelsDD[0].id, 0);
    ASSERT_EQ(piModel.modelsDD[1].id, 1);
}

TEST_F(PimodelTest, OperatorTest) {
    auto eval = piModel(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val.size(), 2);
    ASSERT_DOUBLE_EQ(eval.val[0], p0(p0Parameters).val);
    ASSERT_DOUBLE_EQ(eval.val[1], p1(p1Parameters).val);

    // d/dT = dt/dT*d/dt
    ASSERT_FALSE(p0.Dxi(0).isError);
    p0 *= 1 / (TMax - TMin);
    ASSERT_FALSE(p1.Dxi(0).isError);
    p1 *= 1 / (TMax - TMin);

    eval = piModel.GetVelocities(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val.size(), 2);
    ASSERT_DOUBLE_EQ(eval.val[0], p0(p0Parameters).val);
    ASSERT_DOUBLE_EQ(eval.val[1], p1(p1Parameters).val);

    // Test error cases
    // tkc too large
    TKC = std::vector<double>{1.1, 1.0, 2.0, 3.0};
    eval = piModel(&TKC);
    ASSERT_TRUE(eval.isError);

    // tkc too small
    TKC = std::vector<double>{1.1, 1.0};
    eval = piModel(&TKC);
    ASSERT_TRUE(eval.isError);
}

TEST_F(PimodelTest, nParametersTest) {
    ASSERT_EQ(piModel.nParameters(), p0.nMonomials() + p1.nMonomials());
}

TEST_F(PimodelTest, SetParametersTest) {
    double p = 99;
    std::vector<double> params = std::vector<double>(parametersPerModel * 2, p);
    auto r = piModel.SetParameters(&params);
    ASSERT_FALSE(r.isError);

    // piModel->modelsCoefficients[mass][coeff]
    for (int i = 0; i < parametersPerModel; i++) {
        ASSERT_EQ(piModel.modelsCoefficients[0][i], p);
    }
    for (int i = 0; i < parametersPerModel; i++) {
        ASSERT_EQ(piModel.modelsCoefficients[1][i], p);
    }

    // Test error cases
    params = std::vector<double>(parametersPerModel * 2 - 1);
    r = piModel.SetParameters(&params);
    ASSERT_TRUE(r.isError);
    params = std::vector<double>(parametersPerModel * 2 + 1);
    r = piModel.SetParameters(&params);
    ASSERT_TRUE(r.isError);
}

TEST_F(PimodelTest, GetParametersTest) {
    std::vector<double> params = std::vector<double>(parametersPerModel * 2);
    ASSERT_FALSE(piModel.GetParameters(&params).isError);

    int paramIndex = 0;
    for (int m = 0; m < 2; m++) {
        for (int i = 0; i < parametersPerModel; i++) {
            ASSERT_EQ(piModel.modelsCoefficients[m][i], params[paramIndex]);
            paramIndex += 1;
        }
    }

    // Test error cases
    params = std::vector<double>(parametersPerModel * 2 - 1);
    ASSERT_TRUE(piModel.GetParameters(&params).isError);
    params = std::vector<double>(parametersPerModel * 2 + 1);
    ASSERT_TRUE(piModel.GetParameters(&params).isError);
}

TEST_F(PimodelTest, ProblemFromTkcTest) {
    Problem p = piModel.problemFromTkc(&tkc);

    ASSERT_EQ(p.masses.size(), 2);
    ASSERT_EQ(p.masses[0].m, m);
    ASSERT_EQ(p.masses[1].m, m);

    ASSERT_EQ(p.springs.size(), 1);
    ASSERT_EQ(p.springs[0].Get_k(), Unnormalize(tkc[1], KMin, KMax));
    ASSERT_EQ(p.dampers.size(), 1);
    ASSERT_EQ(p.dampers[0].Get_c(), Unnormalize(tkc[2], CMin, CMax));
}

TEST_F(PimodelTest, getXModelTest) {
    std::vector<double> XModel = piModel.getXModel(&tkc);

    ASSERT_EQ(XModel.size(), 4);
    ASSERT_EQ(XModel[0], p0(p0Parameters).val);
    ASSERT_EQ(XModel[1], p1(p1Parameters).val);

    // d/dT = dt/dT*d/dt
    ASSERT_FALSE(p0.Dxi(0).isError);
    p0 *= 1 / (TMax - TMin);
    ASSERT_FALSE(p1.Dxi(0).isError);
    p1 *= 1 / (TMax - TMin);

    ASSERT_EQ(XModel[2], p0(p0Parameters).val);  // dx0/dT = 1
    ASSERT_EQ(XModel[3], p1(p1Parameters).val);  // dx1/dT = 5
}

TEST_F(PimodelTest, getAccelsFromDiffEqTest) {
    Maybe<Problem> mProblem =
        this->pd.BuildFromVector(std::vector<double>{KMin, CMin});
    ASSERT_FALSE(mProblem.isError);
    Problem problem = mProblem.val;

    auto A = piModel.getAccelsFromDiffEq(&problem, tkc);
    ASSERT_EQ(A.size1(), 2);
    ASSERT_EQ(A.size2(), 1);

    Poly dp0dt = p0;
    ASSERT_FALSE(dp0dt.Dxi(0).isError);
    dp0dt *= 1 / (TMax - TMin);

    Poly dp1dt = p1;
    ASSERT_FALSE(dp1dt.Dxi(0).isError);
    dp1dt *= 1 / (TMax - TMin);

    // expectedAccel0(t,k,c) =
    // 1 / m * (-kMin * p0 + kMin * p1 - cMin * dp0dt + cMin * dp1dt);
    // expectedAccel1(t,k,c) = -expectedAccel0(t,k,c)
    Polys expectedPolys =
        1 / m * (-KMin * p0 + KMin * p1 + (-CMin) * dp0dt + CMin * dp1dt);

    ASSERT_EQ(A(0, 0), expectedPolys);
    ASSERT_EQ(A(1, 0), -1 * expectedPolys);
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
    auto model = Pimodel(pd, TMin, TMax, 0, 0, 1);
    model.SetResidues(true, false);
    ASSERT_EQ(model.initialConditionsResiduesTkc.size(), 1);
    ASSERT_DOUBLE_EQ(model.initialConditionsResiduesTkc[0][0].Get(), 0);
    ASSERT_DOUBLE_EQ(model.initialConditionsResiduesTkc[0][1].Get(), 0.5);
    ASSERT_DOUBLE_EQ(model.initialConditionsResiduesTkc[0][2].Get(), 0.5);

    double t = 0;
    double k = 0;
    double c = 0;

    int i = 0;
    while (k <= 1) {
        c = 0;
        while (c <= 1) {
            ASSERT_DOUBLE_EQ(piModel.initialConditionsResiduesTkc[i][0].Get(),
                             t);
            ASSERT_DOUBLE_EQ(piModel.initialConditionsResiduesTkc[i][1].Get(),
                             k);
            ASSERT_DOUBLE_EQ(piModel.initialConditionsResiduesTkc[i][2].Get(),
                             c);
            i++;
            c += 1 / (kcDiscretization + 0.0);
        }
        k += 1 / (kcDiscretization + 0.0);
    }
}

TEST_F(PimodelTest, PhysicsResiduesTkcTest) {
    auto model = Pimodel(pd, TMin, TMax, 0, 0, 1);
    model.SetResidues(false, true);
    ASSERT_EQ(model.physicsResiduesTkc.size(), 1);
    ASSERT_DOUBLE_EQ(model.physicsResiduesTkc[0][0].Get(), 0.5);
    ASSERT_DOUBLE_EQ(model.physicsResiduesTkc[0][1].Get(), 0.5);
    ASSERT_DOUBLE_EQ(model.physicsResiduesTkc[0][2].Get(), 0.5);

    double t = 0;
    double k = 0;
    double c = 0;

    int i = 0;
    while (t <= 1) {
        k = 0;
        while (k <= 1) {
            c = 0;
            while (c <= 1) {
                ASSERT_DOUBLE_EQ(piModel.physicsResiduesTkc[i][0].Get(), t);
                ASSERT_DOUBLE_EQ(piModel.physicsResiduesTkc[i][1].Get(), k);
                ASSERT_DOUBLE_EQ(piModel.physicsResiduesTkc[i][2].Get(), c);
                i++;
                c += 1 / (kcDiscretization + 0.0);
            }
            k += 1 / (kcDiscretization + 0.0);
        }
        t += 1 / (kcDiscretization + 0.0);
    }
}

TEST_F(PimodelTest, InitialConditionsResiduesTest) {
    // For each mass, one residue for each tkc value
    // Ex:
    // 2 masses
    // timeDiscretization = 1
    // kcDiscretization = 0
    // -> tkc values = [0,0.5,0.5], [1,0.5,0.5]
    // initialDisplacementResidues:
    //  (x0-x0Model(t = 0, k = 0.5, c = 0.5))
    //  (x0-x0Model(t = 1, k = 0.5, c = 0.5))
    //  (x1-x1Model(t = 0, k = 0.5, c = 0.5))
    //  (x1-x1Model(t = 1, k = 0.5, c = 0.5))
    ASSERT_EQ(piModel.initialDispResidues.size(),
              piModel.initialConditionsResiduesTkc.size() * 2);  // 2 masses
    ASSERT_EQ(piModel.initialVelResidues.size(),
              piModel.initialConditionsResiduesTkc.size() * 2);  // 2 masses

    // Verify all the X values are set correctly in the models
    for (int resId = 0;
         resId < int(piModel.initialConditionsResiduesTkc.size()); resId++) {
        for (int massId = 0; massId < 2; massId++) {
            for (int p = 0;
                 p < int(piModel.initialDispResidues[resId].polys.size());
                 p++) {
                auto tkc = piModel.initialDispResidues[resId].polys[p].GetX();
                ASSERT_DOUBLE_EQ(
                    tkc[0],  // t
                    piModel.initialConditionsResiduesTkc[resId][0].Get());
                ASSERT_DOUBLE_EQ(
                    tkc[1],  // k
                    piModel.initialConditionsResiduesTkc[resId][1].Get());
                ASSERT_DOUBLE_EQ(
                    tkc[2],  // c
                    piModel.initialConditionsResiduesTkc[resId][2].Get());

                tkc = piModel.initialVelResidues[resId].polys[p].GetX();
                ASSERT_DOUBLE_EQ(
                    tkc[0],  // t
                    piModel.initialConditionsResiduesTkc[resId][0].Get());
                ASSERT_DOUBLE_EQ(
                    tkc[1],  // k
                    piModel.initialConditionsResiduesTkc[resId][1].Get());
                ASSERT_DOUBLE_EQ(
                    tkc[2],  // c
                    piModel.initialConditionsResiduesTkc[resId][2].Get());
            }
        }
    }

    // Check initial displacement residues
    std::vector<Polys> residues = std::vector<Polys>();
    std::vector<double> initial = std::vector<double>{0.0, initialDisplacement};
    for (int massId = 0; massId < 2; massId++) {
        for (int i = 0; i < int(piModel.initialConditionsResiduesTkc.size());
             i++) {
            Poly p = piModel.models(massId, 0);
            p.SetX(Bounded::Get(piModel.initialConditionsResiduesTkc[i]));

            residues.push_back(p + (-initial[massId]));
        }
    }
    ASSERT_EQ(residues.size(), piModel.initialDispResidues.size());
    for (int i = 0; i < int(residues.size()); i++) {
        ASSERT_EQ(residues[i], piModel.initialDispResidues[i]);
    }

    // Check initial velocity residues
    residues = std::vector<Polys>();
    initial = std::vector<double>{0.0, initialVelocity};
    for (int massId = 0; massId < 2; massId++) {
        for (int i = 0; i < int(piModel.initialConditionsResiduesTkc.size());
             i++) {
            Poly p = piModel.models(massId, 0);
            p.SetX(Bounded::Get(piModel.initialConditionsResiduesTkc[i]));

            ASSERT_FALSE(p.Dxi(0).isError);
            p *= 1 / (TMax - TMin);

            residues.push_back(p + (-initial[massId]));
        }
    }
    ASSERT_EQ(residues.size(), piModel.initialVelResidues.size());
    for (int i = 0; i < int(residues.size()); i++) {
        ASSERT_EQ(residues[i], piModel.initialVelResidues[i]);
    }
}

TEST_F(PimodelTest, PhysicsResiduesTest) {
    // For each mass, one residue for each tkc value
    // Ex:
    // 2 masses
    // timeDiscretization = 1
    // kcDiscretization = 0
    // -> tkc values = [0,0.5,0.5], [1,0.5,0.5]
    // physicsResidues:
    //  (x0ModelDotDot-x0DotDot(x0Model(t = 0, k = 0.5, c = 0.5),
    //                          x0ModelDot(t = 0, k = 0.5, c = 0.5))
    //
    //  (x0ModelDotDot-x0DotDot(x0Model(t = 1, k = 0.5, c = 0.5),
    //                          x0ModelDot(t = 1, k = 0.5, c = 0.5))
    //
    //  (x1ModelDotDot-x1DotDot(x1Model(t = 0, k = 0.5, c = 0.5),
    //                          x1ModelDot(t = 0, k = 0.5, c = 0.5))
    //
    //  (x1ModelDotDot-x1DotDot(x1Model(t = 1, k = 0.5, c = 0.5),
    //                          x1ModelDot(t = 1, k = 0.5, c = 0.5))

    ASSERT_EQ(piModel.physicsResidues.size(),
              piModel.physicsResiduesTkc.size() * 2);  // 2 masses

    // Verify all the X values are set correctly in the models
    for (int resId = 0; resId < int(piModel.physicsResiduesTkc.size());
         resId++) {
        for (int massId = 0; massId < 2; massId++) {
            for (int p = 0;
                 p < int(piModel.physicsResidues[resId].polys.size()); p++) {
                auto tkc = piModel.physicsResidues[resId].polys[p].GetX();
                ASSERT_DOUBLE_EQ(tkc[0],  // t
                                 piModel.physicsResiduesTkc[resId][0].Get());
                ASSERT_DOUBLE_EQ(tkc[1],  // k
                                 piModel.physicsResiduesTkc[resId][1].Get());
                ASSERT_DOUBLE_EQ(tkc[2],  // c
                                 piModel.physicsResiduesTkc[resId][2].Get());
            }
        }
    }

    std::vector<Polys> residues = std::vector<Polys>();

    for (int i = 0; i < int(piModel.physicsResiduesTkc.size()); i++) {
        Poly x0Model = piModel.models(0, 0);

        Poly x0ModelDot = piModel.models(0, 0);
        ASSERT_FALSE(x0ModelDot.Dxi(0).isError);
        x0ModelDot *= 1 / (TMax - TMin);

        Poly x0ModelDotDot = piModel.models(0, 0);
        ASSERT_FALSE(x0ModelDotDot.Dxi(0).isError);
        ASSERT_FALSE(x0ModelDotDot.Dxi(0).isError);
        x0ModelDotDot *= 1 / (TMax - TMin);
        x0ModelDotDot *= 1 / (TMax - TMin);

        x0Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDotDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        residues.push_back(x0ModelDotDot + (-0));  // Initial mass is fixed
    }

    Bounded k;  // normalized k (from 0 to 1)
    Bounded c;  // normalized c (from 0 to 1)
    double K;   // actual value of K
    double C;   // actual value of C
    for (int i = 0; i < int(piModel.physicsResiduesTkc.size()); i++) {
        Poly x0Model = piModel.models(0, 0);

        Poly x0ModelDot = piModel.models(0, 0);
        ASSERT_FALSE(x0ModelDot.Dxi(0).isError);
        x0ModelDot *= 1 / (TMax - TMin);

        x0Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        Poly x1Model = piModel.models(1, 0);

        Poly x1ModelDot = piModel.models(1, 0);
        ASSERT_FALSE(x1ModelDot.Dxi(0).isError);
        x1ModelDot *= 1 / (TMax - TMin);

        Poly x1ModelDotDot = piModel.models(1, 0);
        ASSERT_FALSE(x1ModelDotDot.Dxi(0).isError);
        ASSERT_FALSE(x1ModelDotDot.Dxi(0).isError);
        x1ModelDotDot *= 1 / (TMax - TMin);
        x1ModelDotDot *= 1 / (TMax - TMin);

        x1Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x1ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x1ModelDotDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        // x1DotDot(x0Model, x1Model) =
        //  1/m1*
        //      ( K*x0Model -  K*x1Model
        //       + C*x0DotModel -  C*x1DotModel)
        k = piModel.physicsResiduesTkc[i][1];
        c = piModel.physicsResiduesTkc[i][2];
        K = Unnormalize(k, KMin, KMax);
        C = Unnormalize(c, CMin, CMax);
        residues.push_back(Polys(x1ModelDotDot) +
                           (-1 / m) * (K * x0Model + (-K) * x1Model +
                                       C * x0ModelDot + (-C) * x1ModelDot));
    }

    ASSERT_EQ(residues.size(), piModel.physicsResidues.size());
    for (int i = 0; i < int(residues.size()); i++) {
        ASSERT_EQ(piModel.physicsResidues[i], residues[i]);
    }
}

TEST_F(PimodelTest, nResiduesTest) {
    auto pd_ = ProblemDescription();
    pd_.AddMass(m, 0.0, 0.0);
    pd_.AddMass(m, 1.0, 0.0);
    pd_.AddSpring(0, 1, KMin, KMax);
    pd_.AddDamper(0, 1, CMin, CMax);

    auto model_ = Pimodel(pd, TMin, TMax, 1, 1, 1);
    // Before calling AddResidues, the function returns 0
    ASSERT_EQ(model_.nResidues(), 0);

    model_.SetResidues(true, true);
    // Residues:
    // Initial conditions:
    // KcDisc. = 1 -> residues will be evaluated for [k=0,k=1]X[c=0,c=1]
    // for k in (0,1):
    //   for c in (0,1):
    //     x0Model(t=0,k,c) - x0_t=0
    //     x1Model(t=0,k,c) - x1_t=0
    //     x0ModelDot(t=0,k,c) - x0Dot_t=0
    //     x1ModelDot(t=0,k,c) - x1Dot_t=0
    // -> 2*2*4 = 16 residues
    //
    // TimeDiscretization = 1 -> residues will be evaluated for t = 0 and t=1
    // KcDisc. = 1 -> residues will be evaluated for [k=0,k=1]X[c=0,c=1]
    // Physics:
    // for t in (0, 1):
    //   for k in (0, 1):
    //     for c in (0, 1):
    //       x0ModelDotDot(t,k,c) - x0DotDot(x0Model(t,k,c),
    //       x0ModelDot(t,k,c))
    //       x1ModelDotDot(t,k,c) - x1DotDot(x1Model(t,k,c),
    //       x1ModelDot(t,k,c))
    // -> 2*2*2*2 = 16 residues
    ASSERT_EQ(model_.nResidues(), 16 + 16);

    model_.SetResidues(true, false);
    ASSERT_EQ(model_.nResidues(), 16);

    model_.SetResidues(false, true);
    ASSERT_EQ(model_.nResidues(), 16);
}

TEST_F(PimodelTest, LossTest) {
    // Basically a combination of InitialConditionsResiduesTest
    // and PhysicsResiduesTest

    int nInitialDispResidues = 0;
    int nInitialVelResidues = 0;
    double k = 0;
    double c = 0;
    while (k <= 1) {
        c = 0;
        while (c <= 1) {
            nInitialDispResidues++;
            nInitialVelResidues++;

            c += 1 / (kcDiscretization + 0.0);
        }
        k += 1 / (kcDiscretization + 0.0);
    }
    int nInitialCondResidues = nInitialDispResidues + nInitialVelResidues;

    int nPhysicsResidues = 0;
    double t = 0;
    while (t <= 1) {
        k = 0;
        while (k <= 1) {
            c = 0;
            while (c <= 1) {
                nPhysicsResidues++;
                c += 1 / (kcDiscretization + 0.0);
            }
            k += 1 / (kcDiscretization + 0.0);
        }
        t += 1 / (kcDiscretization + 0.0);
    }
    double totalResiduesTerms =
        nInitialDispResidues + nInitialDispResidues + nPhysicsResidues;
    double initialConditionsResiduesWeight =
        nPhysicsResidues / totalResiduesTerms;
    double physicsResiduesWeight = nInitialCondResidues / totalResiduesTerms;

    double expectedLoss = 0.0;
    std::vector<double> initial = std::vector<double>{0.0, initialDisplacement};
    for (int massId = 0; massId < 2; massId++) {
        for (int i = 0; i < int(piModel.initialConditionsResiduesTkc.size());
             i++) {
            Poly xModel = piModel.models(massId, 0);
            xModel.SetX(Bounded::Get(piModel.initialConditionsResiduesTkc[i]));

            Polys residue = xModel + (-initial[massId]);
            expectedLoss += initialConditionsResiduesWeight *
                            pow(residue(piModel.modelsCoefficients).val, 2);
        }
    }
    initial = std::vector<double>{0.0, initialVelocity};
    for (int massId = 0; massId < 2; massId++) {
        for (int i = 0; i < int(piModel.initialConditionsResiduesTkc.size());
             i++) {
            Poly xModelDot = piModel.models(massId, 0);
            ASSERT_FALSE(xModelDot.Dxi(0).isError);
            xModelDot *= 1 / (TMax - TMin);
            xModelDot.SetX(
                Bounded::Get(piModel.initialConditionsResiduesTkc[i]));

            Polys residue = (xModelDot + (-initial[massId]));
            expectedLoss += initialConditionsResiduesWeight *
                            pow(residue(piModel.modelsCoefficients).val, 2);
        }
    }
    for (int i = 0; i < int(piModel.physicsResiduesTkc.size()); i++) {
        Poly x0Model = piModel.models(0, 0);
        Poly x0ModelDot = piModel.models(0, 0);

        ASSERT_FALSE(x0ModelDot.Dxi(0).isError);
        x0ModelDot *= 1 / (TMax - TMin);

        Poly x0ModelDotDot = piModel.models(0, 0);
        ASSERT_FALSE(x0ModelDotDot.Dxi(0).isError);
        ASSERT_FALSE(x0ModelDotDot.Dxi(0).isError);
        x0ModelDotDot *= 1 / (TMax - TMin);
        x0ModelDotDot *= 1 / (TMax - TMin);

        x0Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDotDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        Polys residue = (x0ModelDotDot + (-0));  // Initial mass is fixed
        expectedLoss += physicsResiduesWeight *
                        pow(residue(piModel.modelsCoefficients).val, 2);
    }
    Bounded kb;  // normalized k (from 0 to 1)
    Bounded cb;  // normalized c (from 0 to 1)
    double K;    // actual value of K
    double C;    // actual value of C
    for (int i = 0; i < int(piModel.physicsResiduesTkc.size()); i++) {
        Poly x0Model = piModel.models(0, 0);
        Poly x0ModelDot = piModel.models(0, 0);

        ASSERT_FALSE(x0ModelDot.Dxi(0).isError);
        x0ModelDot *= 1 / (TMax - TMin);

        x0Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x0ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        Poly x1Model = piModel.models(1, 0);
        Poly x1ModelDot = piModel.models(1, 0);
        ASSERT_FALSE(x1ModelDot.Dxi(0).isError);
        x1ModelDot *= 1 / (TMax - TMin);

        Poly x1ModelDotDot = piModel.models(1, 0);
        ASSERT_FALSE(x1ModelDotDot.Dxi(0).isError);
        ASSERT_FALSE(x1ModelDotDot.Dxi(0).isError);
        x1ModelDotDot *= 1 / (TMax - TMin);
        x1ModelDotDot *= 1 / (TMax - TMin);

        x1Model.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x1ModelDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));
        x1ModelDotDot.SetX(Bounded::Get(piModel.physicsResiduesTkc[i]));

        // x1DotDot(x0Model, x1Model) =
        //  1/m1*
        //      ( K*x0Model -  K*x1Model
        //       + C*x0DotModel -  C*x1DotModel)
        kb = piModel.physicsResiduesTkc[i][1];
        cb = piModel.physicsResiduesTkc[i][2];
        K = Unnormalize(kb, KMin, KMax);
        C = Unnormalize(cb, CMin, CMax);

        Polys residue = (Polys(x1ModelDotDot) +
                         (-1 / m) * (K * x0Model + (-K) * x1Model +
                                     C * x0ModelDot + (-C) * x1ModelDot));
        expectedLoss += physicsResiduesWeight *
                        pow(residue(piModel.modelsCoefficients).val, 2);
    }

    ASSERT_DOUBLE_EQ(piModel.Loss(), expectedLoss);
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
    int batchSize;
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
        this->learningRate = 0.1;
        this->batchSize = 1000;  // huge batch size -> single batch is used
        this->maxSteps = 100;
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
        this->model.SetResidues(true, true);

        double initialLoss = this->model.Loss();
        this->model.Train(this->learningRate, this->batchSize, this->maxSteps,
                          this->log);
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
    ASSERT_DEATH({ Pimodels(this->pd, 1.0, 1, 1, -1, 1); }, "");
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
    double totalT = 100 + Random();
    double K = Random(KMin, KMax);
    double C = Random(CMin, CMax);
    std::vector<double> TKC = std::vector<double>{0, K, C};

    // Single Pimodel
    Pimodels pimodels = Pimodels(this->pd, totalT, 1, 1, 1, 1);
    ASSERT_DEATH({ pimodels.setContinuity(0, TKC); }, "");

    // Two Pimodels
    // -> t in [0, totalT/2)
    // -> t in [totalT/2, totalT]
    pimodels = Pimodels(this->pd, totalT, 2, 1, 1, 1);
    std::vector<double> params0 =
        std::vector<double>(pimodels.pimodels[0].nParameters());
    for (int i = 0; i < int(params0.size()); i++) {
        params0[i] = Random();
    }
    std::vector<double> params1 =
        std::vector<double>(pimodels.pimodels[1].nParameters());
    for (int i = 0; i < int(params1.size()); i++) {
        params1[i] = Random();
    }
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&params0).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&params1).isError);
    ASSERT_DEATH({ pimodels.setContinuity(0, TKC); }, "");
    pimodels.setContinuity(1, TKC);

    // No initial conditions are added to the first pimodel
    ASSERT_EQ(pimodels.pimodels[0].p.initialDisps.size(), 0);
    ASSERT_EQ(pimodels.pimodels[0].p.initialVels.size(), 0);

    // One initial condition (disp and vel) is added for each mass
    // for the second pimodel
    TKC[0] = totalT / 2;
    ASSERT_EQ(pimodels.pimodels[1].p.initialDisps.size(), 2);
    ASSERT_EQ(pimodels.pimodels[1].p.initialVels.size(), 2);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].val,
                     pimodels.pimodels[0](&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].val,
                     pimodels.pimodels[0](&TKC).val[1]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].val,
                     pimodels.pimodels[0].GetVelocities(&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].val,
                     pimodels.pimodels[0].GetVelocities(&TKC).val[1]);

    // Three Pimodels
    // t in [0, totalT/3)
    // t in [totalT/3, 2*totalT/3)
    // t in [2*totalT/3, totalT]

    pimodels = Pimodels(this->pd, totalT, 3, 1, 1, 1);
    params0 = std::vector<double>(pimodels.pimodels[0].nParameters());
    for (int i = 0; i < int(params0.size()); i++) {
        params0[i] = Random();
    }
    params1 = std::vector<double>(pimodels.pimodels[1].nParameters());
    for (int i = 0; i < int(params1.size()); i++) {
        params1[i] = Random();
    }
    auto params2 = std::vector<double>(pimodels.pimodels[2].nParameters());
    for (int i = 0; i < int(params1.size()); i++) {
        params2[i] = Random();
    }
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&params0).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&params1).isError);
    ASSERT_FALSE(pimodels.pimodels[2].SetParameters(&params2).isError);
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
    TKC[0] = totalT / 3;
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[0].val,
                     pimodels.pimodels[0](&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialDisps[1].val,
                     pimodels.pimodels[0](&TKC).val[1]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[0].val,
                     pimodels.pimodels[0].GetVelocities(&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[1].p.initialVels[1].val,
                     pimodels.pimodels[0].GetVelocities(&TKC).val[1]);

    // Pimodel 2
    TKC[0] = 2 * totalT / 3;
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[0].val,
                     pimodels.pimodels[1](&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialDisps[1].val,
                     pimodels.pimodels[1](&TKC).val[1]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[0].massId, 0);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[0].val,
                     pimodels.pimodels[1].GetVelocities(&TKC).val[0]);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[1].massId, 1);
    ASSERT_DOUBLE_EQ(pimodels.pimodels[2].p.initialVels[1].val,
                     pimodels.pimodels[1].GetVelocities(&TKC).val[1]);
}

TEST_F(PimodelsTest, OperatorTest) {
    double totalT = 100 + Random();
    double K = Random(KMin, KMax);
    double C = Random(CMin, CMax);
    std::vector<double> TKC = std::vector<double>{0, K, C};

    // Three Pimodels
    // t in [0, tMax/3):
    //  x0(t,k,c) = p0_0
    //  x1(t,k,c) = p1_0
    // t in [tMax/3, 2*tMax/3):
    //  x0(t,k,c) = p0_1
    //  x1(t,k,c) = p1_1
    // t in [2*tMax/3, tMax]:
    //  x0(t,k,c) = p0_2
    //  x1(t,k,c) = p1_2
    Pimodels pimodels = Pimodels(this->pd, totalT, 3, 1, 1, 1);
    auto params0 = std::vector<double>(pimodels.pimodels[0].nParameters());
    for (int i = 0; i < int(params0.size()); i++) {
        params0[i] = Random();
    }
    auto params1 = std::vector<double>(pimodels.pimodels[1].nParameters());
    for (int i = 0; i < int(params1.size()); i++) {
        params1[i] = Random();
    }
    auto params2 = std::vector<double>(pimodels.pimodels[2].nParameters());
    for (int i = 0; i < int(params1.size()); i++) {
        params2[i] = Random();
    }
    ASSERT_FALSE(pimodels.pimodels[0].SetParameters(&params0).isError);
    ASSERT_FALSE(pimodels.pimodels[1].SetParameters(&params1).isError);
    ASSERT_FALSE(pimodels.pimodels[2].SetParameters(&params2).isError);

    double T;
    // Pimodel 0
    T = 0;
    TKC[0] = T;
    auto eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, pimodels.pimodels[0](&TKC).val);

    // Pimodel 1
    T = totalT / 3;
    TKC[0] = T;
    eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, pimodels.pimodels[1](&TKC).val);

    // Pimodel 2
    T = 2 * totalT / 3;
    TKC[0] = T;
    eval = pimodels(&TKC);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, pimodels.pimodels[2](&TKC).val);
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
    double learningRate = 1.0;
    int batchSize = 1000;  // huge batch size -> single batch is used
    int maxSteps = 300;
    bool log = false;

    // Train all models
    Pimodels models = Pimodels(pd, tMax, nModels, timeDiscretization,
                               kcDiscretization, order);
    models.Train(learningRate, batchSize, maxSteps, log);

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

TEST(Debug, debugtest) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 1.0, 0.0);  // m2
    pd.AddMass(150, 1.0, 3.0);  // m3
    pd.AddMass(700, 2.0, 0.0);  // m4
    pd.AddMass(80, 3.0, 0.0);   // m5

    double min = 100.0;
    double max = 100000;
    pd.AddSpring(0, 1, min, max);  // k01
    pd.AddSpring(1, 2, min, max);  // k12
    pd.AddSpring(1, 3, min, max);  // k13
    pd.AddSpring(1, 4, min, max);  // k14
    pd.AddDamper(1, 4, min, max);  // c14
    pd.AddSpring(0, 2, min, max);  // k02
    pd.AddDamper(0, 2, min, max);  // c02
    pd.AddSpring(2, 4, min, max);  // k24
    pd.AddDamper(2, 4, min, max);  // c24
    pd.AddSpring(0, 3, min, max);  // k03
    pd.AddDamper(0, 3, min, max);  // c03
    pd.AddSpring(3, 4, min, max);  // k34
    pd.AddDamper(3, 4, min, max);  // c34
    pd.AddSpring(4, 5, min, max);  // k45
    pd.AddDamper(4, 5, min, max);  // c45
    pd.SetFixedMass(0);
    pd.AddInitialVel(200.0);  // initial speed
    assert(pd.IsOk());

    // Learning parameters
    double finalT = 0.05;
    int nModels = 3;
    int timeDiscretization = 1;
    int kcDiscretization = 0;
    int order = 3;
    double learningRate = 0.1;
    int batchSize = 1000;
    int maxSteps = 10;
    bool log = true;

    // Train all models
    Pimodels models = Pimodels(pd, finalT, nModels, timeDiscretization,
                               kcDiscretization, order);
    assert(!models.Train(learningRate, batchSize, maxSteps, log).isError);

    // Get problem using intermediate value for k and c, and integrate it.
    double mean = (min + max) / 2;
    Maybe<Problem> mP = pd.BuildFromVector(
        std::vector<double>{mean, mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean});
    assert(!mP.isError);
    Problem p = mP.val;
    assert(!p.Integrate(finalT).isError);

    std::vector<double> tkc =
        std::vector<double>{0.0,  mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean, mean};
    Maybe<std::vector<double>> X;
    Maybe<std::vector<double>> XDot;
    std::cout << "t,x5,x5Dot,modelX5,modelX5Dot" << std::endl;
    for (int i = 0; i < int(p.t.size()); i += 1) {
        tkc[0] = p.t[i];
        X = models(&tkc);
        XDot = models.GetVelocities(&tkc);

        std::cout << p.t[i] << ",";                                // t,
        std::cout << p.XHistory[i][p.GetMassDispIndex(5)] << ",";  // x5,
        std::cout << p.XHistory[i][p.GetMassVelIndex(5)] << ",";   // x5Dot,
        std::cout << X.val[5] << ",";                              // modelX5,
        std::cout << XDot.val[5] << std::endl;                     // modelX5Dot
    }
}