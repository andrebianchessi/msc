#include "pimodel.h"

#include <gtest/gtest.h>

#include <vector>

#include "bounded.h"
#include "utils.h"

class PimodelTest : public testing::Test {
   public:
    ProblemDescription pd;

    void SetUp() {
        // Called before every TEST_F
        // Create problem description of two masses,
        // connected with a spring and damper. One is fixed
        // and the other has initial displacement.
        this->pd = ProblemDescription();
        this->pd.AddMass(1.0, 0.0, 0.0);
        this->pd.AddMass(1.0, 1.0, 0.0);
        this->pd.AddSpring(0, 1, 0.5, 1.0);
        this->pd.AddDamper(0, 1, 10, 50);
        this->pd.SetFixedMass(0);
        this->pd.AddInitialDisp(1, 10.0);
        std::vector<Bounded> dna = std::vector<Bounded>(2);
        auto e0 = this->pd.BuildFromDNA(dna);
        ASSERT_FALSE(e0.isError);
    }
};

TEST_F(PimodelTest, ConstructorTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 2);

    // There should be one polynomial for each mass
    ASSERT_EQ(model.polys.size(), 2);
}

TEST_F(PimodelTest, OperatorTest) {
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 2);

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
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 1);

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
    Pimodel model = Pimodel(&this->pd, 1.0, 10, 1);

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

TEST_F(PimodelTest, LossTest) {
    // The polynomial of each mass will be of order 1:
    // x0(t,k,c) = a0*t + a1*k + a2*c + a3*1
    // x1(t,k,c) = b0*t + b1*k + b2*c + b3*1
    Pimodel model = Pimodel(&this->pd, 1.0, 2, 1);

    double loss = model.Loss();
}