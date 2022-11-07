#include "linear_regression.h"

#include <gtest/gtest.h>

#include "maybe.h"

TEST(LinRegTest, ConstructorTest) {
    // Invalid values
    Maybe<LinReg> lr = LinReg::NewLinReg(-1, 0);
    ASSERT_TRUE(lr.isError);

    lr = LinReg::NewLinReg(0, -1);
    ASSERT_TRUE(lr.isError);

    lr = LinReg::NewLinReg(-10, -10);
    ASSERT_TRUE(lr.isError);

    // y() = a
    lr = LinReg::NewLinReg(0, 0);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms, 1);

    // y(x1, ... , x10) = k
    lr = LinReg::NewLinReg(10, 0);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms, 1);

    // y(x1, ... , x10) = a*x1 + ... + j*x10 + k
    lr = LinReg::NewLinReg(10, 1);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms, 11);

    // y(x1,x2) = x1^2 + x1*x2 + x2^2 + x1 + x2 + k
    lr = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms, 6);

    // y(x1,x2) = x1^3 + x1^2*x2 + x1*x2^2 + x2^3 + x1^2 + x1*x2 + x2^2 + x1 +
    // x2 + k
    lr = LinReg::NewLinReg(2, 3);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms, 10);
}

TEST(LinRegTest, coefficientsTest) {
    // y(x1,x2) = x1^3 + x1^2*x2 + x1*x2^2 + x2^3 + x1^2 + x1*x2 + x2^2 + x1 +
    // x2 + k

    // Root
    // children: x1^3 x1^2 x1^1 x1^0

    // x1^3
    // children: x2^0

    // x1^2
    // children: x2^1 x2^0

    // x1
    // children: x2^2 x2^1 x2^0

    // x1^0
    // children: x2^3 x2^2 x2^1 x2^0

    Maybe<LinReg> m = LinReg::NewLinReg(2, 3);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    ASSERT_EQ(lr.coefficients->children.size(), 4);

    ASSERT_EQ(lr.coefficients->children[0]->e, 3);
    ASSERT_EQ(lr.coefficients->children[1]->e, 2);
    ASSERT_EQ(lr.coefficients->children[2]->e, 1);
    ASSERT_EQ(lr.coefficients->children[3]->e, 0);

    ASSERT_EQ(lr.coefficients->children[0]->children.size(), 1);
    ASSERT_EQ(lr.coefficients->children[0]->children[0]->e, 0);
    ASSERT_EQ(lr.coefficients->children[0]->children[0]->a, 0.0);

    ASSERT_EQ(lr.coefficients->children[1]->children.size(), 2);
    ASSERT_EQ(lr.coefficients->children[1]->children[0]->e, 1);
    ASSERT_EQ(lr.coefficients->children[1]->children[0]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[1]->children[1]->e, 0);
    ASSERT_EQ(lr.coefficients->children[1]->children[1]->a, 0.0);

    ASSERT_EQ(lr.coefficients->children[2]->children.size(), 3);
    ASSERT_EQ(lr.coefficients->children[2]->children[0]->e, 2);
    ASSERT_EQ(lr.coefficients->children[2]->children[0]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[2]->children[1]->e, 1);
    ASSERT_EQ(lr.coefficients->children[2]->children[1]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[2]->children[2]->e, 0);
    ASSERT_EQ(lr.coefficients->children[2]->children[2]->a, 0.0);

    ASSERT_EQ(lr.coefficients->children[3]->children.size(), 4);
    ASSERT_EQ(lr.coefficients->children[3]->children[0]->e, 3);
    ASSERT_EQ(lr.coefficients->children[3]->children[0]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[3]->children[1]->e, 2);
    ASSERT_EQ(lr.coefficients->children[3]->children[1]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[3]->children[2]->e, 1);
    ASSERT_EQ(lr.coefficients->children[3]->children[2]->a, 0.0);
    ASSERT_EQ(lr.coefficients->children[3]->children[3]->e, 0);
    ASSERT_EQ(lr.coefficients->children[3]->children[3]->a, 0.0);
}

TEST(LinRegTest, operatorTest) {
    Maybe<LinReg> m = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = lr(&X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    // x1^2*x2^0
    lr.coefficients->children[0]->children[0]->a = 1.0;

    // x1^1*x2^1
    lr.coefficients->children[1]->children[0]->a = 2.0;

    // x1^1*x2^0
    lr.coefficients->children[1]->children[1]->a = 3.0;

    // x1^0*x2^2
    lr.coefficients->children[2]->children[0]->a = 4.0;

    // x1^0*x2^1
    lr.coefficients->children[2]->children[1]->a = 5.0;

    // x1^0*x2^0
    lr.coefficients->children[2]->children[2]->a = 6.0;

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    auto dfsOutput = dfs(lr.coefficients->children[0]->children[0], 1, &X);
    ASSERT_DOUBLE_EQ(dfsOutput, 1.0);

    dfsOutput = dfs(lr.coefficients->children[0], 0, &X);
    ASSERT_DOUBLE_EQ(dfsOutput, x1 * x1);

    v = lr(&X);
    ASSERT_FALSE(v.isError);

    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_DOUBLE_EQ(v.val, expected);
}

TEST(LinRegTest, coefficientAtTest) {
    Maybe<LinReg> m = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    std::vector<int> powers = {0, 0, 0};
    ASSERT_TRUE(lr.CoefficientAt(powers).isError);

    powers = {-1, 0};
    ASSERT_TRUE(lr.CoefficientAt(powers).isError);
    powers = {0, -1};
    ASSERT_TRUE(lr.CoefficientAt(powers).isError);

    powers = {0, 0};
    Maybe<int> c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 0.0);

    // x1^2*x2^0
    lr.coefficients->children[0]->children[0]->a = 1.0;

    // x1^1*x2^1
    lr.coefficients->children[1]->children[0]->a = 2.0;

    // x1^1*x2^0
    lr.coefficients->children[1]->children[1]->a = 3.0;

    // x1^0*x2^2
    lr.coefficients->children[2]->children[0]->a = 4.0;

    // x1^0*x2^1
    lr.coefficients->children[2]->children[1]->a = 5.0;

    // x1^0*x2^0
    lr.coefficients->children[2]->children[2]->a = 6.0;

    powers = {2, 0};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 1.0);

    powers = {1, 1};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 2.0);

    powers = {1, 0};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 3.0);

    powers = {0, 2};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 4.0);

    powers = {0, 1};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 5.0);

    powers = {0, 0};
    c = lr.CoefficientAt(powers);
    ASSERT_FALSE(c.isError);
    ASSERT_DOUBLE_EQ(c.val, 6.0);
}