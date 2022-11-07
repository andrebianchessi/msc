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
    ASSERT_EQ(lr.val.nTerms(), 1);

    // y(x1, ... , x10) = k
    lr = LinReg::NewLinReg(10, 0);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms(), 1);

    // y(x1, ... , x10) = a*x1 + ... + j*x10 + k
    lr = LinReg::NewLinReg(10, 1);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms(), 11);

    // y(x1,x2) = x1^2 + x1*x2 + x2^2 + x1 + x2 + k
    lr = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms(), 6);

    // y(x1,x2) = x1^3 + x1^2*x2 + x1*x2^2 + x2^3 + x1^2 + x1*x2 + x2^2 + x1 +
    // x2 + k
    lr = LinReg::NewLinReg(2, 3);
    ASSERT_FALSE(lr.isError);
    ASSERT_EQ(lr.val.nTerms(), 10);
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

TEST(LinRegTest, DTest) {
    Maybe<LinReg> m = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = lr.D(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);
    v = lr.D(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(lr.D(2, &X).isError);
    ASSERT_TRUE(lr.D(-1, &X).isError);

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

    v = lr.D(0, &X);
    ASSERT_FALSE(v.isError);
    double d1Expected = 1.0 * 2 * x1 + 2.0 * x2 + 3.0;
    ASSERT_DOUBLE_EQ(v.val, d1Expected);

    v = lr.D(1, &X);
    ASSERT_FALSE(v.isError);
    double d2Expected = 2.0 * x1 + 4.0 * 2 * x2 + 5.0;
    ASSERT_DOUBLE_EQ(v.val, d2Expected);
}

TEST(LinRegTest, D2Test) {
    Maybe<LinReg> m = LinReg::NewLinReg(2, 2);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = lr.D2(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);
    v = lr.D2(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(lr.D2(2, &X).isError);
    ASSERT_TRUE(lr.D2(-1, &X).isError);

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

    v = lr.D2(0, &X);
    ASSERT_FALSE(v.isError);
    double dd1Expected = 1.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd1Expected);

    v = lr.D2(1, &X);
    ASSERT_FALSE(v.isError);
    double dd2Expected = 4.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd2Expected);
}

TEST(LinRegTest, DAndD2Test) {
    // Similar to DTest and D2Test, but with a more complex polynomial
    Maybe<LinReg> m = LinReg::NewLinReg(3, 3);
    ASSERT_FALSE(m.isError);
    LinReg lr = m.val;

    std::vector<double> X = {1.0, 2.0, 3.0};

    // I only setup a few parameters, for simplicity

    // x1^3
    lr.coefficients->children[0]->children[0]->children[0]->a = 5.0;

    // x2^3
    lr.coefficients->children[3]->children[0]->children[0]->a = 7.0;

    // x3^3
    lr.coefficients->children[3]->children[3]->children[0]->a = 9.0;

    // x1^2*x2
    lr.coefficients->children[1]->children[0]->children[0]->a = 13.0;

    // x2^2*x3
    lr.coefficients->children[3]->children[1]->children[0]->a = 17.0;

    // x1*x3^2
    lr.coefficients->children[2]->children[2]->children[0]->a = 19.0;

    // y(x1,x2,x3) = 5*x1^3 + 7*x2^3 + 9*x3^3 + 13*x1^2*x2 + 17*x2^2*x3 +
    // 19*x1*x3^2
    double x1 = 9.0;
    double x2 = 17.0;
    double x3 = 21.0;
    X = {x1, x2, x3};

    auto d = lr.D(0, &X);
    ASSERT_FALSE(d.isError);
    double d1Expected = 3 * 5 * x1 * x1 + 2 * 13 * x1 * x2 + 19 * x3 * x3;
    ASSERT_DOUBLE_EQ(d.val, d1Expected);

    d = lr.D(1, &X);
    ASSERT_FALSE(d.isError);
    double d2Expected = 3 * 7 * x2 * x2 + 13 * x1 * x1 + 2 * 17 * x2 * x3;
    ASSERT_DOUBLE_EQ(d.val, d2Expected);

    d = lr.D(2, &X);
    ASSERT_FALSE(d.isError);
    double d3Expected = 3 * 9 * x3 * x3 + 17 * x2 * x2 + 2 * 19 * x1 * x3;
    ASSERT_DOUBLE_EQ(d.val, d3Expected);

    auto dd = lr.D2(0, &X);
    ASSERT_FALSE(dd.isError);
    double dd1Expected = 2 * 3 * 5 * x1 + 2 * 13 * x2;
    ASSERT_DOUBLE_EQ(dd.val, dd1Expected);

    dd = lr.D2(1, &X);
    ASSERT_FALSE(dd.isError);
    double dd2Expected = 2 * 3 * 7 * x2 + 2 * 17 * x3;
    ASSERT_DOUBLE_EQ(dd.val, dd2Expected);

    dd = lr.D2(2, &X);
    ASSERT_FALSE(dd.isError);
    double dd3Expected = 2 * 3 * 9 * x3 + 2 * 19 * x1;
    ASSERT_DOUBLE_EQ(dd.val, dd3Expected);
}