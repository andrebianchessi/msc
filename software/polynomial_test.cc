#include "polynomial.h"

#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operations.hpp>

#include "maybe.h"
#include "utils.h"

TEST(MonomialTest, constructorTest) {
    // Empty monomial
    Monomial m = Monomial(std::vector<int>{});

    EXPECT_DOUBLE_EQ(m.a, 0);
    EXPECT_EQ(m.exps.size(), 0);

    // 0*x0^1*x1^0*x2^3
    m = Monomial(std::vector<int>{1, 0, 3});

    EXPECT_DOUBLE_EQ(m.a, 0);
    ASSERT_EQ(m.exps.size(), 3);
    EXPECT_EQ(m.exps[0], 1);
    EXPECT_EQ(m.exps[1], 0);
    EXPECT_EQ(m.exps[2], 3);
}

TEST(MonomialTest, evalTest) {
    // 1*1*x0^1*x1^0*x2^3
    Monomial m = Monomial(std::vector<int>{1, 0, 3});
    double a = 932847.1;
    m.a = a;

    double x0 = 1234;
    double x1 = 56234;
    double x2 = 323.2;
    auto X = std::vector<double>{x0, x1, x2};
    Maybe<double> eval = m(X);
    ASSERT_FALSE(eval.isError);
    EXPECT_DOUBLE_EQ(eval.val, a * x0 * 1 * x2 * x2 * x2);
}

TEST(PolyConstructorTest, constructorTest) {
    // Invalid values
    Poly p;
    Maybe<Void> r = p.Build(-1, 0);
    ASSERT_TRUE(r.isError);

    r = p.Build(0, -1);
    ASSERT_TRUE(r.isError);

    r = p.Build(-10, -10);
    ASSERT_TRUE(r.isError);

    // y() = 0*1
    r = p.Build(0, 0);
    ASSERT_FALSE(r.isError);
    ASSERT_EQ(p.monomials.size(), 1);
    ASSERT_EQ(p.monomials[0].a, 0);
    ASSERT_EQ(p.monomials[0].exps.size(), 0);

    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0*1
    r = p.Build(2, 2);
    ASSERT_FALSE(r.isError);
    ASSERT_EQ(p.monomials.size(), 6);
    // 0*x^2
    ASSERT_EQ(p.monomials[0].a, 0);
    ASSERT_EQ(p.monomials[0].exps.size(), 2);
    ASSERT_EQ(p.monomials[0].exps[0], 2);
    ASSERT_EQ(p.monomials[0].exps[1], 0);
    // 0*xy
    ASSERT_EQ(p.monomials[1].a, 0);
    ASSERT_EQ(p.monomials[1].exps.size(), 2);
    ASSERT_EQ(p.monomials[1].exps[0], 1);
    ASSERT_EQ(p.monomials[1].exps[1], 1);
    // 0*x
    ASSERT_EQ(p.monomials[2].a, 0);
    ASSERT_EQ(p.monomials[2].exps.size(), 2);
    ASSERT_EQ(p.monomials[2].exps[0], 1);
    ASSERT_EQ(p.monomials[2].exps[1], 0);
    // 0*y^2
    ASSERT_EQ(p.monomials[3].a, 0);
    ASSERT_EQ(p.monomials[3].exps.size(), 2);
    ASSERT_EQ(p.monomials[3].exps[0], 0);
    ASSERT_EQ(p.monomials[3].exps[1], 2);
    // 0*y
    ASSERT_EQ(p.monomials[4].a, 0);
    ASSERT_EQ(p.monomials[4].exps.size(), 2);
    ASSERT_EQ(p.monomials[4].exps[0], 0);
    ASSERT_EQ(p.monomials[4].exps[1], 1);
    // 0*1
    ASSERT_EQ(p.monomials[5].a, 0);
    ASSERT_EQ(p.monomials[5].exps.size(), 2);
    ASSERT_EQ(p.monomials[5].exps[0], 0);
    ASSERT_EQ(p.monomials[5].exps[1], 0);
}

class PolyTest : public testing::Test {
   public:
    // Polynomial of two vars or order two with all coeff. = 0
    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0
    Poly n2o2Zeros;
    // Polynomial of two vars or order two with all coeff. = 2.0
    // P = 2*x^2 + 2*xy + 2*x + 2*y^2 + 2*y + 2
    Poly n2o2Twos;

    // Polynomial of two vars or order two
    // P = 1x^2 + 2xy + 3x + 4y^2 + 5y + 6
    Poly n2o2;

    void SetUp() {
        // Called before every TEST_F
        Maybe<Void> r = n2o2Zeros.Build(2, 2, 0);
        ASSERT_FALSE(r.isError);
        n2o2Zeros.monomials[0].a = 0.0;
        n2o2Zeros.monomials[1].a = 0.0;
        n2o2Zeros.monomials[2].a = 0.0;
        n2o2Zeros.monomials[3].a = 0.0;
        n2o2Zeros.monomials[4].a = 0.0;
        n2o2Zeros.monomials[5].a = 0.0;

        r = n2o2Twos.Build(2, 2, 1);
        ASSERT_FALSE(r.isError);
        n2o2Twos.monomials[0].a = 2.0;
        n2o2Twos.monomials[1].a = 2.0;
        n2o2Twos.monomials[2].a = 2.0;
        n2o2Twos.monomials[3].a = 2.0;
        n2o2Twos.monomials[4].a = 2.0;
        n2o2Twos.monomials[5].a = 2.0;

        r = n2o2.Build(2, 2, 2);
        ASSERT_FALSE(r.isError);
        n2o2.monomials[0].a = 1.0;
        n2o2.monomials[1].a = 2.0;
        n2o2.monomials[2].a = 3.0;
        n2o2.monomials[3].a = 4.0;
        n2o2.monomials[4].a = 5.0;
        n2o2.monomials[5].a = 6.0;
    }
};

TEST_F(PolyTest, operatorTest) {
    std::vector<double> X;
    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    auto eval = n2o2(X);
    ASSERT_FALSE(eval.isError);

    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_DOUBLE_EQ(eval.val, expected);

    eval = n2o2(X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, expected);
}

TEST_F(PolyTest, copyConstructorAndAssignment) {
    Poly copyConstructed = n2o2;
    n2o2.monomials[0].a = 99;
    ASSERT_EQ(copyConstructed.monomials[0].a, 1.0);
    ASSERT_EQ(copyConstructed.monomials[1].a, 2.0);
    ASSERT_EQ(copyConstructed.monomials[2].a, 3.0);
    ASSERT_EQ(copyConstructed.monomials[3].a, 4.0);
    ASSERT_EQ(copyConstructed.monomials[4].a, 5.0);
    ASSERT_EQ(copyConstructed.monomials[5].a, 6.0);

    Poly assigned;
    assigned.Build(2, 3);  // Shouldn't change anything
    assigned = n2o2Twos;
    n2o2Twos.monomials[0].a = 99;
    ASSERT_EQ(assigned.monomials[0].a, 2);
    ASSERT_EQ(assigned.monomials[1].a, 2);
    ASSERT_EQ(assigned.monomials[2].a, 2);
    ASSERT_EQ(assigned.monomials[3].a, 2);
    ASSERT_EQ(assigned.monomials[4].a, 2);
    ASSERT_EQ(assigned.monomials[5].a, 2);
}

TEST_F(PolyTest, DxiTest) {
    double x = 987;
    double y = 91;
    std::vector<double> xy = {x, y};
    Maybe<Void> err;
    Maybe<double> eval;
    Poly p;

    p = n2o2Zeros;
    ASSERT_FALSE(p.Dxi(0).isError);
    eval = p(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    p = n2o2Zeros;
    ASSERT_FALSE(p.Dxi(1).isError);
    eval = p(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    p = n2o2Zeros;
    ASSERT_TRUE(p.Dxi(2).isError);
    ASSERT_TRUE(p.Dxi(-1).isError);

    p = n2o2;
    // P = 1x^2 + 2xy + 3x + 4y^2 + 5y + 6
    ASSERT_FALSE(p.Dxi(0).isError);
    eval = p(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 1.0 * 2 * x + 2.0 * y + 3.0);

    p = n2o2;
    ASSERT_FALSE(p.Dxi(1).isError);
    eval = p(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 2.0 * x + 4.0 * 2 * y + 5.0);
}

TEST_F(PolyTest, DxiTest2) {
    Maybe<Void> err;
    double x = 245;
    double y = 5763;
    auto xy = std::vector<double>{x, y};

    Poly p;
    p = this->n2o2;

    ASSERT_TRUE(p.Dxi(-1).isError);
    ASSERT_TRUE(p.Dxi(2).isError);

    // n2o2 = 1x^2 + 2xy + 3x + 4y^2 + 5y + 6
    // dn2o2/dx = 2x + 2y + 3 + 0 + 0 + 0
    err = p.Dxi(0);
    ASSERT_FALSE(err.isError);
    ASSERT_EQ(p.monomials[0].a, 2);
    ASSERT_EQ(p.monomials[1].a, 2);
    ASSERT_EQ(p.monomials[2].a, 3);
    ASSERT_EQ(p.monomials[3].a, 0);
    ASSERT_EQ(p.monomials[4].a, 0);
    ASSERT_EQ(p.monomials[5].a, 0);
    ASSERT_DOUBLE_EQ(p(xy).val, 2 * x + 2 * y + 3);

    // p = 2x + 2y + 3 + 0 + 0 + 0
    // dp/dy = 0 + 2 + 0 + 0 + 0 + 0
    err = p.Dxi(1);
    ASSERT_FALSE(err.isError);
    ASSERT_EQ(p.monomials[0].a, 0);
    ASSERT_EQ(p.monomials[1].a, 2);
    ASSERT_EQ(p.monomials[2].a, 0);
    ASSERT_EQ(p.monomials[3].a, 0);
    ASSERT_EQ(p.monomials[4].a, 0);
    ASSERT_EQ(p.monomials[5].a, 0);
    ASSERT_DOUBLE_EQ(p(xy).val, 2);

    // p = 0 + 2 + 0 + 0 + 0 + 0
    // dp/dx = 0 + 0 + 0 + 0 + 0 + 0
    err = p.Dxi(0);
    ASSERT_FALSE(err.isError);
    ASSERT_EQ(p.monomials[0].a, 0);
    ASSERT_EQ(p.monomials[1].a, 0);
    ASSERT_EQ(p.monomials[2].a, 0);
    ASSERT_EQ(p.monomials[3].a, 0);
    ASSERT_EQ(p.monomials[4].a, 0);
    ASSERT_EQ(p.monomials[5].a, 0);
    ASSERT_DOUBLE_EQ(p(xy).val, 0.0);
}

TEST_F(PolyTest, DaTest) {
    std::vector<double> X = {1.0, 2.0};
    auto coefs = std::vector<double>(1);
    auto r = n2o2Zeros.Da(X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.Da(X, &coefs);
    ASSERT_TRUE(r.isError);

    double x = 17.0;
    double y = 29.0;
    X = {x, y};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.Da(X, &coefs);
    ASSERT_FALSE(r.isError);
    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0
    ASSERT_DOUBLE_EQ(coefs[0], x * x);
    ASSERT_DOUBLE_EQ(coefs[1], x * y);
    ASSERT_DOUBLE_EQ(coefs[2], x);
    ASSERT_DOUBLE_EQ(coefs[3], y * y);
    ASSERT_DOUBLE_EQ(coefs[4], y);
    ASSERT_DOUBLE_EQ(coefs[5], 1);
}

TEST_F(PolyTest, PolysConstructorTest) {
    Polys p = Polys();
    ASSERT_EQ(p.k.size(), 0);
    ASSERT_EQ(p.polys.size(), 0);

    p = Polys(n2o2);
    ASSERT_EQ(p.k.size(), 1);
    ASSERT_EQ(p.k[0], 1.0);

    ASSERT_EQ(p.polys.size(), 1);
    ASSERT_EQ(p.polys[0].monomials[0].exps, n2o2.monomials[0].exps);
    ASSERT_EQ(p.polys[0].monomials[0].a, n2o2.monomials[0].a);
    ASSERT_EQ(p.polys[0].monomials[1].exps, n2o2.monomials[1].exps);
    ASSERT_EQ(p.polys[0].monomials[1].a, n2o2.monomials[1].a);
    ASSERT_EQ(p.polys[0].monomials[2].exps, n2o2.monomials[2].exps);
    ASSERT_EQ(p.polys[0].monomials[2].a, n2o2.monomials[2].a);
    ASSERT_EQ(p.polys[0].monomials[3].exps, n2o2.monomials[3].exps);
    ASSERT_EQ(p.polys[0].monomials[3].a, n2o2.monomials[3].a);
    ASSERT_EQ(p.polys[0].monomials[4].exps, n2o2.monomials[4].exps);
    ASSERT_EQ(p.polys[0].monomials[4].a, n2o2.monomials[4].a);
    ASSERT_EQ(p.polys[0].monomials[5].exps, n2o2.monomials[5].exps);
    ASSERT_EQ(p.polys[0].monomials[5].a, n2o2.monomials[5].a);
}

TEST_F(PolyTest, MatrixMultiplicationTest) {
    // Tests that the class works correctly with boost matrix multiplication

    Poly p0;
    p0.Build(2, 1, 0);
    p0.monomials[0].a = 1.0;
    p0.monomials[1].a = 2.0;
    p0.monomials[2].a = 3.0;

    Poly p1;
    p1.Build(2, 1, 1);
    p1.monomials[0].a = 4.0;
    p1.monomials[1].a = 5.0;
    p1.monomials[2].a = 6.0;

    using namespace boost::numeric::ublas;

    double k = 5.0;
    matrix<double> K(2, 2);
    K(0, 0) = -k;
    K(0, 1) = k;
    K(1, 0) = k;
    K(1, 1) = -k;

    matrix<Poly> X(2, 1);
    X(0, 0) = p0;
    X(1, 0) = p1;

    // matrix<Poly> mult(2, 1);
    auto mult = prod(K, X);

    // p0 = 1x + 2y + 3
    // p1 = 4x + 5y + 6
    //|-k k|*|p0| = |-k*p0 + k*p1|
    //|k -k| |p1|   |k*p0 - k*p1|
    ASSERT_EQ(mult(0, 0).k.size(), 2);
    int minus_k_p0_index;
    int k_p1_index;
    ASSERT_TRUE(mult(0, 0).k[0] == -k || mult(0, 0).k[0] == k);
    if (mult(0, 0).k[0] == -k) {
        minus_k_p0_index = 0;
        k_p1_index = 1;
    } else {
        minus_k_p0_index = 1;
        k_p1_index = 0;
    }
    ASSERT_EQ(mult(0, 0).k[k_p1_index], k);
    ASSERT_EQ(mult(0, 0).k[minus_k_p0_index], -k);
    ASSERT_EQ(mult(0, 0).polys.size(), 2);

    // -k*p0 = -k*(1x + 2y + 3)
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials.size(), 3);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[0].exps[0], 1);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[0].exps[1], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[0].a, 1);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[1].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[1].exps[1], 1);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[1].a, 2);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[2].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[2].exps[1], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[2].a, 3);
    // k*p1 = k*(4x + 5y + 6)
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials.size(), 3);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[0].exps[0], 1);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[0].exps[1], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[0].a, 4);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[1].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[1].exps[1], 1);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[1].a, 5);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[2].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[2].exps[1], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[2].a, 6);

    double x = 1;
    double y = 2;
    double p0Val = 1 * x + 2 * y + 3;
    double p1Val = 4 * x + 5 * y + 6;
    std::vector<double> xy = {x, y};

    auto eval = mult(0, 0)(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, -k * p0Val + k * p1Val);

    eval = mult(1, 0)(xy);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, k * p0Val - k * p1Val);
}

// TEST_F(PolyTest, GetCoefficientsTest) {
//     auto coefs = std::vector<double>(1);
//     auto r = n2o2Zeros.GetCoefficients(&coefs);
//     ASSERT_TRUE(r.isError);

//     coefs = std::vector<double>(7);
//     r = n2o2Zeros.GetCoefficients(&coefs);
//     ASSERT_TRUE(r.isError);

//     coefs = std::vector<double>(6);

//     r = n2o2Zeros.GetCoefficients(&coefs);
//     ASSERT_FALSE(r.isError);
//     // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0
//     ASSERT_DOUBLE_EQ(coefs[0], 0);
//     ASSERT_DOUBLE_EQ(coefs[1], 0);
//     ASSERT_DOUBLE_EQ(coefs[2], 0);
//     ASSERT_DOUBLE_EQ(coefs[3], 0);
//     ASSERT_DOUBLE_EQ(coefs[4], 0);
//     ASSERT_DOUBLE_EQ(coefs[5], 0);

//     r = n2o2Twos.GetCoefficients(&coefs);
//     ASSERT_FALSE(r.isError);
//     // P = 2*x^2 + 2*xy + 2*x + 2*y^2 + 2*y + 2
//     ASSERT_DOUBLE_EQ(coefs[0], 2);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 2);
//     ASSERT_DOUBLE_EQ(coefs[3], 2);
//     ASSERT_DOUBLE_EQ(coefs[4], 2);
//     ASSERT_DOUBLE_EQ(coefs[5], 2);

//     double x = 3;
//     double y = 5;
//     auto X = std::vector<double>{x, y};
//     ASSERT_DOUBLE_EQ(n2o2Twos(&X).val,
//                      (2 * x * x + 2 * x * y + 2 * x + 2 * y * y + 2 * y +
//                      2));

//     r = n2o2.GetCoefficients(&coefs);
//     ASSERT_FALSE(r.isError);
//     // P = 1*x^2 + 2*xy + 3*x + 4*y^2 + 5*y + 6
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(n2o2(&X).val,
//                      1 * x * x + 2 * x * y + 3 * x + 4 * y * y + 5 * y + 6);
// }

// TEST_F(PolyTest, MultiplicationOperatorTest) {
//     auto coefs = std::vector<double>(6);

//     Poly p = this->n2o2 * 7;
//     auto r = p.GetCoefficients(&coefs);
//     ASSERT_FALSE(r.isError);
//     ASSERT_EQ(p.n, 2);
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(p.monomials[0].k, 7);
//     ASSERT_DOUBLE_EQ(p.monomials[1].k, 7);
//     ASSERT_DOUBLE_EQ(p.monomials[2].k, 7);
//     ASSERT_DOUBLE_EQ(p.monomials[3].k, 7);
//     ASSERT_DOUBLE_EQ(p.monomials[4].k, 7);
//     ASSERT_DOUBLE_EQ(p.monomials[5].k, 7);

//     p = 9 * p;
//     ASSERT_EQ(p.n, 2);
//     ASSERT_DOUBLE_EQ(p.monomials[0].k, 7 * 9);
//     ASSERT_DOUBLE_EQ(p.monomials[1].k, 7 * 9);
//     ASSERT_DOUBLE_EQ(p.monomials[2].k, 7 * 9);
//     ASSERT_DOUBLE_EQ(p.monomials[3].k, 7 * 9);
//     ASSERT_DOUBLE_EQ(p.monomials[4].k, 7 * 9);
//     ASSERT_DOUBLE_EQ(p.monomials[5].k, 7 * 9);

//     // Assert original is not changed
//     r = this->n2o2.GetCoefficients(&coefs);
//     ASSERT_FALSE(r.isError);
//     ASSERT_EQ(this->n2o2.n, 2);
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[0].k, 1);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[1].k, 1);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[2].k, 1);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[3].k, 1);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[4].k, 1);
//     ASSERT_DOUBLE_EQ(n2o2.monomials[5].k, 1);
// }

// TEST_F(PolyTest, PlusOperatorTest) {
//     Poly p0;
//     p0.Build(2, 2, 0.0, 0);
//     auto coefsToSet = std::vector<double>{1, 2, 3, 4, 5, 6};
//     ASSERT_FALSE(p0.SetCoefficients(&coefsToSet).isError);

//     Poly p1;
//     p1.Build(2, 2, 0.0, 6);
//     coefsToSet = std::vector<double>{10, 20, 30, 40, 50, 60};
//     ASSERT_FALSE(p1.SetCoefficients(&coefsToSet).isError);

//     Poly pSum;
//     pSum = p0 + p1;

//     ASSERT_EQ(pSum.n, 2);
//     auto coefs = std::vector<double>(12);
//     ASSERT_FALSE(pSum.GetCoefficients(&coefs).isError);
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(coefs[6], 10);
//     ASSERT_DOUBLE_EQ(coefs[7], 20);
//     ASSERT_DOUBLE_EQ(coefs[8], 30);
//     ASSERT_DOUBLE_EQ(coefs[9], 40);
//     ASSERT_DOUBLE_EQ(coefs[10], 50);
//     ASSERT_DOUBLE_EQ(coefs[11], 60);

//     pSum = pSum + 99.0;
//     coefs = std::vector<double>(12);
//     ASSERT_FALSE(pSum.GetCoefficients(&coefs).isError);
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(coefs[6], 10);
//     ASSERT_DOUBLE_EQ(coefs[7], 20);
//     ASSERT_DOUBLE_EQ(coefs[8], 30);
//     ASSERT_DOUBLE_EQ(coefs[9], 40);
//     ASSERT_DOUBLE_EQ(coefs[10], 50);
//     ASSERT_DOUBLE_EQ(coefs[11], 60 + 99);
// }

// TEST_F(PolyTest, EqualityOperatorTest) {
//     Poly p0;
//     p0.Build(2, 2, 0.0, 0);
//     auto coefsToSet = std::vector<double>{1, 2, 3, 4, 5, 6};
//     ASSERT_FALSE(p0.SetCoefficients(&coefsToSet).isError);

//     Poly p1;
//     p1.Build(2, 2, 0.0, 0);
//     coefsToSet = std::vector<double>{10, -20, 30, 40, -50, 60};
//     ASSERT_FALSE(p1.SetCoefficients(&coefsToSet).isError);

//     Poly p2;
//     p2.Build(2, 2, 0.0, 0);
//     coefsToSet = std::vector<double>{1, 2, 3, 4, 5, 6};
//     ASSERT_FALSE(p2.SetCoefficients(&coefsToSet).isError);

//     ASSERT_FALSE(p0 == p1);
//     ASSERT_TRUE(p0 != p1);
//     ASSERT_FALSE(p1 == p2);
//     ASSERT_TRUE(p1 != p2);
//     ASSERT_TRUE(p0 == p2);
//     ASSERT_FALSE(p0 != p2);

//     Poly p3;
//     Poly p4;
//     ASSERT_TRUE(p3 == p4);
//     ASSERT_FALSE(p3 != p4);

//     p3.Build(2, 3, 0.0, 0);
//     p4.Build(2, 2, 0.0, 0);
//     ASSERT_FALSE(p3 == p4);
//     ASSERT_TRUE(p3 != p4);
// }

// TEST_F(PolyTest, SetCoefficientsTest) {
//     auto coefsToSet = std::vector<double>(1);
//     auto r = n2o2Zeros.GetCoefficients(&coefsToSet);
//     ASSERT_TRUE(r.isError);

//     coefsToSet = std::vector<double>(7);
//     r = n2o2Zeros.GetCoefficients(&coefsToSet);
//     ASSERT_TRUE(r.isError);

//     coefsToSet = std::vector<double>{9, 10, 11, 12, 13, 14};
//     r = n2o2Zeros.SetCoefficients(&coefsToSet);
//     ASSERT_FALSE(r.isError);

//     auto coefsSet = std::vector<double>(6);
//     r = n2o2Zeros.GetCoefficients(&coefsSet);
//     ASSERT_FALSE(r.isError);
//     ASSERT_DOUBLE_EQ(coefsSet[0], 9);
//     ASSERT_DOUBLE_EQ(coefsSet[1], 10);
//     ASSERT_DOUBLE_EQ(coefsSet[2], 11);
//     ASSERT_DOUBLE_EQ(coefsSet[3], 12);
//     ASSERT_DOUBLE_EQ(coefsSet[4], 13);
//     ASSERT_DOUBLE_EQ(coefsSet[5], 14);
// }

// TEST_F(PolyTest, PlusAndMultiplicationTest) {
//     auto coefs = std::vector<double>(4);

//     Poly p0;
//     p0.Build(1, 3, 0.0, 0);
//     coefs = {1, 2, 3, 4};
//     ASSERT_FALSE(p0.SetCoefficients(&coefs).isError);
//     // p0 = x + 2y + 3x + 4

//     Poly p1;
//     p1.Build(1, 3, 0.0, 4);
//     coefs = {5, 6, 7, 8};
//     ASSERT_FALSE(p1.SetCoefficients(&coefs).isError);
//     // p0 = 5x + 6y + 7z + 8

//     Poly p2;
//     p2 = 2 * p0 + (-3) * p1;
//     ASSERT_EQ(p2.n, 1);
//     // p2 = 2*(x + 2y + 3x + 4) -3*(5x + 6y + 7z + 8)
//     coefs = std::vector<double>(8);
//     ASSERT_FALSE(p2.GetCoefficients(&coefs).isError);
//     ASSERT_DOUBLE_EQ(coefs[0], 1);
//     ASSERT_DOUBLE_EQ(coefs[1], 2);
//     ASSERT_DOUBLE_EQ(coefs[2], 3);
//     ASSERT_DOUBLE_EQ(coefs[3], 4);
//     ASSERT_DOUBLE_EQ(coefs[4], 5);
//     ASSERT_DOUBLE_EQ(coefs[5], 6);
//     ASSERT_DOUBLE_EQ(coefs[6], 7);
//     ASSERT_DOUBLE_EQ(coefs[7], 8);
//     ASSERT_DOUBLE_EQ(p2.monomials[0].k, 2);
//     ASSERT_DOUBLE_EQ(p2.monomials[1].k, 2);
//     ASSERT_DOUBLE_EQ(p2.monomials[2].k, 2);
//     ASSERT_DOUBLE_EQ(p2.monomials[3].k, 2);
//     ASSERT_DOUBLE_EQ(p2.monomials[4].k, -3);
//     ASSERT_DOUBLE_EQ(p2.monomials[5].k, -3);
//     ASSERT_DOUBLE_EQ(p2.monomials[6].k, -3);
//     ASSERT_DOUBLE_EQ(p2.monomials[7].k, -3);
// }
