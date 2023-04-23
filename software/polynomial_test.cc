#include "polynomial.h"

#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operations.hpp>

#include "maybe.h"
#include "utils.h"

TEST(MonomialTest, constructorTest) {
    // Empty monomial
    Monomial m = Monomial(std::vector<int>{});

    EXPECT_EQ(m.exps.size(), 0);

    // 0*x0^1*x1^0*x2^3
    m = Monomial(std::vector<int>{1, 0, 3});

    ASSERT_EQ(m.exps.size(), 3);
    EXPECT_EQ(m.exps[0], 1);
    EXPECT_EQ(m.exps[1], 0);
    EXPECT_EQ(m.exps[2], 3);
}

TEST(MonomialTest, evalTest) {
    // 1*1*x0^1*x1^0*x2^3
    Monomial m = Monomial(std::vector<int>{1, 0, 3});
    double a = 932847.1;

    double x0 = 1234;
    double x1 = 56234;
    double x2 = 323.2;
    auto X = std::vector<double>{x0, x1, x2};
    Maybe<double> eval = m(a, X);
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
    ASSERT_EQ(p.monomials[0].k, 1);
    ASSERT_EQ(p.monomials[0].exps.size(), 0);

    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0*1
    r = p.Build(2, 2);
    ASSERT_FALSE(r.isError);
    ASSERT_EQ(p.monomials.size(), 6);
    // 0*x^2
    ASSERT_EQ(p.monomials[0].k, 1);
    ASSERT_EQ(p.monomials[0].exps.size(), 2);
    ASSERT_EQ(p.monomials[0].exps[0], 2);
    ASSERT_EQ(p.monomials[0].exps[1], 0);
    // 0*xy;
    ASSERT_EQ(p.monomials[1].k, 1);
    ASSERT_EQ(p.monomials[1].exps.size(), 2);
    ASSERT_EQ(p.monomials[1].exps[0], 1);
    ASSERT_EQ(p.monomials[1].exps[1], 1);
    // 0*x
    ASSERT_EQ(p.monomials[2].k, 1);
    ASSERT_EQ(p.monomials[2].exps.size(), 2);
    ASSERT_EQ(p.monomials[2].exps[0], 1);
    ASSERT_EQ(p.monomials[2].exps[1], 0);
    // 0*y^2
    ASSERT_EQ(p.monomials[3].k, 1);
    ASSERT_EQ(p.monomials[3].exps.size(), 2);
    ASSERT_EQ(p.monomials[3].exps[0], 0);
    ASSERT_EQ(p.monomials[3].exps[1], 2);
    // 0*y
    ASSERT_EQ(p.monomials[4].k, 1);
    ASSERT_EQ(p.monomials[4].exps.size(), 2);
    ASSERT_EQ(p.monomials[4].exps[0], 0);
    ASSERT_EQ(p.monomials[4].exps[1], 1);
    // 0*1
    ASSERT_EQ(p.monomials[5].k, 1);
    ASSERT_EQ(p.monomials[5].exps.size(), 2);
    ASSERT_EQ(p.monomials[5].exps[0], 0);
    ASSERT_EQ(p.monomials[5].exps[1], 0);

    r = p.Build(3, 3, false, 0);
    // t^3
    // t^2x + t^2y + t^2
    // tx^2 + txy + tx + ty^2 + ty + t
    // x^3
    // x^2y + x^2
    // xy^2 + x*y + x
    // y^3 + y^2 + y + 1
    ASSERT_EQ(p.monomials.size(), 20);
    std::vector<std::vector<int>> exps = std::vector<std::vector<int>>{
        std::vector<int>{3, 0, 0},

        std::vector<int>{2, 1, 0}, std::vector<int>{2, 0, 1},
        std::vector<int>{2, 0, 0},

        std::vector<int>{1, 2, 0}, std::vector<int>{1, 1, 1},
        std::vector<int>{1, 1, 0}, std::vector<int>{1, 0, 2},
        std::vector<int>{1, 0, 1}, std::vector<int>{1, 0, 0},

        std::vector<int>{0, 3, 0}, std::vector<int>{0, 2, 1},
        std::vector<int>{0, 2, 0},

        std::vector<int>{0, 1, 2}, std::vector<int>{0, 1, 1},
        std::vector<int>{0, 1, 0},

        std::vector<int>{0, 0, 3}, std::vector<int>{0, 0, 2},
        std::vector<int>{0, 0, 1}, std::vector<int>{0, 0, 0},
    };
    for (int i = 0; i < 20; i++) {
        ASSERT_EQ(p.monomials[i].exps, exps[i]);
    }

    r = p.Build(3, 3, true, 0);
    // t^3
    // t^2x + t^2y + t^2
    // tx + ty + t
    // 1
    ASSERT_EQ(p.monomials.size(), 8);
    exps = std::vector<std::vector<int>>{
        std::vector<int>{3, 0, 0},

        std::vector<int>{2, 1, 0}, std::vector<int>{2, 0, 1},
        std::vector<int>{2, 0, 0},

        std::vector<int>{1, 1, 0}, std::vector<int>{1, 0, 1},
        std::vector<int>{1, 0, 0},

        std::vector<int>{0, 0, 0},
    };
    for (int i = 0; i < 8; i++) {
        ASSERT_EQ(p.monomials[i].exps, exps[i]);
    }
}

class PolyTest : public testing::Test {
   public:
    // Polynomial of two vars or order two. Monomials are:
    // x^2, xy, x, y^2, y, 1
    Poly n2o2;

    // Other instances of Poly of two vars or order two, but with ids 1 and 2.
    Poly n2o2_1;
    Poly n2o2_2;

    // Coefficients that can be used for n2o2
    std::vector<double> zeros;
    std::vector<double> twos;
    std::vector<double> all;

    void SetUp() {
        // Called before every TEST_F
        Maybe<Void> r = n2o2.Build(2, 2, false, 0);
        ASSERT_FALSE(r.isError);
        r = n2o2_1.Build(2, 2, false, 1);
        ASSERT_FALSE(r.isError);
        r = n2o2_2.Build(2, 2, false, 2);
        ASSERT_FALSE(r.isError);

        zeros = std::vector<double>{0, 0, 0, 0, 0, 0};
        twos = std::vector<double>{2, 2, 2, 2, 2, 2};
        all = std::vector<double>{1, 2, 3, 4, 5, 6};
    }
};

TEST_F(PolyTest, GetSetXTest) {
    auto xGet = n2o2.GetX();
    ASSERT_EQ(xGet.size(), 2);
    ASSERT_DOUBLE_EQ(xGet[0], 0);
    ASSERT_DOUBLE_EQ(xGet[1], 0);

    std::vector<double> X;
    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};
    n2o2.SetX(X);
    ASSERT_EQ(n2o2.X.size(), 2);
    ASSERT_DOUBLE_EQ(n2o2.X[0], x1);
    ASSERT_DOUBLE_EQ(n2o2.X[1], x2);

    xGet = n2o2.GetX();
    ASSERT_EQ(xGet.size(), 2);
    ASSERT_DOUBLE_EQ(xGet[0], x1);
    ASSERT_DOUBLE_EQ(xGet[1], x2);

    X = {};  // too short
    ASSERT_DEATH({ n2o2.SetX(X); }, "");

    X = {1.0, 1.0, 2.0};  // too long
    ASSERT_DEATH({ n2o2.SetX(X); }, "");
}

TEST_F(PolyTest, scalarMultiplicationTest) {
    for (int i = 0; i < n2o2.monomials.size(); i++) {
        n2o2.monomials[i].k = i;
    }

    n2o2 *= 10;
    for (int i = 0; i < n2o2.monomials.size(); i++) {
        ASSERT_DOUBLE_EQ(n2o2.monomials[i].k, i * 10);
    }
    n2o2 *= 0;
    for (int i = 0; i < n2o2.monomials.size(); i++) {
        ASSERT_DOUBLE_EQ(n2o2.monomials[i].k, i * 0);
    }
}

TEST_F(PolyTest, operatorTest) {
    std::vector<double> X;
    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    n2o2.SetX(X);
    auto eval = n2o2(all);
    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, expected);
}

TEST_F(PolyTest, copyConstructorAndAssignment) {
    Poly copyConstructed = n2o2;
    n2o2.monomials[0].k = 99;
    ASSERT_EQ(copyConstructed.monomials[0].k, 1.0);
    ASSERT_EQ(copyConstructed.monomials[1].k, 1.0);
    ASSERT_EQ(copyConstructed.monomials[2].k, 1.0);
    ASSERT_EQ(copyConstructed.monomials[3].k, 1.0);
    ASSERT_EQ(copyConstructed.monomials[4].k, 1.0);
    ASSERT_EQ(copyConstructed.monomials[5].k, 1.0);

    Poly assigned;
    assigned.Build(2, 3);  // Shouldn't change anything
    assigned = n2o2_1;
    n2o2_1.monomials[0].k = 99;
    ASSERT_EQ(assigned.monomials[0].k, 1);
    ASSERT_EQ(assigned.monomials[1].k, 1);
    ASSERT_EQ(assigned.monomials[2].k, 1);
    ASSERT_EQ(assigned.monomials[3].k, 1);
    ASSERT_EQ(assigned.monomials[4].k, 1);
    ASSERT_EQ(assigned.monomials[5].k, 1);
}

TEST_F(PolyTest, DxiTest) {
    double x = 987;
    double y = 91;
    std::vector<double> xy = {x, y};
    Maybe<Void> err;
    Maybe<double> eval;
    Poly p;

    p = n2o2;
    p.SetX(xy);
    ASSERT_FALSE(p.Dxi(0).isError);
    eval = p(zeros);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    p = n2o2;
    p.SetX(xy);
    ASSERT_FALSE(p.Dxi(1).isError);
    eval = p(zeros);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    p = n2o2;
    p.SetX(xy);
    ASSERT_TRUE(p.Dxi(2).isError);
    ASSERT_TRUE(p.Dxi(-1).isError);

    p = n2o2;
    p.SetX(xy);
    // P = a0x^2 + a1*xy + a2*x + a3*y^2 + a4*y + a5*1
    ASSERT_FALSE(p.Dxi(0).isError);
    // P = 2*a0x + a1*y + a2 + 0*a3*y^2 + 0*a4*y + 0*a5*1
    double a0 = 1;
    double a1 = 2;
    double a2 = 3;
    eval = p(all);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 2 * a0 * x + 1 * a1 * y + a2);

    p = n2o2;
    p.SetX(xy);
    // P = a0x^2 + a1*xy + a2*x + a3*y^2 + a4*y + a5*1
    ASSERT_FALSE(p.Dxi(1).isError);
    // P = a1*x + 2*a3*y + a4
    eval = p(all);
    double a3 = 4;
    double a4 = 5;
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, a1 * x + 2 * a3 * y + a4);
}

TEST_F(PolyTest, DxiTest2) {
    Maybe<Void> err;
    double x = 245;
    double y = 5763;
    auto xy = std::vector<double>{x, y};

    Poly p;
    p = this->n2o2;
    p.SetX(xy);

    ASSERT_TRUE(p.Dxi(-1).isError);
    ASSERT_TRUE(p.Dxi(2).isError);

    // a1 = 1, a2 = 2 ... a6 = 6
    // n2o2 = a1x^2 + a2xy + a3x + a4y^2 + a5y + a6
    // dn2o2/dx = 2*a1*x + 1*a2*y + 1*a3*1 + 0 + 0 + 0
    err = p.Dxi(0);
    ASSERT_FALSE(err.isError);
    ASSERT_EQ(p.monomials[0].k, 2);
    ASSERT_EQ(p.monomials[1].k, 1);
    ASSERT_EQ(p.monomials[2].k, 1);
    ASSERT_EQ(p.monomials[3].k, 0);
    ASSERT_EQ(p.monomials[4].k, 0);
    ASSERT_EQ(p.monomials[5].k, 0);
    ASSERT_DOUBLE_EQ(p(all).val, 2 * x + 2 * y + 3);

    // p = 2*a1*x + 1*a2*y + 1*a3*1 + 0 + 0 + 0
    // dp/dy = 0 + 1*a2*1 + 0 + 0 + 0 + 0
    err = p.Dxi(1);
    ASSERT_FALSE(err.isError);
    ASSERT_EQ(p.monomials[0].k, 0);
    ASSERT_EQ(p.monomials[1].k, 1);
    ASSERT_EQ(p.monomials[2].k, 0);
    ASSERT_EQ(p.monomials[3].k, 0);
    ASSERT_EQ(p.monomials[4].k, 0);
    ASSERT_EQ(p.monomials[5].k, 0);
    ASSERT_DOUBLE_EQ(p(all).val, 2);
}

TEST_F(PolyTest, DaTest) {
    auto r = n2o2.Dai(-1);
    ASSERT_TRUE(r.isError);

    std::vector<double> X;
    double x = 17.0;
    double y = 29.0;
    X = {x, y};
    n2o2.SetX(X);

    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0
    r = n2o2.Dai(0);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, x * x);

    r = n2o2.Dai(1);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, x * y);

    r = n2o2.Dai(2);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, x);

    r = n2o2.Dai(3);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, y * y);

    r = n2o2.Dai(4);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, y);

    r = n2o2.Dai(5);
    ASSERT_FALSE(r.isError);
    ASSERT_DOUBLE_EQ(r.val, 1);
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
    ASSERT_EQ(p.polys[0].monomials[1].exps, n2o2.monomials[1].exps);
    ASSERT_EQ(p.polys[0].monomials[2].exps, n2o2.monomials[2].exps);
    ASSERT_EQ(p.polys[0].monomials[3].exps, n2o2.monomials[3].exps);
    ASSERT_EQ(p.polys[0].monomials[4].exps, n2o2.monomials[4].exps);
    ASSERT_EQ(p.polys[0].monomials[5].exps, n2o2.monomials[5].exps);
}

TEST_F(PolyTest, MatrixMultiplicationTest) {
    // Tests that the class works correctly with boost matrix multiplication
    double x = 1;
    double y = 2;
    std::vector<double> xy = {x, y};

    Poly p0;
    p0.Build(2, 1, false, 0);
    p0.SetX(xy);

    Poly p1;
    p1.Build(2, 1, false, 1);
    p1.SetX(xy);

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
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[1].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[1].exps[1], 1);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[2].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[minus_k_p0_index].monomials[2].exps[1], 0);
    // k*p1 = k*(4x + 5y + 6)
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials.size(), 3);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[0].exps[0], 1);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[0].exps[1], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[1].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[1].exps[1], 1);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[2].exps[0], 0);
    ASSERT_EQ(mult(0, 0).polys[k_p1_index].monomials[2].exps[1], 0);

    std::vector<double> a0 = std::vector<double>{1, 2, 3};
    std::vector<double> a1 = std::vector<double>{4, 5, 6};
    std::vector<std::vector<double>> a =
        std::vector<std::vector<double>>{a0, a1};
    double p0Val = 1 * x + 2 * y + 3;
    double p1Val = 4 * x + 5 * y + 6;

    auto eval = mult(0, 0)(a);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, -k * p0Val + k * p1Val);

    eval = mult(1, 0)(a);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, k * p0Val - k * p1Val);

    // Test error: a is missing some coefficients
    a = std::vector<std::vector<double>>{a1};
    eval = mult(1, 0)(a);
    ASSERT_TRUE(eval.isError);
    a = std::vector<std::vector<double>>{a0};
    eval = mult(1, 0)(a);
    ASSERT_TRUE(eval.isError);
}

TEST_F(PolyTest, DxiPolysTest) {
    double x = 987;
    double y = 91;
    std::vector<double> xy = {x, y};
    n2o2.SetX(xy);

    std::vector<std::vector<double>> coefs;
    Maybe<Void> err;
    Maybe<double> eval;
    Polys ps;

    ps = Polys(n2o2);
    ASSERT_FALSE(ps.Dxi(0).isError);
    coefs = {zeros};
    eval = ps(coefs);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    ps = Polys(n2o2);
    ASSERT_FALSE(ps.Dxi(1).isError);
    coefs = {zeros};
    eval = ps(coefs);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 0);

    ps = Polys(n2o2);
    ASSERT_TRUE(ps.Dxi(2).isError);
    ASSERT_TRUE(ps.Dxi(-1).isError);

    ps = Polys(n2o2);
    // P = 1x^2 + 2xy + 3x + 4y^2 + 5y + 6
    ASSERT_FALSE(ps.Dxi(0).isError);
    coefs = {all};
    eval = ps(coefs);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 1.0 * 2 * x + 2.0 * y + 3.0);

    ps = Polys(n2o2);
    ASSERT_FALSE(ps.Dxi(1).isError);
    coefs = {all};
    eval = ps(coefs);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 2.0 * x + 4.0 * 2 * y + 5.0);

    ps = Polys(n2o2);
    ps *= 7;
    ps += n2o2;
    // ps = 8*(1x^2 + 2xy + 3x + 4y^2 + 5y + 6)
    ASSERT_FALSE(ps.Dxi(0).isError);
    coefs = {all};
    eval = ps(coefs);
    ASSERT_FALSE(eval.isError);
    ASSERT_EQ(eval.val, 8 * (1.0 * 2 * x + 2.0 * y + 3.0));
}

TEST_F(PolyTest, DaPolysTest) {
    double x = 7;
    double y = 8;
    std::vector<double> xy = {x, y};
    std::map<std::tuple<int, int>, double> map;
    Polys ps;

    n2o2.SetX(xy);
    n2o2_1.SetX(xy);
    n2o2_2.SetX(xy);

    ps = n2o2 + n2o2_1 + 7.0 * n2o2_2;
    // ps = p0 + p1 + 7*p2
    // ps = x^2 + xy + x + y^2 + y + 1
    //      + x^2 + xy + x + y^2 + y + 1
    //      + 7*(x^2 + xy + x + y^2 + y + 1)

    map = ps.Da();
    ASSERT_EQ(map.size(), 6 * 3);

    // p0
    auto val = map[std::tuple<int, int>{0, 0}];
    ASSERT_DOUBLE_EQ(val, x * x);
    val = map[std::tuple<int, int>{0, 1}];
    ASSERT_DOUBLE_EQ(val, x * y);
    val = map[std::tuple<int, int>{0, 2}];
    ASSERT_DOUBLE_EQ(val, x);
    val = map[std::tuple<int, int>{0, 3}];
    ASSERT_DOUBLE_EQ(val, y * y);
    val = map[std::tuple<int, int>{0, 4}];
    ASSERT_DOUBLE_EQ(val, y);
    val = map[std::tuple<int, int>{0, 5}];
    ASSERT_DOUBLE_EQ(val, 1);

    // p1
    val = map[std::tuple<int, int>{1, 0}];
    ASSERT_DOUBLE_EQ(val, x * x);
    val = map[std::tuple<int, int>{1, 1}];
    ASSERT_DOUBLE_EQ(val, x * y);
    val = map[std::tuple<int, int>{1, 2}];
    ASSERT_DOUBLE_EQ(val, x);
    val = map[std::tuple<int, int>{1, 3}];
    ASSERT_DOUBLE_EQ(val, y * y);
    val = map[std::tuple<int, int>{1, 4}];
    ASSERT_DOUBLE_EQ(val, y);
    val = map[std::tuple<int, int>{1, 5}];
    ASSERT_DOUBLE_EQ(val, 1);

    // p2
    val = map[std::tuple<int, int>{2, 0}];
    ASSERT_DOUBLE_EQ(val, 7 * x * x);
    val = map[std::tuple<int, int>{2, 1}];
    ASSERT_DOUBLE_EQ(val, 7 * x * y);
    val = map[std::tuple<int, int>{2, 2}];
    ASSERT_DOUBLE_EQ(val, 7 * x);
    val = map[std::tuple<int, int>{2, 3}];
    ASSERT_DOUBLE_EQ(val, 7 * y * y);
    val = map[std::tuple<int, int>{2, 4}];
    ASSERT_DOUBLE_EQ(val, 7 * y);
    val = map[std::tuple<int, int>{2, 5}];
    ASSERT_DOUBLE_EQ(val, 7 * 1);

    ps = n2o2;
    // a = 2
    // ps = x^2 + xy + x + y^2 + y + 1
    ps.Dxi(0);
    // ps = 2*x + 1*y + 1 + 0*y^2 + 0*y + 0*1
    ps += 5 * n2o2;
    // ps = (2*x + 1*y + 1 + 0*y^2 + 0*y + 0*1)
    //     +5*(x^2 + xy + x + y^2 + y + 1)
    ASSERT_EQ(ps.polys.size(), 2);

    // dn2o2/dx
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials.size(), 6);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[0].k, 2);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[1].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[2].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[3].k, 0);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[4].k, 0);
    ASSERT_DOUBLE_EQ(ps.polys[0].monomials[5].k, 0);
    ASSERT_EQ(ps.polys[0].monomials[0].exps[0], 1);
    ASSERT_EQ(ps.polys[0].monomials[0].exps[1], 0);
    ASSERT_EQ(ps.polys[0].monomials[1].exps[0], 0);
    ASSERT_EQ(ps.polys[0].monomials[1].exps[1], 1);
    ASSERT_EQ(ps.polys[0].monomials[2].exps[0], 0);
    ASSERT_EQ(ps.polys[0].monomials[2].exps[1], 0);
    ASSERT_EQ(ps.polys[0].monomials[3].exps[0], 0);
    ASSERT_EQ(ps.polys[0].monomials[3].exps[1], 2);
    ASSERT_EQ(ps.polys[0].monomials[4].exps[0], 0);
    ASSERT_EQ(ps.polys[0].monomials[4].exps[1], 1);
    ASSERT_EQ(ps.polys[0].monomials[5].exps[0], 0);
    ASSERT_EQ(ps.polys[0].monomials[5].exps[1], 0);

    // 5*n2o2
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials.size(), 6);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[0].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[1].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[2].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[3].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[4].k, 1);
    ASSERT_DOUBLE_EQ(ps.polys[1].monomials[5].k, 1);
    ASSERT_EQ(ps.polys[1].monomials[0].exps[0], 2);
    ASSERT_EQ(ps.polys[1].monomials[0].exps[1], 0);
    ASSERT_EQ(ps.polys[1].monomials[1].exps[0], 1);
    ASSERT_EQ(ps.polys[1].monomials[1].exps[1], 1);
    ASSERT_EQ(ps.polys[1].monomials[2].exps[0], 1);
    ASSERT_EQ(ps.polys[1].monomials[2].exps[1], 0);
    ASSERT_EQ(ps.polys[1].monomials[3].exps[0], 0);
    ASSERT_EQ(ps.polys[1].monomials[3].exps[1], 2);
    ASSERT_EQ(ps.polys[1].monomials[4].exps[0], 0);
    ASSERT_EQ(ps.polys[1].monomials[4].exps[1], 1);
    ASSERT_EQ(ps.polys[1].monomials[5].exps[0], 0);
    ASSERT_EQ(ps.polys[1].monomials[5].exps[1], 0);

    map = ps.Da();
    // dps/da0 = 2*x + 5*x^2
    val = map[std::tuple<int, int>{n2o2.id, 0}];
    ASSERT_DOUBLE_EQ(val, 2 * x + 5 * x * x);
}

TEST_F(PolyTest, EqualityOperatorTest) {
    ASSERT_EQ(Polys(n2o2), Polys(n2o2));
    ASSERT_EQ(Polys(n2o2) + Polys(n2o2_1), Polys(n2o2) + Polys(n2o2_1));

    ASSERT_NE(Polys(n2o2_1), Polys(n2o2));
    ASSERT_NE(2 * Polys(n2o2), Polys(n2o2));
    ASSERT_NE(Polys(n2o2) + Polys(n2o2_1), Polys(n2o2_1) + Polys(n2o2));
}

TEST_F(PolyTest, PolysNMonomialsTest) {
    ASSERT_EQ(Polys(n2o2).nMonomials(), 6);
    ASSERT_EQ((Polys(n2o2) + Polys(n2o2)).nMonomials(), 12);
}