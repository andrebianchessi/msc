#include "polynomial.h"

#include <gtest/gtest.h>

#include "maybe.h"
#include "utils.h"

TEST(NodeTest, AddChildTest) {
    // P = x^2 + xy + x + y^2 + y + 1
    auto root = Node(3);
    ASSERT_EQ(root.children.size(), 3);

    root.AddChild(2, 0, 1);
    root.AddChild(1, 0, 2);
    root.AddChild(0, 0, 3);

    // x^2*y^0
    root.children[0]->AddChild(0, 2, 0);

    // xy
    root.children[1]->AddChild(1, 1, 0);

    // xy^0
    root.children[1]->AddChild(0, 1, 0);

    // x^0y^2
    root.children[2]->AddChild(2, 0, 0);

    // x^0y^1
    root.children[2]->AddChild(1, 0, 0);

    // x^0y^1
    root.children[2]->AddChild(0, 0, 0);

    ASSERT_EQ(root.children.size(), 3);
    ASSERT_EQ(root.children[0]->exp, 2);
    ASSERT_EQ(root.children[1]->exp, 1);
    ASSERT_EQ(root.children[2]->exp, 0);

    ASSERT_EQ(root.children[0]->children.size(), 1);
    ASSERT_EQ(root.children[0]->children[0]->exp, 0);

    ASSERT_EQ(root.children[1]->children.size(), 2);
    ASSERT_EQ(root.children[1]->children[0]->exp, 1);
    ASSERT_EQ(root.children[1]->children[1]->exp, 0);

    ASSERT_EQ(root.children[2]->children.size(), 3);
    ASSERT_EQ(root.children[2]->children[0]->exp, 2);
    ASSERT_EQ(root.children[2]->children[1]->exp, 1);
    ASSERT_EQ(root.children[2]->children[2]->exp, 0);
}

TEST(PolyConstructorTest, constructorTest) {
    // Invalid values
    Maybe<Poly> r = Poly::NewPoly(-1, 0);
    ASSERT_TRUE(r.isError);

    r = Poly::NewPoly(0, -1);
    ASSERT_TRUE(r.isError);

    r = Poly::NewPoly(-10, -10);
    ASSERT_TRUE(r.isError);

    // y() = a
    r = Poly::NewPoly(0, 0);
    ASSERT_FALSE(r.isError);
    Poly p = r.val;
    ASSERT_EQ(p.coefficients->children.size(), 1);

    r = Poly::NewPoly(2, 2, 3.0);
    ASSERT_FALSE(r.isError);
    p = r.val;
    // Same tree as in NodeTest
    // P = x^2 + xy + x + y^2 + y + 1
    ASSERT_EQ(p.coefficients->children.size(), 3);
    ASSERT_EQ(p.coefficients->children[0]->exp, 2);
    ASSERT_EQ(p.coefficients->children[1]->exp, 1);
    ASSERT_EQ(p.coefficients->children[2]->exp, 0);
    ASSERT_EQ(p.coefficients->children[0]->children.size(), 1);
    ASSERT_EQ(p.coefficients->children[0]->children[0]->exp, 0);
    ASSERT_EQ(p.coefficients->children[1]->children.size(), 2);
    ASSERT_EQ(p.coefficients->children[1]->children[0]->exp, 1);
    ASSERT_EQ(p.coefficients->children[1]->children[1]->exp, 0);
    ASSERT_EQ(p.coefficients->children[2]->children.size(), 3);
    ASSERT_EQ(p.coefficients->children[2]->children[0]->exp, 2);
    ASSERT_EQ(p.coefficients->children[2]->children[1]->exp, 1);
    ASSERT_EQ(p.coefficients->children[2]->children[2]->exp, 0);

    ASSERT_DOUBLE_EQ(p.coefficients->children[0]->children[0]->a, 3.0);
    ASSERT_DOUBLE_EQ(p.coefficients->children[1]->children[0]->a, 3.0);
    ASSERT_DOUBLE_EQ(p.coefficients->children[1]->children[1]->a, 3.0);
    ASSERT_DOUBLE_EQ(p.coefficients->children[2]->children[0]->a, 3.0);
    ASSERT_DOUBLE_EQ(p.coefficients->children[2]->children[1]->a, 3.0);
    ASSERT_DOUBLE_EQ(p.coefficients->children[2]->children[2]->a, 3.0);

    ASSERT_EQ(p.coefficients->children[0]->children[0]->children.size(), 0);
    ASSERT_EQ(p.coefficients->children[1]->children[0]->children.size(), 0);
    ASSERT_EQ(p.coefficients->children[1]->children[1]->children.size(), 0);
    ASSERT_EQ(p.coefficients->children[2]->children[0]->children.size(), 0);
    ASSERT_EQ(p.coefficients->children[2]->children[1]->children.size(), 0);
    ASSERT_EQ(p.coefficients->children[2]->children[2]->children.size(), 0);
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
        Maybe<Poly> r = Poly::NewPoly(2, 2);
        ASSERT_FALSE(r.isError);
        Poly p = r.val;
        this->n2o2Zeros = p;

        r = Poly::NewPoly(2, 2, 2.0);
        ASSERT_FALSE(r.isError);
        p = r.val;
        this->n2o2Twos = p;

        r = Poly::NewPoly(2, 2);
        ASSERT_FALSE(r.isError);
        p = r.val;
        p.coefficients->children[0]->children[0]->a = 1.0;
        p.coefficients->children[1]->children[0]->a = 2.0;
        p.coefficients->children[1]->children[1]->a = 3.0;
        p.coefficients->children[2]->children[0]->a = 4.0;
        p.coefficients->children[2]->children[1]->a = 5.0;
        p.coefficients->children[2]->children[2]->a = 6.0;
        this->n2o2 = p;
    }
};

TEST_F(PolyTest, dfsTest) {
    double a = 2.0;  // coefficients
    double x = 2.0;
    double y = 3.0;
    std::vector<double> X = {x, y};

    double a_y_0 =
        dfs(this->n2o2Twos.coefficients->children[0]->children[0], 2, &X);
    ASSERT_DOUBLE_EQ(a_y_0, a * 1);
    double x_2_a_y_0 = dfs(this->n2o2Twos.coefficients->children[0], 1, &X);
    ASSERT_DOUBLE_EQ(x_2_a_y_0, x * x * a * 1);
}

TEST_F(PolyTest, operatorTest) {
    std::vector<double> X;
    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    auto eval = n2o2(&X);
    ASSERT_FALSE(eval.isError);

    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_DOUBLE_EQ(eval.val, expected);

    n2o2.Multiply(5.0);
    eval = n2o2(&X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, expected * 5);
}

TEST_F(PolyTest, DxiTest) {
    std::vector<double> X = {1.0, 2.0};
    auto eval = n2o2Zeros.Dxi(0, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, 0.0);

    eval = n2o2Zeros.Dxi(1, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, 0.0);

    ASSERT_TRUE(n2o2Zeros.Dxi(2, &X).isError);
    ASSERT_TRUE(n2o2Zeros.Dxi(-1, &X).isError);

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    eval = n2o2.Dxi(0, &X);
    ASSERT_FALSE(eval.isError);
    double d1Expected = 1.0 * 2 * x1 + 2.0 * x2 + 3.0;
    ASSERT_DOUBLE_EQ(eval.val, d1Expected);

    eval = n2o2.Dxi(1, &X);
    ASSERT_FALSE(eval.isError);
    double d2Expected = 2.0 * x1 + 4.0 * 2 * x2 + 5.0;
    ASSERT_DOUBLE_EQ(eval.val, d2Expected);

    n2o2.Multiply(5.0);
    eval = n2o2.Dxi(0, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, d1Expected * 5);
    eval = n2o2.Dxi(1, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, d2Expected * 5);
}

TEST_F(PolyTest, D2xiTest) {
    std::vector<double> X = {1.0, 2.0};

    auto eval = n2o2Zeros.D2xi(0, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, 0.0);
    eval = n2o2Zeros.D2xi(1, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, 0.0);

    ASSERT_TRUE(n2o2Zeros.D2xi(2, &X).isError);
    ASSERT_TRUE(n2o2Zeros.D2xi(-1, &X).isError);

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    eval = n2o2.D2xi(0, &X);
    ASSERT_FALSE(eval.isError);
    double dd1Expected = 1.0 * 2;
    ASSERT_DOUBLE_EQ(eval.val, dd1Expected);

    eval = n2o2.D2xi(1, &X);
    ASSERT_FALSE(eval.isError);
    double dd2Expected = 4.0 * 2;
    ASSERT_DOUBLE_EQ(eval.val, dd2Expected);

    n2o2.Multiply(5.0);
    eval = n2o2.D2xi(0, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, dd1Expected * 5);
    eval = n2o2.D2xi(1, &X);
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(eval.val, dd2Expected * 5);
}

TEST_F(PolyTest, DaTest) {
    std::vector<double> X = {1.0, 2.0};
    auto coefs = std::vector<double>(1);
    auto r = n2o2Zeros.Da(&X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.Da(&X, &coefs);
    ASSERT_TRUE(r.isError);

    double x = 17.0;
    double y = 29.0;
    X = {x, y};
    coefs = std::vector<double>(6);
    n2o2Zeros.Multiply(10);  // Should not change anything
    r = n2o2Zeros.Da(&X, &coefs);
    ASSERT_FALSE(r.isError);
    // P = 10*(0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0)
    ASSERT_DOUBLE_EQ(coefs[0], x * x);
    ASSERT_DOUBLE_EQ(coefs[1], x * y);
    ASSERT_DOUBLE_EQ(coefs[2], x);
    ASSERT_DOUBLE_EQ(coefs[3], y * y);
    ASSERT_DOUBLE_EQ(coefs[4], y);
    ASSERT_DOUBLE_EQ(coefs[5], 1);
}

TEST_F(PolyTest, DaDxiTest) {
    std::vector<double> X = {1.0, 2.0};
    auto coefs = std::vector<double>(1);

    auto r = n2o2Zeros.DaDxi(0, &X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.DaDxi(0, &X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0, 2.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.DaDxi(-1, &X, &coefs);
    ASSERT_TRUE(r.isError);
    r = n2o2Zeros.DaDxi(2, &X, &coefs);
    ASSERT_TRUE(r.isError);

    // Should not change any of the results,
    // as we're testing the derivatives with respect to the coeffs
    n2o2.Multiply(5);
    double x = 9.0;
    double y = 17.0;
    X = {x, y};

    auto eval = n2o2.DaDxi(0, &X, &coefs);
    ASSERT_FALSE(eval.isError);
    // P = 1 * x^2 + 2 * xy + 3 * x + 4 * y^2 + 5 * y + 6 * 1
    // Dx0 = d/dx(P) = 1 * 2x + 2 * y + 3 * 1
    // DaDx0 = [2*x, y, 1]
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(coefs[0], 2.0 * x);
    ASSERT_DOUBLE_EQ(coefs[1], y);
    ASSERT_DOUBLE_EQ(coefs[2], 1.0);
    ASSERT_DOUBLE_EQ(coefs[3], 0);
    ASSERT_DOUBLE_EQ(coefs[4], 0);
    ASSERT_DOUBLE_EQ(coefs[5], 0);

    eval = n2o2.DaDxi(1, &X, &coefs);
    ASSERT_FALSE(eval.isError);
    // P = 1 * x^2 + 2 * xy + 3 * x + 4 * y^2 + 5 * y + 6 * 1
    // Dx1 = d/dy(P) = 1 * 0 + 2 * x + 3 * 0 + 4 * 2y + 5 * 1 + 6 * 0
    // DaDx0 = [0, x, 0, 2y, 1, 0]
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(coefs[0], 0);
    ASSERT_DOUBLE_EQ(coefs[1], x);
    ASSERT_DOUBLE_EQ(coefs[2], 0);
    ASSERT_DOUBLE_EQ(coefs[3], 2 * y);
    ASSERT_DOUBLE_EQ(coefs[4], 1);
    ASSERT_DOUBLE_EQ(coefs[5], 0);
}

TEST_F(PolyTest, DaD2xiTest) {
    std::vector<double> X = {1.0, 2.0};
    auto coefs = std::vector<double>(1);

    auto r = n2o2Zeros.DaD2xi(0, &X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.DaD2xi(0, &X, &coefs);
    ASSERT_TRUE(r.isError);

    X = {1.0, 2.0};
    coefs = std::vector<double>(6);
    r = n2o2Zeros.DaD2xi(-1, &X, &coefs);
    ASSERT_TRUE(r.isError);
    r = n2o2Zeros.DaD2xi(2, &X, &coefs);
    ASSERT_TRUE(r.isError);

    // Should not change any of the results,
    // as we're testing the derivatives with respect to the coeffs
    n2o2.Multiply(5);
    double x = 9.0;
    double y = 17.0;
    X = {x, y};

    auto eval = n2o2.DaD2xi(0, &X, &coefs);
    ASSERT_FALSE(eval.isError);
    // P = 1 * x^2 + 2 * xy + 3 * x + 4 * y^2 + 5 * y + 6 * 1
    // Dx0 = d/dx(P) = 1 * 2x + 2 * y + 3 * 1
    // D2x0 = d2/dx2(P) = 1 * 2
    // DaDx0 = [2, 0, 0, 0, 0, 0]
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(coefs[0], 2.0);
    ASSERT_DOUBLE_EQ(coefs[1], 0);
    ASSERT_DOUBLE_EQ(coefs[2], 0);
    ASSERT_DOUBLE_EQ(coefs[3], 0);
    ASSERT_DOUBLE_EQ(coefs[4], 0);
    ASSERT_DOUBLE_EQ(coefs[5], 0);

    eval = n2o2.DaD2xi(1, &X, &coefs);
    ASSERT_FALSE(eval.isError);
    // P = 1 * x^2 + 2 * xy + 3 * x + 4 * y^2 + 5 * y + 6 * 1
    // Dx1 = d/dy(P) = 1 * 0 + 2 * x + 3 * 0 + 4 * 2y + 5 * 1 + 6 * 0
    // D2x1 = d2/dy2(P) = 1 * 0 + 2 * 0 + 3 * 0 + 4 * 2 + 5 * 0 + 6 * 0
    // DaDx0 = [0, x, 0, 2y, 1, 0]
    ASSERT_FALSE(eval.isError);
    ASSERT_DOUBLE_EQ(coefs[0], 0);
    ASSERT_DOUBLE_EQ(coefs[1], 0);
    ASSERT_DOUBLE_EQ(coefs[2], 0);
    ASSERT_DOUBLE_EQ(coefs[3], 2);
    ASSERT_DOUBLE_EQ(coefs[4], 0);
    ASSERT_DOUBLE_EQ(coefs[5], 0);
}
TEST_F(PolyTest, GetCoefficientsTest) {
    auto coefs = std::vector<double>(1);
    auto r = n2o2Zeros.GetCoefficients(&coefs);
    ASSERT_TRUE(r.isError);

    coefs = std::vector<double>(7);
    r = n2o2Zeros.GetCoefficients(&coefs);
    ASSERT_TRUE(r.isError);

    coefs = std::vector<double>(6);

    n2o2Zeros.Multiply(10);
    r = n2o2Zeros.GetCoefficients(&coefs);
    ASSERT_FALSE(r.isError);
    // P = 0*x^2 + 0*xy + 0*x + 0*y^2 + 0*y + 0
    ASSERT_DOUBLE_EQ(coefs[0], 0);
    ASSERT_DOUBLE_EQ(coefs[1], 0);
    ASSERT_DOUBLE_EQ(coefs[2], 0);
    ASSERT_DOUBLE_EQ(coefs[3], 0);
    ASSERT_DOUBLE_EQ(coefs[4], 0);
    ASSERT_DOUBLE_EQ(coefs[5], 0);

    n2o2Twos.Multiply(17);
    r = n2o2Twos.GetCoefficients(&coefs);
    ASSERT_FALSE(r.isError);
    // P = 17*(2*x^2 + 2*xy + 2*x + 2*y^2 + 2*y + 2)
    ASSERT_DOUBLE_EQ(coefs[0], 2 * 17);
    ASSERT_DOUBLE_EQ(coefs[1], 2 * 17);
    ASSERT_DOUBLE_EQ(coefs[2], 2 * 17);
    ASSERT_DOUBLE_EQ(coefs[3], 2 * 17);
    ASSERT_DOUBLE_EQ(coefs[4], 2 * 17);
    ASSERT_DOUBLE_EQ(coefs[5], 2 * 17);

    n2o2.Multiply(19);
    r = n2o2.GetCoefficients(&coefs);
    ASSERT_FALSE(r.isError);
    // P = 19*(1*x^2 + 2*xy + 3*x + 3*y^2 + 5*y + 6)
    ASSERT_DOUBLE_EQ(coefs[0], 1 * 19);
    ASSERT_DOUBLE_EQ(coefs[1], 2 * 19);
    ASSERT_DOUBLE_EQ(coefs[2], 3 * 19);
    ASSERT_DOUBLE_EQ(coefs[3], 4 * 19);
    ASSERT_DOUBLE_EQ(coefs[4], 5 * 19);
    ASSERT_DOUBLE_EQ(coefs[5], 6 * 19);
}