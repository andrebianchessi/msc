#include "polynomial.h"

#include <gtest/gtest.h>

#include "maybe.h"

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

    auto v = n2o2(&X);
    ASSERT_FALSE(v.isError);

    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_DOUBLE_EQ(v.val, expected);
}

TEST_F(PolyTest, DxiTest) {
    std::vector<double> X = {1.0, 2.0};
    auto v = n2o2Zeros.Dxi(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    v = n2o2Zeros.Dxi(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(n2o2Zeros.Dxi(2, &X).isError);
    ASSERT_TRUE(n2o2Zeros.Dxi(-1, &X).isError);

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    v = n2o2.Dxi(0, &X);
    ASSERT_FALSE(v.isError);
    double d1Expected = 1.0 * 2 * x1 + 2.0 * x2 + 3.0;
    ASSERT_DOUBLE_EQ(v.val, d1Expected);

    v = n2o2.Dxi(1, &X);
    ASSERT_FALSE(v.isError);
    double d2Expected = 2.0 * x1 + 4.0 * 2 * x2 + 5.0;
    ASSERT_DOUBLE_EQ(v.val, d2Expected);
}

TEST_F(PolyTest, DDxiTest) {
    std::vector<double> X = {1.0, 2.0};

    auto v = n2o2Zeros.DDxi(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);
    v = n2o2Zeros.DDxi(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(n2o2Zeros.DDxi(2, &X).isError);
    ASSERT_TRUE(n2o2Zeros.DDxi(-1, &X).isError);

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    v = n2o2.DDxi(0, &X);
    ASSERT_FALSE(v.isError);
    double dd1Expected = 1.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd1Expected);

    v = n2o2.DDxi(1, &X);
    ASSERT_FALSE(v.isError);
    double dd2Expected = 4.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd2Expected);
}
