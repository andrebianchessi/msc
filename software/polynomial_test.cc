#include "polynomial.h"

#include <gtest/gtest.h>

#include "maybe.h"

TEST(PolyTest, NodeTest) {
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

TEST(PolyTest, constructorTest) {
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

TEST(PolyTest, dfsTest) {
    double x = 2.0;
    double y = 3.0;
    double a = 5.0;

    Maybe<Poly> m = Poly::NewPoly(2, 2, a);
    ASSERT_FALSE(m.isError);
    Poly p = m.val;

    std::vector<double> X = {x, y};
    double a_y_0 = dfs(p.coefficients->children[0]->children[0], 2, &X);
    ASSERT_DOUBLE_EQ(a_y_0, a * 1);
    double x_2_a_y_0 = dfs(p.coefficients->children[0], 1, &X);
    ASSERT_DOUBLE_EQ(x_2_a_y_0, x * x * a * 1);
}

TEST(PolyTest, operatorTest) {
    Maybe<Poly> m = Poly::NewPoly(2, 2);
    ASSERT_FALSE(m.isError);
    Poly p = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = p(&X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    // x1^2*x2^0
    p.coefficients->children[0]->children[0]->a = 1.0;

    // x1^1*x2^1
    p.coefficients->children[1]->children[0]->a = 2.0;

    // x1^1*x2^0
    p.coefficients->children[1]->children[1]->a = 3.0;

    // x1^0*x2^2
    p.coefficients->children[2]->children[0]->a = 4.0;

    // x1^0*x2^1
    p.coefficients->children[2]->children[1]->a = 5.0;

    // x1^0*x2^0
    p.coefficients->children[2]->children[2]->a = 6.0;

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    v = p(&X);
    ASSERT_FALSE(v.isError);

    double expected = 1.0 * x1 * x1 + 2.0 * x1 * x2 + 3.0 * x1 + 4.0 * x2 * x2 +
                      5.0 * x2 + 6.0;
    ASSERT_DOUBLE_EQ(v.val, expected);
}

TEST(PolyTest, DxiTest) {
    Maybe<Poly> m = Poly::NewPoly(2, 2);
    ASSERT_FALSE(m.isError);
    Poly p = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = p.Dxi(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);
    v = p.Dxi(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(p.Dxi(2, &X).isError);
    ASSERT_TRUE(p.Dxi(-1, &X).isError);

    // x1^2*x2^0
    p.coefficients->children[0]->children[0]->a = 1.0;

    // x1^1*x2^1
    p.coefficients->children[1]->children[0]->a = 2.0;

    // x1^1*x2^0
    p.coefficients->children[1]->children[1]->a = 3.0;

    // x1^0*x2^2
    p.coefficients->children[2]->children[0]->a = 4.0;

    // x1^0*x2^1
    p.coefficients->children[2]->children[1]->a = 5.0;

    // x1^0*x2^0
    p.coefficients->children[2]->children[2]->a = 6.0;

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    v = p.Dxi(0, &X);
    ASSERT_FALSE(v.isError);
    double d1Expected = 1.0 * 2 * x1 + 2.0 * x2 + 3.0;
    ASSERT_DOUBLE_EQ(v.val, d1Expected);

    v = p.Dxi(1, &X);
    ASSERT_FALSE(v.isError);
    double d2Expected = 2.0 * x1 + 4.0 * 2 * x2 + 5.0;
    ASSERT_DOUBLE_EQ(v.val, d2Expected);
}

TEST(PolyTest, DDxiTest) {
    Maybe<Poly> m = Poly::NewPoly(2, 2);
    ASSERT_FALSE(m.isError);
    Poly p = m.val;

    std::vector<double> X = {1.0, 2.0};
    auto v = p.DDxi(0, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);
    v = p.DDxi(1, &X);
    ASSERT_FALSE(v.isError);
    ASSERT_DOUBLE_EQ(v.val, 0.0);

    ASSERT_TRUE(p.DDxi(2, &X).isError);
    ASSERT_TRUE(p.DDxi(-1, &X).isError);

    // x1^2*x2^0
    p.coefficients->children[0]->children[0]->a = 1.0;

    // x1^1*x2^1
    p.coefficients->children[1]->children[0]->a = 2.0;

    // x1^1*x2^0
    p.coefficients->children[1]->children[1]->a = 3.0;

    // x1^0*x2^2
    p.coefficients->children[2]->children[0]->a = 4.0;

    // x1^0*x2^1
    p.coefficients->children[2]->children[1]->a = 5.0;

    // x1^0*x2^0
    p.coefficients->children[2]->children[2]->a = 6.0;

    double x1 = 9.0;
    double x2 = 17.0;
    X = {x1, x2};

    v = p.DDxi(0, &X);
    ASSERT_FALSE(v.isError);
    double dd1Expected = 1.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd1Expected);

    v = p.DDxi(1, &X);
    ASSERT_FALSE(v.isError);
    double dd2Expected = 4.0 * 2;
    ASSERT_DOUBLE_EQ(v.val, dd2Expected);
}
