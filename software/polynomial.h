#pragma once

#include <gtest/gtest.h>
#include <math.h>

#include <vector>

#include "maybe.h"

class Node {
    friend class Poly;
    /*     Node of tree used to represent coefficients of multivariate
       polynomial The height of each node represents which variable this node
       refers to.

        Ex:
        root
        |   \
        N    N -> these nodes represent powers of the first variable
        |\   |\
        N N  N N -> these nodes represent powers of the second variable

        The terms are obtained by traversing the tree
        Ex the representation of P = 1*x^2 + 2*x*y + 3*x + 4*y^2 + 5*y + 6:
            root
        /     |     \
        N0   N1      N2
        |    |  \   | \  \
        N3   N4  N5 N6 N7 N8

        N0.exp = 2
        N3.exp = 0, N3.a = 1

        N1.exp = 1
        N4.exp = 1, N4.a = 2
        N5.exp = 0, N5.a = 3

        N2.exp = 0
        N6.exp = 2, N6.a = 4
        N7.exp = 1, N7.a = 5
        N8.exp = 0, N8.a = 6

        x^2: root.children[0].children[0]
        xy: root.children[1].children[0]
        x: root.children[1].children[1]
        y^2: root.children[2].children[0]
        y: root.children[2].children[1]
        1: root.children[2].children[2] */

   public:
    int exp;            // exponent of this node
    int parentsExpSum;  // sum of the exps of parents

    double a;  // coefficient of this monomial; only present at leaf node

    // child nodes
    // has constant size after Node construction
    std::vector<Node*> children;

    // Empty constructor
    Node(){};

    // Root node constructor
    Node(int nChildren);

    // Non-leaf node constructor
    Node(int exp, int parentsExpSum, int nChildren);

    // Leaf node constructor
    Node(int exp, int parentsExpSum, double coefficient);

    bool IsLeaf() { return this->children.size() == 0; }

   private:
    FRIEND_TEST(NodeTest, AddChildTest);
    int childrenI;  // index of first non empty index at children vector
    Node* AddChild(int exp, int parentsExpSum, int nChildren);
};

class Poly {
    // Class that represents a multivariate polynomial
    // These polynomials are sum of terms in which the sum of the exponents
    // is <= than a certain number - the order of the polynomial.
    //
    // For example, consider a polynomial of 2 variables, x and y
    // Let a,b,c,d,e,f be constants, the polynomials of order 0, 1 and 2 are:
    // P_0 = a
    // P_1 = a*x + b*y + c
    // P_2 = a*x^2 + b*x*y + c*x + d*y^2 + e*y + f
    //
    // Through all the comments related to this class, the constants
    // (a,b,c...) are considered 1 for simplicity.
    // With this simplification, the P_i above become:
    // P_0 = 1
    // P_1 = x + y + 1
    // P_2 = x^2 + x*y + x + y^2 + y + 1
   public:
    // Number of variables
    // P = x^2 + xy + x + y^2 + y + 1 -> n = 2
    // P = x + y + z + 1 -> n = 3
    int n;

    // Highest order of the polynomial
    // P = x^2 + xy + x + y^2 + y + 1 -> order = 2
    // P = x + y + z + 1 -> oder = 1
    int order;
    int nTerms;  // number of terms the polynomial has

    // Constructor
    // Arguments are the number of variables, the order of the polynomial
    // and the value for the coefficients.
    // Ex:
    // NewPoly(2,2,2.0) returns an instance that represents:
    // P = 2.0*x^2 + 2.0*xy + 2.0*x + 2.0*y^2 + 2.0*y + 2.0
    static Maybe<Poly> NewPoly(int n, int order, double coefficients);
    static Maybe<Poly> NewPoly(int n, int order) {
        return NewPoly(n, order, 0.0);
    };

    // Evaluates the polynomial for the given values of the variables
    Maybe<double> operator()(std::vector<double>* X);

    // Returns the derivative with respect to i-th variable
    // Ex:
    // P = x^2 + xy + x + y^2 + y + 1
    // Dxi(0, {3.0, 0.0}) = 2*x + y + 1 = 2*3.0 + 0.0 + 1
    Maybe<double> Dxi(int i, std::vector<double>* X);

    // Returns the second derivative with respect to the i-th variable,
    // similar to Poly::Dxi
    Maybe<double> D2xi(int i, std::vector<double>* X);

    // Get the derivatives with respect to each coefficient
    // Ex:
    // X = {7,8}
    // P = 1.0x^2 + 2.0xy + 3.0x + 4.0y^2 + 5.0y + 6.0
    // target will be set to {7^2, 7*8, 7, 8^2, 8, 1.0}
    Maybe<Void> Da(std::vector<double>* X, std::vector<double>* target);

    // Get the derivatives with respect to each coefficient of the derivative
    // with respect to the i-th variable
    // Ex:
    // P = x^2 + xy + x + y^2 + y + 1
    // Dx0 = 2*x + y + 1
    // DaDx0 = [x, y, 1, 0, 0, 0]
    Maybe<Void> DaDxi(int i, std::vector<double>* X,
                      std::vector<double>* target);

    // Get the derivatives with respect to each coefficient of the second
    // derivative with respect to the i-th variable
    // Ex: P = x^2 + xy + x + y^2 + y + 1
    // Dx0 = 2*x + y + 1
    // D2x0 = 2
    // DaD2xi = [2, 0, 0, 0, 0, 0]
    Maybe<Void> DaD2xi(int i, std::vector<double>* X,
                       std::vector<double>* target);

    // Sets the coefficients at target vector
    // Ex:
    // P = 1.0x^2 + 2.0xy + 3.0x + 4.0y^2 + 5.0y + 6.0
    // target will be set to {1,2,3,4,5,6}
    Maybe<Void> GetCoefficients(std::vector<double>* target);

    // Sets the coefficients of the polynomial with the values provided
    // Ex:
    // coefficients = {5,6,7,8,9,10}
    // P will be set to:
    // P = 5x^2 + 6xy + 7x + 8y^2 + 9y + 10*1
    Maybe<Void> SetCoefficients(std::vector<double>* coefficients);

    // Multiplies all coefficients by k
    // Ex:
    // P = 1*x^2 + 1*x + 1*1
    // P.Multiply(9)
    // P = 9*x^2 + 9*x + 9*1
    void Multiply(double k);

   private:
    friend class PolyTest;
    FRIEND_TEST(PolyConstructorTest, constructorTest);
    FRIEND_TEST(PolyTest, dfsTest);
    Node* coefficients;

    // vector to quickly access the leaf nodes (which contain the coefficients)
    std::vector<Node*> leafNodes;

    // Number that multiplies this polynomial
    // Note that we consider that this is part of each coefficient
    // Ex:
    // P = 7*(3x^2 + 4x + 5)
    // The coefficients are: 7*3, 7*4, 7*5
    double k;
};

double dfs(Node* node, int nodeTreeDepth, std::vector<double>* X);