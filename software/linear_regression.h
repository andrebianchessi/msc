#pragma once

#include <gtest/gtest.h>
#include <math.h>

#include <vector>

#include "maybe.h"

class Node {
    /*     Tree node to represent polynomial terms
        Ex:
        y(x1,x2) = 1.0*x1^2 + 2.0*x1*x2 + 3.0*x2^2 + 4.0*x1 + 5.0*x2 + 6.0

        Represented as:

           Root
        /   \    \
        a0   a1   a2
        |    | \   |  \ \
        b0   b1 b2 b3 b4 b5

        a0={e:2, l:0}
        a1={e:1, l:1}
        a2={e:0, l:2}

        b0={e:0,l:0,a:1.0}
        b1={e:1,l:0,a:2.0}
        b2={e:0,l:1,a:4.0}
        b3={e:2,l:0,a:3.0}
        b4={e:1,l:1,a:5.0}
        b5={e:0,l:2,a:6.0} */
   public:
    int e;     // polynomial order of this node
    int l;     // polynomial order left for children of this node
    double a;  // coefficient of this monomial; only present at leaf node
    std::vector<Node*> children;

    bool IsLeaf() { return this->children.size() == 0; }
};

// Class that represents a linear regression
class LinReg {
   public:
    // Parameters that define the regression
    // Ex:
    // y(x1,x2) = 1.0*x1 + 3.0*x2 + 0.1
    //   XSize = 2
    //   order = 1
    //   nTerms = 3

    // y(x1,x2,x3) = 1.0*x1 + 3.0*x2 + 0.5*x3 + 7.0
    //   XSize = 3
    //   order = 1
    //   nTerms = 4

    // y(x1,x2) = 1.0*x1^2 + 0.1*x1*x2 + 1.2*x2^2 + 1.0*x1 + 3.0*x2 + 1.1
    //   XSize = 2
    //   order = 2
    //   nTerms = 6

    // y() = 7.0
    //  XSize = 0
    //  order = 0
    //  nTerms = 1
    int XSize;
    int order;
    int nTerms();

    // Initializes regression considering all coefficients = 0
    // Ex:
    // LinReg(2, 1) -> y(x1,x2) = 0*x1 + 0*x2 + 0
    static Maybe<LinReg> NewLinReg(int XSize, int order);
    Maybe<double> operator()(std::vector<double>* X);

    // Returns the component i of the gradient at X
    // Ex:
    // y(x1,x2) = 1.0*x1^2 + 3.0*x2^2 + 0.1
    // Grad({1,1}, 0) -> 1.0*2*x1 = 1*2*1 = 2
    // Grad({1,1}, 1) -> 3.0*2*x2 = 3*2*1 = 6
    Maybe<double> Grad(int di, std::vector<double>* X);

    // Returns the coefficient for the powers of the inputs
    // Ex:
    // y(x1,x2) = 1.0*x1^2 + 2.0*x1*x2 + 3.0*x2^2 + 4.0*2x1 + 5.0*x2 + 6.0
    // CoefficientAt({0,0}) -> 6.0
    // CoefficientAt({2,0}) -> 1.0
    // CoefficientAt({1,1}) -> 2.0
    Maybe<int> CoefficientAt(std::vector<int> powers);

   private:
    FRIEND_TEST(LinRegTest, coefficientsTest);
    FRIEND_TEST(LinRegTest, coefficientAtTest);
    FRIEND_TEST(LinRegTest, operatorTest);
    FRIEND_TEST(LinRegTest, gradTest);
    Node* coefficients;
};

double dfs(Node* n, int index, std::vector<double>* X);