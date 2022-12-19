#pragma once

#include <gtest/gtest.h>
#include <math.h>

#include <cassert>
#include <memory>
#include <vector>

#include "maybe.h"

class Monomial {
    // Represents k*(a*x0^e0*x1^e1*)...
    // Note: k is the scalar that changes when we multiply this monomial by a
    // scalar
   public:
    // [e0, e1, ...]
    std::vector<int> exponents;
    double a;
    double k;

    Monomial(double a, std::vector<int> exponents);
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
    int nMonomials();  // number of monomials the polynomial has

    // Constructor
    Poly();
    // Same as empty constructor, only exists for compatibility reasons
    Poly(int);

    // Build method to be called after the constructor
    // Arguments are the number of variables, the order of the polynomial
    // and the value for the coefficients.
    // Ex:
    // Build(2,2,2.0) will set:
    // P = 2.0*x^2 + 2.0*xy + 2.0*x + 2.0*y^2 + 2.0*y + 2.0
    Maybe<Void> Build(int n, int order, double coefficients);
    Maybe<Void> Build(int n, int order) { return Build(n, order, 0.0); };

    // Evaluates the polynomial for the given values of the variables
    Maybe<double> operator()(std::vector<double>* X);

    // Sets this polynomial as the derivative of this polynomial
    // with respect to i-th variable.
    // Ex:
    // P = x^2 + xy + x + y^2 + y + 1
    // Dxi(0) -> P will be set to: 2x + y + 1 + 0 + 0 + 0
    Maybe<Void> Dxi(int i);

    // Get the derivatives with respect to each coefficient
    // Ex:
    // X = {7,8}
    // P = 1.0x^2 + 2.0xy + 3.0x + 4.0y^2 + 5.0y + 6.0
    // target will be set to {7^2, 7*8, 7, 8^2, 8, 1.0}
    Maybe<Void> Da(std::vector<double>* X, std::vector<double>* target);

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

    friend void buildDfs(std::vector<int>* exponents, int exponentsSum,
                         Poly& p);

    friend Poly operator+(Poly const& left, Poly const& right);
    friend Poly operator*(double x, const Poly& p);
    friend Poly operator*(const Poly& p, double x);
    friend Poly operator+(double x, const Poly& p);
    friend Poly operator+(const Poly& p, double x);
    Poly& operator+=(const Poly& right);
    bool operator==(Poly const& right);
    bool operator!=(Poly const& right);

   private:
    // Auxiliary function used in Build method.
    void buildDfs(std::vector<int>& exponents, int exponentsSum,
                  int indexAtExponents, double coefficientToSet);

    friend class PolyTest;
    FRIEND_TEST(PolyConstructorTest, constructorTest);
    FRIEND_TEST(PolyTest, copyConstructorAndAssignment);
    FRIEND_TEST(PolyTest, MultiplicationOperatorTest);
    FRIEND_TEST(PolyTest, PlusAndMultiplicationTest);

    // Indicates Build method hasn't been called on this entity
    bool isZero;

    std::vector<Monomial> monomials;
};

// double dfs(std::shared_ptr<Node> node, int nodeTreeDepth,
//            std::vector<double>* X);
Poly operator*(double x, const Poly& p);
Poly operator*(const Poly& p, double x);
Poly operator+(double x, const Poly& p);
Poly operator+(const Poly& p, double x);
