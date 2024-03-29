#pragma once

#include <gtest/gtest.h>
#include <math.h>

#include <cassert>
#include <memory>
#include <tuple>
#include <vector>

#include "maybe.h"

class Monomial {
    // Represents a monomial in the form:
    // a*x0^e0*x1^e1*...*xi^ei

   public:
    double k;  // scalar that multiplies this entity

    // [e0, e1, ...]
    std::vector<int> exps;

    Monomial(std::vector<int> exps);

    Maybe<double> operator()(double a, std::vector<double>& X) const;

    // Differentiate with respect to xi
    Maybe<Void> Dxi(int i);

    // Value of derivative with respect to a
    Maybe<double> Da(std::vector<double>& X) const;
};

class Polys;

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
    // Identifier
    // Used to order polynomial coefficients when calculating gradient of
    // linear combination of Poly instances (see Polys class)
    int id;
    // Used to identify the number of times this instance was differentiated,
    // i.e. how many times Dxi was called. Mostly used for debugging.
    int dxCount;

    // Number of variables
    // P = x^2 + xy + x + y^2 + y + 1 -> n = 2
    // P = x + y + z + 1 -> n = 3
    int n;

    // Highest order of the polynomial
    // P = x^2 + xy + x + y^2 + y + 1 -> order = 2
    // P = x + y + z + 1 -> oder = 1
    int order;
    int nMonomials() const;  // number of monomials the polynomial has
    // Used when we want to reduce the number of therms in the polynomial.
    // See Build method for more details.
    bool tiny;

    // Constructor
    Poly();
    Poly(int);  // Same as default constructor, only exists for compatibility

    // Build method to be called after the constructor
    // Arguments are the number:
    // n: number of inputs
    // order: order of polynomial
    // tiny: parameter used to reduce the number of monomials. When set true,
    // the resulting polynomial will only contain the monomials in which the
    // first input variable has exponent in [order, 0], and the other variables
    // have either exponent 1 or zero. When the first variable has
    // exponent zero, all the other ones also do so. See the constructor tests
    // for more details.
    // id: identifier of polynomial (only used
    // when doing linear combination of polynomials)
    Maybe<Void> Build(int n, int order, bool tiny, int id);
    Maybe<Void> Build(int n, int order) {
        return this->Build(n, order, false, 0);
    };

    // Getter and setter for X, which hold the values of the variables in which
    // this instance is evaluated
    std::vector<double> GetX() const;
    void SetX(std::vector<double> X);

    // Evaluates the polynomial for the given values of its coefficients and
    // X value set. Complexity is O(nMonomials)
    Maybe<double> operator()(std::vector<double>& a) const;

    // Multiply by scalar
    Poly& operator*=(double s);

    // Sets this polynomial as the derivative of this polynomial
    // with respect to i-th variable.
    // Ex:
    // P = x^2 + xy + x + y^2 + y + 1
    // Dxi(0) -> P will be set to: 2x + y + 1 + 0 + 0 + 0
    Maybe<Void> Dxi(int i);

    // Get the derivative with respect to the i-th coefficient
    // Ex:
    // X = {7,8}
    // P = 1.0x^2 + 2.0xy + 3.0x + 4.0y^2 + 5.0y + 6.0
    // Da0 = 7^2
    // Da1 = 7*8
    Maybe<double> Dai(int i) const;

   private:
    std::vector<double> X;

    FRIEND_TEST(PolyConstructorTest, constructorTest);
    FRIEND_TEST(PolyTest, GetSetXTest);
    FRIEND_TEST(PolyTest, copyConstructorAndAssignment);
    FRIEND_TEST(PolyTest, DxiTest2);
    FRIEND_TEST(PolyTest, PolysConstructorTest);
    FRIEND_TEST(PolyTest, MatrixMultiplicationTest);
    FRIEND_TEST(PimodelTest, getAccelsFromDiffEqTest);
    FRIEND_TEST(PolyTest, DaPolysTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTest);
    FRIEND_TEST(PimodelTest, OperatorTest);
    FRIEND_TEST(PolyTest, scalarMultiplicationTest);

    // Auxiliary function used in Build method.
    void buildDfs(std::vector<int>& exponents, int exponentsSum,
                  int exponentsIndex);

    friend class Polys;
    friend class PolyTest;
    friend class Pimodel;
    friend bool operator==(Polys const& right, Polys const& left);
    friend Polys operator+(Poly const& p, double k);
    friend Polys operator+(double k, Poly const& p);

    // Indicates Build still wasn't called yet
    bool isZero;

    std::vector<Monomial> monomials;
};

class Polys {
    // Represents linear combination of Poly instances
    // It's assumed that all polys have the same number of inputs
   public:
    std::vector<Poly> polys;
    std::vector<double> k;
    double plus;  // Scalar that is added to the Polys instance. Defaults to 0
    Polys();
    Polys(const Poly& p);

    int nMonomials() const;

    Polys& operator+=(const Poly& right);
    Polys& operator+=(const Polys& right);
    Polys& operator*=(double k);

    // Operator complexity is O(nMonomials)
    Maybe<double> operator()(std::vector<std::vector<double>>& a) const;
    Maybe<Void> Dxi(int i);

    // Returns a map in which the key (i, j) represents the non-zero derivatives
    // of the polynomial with id=i with respect to the j-th coefficient
    std::map<std::tuple<int, int>, double> Da() const;
};

Polys operator*(double k, const Poly& p);
Polys operator*(const Poly& p, double k);
Polys operator+(Poly const& left, Poly const& right);
Polys operator+(Polys const& left, Poly const& right);
Polys operator+(Poly const& left, Polys const& right);
Polys operator+(Polys const& left, Polys const& right);
Polys operator*(double k, Polys const& right);

// Scalar sum
Polys operator+(Poly const& p, double k);
Polys operator+(double k, Poly const& p);

// Note: order of monomials is considered for equality
bool operator==(Polys const& right, Polys const& left);
bool operator!=(Polys const& right, Polys const& left);
