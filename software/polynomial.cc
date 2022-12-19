#include "polynomial.h"

#include <boost/math/special_functions/binomial.hpp>
#include <cassert>
#include <memory>
#include <queue>

#include "maybe.h"
#include "utils.h"

Monomial::Monomial(double a, std::vector<int> exponents) {
    this->a = a;
    this->exponents = exponents;
    this->k = 1;
};

Poly::Poly() { this->isZero = true; }

Poly::Poly(int k) { this->isZero = true; }

void Poly::buildDfs(std::vector<int>& exponents, int exponentsSum,
                    int indexAtExponents, double coefficientToSet) {
    if (indexAtExponents == int(exponents.size())) {
        this->monomials.push_back(Monomial(coefficientToSet, exponents));
        return;
    }

    int exponentsLeft = this->order - exponentsSum;
    for (int e = exponentsLeft; e >= 0; e--) {
        exponents[indexAtExponents] = e;
        buildDfs(exponents, exponentsSum + e, indexAtExponents + 1,
                 coefficientToSet);
    }
}
Maybe<Void> Poly::Build(int n, int order, double coefficients) {
    Maybe<Void> r;
    if (n < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "n and order must be >=0";
        return r;
    }

    this->n = n;
    this->order = order;
    // Source:
    //
    // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
    int nMonomials = int(boost::math::binomial_coefficient<double>(
        this->n + this->order, this->order));

    this->monomials = std::vector<Monomial>();
    this->monomials.reserve(nMonomials);

    std::vector<int> exponents = std::vector<int>(n);
    this->buildDfs(exponents, 0, 0, coefficients);

    this->isZero = false;
    return r;
};

int Poly::nMonomials() { return this->monomials.size(); }

Maybe<double> Poly::operator()(std::vector<double>* X) {
    Maybe<double> r;
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    double val = 0;
    double monomialVal;
    for (int m = 0; m < int(this->monomials.size()); m++) {
        monomialVal = this->monomials[m].a * this->monomials[m].k;
        for (int i = 0; i < this->n; i++) {
            monomialVal *= std::pow(X->at(i), this->monomials[m].exponents[i]);
        }
        val += monomialVal;
    }
    r.val = val;
    return r;
};

Poly operator*(double x, const Poly& p) {
    Poly newP = p;
    for (int m = 0; m < int(newP.monomials.size()); m++) {
        newP.monomials[m].k *= x;
    }
    return newP;
}
Poly operator*(const Poly& p, double x) { return x * p; }

Poly operator+(Poly const& left, Poly const& right) {
    if (left.isZero) {
        return right;
    }
    if (right.isZero) {
        return left;
    }
    Poly newP = left;
    newP.order = std::max(left.order, right.order);
    newP.monomials.reserve(left.monomials.size() + right.monomials.size());
    for (int i = 0; i < int(right.monomials.size()); i++) {
        newP.monomials.push_back(right.monomials[i]);
    }
    return newP;
};
Poly& Poly::operator+=(const Poly& right) {
    if (this->isZero) {
        *this = right;
        return *this;
    }
    this->order = std::max(this->order, right.order);
    this->monomials.reserve(this->monomials.size() + right.monomials.size());
    for (int i = 0; i < int(right.monomials.size()); i++) {
        this->monomials.push_back(right.monomials[i]);
    }
    return *this;
}
Poly operator+(double x, const Poly& p) {
    assert(!p.isZero);
    Poly newP = p;
    newP.monomials[newP.monomials.size() - 1].a += x;
    return newP;
}
Poly operator+(const Poly& p, double x) {
    assert(!p.isZero);
    Poly newP = p;
    newP.monomials[newP.monomials.size() - 1].a += x;
    return newP;
}

bool Poly::operator==(Poly const& right) {
    if (right.n != this->n || right.order != this->order) {
        return false;
    }
    if (this->monomials.size() != right.monomials.size()) {
        return false;
    }
    for (int i = 0; i < int(this->monomials.size()); i++) {
        if (this->monomials[i].a != right.monomials[i].a) {
            return false;
        }
        if (this->monomials[i].exponents != right.monomials[i].exponents) {
            return false;
        }
    }
    return true;
}

bool Poly::operator!=(Poly const& right) { return !((*this) == right); }

Maybe<Void> Poly::GetCoefficients(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nMonomials()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    for (int i = 0; i < this->nMonomials(); i++) {
        (*target)[i] = this->monomials[i].a;
    }
    return r;
}

Maybe<Void> Poly::SetCoefficients(std::vector<double>* coefficients) {
    Maybe<Void> r;
    if (int(coefficients->size()) != this->nMonomials()) {
        r.isError = true;
        r.errMsg = "coefficients must have same length as the number of terms";
        return r;
    }
    for (int i = 0; i < this->nMonomials(); i++) {
        this->monomials[i].a = (*coefficients)[i];
    }
    return r;
}

Maybe<Void> Poly::Dxi(int i) {
    Maybe<Void> r;
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "Invalid i";
        return r;
    }

    int xiExp;
    for (int m = 0; m < int(this->monomials.size()); m++) {
        xiExp = this->monomials[m].exponents[i];
        // If this monomial doesn't depend on variable i,
        // this monomial becomes 0
        if (xiExp == 0) {
            for (int e = 0; e < int(this->monomials[m].exponents.size()); e++) {
                this->monomials[m].exponents[e] = 0;
            }
            this->monomials[m].a = 0;
            continue;
        }
        this->monomials[m].exponents[i] -= 1;
        this->monomials[m].a *= xiExp;
    }

    return r;
}

Maybe<Void> Poly::Da(std::vector<double>* X, std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nMonomials()) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }
    if (int(X->size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }

    for (int i = 0; i < this->nMonomials(); i++) {
        target->at(i) = 1;
        for (int xi = 0; xi < this->n; xi++) {
            target->at(i) *=
                std::pow(X->at(xi), this->monomials[i].exponents[xi]);
        }
    }
    return r;
}
