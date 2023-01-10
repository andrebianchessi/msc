#include "polynomial.h"

#include <boost/math/special_functions/binomial.hpp>
#include <cassert>
#include <memory>
#include <queue>

#include "maybe.h"
#include "utils.h"

Monomial::Monomial(std::vector<int> exps) {
    for (auto e : exps) {
        assert(e >= 0);
    }
    this->a = 0;
    this->exps = exps;
}

Maybe<double> Monomial::operator()(std::vector<double>& X) const {
    Maybe<double> r;
    if (X.size() != this->exps.size()) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }

    double prod = this->a;
    for (int e = 0; e < int(this->exps.size()); e++) {
        prod *= pow(X.at(e), this->exps[e]);
    }

    r.val = prod;
    return r;
}

Maybe<Void> Monomial::Dxi(int i) {
    Maybe<Void> r;
    if (i < 0 || i >= int(this->exps.size())) {
        r.isError = true;
        r.errMsg = "Invalid i";
        return r;
    }

    if (this->exps.at(i) >= 1) {
        this->a *= this->exps.at(i);
        this->exps.at(i) -= 1;
        return r;
    }
    this->a = 0;
    return r;
}

Maybe<double> Monomial::Da(std::vector<double>& X) const {
    Maybe<double> r;
    if (X.size() != this->exps.size()) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    double prod = 1;
    for (int e = 0; e < int(this->exps.size()); e++) {
        prod *= pow(X.at(e), this->exps[e]);
    }

    r.val = prod;
    return r;
}

Poly::Poly() { this->isZero = true; }

Poly::Poly(int k) { this->isZero = true; }

void Poly::buildDfs(std::vector<int>& exponents, int exponentsSum,
                    int exponentsIndex) {
    if (exponentsIndex == int(exponents.size())) {
        this->monomials.push_back(Monomial(exponents));
        return;
    }

    int exponentsLeft = this->order - exponentsSum;
    for (int e = exponentsLeft; e >= 0; e--) {
        exponents[exponentsIndex] = e;
        buildDfs(exponents, exponentsSum + e, exponentsIndex + 1);
    }
}
Maybe<Void> Poly::Build(int n, int order, int id) {
    Maybe<Void> r;
    if (n < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "n and order must be >=0";
        return r;
    }

    this->id = id;
    this->n = n;
    this->order = order;
    // Source:
    // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
    int nMonomials = int(boost::math::binomial_coefficient<double>(
        this->n + this->order, this->order));

    this->monomials = std::vector<Monomial>();
    this->monomials.reserve(nMonomials);

    std::vector<int> exponents = std::vector<int>(n);
    this->buildDfs(exponents, 0, 0);

    this->isZero = false;
    return r;
};

int Poly::nMonomials() const { return this->monomials.size(); }

Maybe<double> Poly::operator()(std::vector<double>& X) const {
    Maybe<double> r;
    if (int(X.size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    Maybe<double> maybeVal;
    double val = 0;
    for (int m = 0; m < int(this->monomials.size()); m++) {
        maybeVal = this->monomials[m](X);
        if (maybeVal.isError) {
            r.isError = true;
            r.errMsg = maybeVal.errMsg;
            return r;
        }
        val += maybeVal.val;
    }
    r.val = val;
    return r;
};

Polys::Polys() {
    this->polys = std::vector<Poly>(0);
    this->k = std::vector<double>(0);
}

Polys::Polys(const Poly& p) {
    if (p.isZero) {
        this->polys = std::vector<Poly>(0);
        this->k = std::vector<double>(0);
        return;
    }
    this->polys = std::vector<Poly>{p};
    this->k = std::vector<double>{1.0};
}
Polys operator*(double k, const Poly& p) {
    Polys ps = Polys(p);
    ps.k[ps.k.size() - 1] *= k;
    return ps;
};
Polys operator*(const Poly& p, double k) { return k * p; };

Polys operator+(Poly const& left, Poly const& right) {
    assert(left.n == right.n);
    Polys p = Polys(left);
    p.k.push_back(1.0);
    p.polys.push_back(right);
    return p;
};
Polys operator+(Polys const& left, Poly const& right) {
    if (left.polys.size() > 0) {
        assert(left.polys[0].n == right.n);
    }
    Polys p = left;
    p.k.push_back(1.0);
    p.polys.push_back(right);
    return p;
};
Polys operator+(Poly const& left, Polys const& right) {
    if (right.polys.size() > 0) {
        assert(right.polys[0].n == left.n);
    }
    Polys p = right;
    p.k.push_back(1.0);
    p.polys.push_back(left);
    return p;
};
Polys operator+(Polys const& left, Polys const& right) {
    Polys p = left;
    for (int i = 0; i < int(right.polys.size()); i++) {
        p.polys.push_back(right.polys[i]);
        p.k.push_back(right.k[i]);
    }
    return p;
};

Polys& Polys::operator+=(const Poly& right) {
    (*this) = (*this) + right;
    return (*this);
}
Polys& Polys::operator+=(const Polys& right) {
    (*this) = (*this) + right;
    return (*this);
}
Polys& Polys::operator*=(double k) {
    for (int i = 0; i < int(this->k.size()); i++) {
        this->k[i] *= k;
    }
    return (*this);
}

Maybe<double> Polys::operator()(std::vector<double>& X) const {
    Maybe<double> r;
    double val = 0;
    for (int p = 0; p < int(this->polys.size()); p++) {
        r = this->polys[p](X);
        if (r.isError) {
            return r;
        }
        val += this->k[p] * r.val;
    }
    r.val = val;
    return r;
};

Maybe<Void> Poly::Dxi(int i) {
    Maybe<Void> r;
    if (i < 0 || i >= this->n) {
        r.isError = true;
        r.errMsg = "Invalid i";
        return r;
    }

    for (int m = 0; m < int(this->monomials.size()); m++) {
        r = this->monomials[m].Dxi(i);
        if (r.isError) {
            return r;
        }
    }

    return r;
}

Maybe<double> Poly::Dai(int i, std::vector<double>& X) const {
    Maybe<double> r;
    if (int(X.size()) != this->n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    if (i < 0 || i >= this->nMonomials()) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }

    Maybe<double> maybeDa = this->monomials[i].Da(X);
    if (maybeDa.isError) {
        return maybeDa;
    }
    r.val = maybeDa.val;
    return r;
}

Maybe<Void> Polys::Dxi(int i) {
    Maybe<Void> r;
    if (this->polys.size() == 0) {
        return r;
    }
    if (i < 0 || i >= this->polys[0].n) {
        r.isError = true;
        r.errMsg = "Invalid i";
        return r;
    }

    for (int p = 0; p < int(this->polys.size()); p++) {
        r = this->polys[p].Dxi(i);
        if (r.isError) {
            return r;
        }
    }

    return r;
}

Maybe<double> Polys::Dai(int pId, int i, std::vector<double>& X) const {
    Maybe<double> r;
    r.val = 0;

    Maybe<double> maybeVal;
    for (int pi = 0; pi < int(this->polys.size()); pi++) {
        if (this->polys[pi].id == pId) {
            maybeVal = this->polys[pi].Dai(i, X);
            if (maybeVal.isError) {
                return maybeVal;
            }
            r.val += this->k[pi] * maybeVal.val;
        }
    }
    return r;
}
