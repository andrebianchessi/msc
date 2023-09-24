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
    this->k = 1;
    this->exps = exps;
}

Maybe<double> Monomial::operator()(double a, std::vector<double>& X) const {
    Maybe<double> r;
    if (X.size() != this->exps.size()) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }

    double prod = this->k * a;
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
        this->k *= this->exps.at(i);
        this->exps.at(i) -= 1;
        return r;
    }
    this->k = 0;
    return r;
}

Maybe<double> Monomial::Da(std::vector<double>& X) const {
    Maybe<double> r;
    if (X.size() != this->exps.size()) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    double prod = this->k;
    for (int e = 0; e < int(this->exps.size()); e++) {
        if (prod == 0) {
            // small optimization
            break;
        }
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

    int exponentsLeft;
    if (this->tiny) {
        // For the first variable we use all exponents
        if (exponentsIndex == 0) {
            exponentsLeft = this->order - exponentsSum;
        } else {
            // When the first variable has exponent 0, we only take
            // exponent 0 of all the other variables
            if (exponents[0] == 0) {
                exponentsLeft = 0;
            } else {
                // If another variable other than the first one has non-zero exp
                if (exponents[0] != exponentsSum) {
                    exponentsLeft = 0;
                } else {
                    if (exponents[0] < order) {
                        exponentsLeft = 1;
                    } else {
                        exponentsLeft = 0;
                    }
                }
            }
        }
    } else {
        exponentsLeft = this->order - exponentsSum;
    }

    for (int e = exponentsLeft; e >= 0; e--) {
        exponents[exponentsIndex] = e;
        buildDfs(exponents, exponentsSum + e, exponentsIndex + 1);
    }
}
Maybe<Void> Poly::Build(int n, int order, bool tiny, int id) {
    Maybe<Void> r;
    if (n < 0 || order < 0) {
        r.isError = true;
        r.errMsg = "n and order must be >=0";
        return r;
    }

    this->id = id;
    this->dxCount = 0;
    this->n = n;
    this->X = std::vector<double>(n);

    this->order = order;
    this->tiny = tiny;

    int _nMonomials = 0;
    if (!tiny) {
        // Source:
        // https://mathoverflow.net/questions/225953/number-of-polynomial-terms-for-certain-degree-and-certain-number-of-variables/225963#225963?newreg=2a0208ceb740461d8eaa21e304b0e341
        // Note:
        _nMonomials = int(boost::math::binomial_coefficient<double>(
            this->n + this->order, this->order));
    }

    this->monomials = std::vector<Monomial>();
    this->monomials.reserve(_nMonomials);

    std::vector<int> exponents = std::vector<int>(n);
    this->buildDfs(exponents, 0, 0);

    this->isZero = false;
    return r;
};

std::vector<double> Poly::GetX() const { return this->X; }
void Poly::SetX(std::vector<double> X) {
    assert(int(X.size()) == this->n);
    this->X = X;
}

int Poly::nMonomials() const { return this->monomials.size(); }

Maybe<double> Poly::operator()(std::vector<double>& a) const {
    Maybe<double> r;
    Maybe<double> maybeVal;
    std::vector<double> X = this->GetX();
    double val = 0;
    for (int m = 0; m < int(this->monomials.size()); m++) {
        maybeVal = this->monomials[m](a[m], X);
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

Poly& Poly::operator*=(double s) {
    for (int i = 0; i < this->monomials.size(); i++) {
        this->monomials[i].k *= s;
    }
}

Polys::Polys() {
    this->polys = std::vector<Poly>(0);
    this->k = std::vector<double>(0);
    this->plus = 0;
}

Polys::Polys(const Poly& p) {
    if (p.isZero) {
        this->polys = std::vector<Poly>(0);
        this->k = std::vector<double>(0);
        this->plus = 0;
        return;
    }
    this->polys = std::vector<Poly>{p};
    this->k = std::vector<double>{1.0};
    this->plus = 0;
}

int Polys::nMonomials() const {
    int n = 0;
    for (int i = 0; i < int(this->polys.size()); i++) {
        n += this->polys[i].nMonomials();
    }
    return n;
}

Polys operator*(double k, const Poly& p) {
    if (k == 0) {
        Polys ps;
        return ps;
    }
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

Polys operator+(Poly const& p, double k) {
    if (k == 0) {
        return Polys(p);
    }
    Polys copy = Polys(p);
    copy.plus += k;
    return copy;
}
Polys operator+(double k, Poly const& p) { return p + k; }

Polys& Polys::operator+=(const Poly& right) {
    (*this) = (*this) + right;
    return (*this);
}
Polys& Polys::operator+=(const Polys& right) {
    (*this) = (*this) + right;
    return (*this);
}
Polys& Polys::operator*=(double k) {
    if (k == 0) {
        this->polys.clear();
        this->k.clear();
        this->plus = 0;
        return (*this);
    }
    for (int i = 0; i < int(this->k.size()); i++) {
        this->k[i] *= k;
    }
    return (*this);
}

Polys operator*(double k, Polys const& right) {
    if (k == 0) {
        return Polys();
    }
    Polys p = right;
    for (int i = 0; i < int(p.k.size()); i++) {
        p.k[i] *= k;
    }
    return p;
}

bool operator==(Polys const& right, Polys const& left) {
    if (right.polys.size() != left.polys.size()) {
        return false;
    }
    for (int i = 0; i < int(right.polys.size()); i++) {
        if (right.k[i] != left.k[i]) {
            return false;
        }
        if (right.polys[i].id != left.polys[i].id) {
            return false;
        }
        if (right.polys[i].n != left.polys[i].n) {
            return false;
        }
        if (right.polys[i].order != left.polys[i].order) {
            return false;
        }
        if (right.polys[i].nMonomials() != left.polys[i].nMonomials()) {
            return false;
        }
        for (int m = 0; m < right.polys[i].nMonomials(); m++) {
            if (right.polys[i].monomials[m].exps !=
                left.polys[i].monomials[m].exps) {
                return false;
            }
            if (right.polys[i].monomials[m].k != left.polys[i].monomials[m].k) {
                return false;
            }
        }
    }
    return true;
}

bool operator!=(Polys const& right, Polys const& left) {
    return !(right == left);
}

Maybe<double> Polys::operator()(std::vector<std::vector<double>>& a) const {
    Maybe<double> r;
    double val = this->plus;
    for (int p = 0; p < int(this->polys.size()); p++) {
        if (this->polys[p].id >= int(a.size())) {
            r.isError = true;
            r.errMsg = "a[i] must contain the coefficients of poly with id = i";
            return r;
        }
        r = this->polys[p](a[this->polys[p].id]);
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
    this->dxCount += 1;
    return r;
}

Maybe<double> Poly::Dai(int i) const {
    Maybe<double> r;
    if (i < 0 || i >= this->nMonomials()) {
        r.isError = true;
        r.errMsg = "invalid i";
        return r;
    }
    std::vector<double> X = this->GetX();
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

std::map<std::tuple<int, int>, double> Polys::Da() const {
    std::map<std::tuple<int, int>, double> map;
    Maybe<double> eval;
    std::tuple<int, int> t;
    int pId;
    for (int p = 0; p < int(this->polys.size()); p++) {
        for (int m = 0; m < this->polys[p].nMonomials(); m++) {
            pId = this->polys[p].id;
            t = std::tuple<int, int>{pId, m};
            if (map.find(t) == map.end()) {
                map.insert({t, 0.0});
            }
            eval = this->polys[p].Dai(m);
            assert(!eval.isError);
            map[t] += this->k[p] * eval.val;
        }
    }
    return map;
}
