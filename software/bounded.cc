#include "bounded.h"

#include "maybe.h"

Bounded::Bounded() { this->val = 0.0; }

Maybe<Bounded> Bounded::CreateBounded(double x) {
    Maybe<Bounded> r;
    Bounded b = Bounded();
    auto e = b.Set(x);
    r.val = b;
    r.isError = e.isError;
    r.errMsg = e.errMsg;
    return r;
}

double Bounded::Get() { return this->val; }

std::vector<double> Bounded::Get(std::vector<Bounded>& v) {
    std::vector<double> d = std::vector<double>(v.size());
    for (int i = 0; i < d.size(); i++) {
        d[i] = v[i].Get();
    }
    return d;
}

Maybe<Void> Bounded::Set(double val) {
    double precision = 0.00000000000001;
    Maybe<Void> r;
    if (!(Bounded::min <= val + precision && val <= Bounded::max + precision)) {
        r.isError = true;
        r.errMsg = "Tried to set Bounded with value outside [min,max] range";
        return r;
    }
    this->val = val;
    return r;
}