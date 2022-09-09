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

Maybe<Void> Bounded::Set(double val) {
    Maybe<Void> r;
    if (!(Bounded::min <= val && val <= Bounded::max)) {
        r.isError = true;
        r.errMsg = "Tried to set Bounded with value outside [min,max] range";
        return r;
    }
    this->val = val;
    return r;
}