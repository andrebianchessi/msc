#include "problem_description.h"

#include <memory>

#include "maybe.h"
#include "utils.h"

void ProblemDescription::AddMass(double m, double px, double py) {
    this->masses.push_back({m, px, py});
};
void ProblemDescription::AddSpring(int m0, int m1, double kMin, double kMax) {
    this->springs.push_back({m0, m1, kMin, kMax});
};
void ProblemDescription::AddDamper(int m0, int m1, double cMin, double cMax) {
    this->dampers.push_back({m0, m1, cMin, cMax});
};
void ProblemDescription::AddInitialVel(int massId, double value) {
    this->initialVels.push_back({value, massId});
};
void ProblemDescription::AddInitialVel(double value) {
    this->initialVels.push_back({value, -1});
};
void ProblemDescription::AddInitialDisp(int massId, double value) {
    this->initialDisps.push_back({value, massId});
};
void ProblemDescription::AddInitialDisp(double value) {
    this->initialDisps.push_back({value, -1});
};
void ProblemDescription::SetFixedMass(int massId) {
    this->fixedMasses.push_back(massId);
};

int ProblemDescription::NumberOfSpringsAndDampers() const {
    return this->springs.size() + this->dampers.size();
};

int ProblemDescription::NumberOfMasses() { return this->masses.size(); }

bool ProblemDescription::IsOk() {
    std::vector<Bounded> dna =
        std::vector<Bounded>(this->NumberOfSpringsAndDampers());
    auto p = this->BuildFromDNA(dna);
    if (p.isError) {
        return false;
    }
    return true;
}

Maybe<Problem> ProblemDescription::BuildFromDNA(std::vector<Bounded> dna) {
    Maybe<Problem> r;

    if (int(dna.size()) != this->NumberOfSpringsAndDampers()) {
        r.isError = true;
        r.errMsg = "dna must be the same size of NumberOfSpringsAndDampers()";
        return r;
    }

    Problem p;
    for (auto m : this->masses) {
        auto e = p.AddMass(m.m, m.px, m.py);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    int i = 0;
    for (auto s : this->springs) {
        double k = Unnormalize(dna[i], s.kMin, s.kMax);
        auto e = p.AddSpring(s.m0, s.m1, k);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
        i += 1;
    }
    for (auto d : this->dampers) {
        double c = Unnormalize(dna[i], d.cMin, d.cMax);
        auto e = p.AddDamper(d.m0, d.m1, c);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
        i += 1;
    }
    p.Build();

    for (auto v : this->initialVels) {
        Maybe<Void> e;
        if (v.massId != -1) {
            e = p.SetInitialVel(v.massId, v.val);
        } else {
            e = p.SetInitialVel(v.val);
        }
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    for (auto d : this->initialDisps) {
        Maybe<Void> e;
        if (d.massId != -1) {
            e = p.SetInitialDisp(d.massId, d.val);
        } else {
            e = p.SetInitialDisp(d.val);
        }
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }

    for (auto m : this->fixedMasses) {
        Maybe<Void> e;
        e = p.FixMass(m);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    r.val = p;
    return r;
};

Maybe<Problem> ProblemDescription::BuildFromVector(
    std::vector<double> springsAndDampers) {
    Maybe<Problem> r;

    if (int(springsAndDampers.size()) != this->NumberOfSpringsAndDampers()) {
        r.isError = true;
        r.errMsg = "dna must be the same size of NumberOfSpringsAndDampers()";
        return r;
    }

    Problem p;
    for (auto m : this->masses) {
        auto e = p.AddMass(m.m, m.px, m.py);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    int i = 0;
    for (auto s : this->springs) {
        double k = springsAndDampers[i];
        if (k < s.kMin || k > s.kMax) {
            r.isError = true;
            r.errMsg = "Invalid spring value";
            return r;
        }
        auto e = p.AddSpring(s.m0, s.m1, k);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
        i += 1;
    }
    for (auto d : this->dampers) {
        double c = springsAndDampers[i];
        if (c < d.cMin || c > d.cMax) {
            r.isError = true;
            r.errMsg = "Invalid damper value";
            return r;
        }
        auto e = p.AddDamper(d.m0, d.m1, c);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
        i += 1;
    }
    p.Build();

    for (auto v : this->initialVels) {
        Maybe<Void> e;
        if (v.massId != -1) {
            e = p.SetInitialVel(v.massId, v.val);
        } else {
            e = p.SetInitialVel(v.val);
        }
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    for (auto d : this->initialDisps) {
        Maybe<Void> e;
        if (d.massId != -1) {
            e = p.SetInitialDisp(d.massId, d.val);
        } else {
            e = p.SetInitialDisp(d.val);
        }
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }

    for (auto m : this->fixedMasses) {
        Maybe<Void> e;
        e = p.FixMass(m);
        if (e.isError) {
            r.isError = true;
            r.errMsg = e.isError;
            return r;
        }
    }
    r.val = p;
    return r;
};
