#include "pimodel.h"

#include <algorithm>
#include <boost/numeric/ublas/operations.hpp>
#include <cassert>
#include <cmath>
#include <vector>

#include "bounded.h"
#include "maybe.h"
#include "problem.h"
#include "problem_description.h"
#include "utils.h"

namespace bst = boost::numeric::ublas;

Pimodel::Pimodel(ProblemDescription* p, double finalT, int nTimeBuckets,
                 int timeDiscretization, int kcDiscretization, int order) {
    assert(p->IsOk());
    assert(timeDiscretization >= 1);
    assert(kcDiscretization >= 1);
    assert(nTimeBuckets >= 1);
    this->p = p;
    int nSprings = p->springs.size();
    int nDampers = p->dampers.size();

    this->timeBuckets = std::vector<double>(nTimeBuckets + 1);
    double timePerTimeBucket = (finalT - 0) / nTimeBuckets;
    for (int b = 0; b < nTimeBuckets; b++) {
        this->timeBuckets[b] = b * timePerTimeBucket;
    }
    this->timeBuckets[nTimeBuckets] = finalT;

    this->models = std::vector<bst::matrix<Poly>>(nTimeBuckets);
    this->modelsCoefficients =
        std::vector<std::vector<std::vector<double>>>(nTimeBuckets);
    for (int b = 0; b < nTimeBuckets; b++) {
        this->models[b] = bst::matrix<Poly>(p->NumberOfMasses(), 1);
        this->modelsCoefficients[b] =
            std::vector<std::vector<double>>(p->NumberOfMasses());
    }

    Poly poly;
    Maybe<Void> r;
    for (int b = 0; b < nTimeBuckets; b++) {
        for (int i = 0; i < p->NumberOfMasses(); i++) {
            r = poly.Build(nSprings + nDampers + 1, order, i);
            assert(!r.isError);
            this->models[b](i, 0) = poly;
            this->modelsCoefficients[b][i] =
                std::vector<double>(poly.nMonomials());
        }
    }

    this->finalT = finalT;
    this->timeDiscretization = timeDiscretization;
    this->kcDiscretization = kcDiscretization;

    this->initialConditionsResiduesTkc = std::vector<std::vector<double>>();
    this->physicsResiduesTkc = std::vector<std::vector<double>>();
    this->AddResiduesTkc();

    this->initialDispResidues = std::vector<Polys>();
    this->initialVelResidues = std::vector<Polys>();
    AddInitialConditionsResidues();
    this->physicsResidues = std::vector<Polys>();
    AddPhysicsResidues();
}

int Pimodel::timeBucket(double t) {
    if (t == 0) {
        return 0;
    }
    if (t == this->timeBuckets[this->timeBuckets.size() - 1]) {
        return this->timeBuckets.size() - 2;
    }
    assert(t >= 0);
    assert(t <= this->timeBuckets[this->timeBuckets.size() - 1]);
    // Last position in which t could be inserted without changing the order
    auto bound =
        std::upper_bound(this->timeBuckets.begin(), this->timeBuckets.end(), t);
    return (bound - this->timeBuckets.begin()) - 1;
}
int Pimodel::timeBucket(std::vector<double>* tkc) {
    assert(tkc->size() > 0);
    return this->timeBucket(tkc->at(0));
}

int Pimodel::inputSize() {
    return 1 + this->p->springs.size() + this->p->dampers.size();
}

Maybe<std::vector<double>> Pimodel::operator()(std::vector<double>* tkc) {
    Maybe<std::vector<double>> r;
    if (int(tkc->size()) != this->inputSize()) {
        r.errMsg = "Invalid X length";
        r.isError = true;
        return r;
    }
    if (tkc->at(0) < 0 || tkc->at(0) > this->finalT) {
        r.errMsg = "Invalid t";
        r.isError = true;
        return r;
    }
    std::vector<double> positions =
        std::vector<double>(this->p->NumberOfMasses());
    Maybe<double> position;
    for (int i = 0; i < int(positions.size()); i++) {
        this->models[this->timeBucket(tkc)](i, 0).SetX(*tkc);
        auto x = this->models[this->timeBucket(tkc)](i, 0);
        position = this->models[this->timeBucket(tkc)](
            i, 0)(this->modelsCoefficients[this->timeBucket(tkc)][i]);
        assert(!position.isError);
        positions[i] = position.val;
    }
    r.val = positions;
    return r;
}

int Pimodel::nParameters() {
    // Optimization: assume all polynomials have the same order (hence, same
    // numer of monomials), so that this function becomes O(1), and not O(n).
    // This is currently the case, but we might want to change this in the
    // future.
    bool optimize = true;
    if (optimize) {
        return this->models.size() *
               this->models[this->timeBucket(0.0)](0, 0).nMonomials() *
               this->p->NumberOfMasses();
    }
    int rv = 0;
    for (int b = 0; b < this->models.size(); b++) {
        for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
            rv += this->models[b](pIndex, 0).nMonomials();
        }
    }
    return rv;
}

Maybe<Void> Pimodel::SetParameters(std::vector<double>* parameters) {
    Maybe<Void> r;
    if (int(parameters->size()) != this->nParameters()) {
        r.isError = true;
        r.errMsg = "Invalid parameters length";
        return r;
    }

    int parametersI = 0;
    for (int b = 0; b < this->models.size(); b++) {
        for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
            for (int m = 0; m < this->models[b](pIndex, 0).nMonomials(); m++) {
                this->modelsCoefficients[b][pIndex][m] =
                    parameters->at(parametersI);
                parametersI += 1;
            }
        }
    }
    return r;
}

Maybe<Void> Pimodel::GetParameters(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nParameters()) {
        r.isError = true;
        r.errMsg = "Invalid target length";
        return r;
    }

    int targetI = 0;
    for (int b = 0; b < this->models.size(); b++) {
        for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
            for (int m = 0; m < this->models[b](pIndex, 0).nMonomials(); m++) {
                target->at(targetI) = this->modelsCoefficients[b][pIndex][m];
                targetI += 1;
            }
        }
    }
    return r;
}

Problem Pimodel::problemFromTkc(std::vector<double>* tkc) {
    std::vector<double> kc = std::vector<double>(tkc->size() - 1);
    std::copy(tkc->begin() + 1, tkc->end(), kc.begin());
    Maybe<Problem> problem = this->p->BuildFromVector(kc);
    assert(!problem.isError);
    return problem.val;
}

std::vector<double> Pimodel::getXModel(std::vector<double>* tkc) {
    // X: State vector. Displacements followed by velocities.
    auto X = std::vector<double>(this->p->NumberOfMasses() * 2);

    int b = this->timeBucket(tkc);
    Maybe<Void> err;
    Maybe<double> eval;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        // Fill displacements
        this->models[b](i, 0).SetX(*tkc);
        eval = this->models[b](i, 0)(this->modelsCoefficients[b][i]);
        assert(!eval.isError);
        X[i] = eval.val;
    }
    Poly dp_dt;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        this->models[b](i, 0).SetX(*tkc);
        // Fill velocities
        dp_dt = this->models[b](i, 0);
        err = dp_dt.Dxi(0);
        assert(!err.isError);
        eval = dp_dt(this->modelsCoefficients[b][i]);
        assert(!eval.isError);
        X[Problem::GetMassVelIndex(this->p->NumberOfMasses(), i)] = eval.val;
    }
    return X;
}

bst::matrix<Polys> Pimodel::getAccelsFromDiffEq(Problem* problem,
                                                std::vector<double>& tkc) {
    // See SetXDot at problem.h for reference
    int nMasses = problem->masses.size();

    int b = this->timeBucket(&tkc);

    // List of displacements and velocities
    bst::matrix<Poly> Disps = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> Vels = bst::matrix<Poly>(nMasses, 1);
    for (int i = 0; i < nMasses; i++) {
        Disps(i, 0) = this->models[b](i, 0);
        Disps(i, 0).SetX(tkc);

        Vels(i, 0) = this->models[b](i, 0);
        assert(!Vels(i, 0).Dxi(0).isError);
        Vels(i, 0).SetX(tkc);
    }

    matrix<Polys> kx = prod(problem->K, Disps);
    matrix<Polys> cxDot = prod(problem->C, Vels);
    matrix<Polys> Accels = prod(problem->MInv, kx + cxDot);
    return Accels;
}

std::vector<double> Pimodel::getInitialX() {
    int nMasses = this->p->masses.size();
    auto initialX = std::vector<double>(nMasses * 2);
    for (auto v : this->p->initialVels) {
        assert(v.massId >= 0 && v.massId < nMasses);
        initialX[Problem::GetMassVelIndex(nMasses, v.massId)] = v.val;
    }
    for (auto d : this->p->initialDisps) {
        assert(d.massId >= 0 && d.massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, d.massId)] = d.val;
    }
    for (auto massId : this->p->fixedMasses) {
        assert(massId >= 0 && massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, massId)] = 0;
        initialX[Problem::GetMassVelIndex(nMasses, massId)] = 0;
    }
    return initialX;
}

void Pimodel::AddInitialConditionsResiduesTkc(std::vector<double>* tkc,
                                              int tkcIndex) {
    if (tkcIndex == int(tkc->size())) {
        this->initialConditionsResiduesTkc.push_back(*tkc);
        return;
    }

    double min;
    double max;
    if (tkcIndex - 1 < int(this->p->springs.size())) {
        min = this->p->springs[tkcIndex - 1].kMin;
        max = this->p->springs[tkcIndex - 1].kMax;
    } else {
        min = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMin;
        max = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMax;
    }
    for (int i = 0; i <= this->kcDiscretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / this->kcDiscretization * i;
        this->AddInitialConditionsResiduesTkc(tkc, tkcIndex + 1);
    }
}
void Pimodel::AddPhysicsResiduesTkc(std::vector<double>* tkc, int tkcIndex) {
    if (tkcIndex == int(tkc->size())) {
        this->physicsResiduesTkc.push_back(*tkc);
        return;
    }

    double min;
    double max;
    int discretization;
    if (tkcIndex == 0) {
        min = 0;
        max = this->finalT;
        discretization = this->timeDiscretization;
    } else {
        if (tkcIndex - 1 < int(this->p->springs.size())) {
            min = this->p->springs[tkcIndex - 1].kMin;
            max = this->p->springs[tkcIndex - 1].kMax;
        } else {
            min = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMin;
            max = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMax;
        }
        discretization = this->kcDiscretization;
    }
    for (int i = 0; i <= discretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / discretization * i;
        this->AddPhysicsResiduesTkc(tkc, tkcIndex + 1);
    }
}
void Pimodel::AddResiduesTkc() {
    std::vector<double> tkc = std::vector<double>(this->inputSize());

    tkc[0] = 0;  // t = 0
    this->AddInitialConditionsResiduesTkc(&tkc, 1);

    this->AddPhysicsResiduesTkc(&tkc, 0);
}

void Pimodel::AddInitialConditionsResidues() {
    auto initialX = this->getInitialX();

    Poly model;
    int b;
    int nMasses = this->models[0].size1();
    for (int m = 0; m < nMasses; m++) {
        for (auto tkc : this->initialConditionsResiduesTkc) {
            b = this->timeBucket(&tkc);
            model = this->models[b](m, 0);
            model.SetX(tkc);
            this->initialDispResidues.push_back(
                model + (-1) * initialX[Problem::GetMassDispIndex(nMasses, m)]);
        }
    }
    Poly modelDot;
    for (int m = 0; m < nMasses; m++) {
        for (auto tkc : this->initialConditionsResiduesTkc) {
            b = this->timeBucket(&tkc);
            modelDot = this->models[b](m, 0);
            assert(!modelDot.Dxi(0).isError);
            modelDot.SetX(tkc);
            this->initialVelResidues.push_back(
                modelDot +
                (-1) * initialX[Problem::GetMassVelIndex(nMasses, m)]);
        }
    }
}

void Pimodel::AddPhysicsResidues() {
    int nMasses = this->models[0].size1();

    int b;
    Poly modelXDotDot;
    bst::matrix<Polys> AccelsFromDiffEq;
    for (int m = 0; m < nMasses; m++) {
        for (auto tkc : this->physicsResiduesTkc) {
            b = this->timeBucket(&tkc);
            modelXDotDot = this->models[b](m, 0);
            assert(!modelXDotDot.Dxi(0).isError);
            assert(!modelXDotDot.Dxi(0).isError);
            modelXDotDot.SetX(tkc);

            Problem problem = this->problemFromTkc(&tkc);
            if (problem.massIsFixed(m)) {
                this->physicsResidues.push_back(Polys(modelXDotDot));
            } else {
                AccelsFromDiffEq = getAccelsFromDiffEq(&problem, tkc);
                this->physicsResidues.push_back(
                    Polys(modelXDotDot) + (-1.0) * AccelsFromDiffEq(m, 0));
            }
        }
    }
}

double Pimodel::Loss() {
    double rv = 0;
    Maybe<double> maybe;

    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->initialDispResidues.size()); i++) {
            maybe = this->initialDispResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            rv += pow(maybe.val, 2);
        }
    }
    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->initialVelResidues.size()); i++) {
            maybe = this->initialVelResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            rv += pow(maybe.val, 2);
        }
    }
    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->physicsResidues.size()); i++) {
            maybe = this->physicsResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            rv += pow(maybe.val, 2);
        }
    }

    return rv;
}

std::vector<double> Pimodel::LossGradient() {
    std::vector<double> grad = std::vector<double>(this->nParameters());

    Maybe<double> maybe;
    Polys residueD = Polys();  // "derivatives" of the residues

    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->initialDispResidues.size()); i++) {
            maybe = this->initialDispResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            residueD += maybe.val * this->initialDispResidues[i];
        }
    }
    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->initialVelResidues.size()); i++) {
            maybe = this->initialVelResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            residueD += maybe.val * this->initialVelResidues[i];
        }
    }
    for (int b = 0; b < this->models.size(); b++) {
        for (int i = 0; i < int(this->physicsResidues.size()); i++) {
            maybe = this->physicsResidues[i](this->modelsCoefficients[b]);
            assert(!maybe.isError);
            residueD += maybe.val * this->physicsResidues[i];
        }
    }

    std::map<std::tuple<int, int>, double> gradMap = residueD.Da();

    std::tuple<int, int> gradMapKey;
    int gradI = 0;
    for (int b = 0; b < this->models.size(); b++) {
        for (int pId = 0; pId < int(this->models[b].size1()); pId++) {
            for (int i = 0; i < this->models[b](pId, 0).nMonomials(); i++) {
                gradMapKey = {pId, i};
                if (gradMap.find(gradMapKey) != gradMap.end()) {
                    grad[gradI] += gradMap[gradMapKey];
                }
                gradI += 1;
            }
        }
    }

    return grad;
}