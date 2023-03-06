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

Pimodel::Pimodel(ProblemDescription p, double initialT, double finalT,
                 int timeDiscretization, int kcDiscretization, int order) {
    assert(p.IsOk());
    assert(timeDiscretization >= 1);
    assert(kcDiscretization >= 1);
    assert(initialT <= finalT);
    this->p = p;
    this->nMasses = p.masses.size();
    int nSprings = p.springs.size();
    int nDampers = p.dampers.size();

    this->models = bst::matrix<Poly>(p.NumberOfMasses(), 1);
    this->modelsCoefficients =
        std::vector<std::vector<double>>(p.NumberOfMasses());

    Poly poly;
    Maybe<Void> r;
    for (int massId = 0; massId < p.NumberOfMasses(); massId++) {
        r = poly.Build(nSprings + nDampers + 1, order, massId);
        assert(!r.isError);
        this->models(massId, 0) = poly;
        this->modelsCoefficients[massId] =
            std::vector<double>(poly.nMonomials());
    }

    this->t0 = initialT;
    this->t1 = finalT;
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

int Pimodel::inputSize() {
    return 1 + this->p.springs.size() + this->p.dampers.size();
}

Maybe<std::vector<double>> Pimodel::operator()(std::vector<double>* tkc) {
    Maybe<std::vector<double>> r;
    if (int(tkc->size()) != this->inputSize()) {
        r.errMsg = "Invalid X length";
        r.isError = true;
        return r;
    }
    // if (tkc->at(0) < this->t0 || tkc->at(0) > this->t1) {
    //     r.errMsg = "Invalid t";
    //     r.isError = true;
    //     return r;
    // }
    std::vector<double> positions =
        std::vector<double>(this->p.NumberOfMasses());
    Maybe<double> position;
    for (int massId = 0; massId < int(positions.size()); massId++) {
        this->models(massId, 0).SetX(*tkc);
        position = this->models(massId, 0)(this->modelsCoefficients[massId]);
        assert(!position.isError);
        positions[massId] = position.val;
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
        return this->models(0, 0).nMonomials() * this->nMasses;
    }
    int rv = 0;
    for (int m = 0; m < int(this->models.size1()); m++) {
        rv += this->models(m, 0).nMonomials();
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
    for (int massId = 0; massId < this->nMasses; massId++) {
        for (int mon = 0; mon < this->models(massId, 0).nMonomials(); mon++) {
            this->modelsCoefficients[massId][mon] = parameters->at(parametersI);
            parametersI += 1;
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
    for (int massId = 0; massId < this->nMasses; massId++) {
        for (int mon = 0; mon < this->models(massId, 0).nMonomials(); mon++) {
            target->at(targetI) = this->modelsCoefficients[massId][mon];
            targetI += 1;
        }
    }
    return r;
}

Problem Pimodel::problemFromTkc(std::vector<double>* tkc) {
    std::vector<double> kc = std::vector<double>(tkc->size() - 1);
    std::copy(tkc->begin() + 1, tkc->end(), kc.begin());
    Maybe<Problem> problem = this->p.BuildFromVector(kc);
    assert(!problem.isError);
    return problem.val;
}

std::vector<double> Pimodel::getXModel(std::vector<double>* tkc) {
    // X: State vector. Displacements followed by velocities.
    auto X = std::vector<double>(this->p.NumberOfMasses() * 2);

    Maybe<Void> err;
    Maybe<double> eval;
    for (int massId = 0; massId < this->nMasses; massId++) {
        // Fill displacements
        this->models(massId, 0).SetX(*tkc);
        eval = this->models(massId, 0)(this->modelsCoefficients[massId]);
        assert(!eval.isError);
        X[massId] = eval.val;
    }
    Poly dp_dt;
    for (int massId = 0; massId < this->nMasses; massId++) {
        this->models(massId, 0).SetX(*tkc);
        // Fill velocities
        dp_dt = this->models(massId, 0);
        err = dp_dt.Dxi(0);
        assert(!err.isError);
        eval = dp_dt(this->modelsCoefficients[massId]);
        assert(!eval.isError);
        X[Problem::GetMassVelIndex(this->p.NumberOfMasses(), massId)] =
            eval.val;
    }
    return X;
}

bst::matrix<Polys> Pimodel::getAccelsFromDiffEq(Problem* problem,
                                                std::vector<double>& tkc) {
    // List of displacements and velocities
    bst::matrix<Poly> Disps = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> Vels = bst::matrix<Poly>(nMasses, 1);
    for (int massId = 0; massId < nMasses; massId++) {
        Disps(massId, 0) = this->models(massId, 0);
        Disps(massId, 0).SetX(tkc);

        Vels(massId, 0) = this->models(massId, 0);
        assert(!Vels(massId, 0).Dxi(0).isError);
        Vels(massId, 0).SetX(tkc);
    }

    matrix<Polys> kx = prod(problem->K, Disps);
    matrix<Polys> cxDot = prod(problem->C, Vels);
    matrix<Polys> Accels = prod(problem->MInv, kx + cxDot);
    return Accels;
}

std::vector<double> Pimodel::getInitialX() {
    int nMasses = this->p.masses.size();
    auto initialX = std::vector<double>(nMasses * 2);
    for (auto v : this->p.initialVels) {
        assert(v.massId >= 0 && v.massId < nMasses);
        initialX[Problem::GetMassVelIndex(nMasses, v.massId)] = v.val;
    }
    for (auto d : this->p.initialDisps) {
        assert(d.massId >= 0 && d.massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, d.massId)] = d.val;
    }
    for (auto massId : this->p.fixedMasses) {
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
    if (tkcIndex - 1 < int(this->p.springs.size())) {
        min = this->p.springs[tkcIndex - 1].kMin;
        max = this->p.springs[tkcIndex - 1].kMax;
    } else {
        min = this->p.dampers[tkcIndex - 1 - this->p.springs.size()].cMin;
        max = this->p.dampers[tkcIndex - 1 - this->p.springs.size()].cMax;
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
        min = this->t0;
        max = this->t1;
        discretization = this->timeDiscretization;
    } else {
        if (tkcIndex - 1 < int(this->p.springs.size())) {
            min = this->p.springs[tkcIndex - 1].kMin;
            max = this->p.springs[tkcIndex - 1].kMax;
        } else {
            min = this->p.dampers[tkcIndex - 1 - this->p.springs.size()].cMin;
            max = this->p.dampers[tkcIndex - 1 - this->p.springs.size()].cMax;
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

    tkc[0] = this->t0;
    this->AddInitialConditionsResiduesTkc(&tkc, 1);

    this->AddPhysicsResiduesTkc(&tkc, 0);
}

void Pimodel::AddInitialConditionsResidues() {
    auto initialX = this->getInitialX();

    Poly model;
    for (int massId = 0; massId < this->nMasses; massId++) {
        for (auto tkc : this->initialConditionsResiduesTkc) {
            model = this->models(massId, 0);
            model.SetX(tkc);
            this->initialDispResidues.push_back(
                model +
                (-1) * initialX[Problem::GetMassDispIndex(nMasses, massId)]);
        }
    }
    Poly modelDot;
    for (int massId = 0; massId < this->nMasses; massId++) {
        for (auto tkc : this->initialConditionsResiduesTkc) {
            modelDot = this->models(massId, 0);
            assert(!modelDot.Dxi(0).isError);
            modelDot.SetX(tkc);
            this->initialVelResidues.push_back(
                modelDot +
                (-1) * initialX[Problem::GetMassVelIndex(nMasses, massId)]);
        }
    }
}

void Pimodel::AddPhysicsResidues() {
    Poly modelXDotDot;
    bst::matrix<Polys> AccelsFromDiffEq;
    for (int m = 0; m < nMasses; m++) {
        for (auto tkc : this->physicsResiduesTkc) {
            modelXDotDot = this->models(m, 0);
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

#define LOSS_TERMS                                                            \
    int nInitialConditionsLossTerms = int(this->initialDispResidues.size()) + \
                                      int(this->initialVelResidues.size());   \
    int nPhysicsLossTerms = int(this->physicsResidues.size());                \
    double nTotalLossTerms = nInitialConditionsLossTerms + nPhysicsLossTerms;
double Pimodel::InitialConditionsWeight() {
    LOSS_TERMS
    return nPhysicsLossTerms / nTotalLossTerms;
};
double Pimodel::PhysicsWeight() {
    LOSS_TERMS
    return nInitialConditionsLossTerms / nTotalLossTerms;
};

double Pimodel::Loss() {
    double rv = 0;
    Maybe<double> maybe;

    double initialConditionsWeight = this->InitialConditionsWeight();
    double physicsWeight = this->PhysicsWeight();

    for (int i = 0; i < int(this->initialDispResidues.size()); i++) {
        maybe = this->initialDispResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        rv += pow(maybe.val, 2);
    }
    for (int i = 0; i < int(this->initialVelResidues.size()); i++) {
        maybe = this->initialVelResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        rv += pow(maybe.val, 2);
    }

    rv *= initialConditionsWeight;

    for (int i = 0; i < int(this->physicsResidues.size()); i++) {
        maybe = this->physicsResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        rv += physicsWeight * pow(maybe.val, 2);
    }

    return rv;
}

std::vector<double> Pimodel::LossGradient() {
    std::vector<double> grad = std::vector<double>(this->nParameters());

    Maybe<double> maybe;
    Polys residueD = Polys();  // "derivatives" of the residues

    double initialConditionsWeight = this->InitialConditionsWeight();
    double physicsWeight = this->PhysicsWeight();

    for (int i = 0; i < int(this->initialDispResidues.size()); i++) {
        maybe = this->initialDispResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        residueD +=
            initialConditionsWeight * maybe.val * this->initialDispResidues[i];
    }

    for (int i = 0; i < int(this->initialVelResidues.size()); i++) {
        maybe = this->initialVelResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        residueD +=
            initialConditionsWeight * maybe.val * this->initialVelResidues[i];
    }

    for (int i = 0; i < int(this->physicsResidues.size()); i++) {
        maybe = this->physicsResidues[i](this->modelsCoefficients);
        assert(!maybe.isError);
        residueD += physicsWeight * maybe.val * this->physicsResidues[i];
    }

    std::map<std::tuple<int, int>, double> gradMap = residueD.Da();
    std::tuple<int, int> gradMapKey;
    int gradI = 0;
    for (int massId = 0; massId < int(this->models.size1()); massId++) {
        for (int mon = 0; mon < this->models(massId, 0).nMonomials(); mon++) {
            gradMapKey = {massId, mon};
            if (gradMap.find(gradMapKey) != gradMap.end()) {
                grad[gradI] += gradMap[gradMapKey];
            }
            gradI += 1;
        }
    }

    return grad;
}

// Pimodels::Pimodels(ProblemDescription* p, double finalT, int timeBuckets,
//                    int timeDiscretization, int kcDiscretization, int order) {
//     this->timeBuckets = std::vector<double>(timeBuckets + 1);
//     double timePerTimeBucket = (finalT - 0) / timeBuckets;
//     for (int b = 0; b < timeBuckets; b++) {
//         this->timeBuckets[b] = b * timePerTimeBucket;
//     }
//     this->timeBuckets[timeBuckets] = finalT;

//     double t0 = 0;
//     this->pimodels = std::vector<Pimodel>(timeBuckets);
//     for (int b = 0; b < timeBuckets; b++) {
//         this->pimodels[b] =
//             Pimodel(p, t0, t0 + timePerTimeBucket, timeDiscretization,
//                     kcDiscretization, order);
//         t0 += timePerTimeBucket;
//     }
// }

// Maybe<double> Pimodels::Train(double learningRate, int maxSteps, bool log) {
//     Maybe<double> r;
//     r = this->pimodels[0].Train(learningRate, maxSteps, log);

//     ProblemDescription* p = this->pimodels[0].p;
//     std::vector<double> tkc =
//         std::vector<double>(1 + p->NumberOfSpringsAndDampers());
//     int tkcI = 1;
//     for (int s = 0; s < p->springs.size(); s++) {
//         tkc[tkcI] = (p->springs[s].kMin + p->springs[s].kMax) / 2;
//         tkcI++;
//     }
//     for (int d = 0; d < p->dampers.size(); d++) {
//         tkc[tkcI] = (p->dampers[d].cMin + p->dampers[d].cMax) / 2;
//         tkcI++;
//     }

//     double t;
//     int nMasses = this->pimodels[0].nMasses;
//     Pimodel* previousPiModel;
//     double modelDisp;
//     double modelVel;
//     Maybe<std::vector<double>> X;
//     for (int b = 1; b < this->pimodels.size() - 1; b++) {
//         t = this->timeBuckets[b];
//         tkc[0] = t;
//         this->pimodels[b].p->initialDisps.clear();
//         this->pimodels[b].p->initialVels.clear();
//         previousPiModel = &(this->pimodels[b - 1]);
//         X = (*previousPiModel)(&tkc);
//         assert(!X.isError);
//         for (int m = 0; m < nMasses; m++) {
//             this->pimodels[b].p->AddInitialDisp(m, X.val[m]);
//         }
//         this->pimodels[b].Train(learningRate, maxSteps, log);
//     }
//     return r;
// };

// int Pimodels::getTimeBucket(double t) const {
//     if (t == 0) {
//         return 0;
//     }
//     if (t == this->timeBuckets[this->timeBuckets.size() - 1]) {
//         return this->timeBuckets.size() - 2;
//     }
//     assert(t >= 0);
//     assert(t <= this->timeBuckets[this->timeBuckets.size() - 1]);
//     // Last position in which t could be inserted without changing the order
//     auto bound =
//         std::upper_bound(this->timeBuckets.begin(), this->timeBuckets.end(),
//         t);
//     return (bound - this->timeBuckets.begin()) - 1;
// }

// Pimodel* Pimodels::GetPimodel(double t) {
//     int b = this->getTimeBucket(t);
//     return &(this->pimodels[b]);
// }