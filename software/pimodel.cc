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
    assert(timeDiscretization >= 0);
    assert(kcDiscretization >= 0);
    assert(initialT <= finalT);
    this->p = p;
    this->nMasses = p.masses.size();
    this->order = order;
    this->t0 = initialT;
    this->t1 = finalT;
    this->timeDiscretization = timeDiscretization;
    this->kcDiscretization = kcDiscretization;
    this->models = bst::matrix<Poly>(p.NumberOfMasses(), 1);
    this->modelsD = std::vector<Poly>(p.NumberOfMasses(), 1);
    this->modelsDD = std::vector<Poly>(p.NumberOfMasses(), 1);
    this->modelsCoefficients =
        std::vector<std::vector<double>>(p.NumberOfMasses());

    int nSprings = p.springs.size();
    int nDampers = p.dampers.size();
    Poly poly;
    Maybe<Void> r;
    for (int massId = 0; massId < p.NumberOfMasses(); massId++) {
        r = poly.Build(nSprings + nDampers + 1, this->order, massId);
        assert(!r.isError);
        this->models(massId, 0) = poly;
        this->modelsCoefficients[massId] =
            std::vector<double>(poly.nMonomials());
        this->modelsD[massId] = poly;
        assert(!this->modelsD[massId].Dxi(0).isError);
        this->modelsDD[massId] = this->modelsD[massId];
        assert(!this->modelsDD[massId].Dxi(0).isError);
    }
}

void Pimodel::SetResidues() {
    this->initialConditionsResiduesTkc = std::vector<std::vector<Bounded>>();
    this->physicsResiduesTkc = std::vector<std::vector<Bounded>>();
    this->AddResiduesTkc();

    this->initialDispResidues = std::vector<Polys>();
    this->initialVelResidues = std::vector<Polys>();
    AddInitialConditionsResidues();
    this->physicsResidues = std::vector<Polys>();
    AddPhysicsResidues();

    this->residueIsCached = std::vector<bool>(this->nResidues());
    this->residueCache = std::vector<double>(this->nResidues());
    for (int i = 0; i < this->nResidues(); i++) {
        this->residueIsCached[i] = false;
    }

    this->initialDispResiduesDa =
        std::vector<std::map<std::tuple<int, int>, double>>(
            this->initialDispResidues.size());
    this->initialVelResiduesDa =
        std::vector<std::map<std::tuple<int, int>, double>>(
            this->initialVelResidues.size());
    this->physicsResiduesDa =
        std::vector<std::map<std::tuple<int, int>, double>>(
            this->physicsResidues.size());
    this->setResiduesDa();
}

int Pimodel::inputSize() {
    return 1 + this->p.springs.size() + this->p.dampers.size();
}

std::vector<double> Pimodel::normalizeTkc(std::vector<double>* TKC) {
    assert(int(TKC->size()) == this->inputSize());

    Maybe<double> err;
    std::vector<double> tkcNormalized = std::vector<double>(TKC->size());
    err = NormalizeToDouble(TKC->at(0), this->t0, this->t1);
    assert(!err.isError);
    tkcNormalized[0] = err.val;

    int i = 1;
    for (int k = 0; k < int(this->p.springs.size()); k++) {
        err = NormalizeToDouble(TKC->at(i), this->p.springs[k].kMin,
                                this->p.springs[k].kMax);
        assert(!err.isError);
        tkcNormalized[i] = err.val;
        i++;
    }
    for (int c = 0; c < int(this->p.dampers.size()); c++) {
        err = NormalizeToDouble(TKC->at(i), this->p.dampers[c].cMin,
                                this->p.dampers[c].cMax);
        assert(!err.isError);
        tkcNormalized[i] = err.val;
        i++;
    }

    return tkcNormalized;
}

Maybe<std::vector<double>> Pimodel::operator()(std::vector<double>* TKC) {
    Maybe<std::vector<double>> r;
    if (int(TKC->size()) != this->inputSize()) {
        r.errMsg = "Invalid X length";
        r.isError = true;
        return r;
    }

    std::vector<double> positions =
        std::vector<double>(this->p.NumberOfMasses());
    Maybe<double> position;
    for (int massId = 0; massId < int(positions.size()); massId++) {
        this->models(massId, 0).SetX(this->normalizeTkc(TKC));
        position = this->models(massId, 0)(this->modelsCoefficients[massId]);
        assert(!position.isError);
        positions[massId] = position.val;
    }
    r.val = positions;
    return r;
}

Maybe<std::vector<double>> Pimodel::GetVelocities(std::vector<double>* TKC) {
    Maybe<std::vector<double>> r;
    if (int(TKC->size()) != this->inputSize()) {
        r.errMsg = "Invalid X length";
        r.isError = true;
        return r;
    }
    std::vector<double> vels = std::vector<double>(this->p.NumberOfMasses());
    Maybe<double> vel;
    for (int massId = 0; massId < int(vels.size()); massId++) {
        this->modelsD[massId].SetX(this->normalizeTkc(TKC));
        vel = this->modelsD[massId](this->modelsCoefficients[massId]);
        assert(!vel.isError);
        vels[massId] = vel.val;
    }
    r.val = vels;
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
    for (int i = 0; i < this->nResidues(); i++) {
        this->residueIsCached[i] = false;
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

Problem Pimodel::problemFromTkc(std::vector<Bounded>* tkc) {
    std::vector<Bounded> kc = std::vector<Bounded>(tkc->size() - 1);
    std::copy(tkc->begin() + 1, tkc->end(), kc.begin());
    Maybe<Problem> problem = this->p.BuildFromDNA(kc);
    assert(!problem.isError);
    return problem.val;
}

std::vector<double> Pimodel::getXModel(std::vector<Bounded>* tkc) {
    // X: State vector. Displacements followed by velocities.
    auto X = std::vector<double>(this->p.NumberOfMasses() * 2);

    Maybe<Void> err;
    Maybe<double> eval;
    for (int massId = 0; massId < this->nMasses; massId++) {
        // Fill displacements
        this->models(massId, 0).SetX(Bounded::Get(*tkc));
        eval = this->models(massId, 0)(this->modelsCoefficients[massId]);
        assert(!eval.isError);
        X[massId] = eval.val;

        // Fill velocities
        this->modelsD[massId].SetX(Bounded::Get(*tkc));
        eval = this->modelsD[massId](this->modelsCoefficients[massId]);
        assert(!eval.isError);
        X[Problem::GetMassVelIndex(this->p.NumberOfMasses(), massId)] =
            eval.val;
    }
    return X;
}

bst::matrix<Polys> Pimodel::getAccelsFromDiffEq(Problem* problem,
                                                std::vector<Bounded>& tkc) {
    // List of displacements and velocities
    bst::matrix<Poly> Disps = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> Vels = bst::matrix<Poly>(nMasses, 1);
    for (int massId = 0; massId < nMasses; massId++) {
        Disps(massId, 0) = this->models(massId, 0);
        Disps(massId, 0).SetX(Bounded::Get(tkc));

        Vels(massId, 0) = this->modelsD[massId];
        Vels(massId, 0).SetX(Bounded::Get(tkc));
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
        assert((v.massId >= 0 && v.massId < nMasses) || v.massId == -1);
        if (v.massId != -1) {
            initialX[Problem::GetMassVelIndex(nMasses, v.massId)] = v.val;
        } else {
            for (int i = 0; i < this->nMasses; i++) {
                initialX[Problem::GetMassVelIndex(nMasses, i)] = v.val;
            }
        }
    }
    for (auto d : this->p.initialDisps) {
        assert((d.massId >= 0 && d.massId < nMasses) || d.massId == -1);
        if (d.massId != -1) {
            initialX[Problem::GetMassDispIndex(nMasses, d.massId)] = d.val;
        } else {
            for (int i = 0; i < this->nMasses; i++) {
                initialX[Problem::GetMassDispIndex(nMasses, i)] = d.val;
            }
        }
    }
    for (auto massId : this->p.fixedMasses) {
        assert(massId >= 0 && massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, massId)] = 0;
        initialX[Problem::GetMassVelIndex(nMasses, massId)] = 0;
    }
    return initialX;
}

void Pimodel::AddInitialConditionsResiduesTkc(std::vector<Bounded>* tkc,
                                              int tkcIndex) {
    if (tkcIndex == int(tkc->size())) {
        this->initialConditionsResiduesTkc.push_back(*tkc);
        return;
    }

    double min = 0.0;
    double max = 1.0;
    if (this->kcDiscretization == 0) {
        tkc->at(tkcIndex) = Bounded();
        assert(!tkc->at(tkcIndex).Set(0.5).isError);
        this->AddInitialConditionsResiduesTkc(tkc, tkcIndex + 1);
        return;
    }
    for (int i = 0; i <= this->kcDiscretization; i++) {
        tkc->at(tkcIndex) = Bounded();
        assert(!tkc->at(tkcIndex)
                    .Set((max - min) / this->kcDiscretization * i)
                    .isError);
        this->AddInitialConditionsResiduesTkc(tkc, tkcIndex + 1);
    }
}
void Pimodel::AddPhysicsResiduesTkc(std::vector<Bounded>* tkc, int tkcIndex) {
    if (tkcIndex == int(tkc->size())) {
        this->physicsResiduesTkc.push_back(*tkc);
        return;
    }

    double min = 0.0;
    double max = 1.0;
    int discretization;
    if (tkcIndex == 0) {
        discretization = this->timeDiscretization;
    } else {
        discretization = this->kcDiscretization;
    }
    if (discretization == 0) {
        tkc->at(tkcIndex) = Bounded();
        assert(!tkc->at(tkcIndex).Set(0.5).isError);
        this->AddPhysicsResiduesTkc(tkc, tkcIndex + 1);
        return;
    }
    for (int i = 0; i <= discretization; i++) {
        tkc->at(tkcIndex) = Bounded();
        assert(
            !tkc->at(tkcIndex).Set((max - min) / discretization * i).isError);
        this->AddPhysicsResiduesTkc(tkc, tkcIndex + 1);
    }
}
void Pimodel::AddResiduesTkc() {
    std::vector<Bounded> tkc = std::vector<Bounded>(this->inputSize());

    tkc[0].Set(0.0);
    this->AddInitialConditionsResiduesTkc(&tkc, 1);

    this->AddPhysicsResiduesTkc(&tkc, 0);
}

void Pimodel::AddInitialConditionsResidues() {
    auto initialX = this->getInitialX();

    Poly model;
    for (int massId = 0; massId < this->nMasses; massId++) {
        for (auto tkc : this->initialConditionsResiduesTkc) {
            model = this->models(massId, 0);
            model.SetX(Bounded::Get(tkc));
            this->initialDispResidues.push_back(
                model +
                (-1) * initialX[Problem::GetMassDispIndex(nMasses, massId)]);

            model = this->modelsD[massId];
            model.SetX(Bounded::Get(tkc));
            this->initialVelResidues.push_back(
                model +
                (-1) * initialX[Problem::GetMassVelIndex(nMasses, massId)]);
        }
    }
}

void Pimodel::AddPhysicsResidues() {
    Poly modelXDotDot;
    bst::matrix<Polys> AccelsFromDiffEq;
    for (int m = 0; m < nMasses; m++) {
        for (auto tkc : this->physicsResiduesTkc) {
            modelXDotDot = this->modelsDD[m];
            modelXDotDot.SetX(Bounded::Get(tkc));

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
    if (nPhysicsLossTerms == 0) {
        return 1;
    }
    return nPhysicsLossTerms / nTotalLossTerms;
};
double Pimodel::PhysicsWeight() {
    LOSS_TERMS
    if (nInitialConditionsLossTerms == 0) {
        return 1;
    }
    return nInitialConditionsLossTerms / nTotalLossTerms;
};

double Pimodel::residueWeight(int i) {
    if (i < this->initialDispResidues.size()) {
        return this->InitialConditionsWeight();
    }
    i -= this->initialDispResidues.size();
    if (i < this->initialVelResidues.size()) {
        return this->InitialConditionsWeight();
    }
    return this->PhysicsWeight();
}

Polys* Pimodel::residueById(int i) {
    if (i < this->initialDispResidues.size()) {
        return &(this->initialDispResidues[i]);
    }
    i -= this->initialDispResidues.size();
    if (i < this->initialVelResidues.size()) {
        return &(this->initialVelResidues[i]);
    }
    i -= this->initialVelResidues.size();
    return &(this->physicsResidues[i]);
}

std::map<std::tuple<int, int>, double>* Pimodel::residueDaById(int i) {
    if (i < this->initialDispResidues.size()) {
        return &(this->initialDispResiduesDa[i]);
    }
    i -= this->initialDispResidues.size();
    if (i < this->initialVelResidues.size()) {
        return &(this->initialVelResiduesDa[i]);
    }
    i -= this->initialVelResidues.size();
    return &(this->physicsResiduesDa[i]);
}

void Pimodel::setResidueCache(int i) {
    Maybe<double> maybe = (*this->residueById(i))(this->modelsCoefficients);
    assert(!maybe.isError);
    this->residueCache[i] = maybe.val;
    this->residueIsCached[i] = true;
}

int Pimodel::nResidues() {
    return this->initialDispResidues.size() + this->initialVelResidues.size() +
           this->physicsResidues.size();
}

double Pimodel::Residue(int i) {
    if (!this->residueIsCached[i]) {
        this->setResidueCache(i);
    }
    return this->residueWeight(i) * pow(this->residueCache[i], 2);
}

void Pimodel::setResiduesDa() {
    for (int i = 0; i < int(this->initialDispResidues.size()); i++) {
        this->initialDispResiduesDa[i] = this->initialDispResidues[i].Da();
    }
    for (int i = 0; i < int(this->initialVelResidues.size()); i++) {
        this->initialVelResiduesDa[i] = this->initialVelResidues[i].Da();
    }

    for (int i = 0; i < int(this->physicsResidues.size()); i++) {
        this->physicsResiduesDa[i] = this->physicsResidues[i].Da();
    }
}

std::vector<double> Pimodel::ResidueGradient(int i) {
    std::vector<double> grad = std::vector<double>(this->nParameters());

    if (!this->residueIsCached[i]) {
        this->setResidueCache(i);
    }
    double weight = this->residueWeight(i);
    double residueVal = this->residueCache[i];
    std::map<std::tuple<int, int>, double>* gradMap = this->residueDaById(i);

    std::tuple<int, int> gradMapKey;
    int gradI = 0;
    for (int massId = 0; massId < int(this->models.size1()); massId++) {
        for (int mon = 0; mon < this->models(massId, 0).nMonomials(); mon++) {
            gradMapKey = {massId, mon};
            if (gradMap->find(gradMapKey) != gradMap->end()) {
                grad[gradI] += weight * residueVal * ((*gradMap)[gradMapKey]);
            }
            gradI += 1;
        }
    }

    return grad;
}

Pimodels::Pimodels(ProblemDescription p, double finalT, int nModels,
                   int timeDiscretization, int kcDiscretization, int order) {
    assert(finalT > 0);
    assert(nModels >= 1);
    assert(timeDiscretization >= 0);
    assert(kcDiscretization >= 0);
    assert(order >= 0);

    this->timeBuckets = std::vector<double>(nModels + 1);
    double timePerTimeBucket = (finalT - 0) / nModels;
    for (int b = 0; b < nModels; b++) {
        this->timeBuckets[b] = b * timePerTimeBucket;
    }
    this->timeBuckets[nModels] = finalT;

    double t0 = 0;
    this->pimodels = std::vector<Pimodel>(nModels);
    for (int b = 0; b < nModels; b++) {
        this->pimodels[b] =
            Pimodel(p, t0, t0 + timePerTimeBucket, timeDiscretization,
                    kcDiscretization, order);
        t0 += timePerTimeBucket;
    }
}

std::vector<double> Pimodels::continuityTkc() const {
    const ProblemDescription& p = this->pimodels[0].p;
    std::vector<double> tkc =
        std::vector<double>(1 + p.NumberOfSpringsAndDampers());

    tkc[0] = 0;  // t = 0

    int tkcI = 1;
    for (int s = 0; s < int(p.springs.size()); s++) {
        tkc[tkcI] = (p.springs[s].kMin + p.springs[s].kMax) / 2;
        tkcI++;
    }
    for (int d = 0; d < int(p.dampers.size()); d++) {
        tkc[tkcI] = (p.dampers[d].cMin + p.dampers[d].cMax) / 2;
        tkcI++;
    }
    return tkc;
}

void Pimodels::setContinuity(int timeBucket, std::vector<double>& TKC) {
    assert(timeBucket >= 1 && timeBucket < int(this->pimodels.size()));
    int& nMasses = this->pimodels[0].nMasses;

    Pimodel& thisPiModel = this->pimodels[timeBucket];
    Pimodel& previousPiModel = this->pimodels[timeBucket - 1];

    // Set t at TKC
    TKC[0] = this->timeBuckets[timeBucket];

    // Clear initial conditions
    thisPiModel.p.initialDisps.clear();
    thisPiModel.p.initialVels.clear();

    // Calculate displacements and velocities of previous model
    auto disps = previousPiModel(&TKC);
    assert(!disps.isError);
    auto vels = previousPiModel.GetVelocities(&TKC);
    assert(!vels.isError);

    // C0 continuity
    for (int massId = 0; massId < nMasses; massId++) {
        thisPiModel.p.AddInitialDisp(massId, disps.val[massId]);
    }
    // C1 continuity
    for (int massId = 0; massId < nMasses; massId++) {
        thisPiModel.p.AddInitialVel(massId, vels.val[massId]);
    }
}

Maybe<double> Pimodels::Train(double learningRate, int maxSteps, bool log) {
    double learningRate0 = learningRate;
    Maybe<double> r;
    this->pimodels[0].SetResidues();

    if (log) {
        std::cout << "## Complexity Stats ##" << std::endl;

        std::cout << "Number of polynomial models (1 for each mass): "
                  << this->pimodels[0].models.size1() << std::endl;
        std::cout << "Number of monomials per polynomial model: "
                  << this->pimodels[0].models(0, 0).nMonomials() << std::endl;
        std::cout << "Number of initial conditions residues: "
                  << this->pimodels[0].initialDispResidues.size() +
                         this->pimodels[0].initialVelResidues.size()
                  << std::endl;
        std::cout << "Number of physics residues: "
                  << this->pimodels[0].physicsResidues.size() << std::endl;
        int n = 0;
        for (int i = 0; i < this->pimodels[0].initialDispResidues.size(); i++) {
            n += this->pimodels[0].initialDispResidues[i].nMonomials();
        }
        std::cout
            << "Avg number of monomials in initial displacement residues: "
            << n / this->pimodels[0].initialDispResidues.size() << std::endl;
        n = 0;
        for (int i = 0; i < this->pimodels[0].initialVelResidues.size(); i++) {
            n += this->pimodels[0].initialVelResidues[i].nMonomials();
        }
        std::cout << "Avg number of monomials in initial velocity residues: "
                  << n / this->pimodels[0].initialVelResidues.size()
                  << std::endl;
        n = 0;
        for (int i = 0; i < this->pimodels[0].physicsResidues.size(); i++) {
            n += this->pimodels[0].physicsResidues[i].nMonomials();
        }
        std::cout << "Avg number of monomials in physics residues: "
                  << n / this->pimodels[0].physicsResidues.size() << std::endl;

        std::cout << "Note: the cost of calculating each residue is "
                     "proportional to the # of monomials"
                  << std::endl;
        std::cout << "Note: the cost of calculating the gradient of each "
                     "residue is proportional to (nMasses*nMonomials)"
                  << std::endl;
        // calculating the derivatives themselves is O(1), because the
        // derivative of each residue is constant and we cache it. However,
        // every time we calculate the gradient we calculate the derivative with
        // respect to each and all parameters. The total number of parameters is
        // nMonomialsPerModel*nModels

        std::cout << "######################" << std::endl;
        std::cout << std::endl;
    }

    r = this->pimodels[0].Train(learningRate, maxSteps, log);

    std::vector<double> tkc = this->continuityTkc();

    for (int b = 1; b < int(this->pimodels.size()); b++) {
        learningRate = learningRate0;
        this->setContinuity(b, tkc);
        this->pimodels[b].SetResidues();
        this->pimodels[b].Train(learningRate, maxSteps, log);
    }
    return r;
};

int Pimodels::getTimeBucket(double t) const {
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

Maybe<std::vector<double>> Pimodels::operator()(std::vector<double>* TKC) {
    assert(TKC->size() > 0);
    int b = this->getTimeBucket(TKC->at(0));
    return this->pimodels[b](TKC);
}

Maybe<std::vector<double>> Pimodels::GetVelocities(std::vector<double>* TKC) {
    assert(TKC->size() > 0);
    int b = this->getTimeBucket(TKC->at(0));
    return this->pimodels[b].GetVelocities(TKC);
}
