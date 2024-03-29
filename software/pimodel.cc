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
                 int initialConditionTrainingPoints, int physicsTrainingPoints,
                 int order) {
    assert(p.IsOk());
    assert(initialConditionTrainingPoints >= 1);
    assert(physicsTrainingPoints >= 2);
    assert(initialT <= finalT);
    this->p = p;
    this->nMasses = p.masses.size();
    this->order = order;
    this->t0 = initialT;
    this->t1 = finalT;
    this->initialConditionTrainingPoints = initialConditionTrainingPoints;
    this->physicsTrainingPoints = physicsTrainingPoints;
    this->models = bst::matrix<Poly>(p.NumberOfMasses(), 1);
    this->modelsD = std::vector<Poly>(p.NumberOfMasses(), 1);
    this->modelsDD = std::vector<Poly>(p.NumberOfMasses(), 1);
    this->modelsCoefficients =
        std::vector<std::vector<double>>(p.NumberOfMasses());

    int nSprings = p.springs.size();
    int nDampers = p.dampers.size();
    Poly poly;
    Maybe<Void> r;

    // T = T0 + (T1-T0)*t
    // e.g.:
    //      t = 0.0 -> T = T0
    //      t = 1.0 -> T = T1
    //      t = 0.5 -> T = (T1+T0)/2
    this->dTdt = finalT - initialT;
    this->dtdT = 1.0 / (finalT - initialT);

    for (int massId = 0; massId < p.NumberOfMasses(); massId++) {
        r = poly.Build(nSprings + nDampers + 1, this->order, true, massId);
        assert(!r.isError);
        this->models(massId, 0) = poly;
        this->modelsCoefficients[massId] =
            std::vector<double>(poly.nMonomials());

        // Calculate first time derivative
        this->modelsD[massId] = poly;
        assert(!this->modelsD[massId].Dxi(0).isError);

        // Calculate second derivative
        this->modelsDD[massId] = this->modelsD[massId];
        assert(!this->modelsDD[massId].Dxi(0).isError);
    }
}

void Pimodel::SetResidues(bool initialConditions, bool physics) {
    this->initialConditionsResiduesTkc = std::vector<std::vector<Bounded>>();
    this->physicsResiduesTkc = std::vector<std::vector<Bounded>>();
    this->AddResiduesTkc();

    this->initialDispResidues = std::vector<Polys>();
    this->initialVelResidues = std::vector<Polys>();
    if (initialConditions) {
        AddInitialConditionsResidues();
    }
    this->physicsResidues = std::vector<Polys>();
    if (physics) {
        AddPhysicsResidues();
    }

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
    // d/dt(model)
    Maybe<double> pDot;
    for (int massId = 0; massId < int(vels.size()); massId++) {
        this->modelsD[massId].SetX(this->normalizeTkc(TKC));
        pDot = this->modelsD[massId](this->modelsCoefficients[massId]);
        assert(!pDot.isError);
        // d/dt * dt/dT = d/dT
        vels[massId] = pDot.val * this->dtdT;
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

bst::matrix<Polys> Pimodel::getAccelsFromDiffEq(Problem* problem,
                                                std::vector<Bounded>& tkc) {
    // List of displacements and velocities
    bst::matrix<Poly> x = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> dxdT = bst::matrix<Poly>(nMasses, 1);
    for (int massId = 0; massId < nMasses; massId++) {
        x(massId, 0) = this->models(massId, 0);
        x(massId, 0).SetX(Bounded::Get(tkc));

        // d/dT(x_i) = (dt/dT)*d/dt(model_i)
        dxdT(massId, 0) = this->modelsD[massId];
        dxdT(massId, 0) *= this->dtdT;
        dxdT(massId, 0).SetX(Bounded::Get(tkc));
    }

    matrix<Polys> k_times_x = prod(problem->K, x);
    matrix<Polys> c_times_dxdT = prod(problem->C, dxdT);
    // d²/dT2²(x)
    matrix<Polys> d2xdT2 = prod(problem->MInv, k_times_x + c_times_dxdT);
    return d2xdT2;
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

void Pimodel::AddInitialConditionsResiduesTkc() {
    while (int(this->initialConditionsResiduesTkc.size()) <
           this->initialConditionTrainingPoints) {
        std::vector<Bounded> tkc = std::vector<Bounded>(this->inputSize());
        assert(!tkc[0].Set(0.0).isError);

        // First we set 0.0 and 1.0. Then, random
        if (int(this->initialConditionsResiduesTkc.size()) == 0) {
            for (int i = 1; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(0.0).isError);
            }
        } else if (int(this->initialConditionsResiduesTkc.size()) == 1) {
            for (int i = 1; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(1.0).isError);
            }
        } else if (int(this->initialConditionsResiduesTkc.size()) >= 1) {
            for (int i = 1; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(Random()).isError);
            }
        }

        this->initialConditionsResiduesTkc.push_back(tkc);
    }
}
void Pimodel::AddPhysicsResiduesTkc() {
    while (int(this->physicsResiduesTkc.size()) < this->physicsTrainingPoints) {
        std::vector<Bounded> tkc = std::vector<Bounded>(this->inputSize());

        // First we set 0.0 and 1.0. Then, random
        if (int(this->physicsResiduesTkc.size()) == 0) {
            for (int i = 0; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(0.0).isError);
            }
        } else if (int(this->physicsResiduesTkc.size()) == 1) {
            for (int i = 0; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(1.0).isError);
            }
        } else if (int(this->physicsResiduesTkc.size()) >= 1) {
            for (int i = 0; i < int(tkc.size()); i++) {
                assert(!tkc[i].Set(Random()).isError);
            }
        }

        this->physicsResiduesTkc.push_back(tkc);
    }
}
void Pimodel::AddResiduesTkc() {
    this->AddInitialConditionsResiduesTkc();
    this->AddPhysicsResiduesTkc();
}

void Pimodel::AddInitialConditionsResidues() {
    auto initialX = this->getInitialX();

    Poly model;
    Poly modelDot;
    for (auto tkc : this->initialConditionsResiduesTkc) {
        for (int massId = 0; massId < this->nMasses; massId++) {
            model = this->models(massId, 0);
            model.SetX(Bounded::Get(tkc));
            this->initialDispResidues.push_back(
                model +
                (-1) * initialX[Problem::GetMassDispIndex(nMasses, massId)]);

            modelDot = this->modelsD[massId];
            modelDot.SetX(Bounded::Get(tkc));
            this->initialVelResidues.push_back(
                modelDot +
                (-this->dTdt) *
                    initialX[Problem::GetMassVelIndex(nMasses, massId)]);
        }
    }
}

void Pimodel::AddPhysicsResidues() {
    // d²/dt²(model)
    Poly modelDotDot;
    // d²/dT²(x)
    bst::matrix<Polys> d2x_dT2;
    for (auto tkc : this->physicsResiduesTkc) {
        Problem problem = this->problemFromTkc(&tkc);
        d2x_dT2 = getAccelsFromDiffEq(&problem, tkc);
        for (int m = 0; m < nMasses; m++) {
            modelDotDot = this->modelsDD[m];
            modelDotDot.SetX(Bounded::Get(tkc));

            if (problem.massIsFixed(m)) {
                this->physicsResidues.push_back(Polys(modelDotDot));
            } else {
                this->physicsResidues.push_back(Polys(modelDotDot) +
                                                (-this->dTdt * this->dTdt) *
                                                    d2x_dT2(m, 0));
            }
        }
    }
}

Polys* Pimodel::residueById(int i) {
    if (i < int(this->initialDispResidues.size())) {
        return &(this->initialDispResidues[i]);
    }
    i -= this->initialDispResidues.size();
    if (i < int(this->initialVelResidues.size())) {
        return &(this->initialVelResidues[i]);
    }
    i -= this->initialVelResidues.size();
    return &(this->physicsResidues[i]);
}

std::map<std::tuple<int, int>, double>* Pimodel::residueDaById(int i) {
    if (i < int(this->initialDispResidues.size())) {
        return &(this->initialDispResiduesDa[i]);
    }
    i -= this->initialDispResidues.size();
    if (i < int(this->initialVelResidues.size())) {
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
    return this->residueCache[i];
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

std::vector<double> Pimodel::LossGradient(int i) {
    std::vector<double> grad = std::vector<double>(this->nParameters());

    double residueVal = this->Residue(i);
    std::map<std::tuple<int, int>, double>* gradMap = this->residueDaById(i);

    std::tuple<int, int> gradMapKey;
    int gradI = 0;
    for (int massId = 0; massId < int(this->models.size1()); massId++) {
        for (int mon = 0; mon < this->models(massId, 0).nMonomials(); mon++) {
            gradMapKey = {massId, mon};
            if (gradMap->find(gradMapKey) != gradMap->end()) {
                grad[gradI] += residueVal * ((*gradMap)[gradMapKey]);
            }
            gradI += 1;
        }
    }

    return grad;
}

Pimodels::Pimodels(ProblemDescription p, double finalT, int nModels,
                   int initialConditionTrainingPoints,
                   int physicsTrainingPoints, int order) {
    assert(finalT > 0);
    assert(nModels >= 1);
    assert(initialConditionTrainingPoints >= 1);
    assert(physicsTrainingPoints >= 2);
    assert(order >= 0);

    this->finalT = finalT;

    this->timeBuckets = std::vector<double>(nModels + 1);
    double timePerTimeBucket = (finalT - 0) / nModels;
    for (int b = 0; b < nModels; b++) {
        this->timeBuckets[b] = b * timePerTimeBucket;
    }
    this->timeBuckets[nModels] = finalT;

    double t0 = 0;
    this->pimodels = std::vector<Pimodel>(nModels);
    for (int b = 0; b < nModels; b++) {
        this->pimodels[b] = Pimodel(p, t0, t0 + timePerTimeBucket,
                                    initialConditionTrainingPoints,
                                    physicsTrainingPoints, order);
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

void Pimodels::logComplexity() {
    // Set all the residues just so that we can print the complexity stats
    std::cout << "## Complexity stats for each time bucket ##" << std::endl;
    this->pimodels[0].SetResidues(true, true);
    std::cout << "Number of models (1 per mass): "
              << this->pimodels[0].models.size1() << std::endl;
    std::cout << "Number of parameters per model: "
              << this->pimodels[0].models(0, 0).nMonomials() << std::endl;
    std::cout << "Total number of parameters: "
              << this->pimodels[0].nParameters() << std::endl;

    std::cout << "Number of initial displacement residues: "
              << this->pimodels[0].initialDispResidues.size() << std::endl;
    std::cout << "Number of initial velocity residues: "
              << this->pimodels[0].initialVelResidues.size() << std::endl;
    std::cout << "Number of physics residues: "
              << this->pimodels[0].physicsResidues.size() << std::endl;

    // Initial displacement residues cost
    int n = 0;
    for (int i = 0; i < int(this->pimodels[0].initialDispResidues.size());
         i++) {
        n += this->pimodels[0].initialDispResidues[i].nMonomials();
    }
    std::cout << "Total cost of calculating initial displacement residues:" << n
              << std::endl;
    // Initial velocities residues cost
    n = 0;
    for (int i = 0; i < int(this->pimodels[0].initialVelResidues.size()); i++) {
        n += this->pimodels[0].initialVelResidues[i].nMonomials();
    }
    std::cout << "Total cost of calculating initial velocity residues:" << n
              << std::endl;
    // Physics residues cost
    n = 0;
    for (int i = 0; i < int(this->pimodels[0].physicsResidues.size()); i++) {
        n += this->pimodels[0].physicsResidues[i].nMonomials();
    }
    std::cout << "Total cost of calculating physics residues:" << n
              << std::endl;

    // std::cout << "Note: The cost of computing the gradient of each residue "
    //              "is proportional to the cost of calculating the residue. The
    //              " "cost of updating the parameters is proportional to the "
    //              "number of parameters."
    //           << std::endl;

    std::cout << "###########################################" << std::endl;
    std::cout << std::endl;
}

Maybe<double> Pimodels::Train(double initialConditionsLearningRate,
                              double physicsLearningRate,
                              double minImprovementToEarlyStop, int maxSteps,
                              bool logComplexity, bool logTraining) {
    Maybe<double> r;
    if (logComplexity) {
        this->logComplexity();
    }

    std::vector<double> tkcCont = this->continuityTkc();
    int tBucket = 0;
    while (tBucket < int(this->pimodels.size())) {
        if (logTraining) {
            std::cout << "## Training model from " << this->timeBuckets[tBucket]
                      << " to " << this->timeBuckets[tBucket + 1] << " ##"
                      << std::endl;
        }
        if (tBucket > 0) {
            this->setContinuity(tBucket, tkcCont);
        }

        if (logTraining) {
            std::cout << "## Training initial conditions ##" << std::endl;
        }
        this->pimodels[tBucket].SetResidues(true, false);
        r = this->pimodels[tBucket].Train(initialConditionsLearningRate,
                                          minImprovementToEarlyStop, maxSteps,
                                          logTraining);
        // std::cout << "## Training physics ##" << std::endl;
        // this->pimodels[tBucket].SetResidues(false, true);
        // this->pimodels[tBucket].Train(physicsLearningRate,
        //                               minImprovementToEarlyStop, maxSteps,
        //                               logTraining);
        if (logTraining) {
            std::cout << "## Training physics and initial conditions ##"
                      << std::endl;
        }
        this->pimodels[tBucket].SetResidues(true, true);
        this->pimodels[tBucket].Train(physicsLearningRate,
                                      minImprovementToEarlyStop, maxSteps,
                                      logTraining);

        tBucket += 1;
        if (logTraining) {
            std::cout << std::endl;
        }
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

double Pimodels::GetMaxAbsAccel(int massId, int timesChecked,
                                std::vector<Bounded>& kc) {
    assert(this->pimodels.size() > 0);
    assert(int(kc.size()) + 1 == this->pimodels[0].inputSize());
    assert(massId < this->pimodels[0].p.NumberOfMasses());
    assert(timesChecked >= 2);

    double timeStep = 1.0 / (timesChecked - 1);

    std::vector<Bounded> tkc = std::vector<Bounded>(kc.size() + 1);
    for (int i = 1; i < int(tkc.size()); i++) {
        tkc[i] = kc[i - 1];
    }

    double maxA = 0;
    for (int b = 0; b < int(this->pimodels.size()); b++) {
        Pimodel& model = this->pimodels[b];
        double t = 0;
        while (t <= 1.0 + 1e-14) {
            assert(!tkc[0].Set(t).isError);
            model.modelsDD[massId].SetX(Bounded::Get(tkc));
            Maybe<double> modelDotDot =
                model.modelsDD[massId](model.modelsCoefficients[massId]);
            assert(!modelDotDot.isError);
            maxA =
                std::max(maxA, abs(modelDotDot.val * model.dtdT * model.dtdT));
            t += timeStep;
        }
    }

    return maxA;
}
