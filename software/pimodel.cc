#include "pimodel.h"

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

Pimodel::Pimodel(ProblemDescription* p, double finalT, int timeDiscretization,
                 int kcDiscretization, int order) {
    assert(p->IsOk());
    assert(timeDiscretization >= 1);
    assert(kcDiscretization >= 1);
    this->p = p;
    int nSprings = p->springs.size();
    int nDampers = p->dampers.size();

    // The model is a set polynomials which represents the dynamic response of
    // the system. The first polynomial represents the position of the first
    // mass , and so on. The input of each polynomial is the values of the
    // springs, the dampers, and time.
    this->polys = bst::matrix<Poly>(p->NumberOfMasses(), 1);
    Poly poly;
    Maybe<Void> r;
    for (int i = 0; i < p->NumberOfMasses(); i++) {
        r = poly.Build(nSprings + nDampers + 1, order);
        assert(!r.isError);
        this->polys(i, 0) = poly;
    }

    this->finalT = finalT;
    this->timeDiscretization = timeDiscretization;
    this->kcDiscretization = kcDiscretization;
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
    std::vector<double> positions =
        std::vector<double>(this->p->NumberOfMasses());
    Maybe<double> position;
    for (int i = 0; i < int(positions.size()); i++) {
        position = this->polys(i, 0)(tkc);
        assert(!position.isError);
        positions[i] = position.val;
    }
    r.val = positions;
    return r;
}

int Pimodel::nParameters() {
    return this->polys(0, 0).nTerms * this->p->NumberOfMasses();
}

Maybe<Void> Pimodel::GetParameters(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nParameters()) {
        r.isError = true;
        r.errMsg = "Invalid target length";
        return r;
    }

    int nTerms = this->polys(0, 0).nTerms;
    std::vector<double> coefs = std::vector<double>(nTerms);
    int targetI = 0;
    for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
        r = this->polys(pIndex, 0).GetCoefficients(&coefs);
        if (r.isError) {
            return r;
        }
        for (int i = 0; i < nTerms; i++) {
            target->at(targetI) = coefs[i];
            targetI += 1;
        }
    }
    return r;
}

Maybe<Void> Pimodel::SetParameters(std::vector<double>* parameters) {
    Maybe<Void> r;
    if (int(parameters->size()) != this->nParameters()) {
        r.isError = true;
        r.errMsg = "Invalid parameters length";
        return r;
    }

    int nTerms = this->polys(0, 0).nTerms;
    std::vector<double> coefs = std::vector<double>(nTerms);
    int parametersI = 0;
    for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
        for (int i = 0; i < nTerms; i++) {
            coefs[i] = parameters->at(parametersI);
            parametersI += 1;
        }
        r = this->polys(pIndex, 0).SetCoefficients(&coefs);
        if (r.isError) {
            return r;
        }
    }
    return r;
}

Problem Pimodel::problemFromTkc(std::vector<double>* tkc) {
    std::vector<double> DNA = std::vector<double>(tkc->size() - 1);
    std::copy(tkc->begin() + 1, tkc->end(), DNA.begin());
    Maybe<Problem> problem = this->p->BuildFromVector(DNA);
    assert(!problem.isError);
    return problem.val;
}

std::vector<double> Pimodel::getXModel(std::vector<double>* tkc) {
    auto X = std::vector<double>(this->p->NumberOfMasses() * 2);

    Maybe<Void> err;
    Maybe<double> eval;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        // Fill positions
        eval = this->polys(i, 0)(tkc);
        assert(!eval.isError);
        X[i] = eval.val;
    }
    Poly dp_dt;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        // Fill velocities
        dp_dt = this->polys(i, 0);
        err = dp_dt.Dxi(0);
        assert(!err.isError);
        eval = dp_dt(tkc);
        assert(!eval.isError);
        X[Problem::GetMassVelIndex(this->p->NumberOfMasses(), i)] = eval.val;
    }
    return X;
}

std::vector<Poly> Pimodel::getAccelsFromModel() {
    std::vector<Poly> A = std::vector<Poly>(this->p->NumberOfMasses());
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        A[i] = this->polys(i, 0);
        assert(!A[i].Dxi(0).isError);
        assert(!A[i].Dxi(0).isError);
    }
    return A;
}

bst::matrix<Poly> Pimodel::getAccelsFromDiffEq(Problem* problem) {
    // See SetXDot at problem.h for reference
    int nMasses = problem->masses.size();
    std::vector<double> coefs = std::vector<double>(4);

    // List of displacements and velocities
    bst::matrix<Poly> Disps = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> Vels = bst::matrix<Poly>(nMasses, 1);
    for (int i = 0; i < nMasses; i++) {
        Disps(i, 0) = this->polys(i, 0);
        Disps(i, 0).GetCoefficients(&coefs);
        Vels(i, 0) = this->polys(i, 0);
        Vels(i, 0).GetCoefficients(&coefs);
        assert(!Vels(i, 0).Dxi(0).isError);
        Vels(i, 0).GetCoefficients(&coefs);
    }

    Disps(0, 0).GetCoefficients(&coefs);
    Disps(1, 0).GetCoefficients(&coefs);
    double k = problem->K(0, 0);
    k = problem->K(0, 1);
    k = problem->K(1, 0);
    k = problem->K(1, 1);
    matrix<Poly> kx = prod(problem->K, Disps);
    kx(0, 0).GetCoefficients(&coefs);
    kx(1, 0).GetCoefficients(&coefs);
    matrix<Poly> cxDot = prod(problem->C, Vels);
    cxDot(0, 0).GetCoefficients(&coefs);
    cxDot(1, 0).GetCoefficients(&coefs);
    matrix<Poly> Accels = prod(problem->MInv, kx + cxDot);
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

void Pimodel::PhysicsLossDfs(std::vector<double>* tkc, int tkcIndex,
                             double* loss) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = this->problemFromTkc(tkc);

        bst::matrix<Poly> AccelsFromDiffEq = getAccelsFromDiffEq(&problem);

        std::vector<Poly> AccelsFromModel = this->getAccelsFromModel();

        Maybe<double> residueEval;
        Poly residue;
        for (int m = 0; m < int(p->masses.size()); m++) {
            residue = AccelsFromModel[m] + (-1) * AccelsFromDiffEq(m, 0);
            residueEval = residue(tkc);
            assert(!residueEval.isError);
            (*loss) = (*loss) + pow(residueEval.val, 2);
        }
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
    for (int i = 0; i <= this->kcDiscretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / this->kcDiscretization * i;
        this->PhysicsLossDfs(tkc, tkcIndex + 1, loss);
    }
}

void Pimodel::InitialConditionsLossDfs(std::vector<double>* tkc, int tkcIndex,
                                       double* loss) {
    if (tkcIndex == int(tkc->size())) {
        std::vector<double> XModel = this->getXModel(tkc);

        // Create a state vector with the true values of initial displacements
        // and velocities
        auto initialX = this->getInitialX();

        // Update the loss value with the errors in initial displacement
        // and velocities
        int nMasses = this->p->masses.size();
        for (int m = 0; m < nMasses; m++) {
            (*loss) = (*loss) +
                      pow(XModel[Problem::GetMassDispIndex(nMasses, m)] -
                              initialX[Problem::GetMassDispIndex(nMasses, m)],
                          2);
        }
        for (int m = 0; m < nMasses; m++) {
            (*loss) = (*loss) +
                      pow(XModel[Problem::GetMassVelIndex(nMasses, m)] -
                              initialX[Problem::GetMassVelIndex(nMasses, m)],
                          2);
        }
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
        this->InitialConditionsLossDfs(tkc, tkcIndex + 1, loss);
    }
}

double Pimodel::Loss() {
    std::vector<double> tkc = std::vector<double>(this->inputSize());
    double loss = 0;

    // Add initial displacements and velocities loss
    tkc[0] = 0;  // t = 0
    this->InitialConditionsLossDfs(&tkc, 1, &loss);

    // Add physics loss
    this->PhysicsLossDfs(&tkc, 0, &loss);
    return loss;
}

Maybe<Void> Da(std::vector<Poly>* polys, std::vector<double>* X,
               std::vector<double>* target) {
    Maybe<Void> r;
    for (int i = 0; i < int(polys->size()) - 1; i++) {
        if ((*polys)[i].n != (*polys)[i + 1].n) {
            r.isError = true;
            r.errMsg = "all polynomials must have the same number of inputs";
            return r;
        }
    }
    if (int(X->size()) != (*polys)[0].n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    int nTermsTotal = 0;
    for (int i = 0; i < int(polys->size()); i++) {
        nTermsTotal += (*polys)[i].nTerms;
    }
    if (int(target->size()) != nTermsTotal) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }

    int gradIndex = 0;

    std::vector<double> piGrad;
    for (int pi = 0; pi < int((*polys).size()); pi++) {
        piGrad = std::vector<double>((*polys)[pi].nTerms);
        (*polys)[pi].Da(X, &piGrad);
        for (int i = 0; i < int(piGrad.size()); i++) {
            target->at(gradIndex) = piGrad[i];
            gradIndex += 1;
        }
    }
    return r;
}

Maybe<Void> Da(boost::numeric::ublas::matrix<Poly>* polys,
               std::vector<double>* X, std::vector<double>* target) {
    Maybe<Void> r;
    if ((*polys).size2() != 1) {
        r.isError = true;
        r.errMsg = "(*polys) must be column matrix";
        return r;
    }
    for (int i = 0; i < int((*polys).size1()) - 1; i++) {
        if ((*polys)(i, 0).n != (*polys)(i + 1, 0).n) {
            r.isError = true;
            r.errMsg = "all polynomials must have the same number of inputs";
            return r;
        }
    }
    if (int(X->size()) != (*polys)(0, 0).n) {
        r.isError = true;
        r.errMsg = "X of invalid length";
        return r;
    }
    int nTermsTotal = 0;
    for (int i = 0; i < int((*polys).size1()); i++) {
        nTermsTotal += (*polys)(i, 0).nTerms;
    }
    if (int(target->size()) != nTermsTotal) {
        r.isError = true;
        r.errMsg = "target must have same length as the number of terms";
        return r;
    }

    int gradIndex = 0;

    std::vector<double> piGrad;
    for (int pi = 0; pi < int((*polys).size1()); pi++) {
        piGrad = std::vector<double>((*polys)(pi, 0).nTerms);
        (*polys)(pi, 0).Da(X, &piGrad);
        for (int i = 0; i < int(piGrad.size()); i++) {
            target->at(gradIndex) = piGrad[i];
            gradIndex += 1;
        }
    }
    return r;
}

void Pimodel::InitialConditionsLossGradientDfs(std::vector<double>* tkc,
                                               int tkcIndex,
                                               std::vector<double>* grad) {
    if (tkcIndex == int(tkc->size())) {
        std::vector<double> XModel = this->getXModel(tkc);

        // Create a state vector with the true values of initial displacements
        // and velocities
        auto initialX = this->getInitialX();

        int nMasses = p->masses.size();

        // Set the gradient
        std::vector<Poly> losses;
        losses.reserve(nMasses);
        for (int m = 0; m < nMasses; m++) {
            losses.push_back(
                2 *
                    (XModel[Problem::GetMassDispIndex(nMasses, m)] -
                     initialX[Problem::GetMassDispIndex(nMasses, m)]) *
                    this->polys(m, 0) +
                2 *
                    (XModel[Problem::GetMassVelIndex(nMasses, m)] -
                     initialX[Problem::GetMassVelIndex(nMasses, m)]) *
                    this->polys(m, 0));
        }

        assert(!Da(&losses, tkc, grad).isError);
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
        this->InitialConditionsLossGradientDfs(tkc, tkcIndex + 1, grad);
    }
}

// THIS IS WRONG AND MUST BE FIXED
// XDotDotFromDiffEq also depends on the parameters of the model,
// thus it can't be treated like a constant as in the initial conditions loss
// gradient
void Pimodel::PhysicsLossGradientDfs(std::vector<double>* tkc, int tkcIndex,
                                     std::vector<double>* grad) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = this->problemFromTkc(tkc);

        bst::matrix<Poly> AccelsFromDiffEq = getAccelsFromDiffEq(&problem);

        std::vector<Poly> AccelsFromModel = this->getAccelsFromModel();

        int nMasses = problem.masses.size();

        auto accels_from_model_m_grad =
            std::vector<double>(this->polys(0, 0).nTerms);
        auto accels_from_diff_eq_m_grad =
            std::vector<double>(this->polys(0, 0).nTerms);
        Poly residue;
        Maybe<double> residueEval;
        int gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            residue = AccelsFromModel[m] + (-1) * AccelsFromDiffEq(m, 0);
            residueEval = residue(tkc);
            assert(!residueEval.isError);

            AccelsFromModel[m].Da(tkc, &accels_from_model_m_grad);
            AccelsFromDiffEq(m, 0).Da(tkc, &accels_from_diff_eq_m_grad);
            for (int i = 0; i < int(accels_from_model_m_grad.size()); i++) {
                grad->at(gradIndex) +=
                    2 * (residueEval.val) * accels_from_model_m_grad[i];
                grad->at(gradIndex) +=
                    2 * (residueEval.val) * accels_from_diff_eq_m_grad[i];
                gradIndex += 1;
            }
        }
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
    for (int i = 0; i <= this->kcDiscretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / this->kcDiscretization * i;
        this->PhysicsLossGradientDfs(tkc, tkcIndex + 1, grad);
    }
}

std::vector<double> Pimodel::LossGradient() {
    std::vector<double> tkc = std::vector<double>(this->inputSize());
    std::vector<double> grad = std::vector<double>(this->nParameters());

    // AInitial displacements and velocities gradient
    tkc[0] = 0;  // t = 0
    this->InitialConditionsLossGradientDfs(&tkc, 1, &grad);

    // Physics loss gradient
    this->PhysicsLossGradientDfs(&tkc, 0, &grad);

    return grad;
}