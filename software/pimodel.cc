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
    this->polys = std::vector<Poly>(p->NumberOfMasses());
    Maybe<Poly> poly;
    for (int i = 0; i < p->NumberOfMasses(); i++) {
        poly = Poly::NewPoly(nSprings + nDampers + 1, order);
        assert(!poly.isError);
        this->polys[i] = poly.val;
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
        position = this->polys[i](tkc);
        assert(!position.isError);
        positions[i] = position.val;
    }
    r.val = positions;
    return r;
}

int Pimodel::nParameters() {
    return this->polys[0].nTerms * this->p->NumberOfMasses();
}

Maybe<Void> Pimodel::GetParameters(std::vector<double>* target) {
    Maybe<Void> r;
    if (int(target->size()) != this->nParameters()) {
        r.isError = true;
        r.errMsg = "Invalid target length";
        return r;
    }

    int nTerms = this->polys[0].nTerms;
    std::vector<double> coefs = std::vector<double>(nTerms);
    int targetI = 0;
    for (int pIndex = 0; pIndex < int(this->polys.size()); pIndex++) {
        r = this->polys[pIndex].GetCoefficients(&coefs);
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

    int nTerms = this->polys[0].nTerms;
    std::vector<double> coefs = std::vector<double>(nTerms);
    int parametersI = 0;
    for (int pIndex = 0; pIndex < int(this->polys.size()); pIndex++) {
        for (int i = 0; i < nTerms; i++) {
            coefs[i] = parameters->at(parametersI);
            parametersI += 1;
        }
        r = this->polys[pIndex].SetCoefficients(&coefs);
        if (r.isError) {
            return r;
        }
    }
    return r;
}

Problem Pimodel::problemFromTkc(Pimodel* model, std::vector<double>* tkc) {
    std::vector<double> DNA = std::vector<double>(tkc->size() - 1);
    std::copy(tkc->begin() + 1, tkc->end(), DNA.begin());
    Maybe<Problem> problem = model->p->BuildFromVector(DNA);
    assert(!problem.isError);
    return problem.val;
}

boost::numeric::ublas::vector<double> Pimodel::getXModel(
    Pimodel* model, std::vector<double>* tkc, Problem* problem) {
    auto X =
        boost::numeric::ublas::vector<double>(model->p->NumberOfMasses() * 2);

    Maybe<double> eval;
    for (int i = 0; i < model->p->NumberOfMasses(); i++) {
        // Fill positions
        eval = model->polys[i](tkc);
        assert(!eval.isError);
        X[i] = eval.val;
    }
    for (int i = 0; i < model->p->NumberOfMasses(); i++) {
        // Fill velocities
        // Note that we calculate the velocities by taking the
        // derivative of the polynomials with respect to time
        eval = model->polys[i].Dxi(0, tkc);
        assert(!eval.isError);
        X[problem->GetMassVelIndex(i)] = eval.val;
    }
    return X;
}

std::vector<double> Pimodel::getXModelDotDot(Pimodel* model,
                                             std::vector<double>* tkc) {
    std::vector<double> XModelDotDot =
        std::vector<double>(model->p->NumberOfMasses());
    Maybe<double> XiModelDotDot;
    for (int i = 0; i < model->p->NumberOfMasses(); i++) {
        XiModelDotDot = model->polys[i].D2xi(0, tkc);
        assert(!XiModelDotDot.isError);
        XModelDotDot[i] = XiModelDotDot.val;
    }
    return XModelDotDot;
}

std::vector<double> Pimodel::getXDotDotFromDiffEq(
    Problem* problem, boost::numeric::ublas::vector<double> X, double t) {
    Maybe<std::vector<double>> XDotDot_XModel = problem->GetAccel(X, t);
    assert(!XDotDot_XModel.isError);
    return XDotDot_XModel.val;
}

std::vector<double> Pimodel::getInitialX(Problem* problem) {
    int nMasses = problem->masses.size();
    auto initialX = std::vector<double>(nMasses * 2);
    for (auto v : p->initialVels) {
        assert(v.massId >= 0 && v.massId < nMasses);
        initialX[Problem::GetMassVelIndex(nMasses, v.massId)] = v.val;
    }
    for (auto d : p->initialDisps) {
        assert(d.massId >= 0 && d.massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, d.massId)] = d.val;
    }
    for (auto massId : p->fixedMasses) {
        assert(massId >= 0 && massId < nMasses);
        initialX[Problem::GetMassDispIndex(nMasses, massId)] = 0;
        initialX[Problem::GetMassVelIndex(nMasses, massId)] = 0;
    }
    return initialX;
}

void Pimodel::PhysicsLossDfs(std::vector<double>* tkc, int tkcIndex,
                             double* loss) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = problemFromTkc(this, tkc);

        boost::numeric::ublas::vector<double> XModel =
            getXModel(this, tkc, &problem);

        std::vector<double> XDotDotFromDiffEq =
            getXDotDotFromDiffEq(&problem, XModel, tkc->at(0));

        std::vector<double> XModelDotDot = getXModelDotDot(this, tkc);

        for (int i = 0; i < this->p->NumberOfMasses(); i++) {
            (*loss) = (*loss) + pow(XModelDotDot[i] - XDotDotFromDiffEq[i], 2);
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
    for (double val = min; val <= max; val += (max - min) / discretization) {
        tkc->at(tkcIndex) = val;
        this->PhysicsLossDfs(tkc, tkcIndex + 1, loss);
    }
}

void Pimodel::InitialConditionsLossDfs(std::vector<double>* tkc, int tkcIndex,
                                       double* loss) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = problemFromTkc(this, tkc);

        boost::numeric::ublas::vector<double> XModel =
            getXModel(this, tkc, &problem);

        // Create a state vector with the true values of initial displacements
        // and velocities
        auto initialX = getInitialX(&problem);

        // Update the loss value with the errors in initial displacement
        // and velocities
        int nMasses = problem.masses.size();
        for (int m = 0; m < nMasses; m++) {
            (*loss) = (*loss) + pow(XModel[problem.GetMassDispIndex(m)] -
                                        initialX[problem.GetMassDispIndex(m)],
                                    2);
        }
        for (int m = 0; m < nMasses; m++) {
            (*loss) = (*loss) + pow(XModel[problem.GetMassVelIndex(m)] -
                                        initialX[problem.GetMassVelIndex(m)],
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
    for (double val = min; val <= max;
         val += (max - min) / this->kcDiscretization) {
        tkc->at(tkcIndex) = val;
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

void Pimodel::InitialConditionsLossGradientDfs(std::vector<double>* tkc,
                                               int tkcIndex,
                                               std::vector<double>* grad) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = problemFromTkc(this, tkc);

        boost::numeric::ublas::vector<double> XModel =
            getXModel(this, tkc, &problem);

        // Create a state vector with the true values of initial displacements
        // and velocities
        auto initialX = getInitialX(&problem);

        int nMasses = p->masses.size();

        // Set the gradient
        std::vector<double> XmModelGrad =
            std::vector<double>(this->polys[0].nTerms);
        double initialXm;
        double XmModel;

        // Initial displacement error
        int gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            initialXm = initialX[Problem::GetMassDispIndex(nMasses, m)];
            XmModel = XModel[Problem::GetMassDispIndex(nMasses, m)];
            this->polys[m].Da(tkc, &XmModelGrad);
            for (int i = 0; i < int(XmModelGrad.size()); i++) {
                grad->at(gradIndex) +=
                    2 * (XmModel - initialXm) * XmModelGrad[i];
                gradIndex += 1;
            }
        }
        // Initial vel error
        gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            initialXm = initialX[Problem::GetMassVelIndex(nMasses, m)];
            XmModel = XModel[Problem::GetMassVelIndex(nMasses, m)];
            this->polys[m].Da(tkc, &XmModelGrad);
            for (int i = 0; i < int(XmModelGrad.size()); i++) {
                grad->at(gradIndex) +=
                    2 * (XmModel - initialXm) * XmModelGrad[i];
                gradIndex += 1;
            }
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
    for (double val = min; val <= max;
         val += (max - min) / this->kcDiscretization) {
        tkc->at(tkcIndex) = val;
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
        Problem problem = problemFromTkc(this, tkc);

        boost::numeric::ublas::vector<double> XModel =
            getXModel(this, tkc, &problem);

        std::vector<double> XDotDotFromDiffEq =
            getXDotDotFromDiffEq(&problem, XModel, tkc->at(0));

        std::vector<double> XModelDotDot = getXModelDotDot(this, tkc);

        int nMasses = problem.masses.size();

        // Set gradient value
        auto XModelDotDotGrad = std::vector<double>(this->polys[0].nTerms);
        int gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            this->polys[m].Da(tkc, &XModelDotDotGrad);
            for (int i = 0; i < int(XModelDotDotGrad.size()); i++) {
                grad->at(gradIndex) +=
                    2 * (XModelDotDot[m] - XDotDotFromDiffEq[m]) *
                    XModelDotDotGrad[i];
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
    for (double val = min; val <= max; val += (max - min) / discretization) {
        tkc->at(tkcIndex) = val;
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