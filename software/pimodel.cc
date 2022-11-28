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

Pimodel::Pimodel(ProblemDescription* p, double finalT, int discretization,
                 int order) {
    assert(discretization >= 1);
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

    this->massIndex = massIndex;
    this->finalT = finalT;
    this->discretization = discretization;
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

void Pimodel::LossDfs(std::vector<double>* tkc, int tkcIndex, double* loss) {
    if (tkcIndex == int(tkc->size())) {
        print("tkc[0]", tkc->at(0));
        print("tkc[1]", tkc->at(1));
        print("tkc[2]", tkc->at(2));
        // Create problem for given tkc (time, springs and dampers)
        std::vector<double> DNA = std::vector<double>(tkc->size() - 1);
        std::copy(tkc->begin() + 1, tkc->end(), DNA.begin());
        Maybe<Problem> problem = this->p->BuildFromVector(DNA);
        assert(!problem.isError);

        // Create state vector and fill it with model predictions
        auto X = boost::numeric::ublas::vector<double>(
            this->p->NumberOfMasses() * 2);
        Maybe<double> eval;
        for (int i = 0; i < this->p->NumberOfMasses(); i++) {
            // Fill positions
            eval = this->polys[i](tkc);
            assert(!eval.isError);
            X[i] = eval.val;
        }
        for (int massIndex = 0; massIndex < this->p->NumberOfMasses();
             massIndex++) {
            // Fill velocities
            // Note that we calculate the velocities by taking the
            // derivative of the polynomials with respect to time
            eval = this->polys[massIndex].Dxi(0, tkc);
            assert(!eval.isError);
            X[problem.val.GetMassVelIndex(massIndex)] = eval.val;
        }

        // Calculate what the true accelerations should be, i.e.
        // plug the values of displacement and velocity the polynomials
        // predict into the discrete element formula
        Maybe<std::vector<double>> trueAccels =
            problem.val.GetAccel(X, tkc->at(0));
        assert(!trueAccels.isError);

        // Calculate the accelerations by taking the second time derivative
        // of the polynomials
        std::vector<double> accels = std::vector<double>(trueAccels.val.size());
        for (int i = 0; i < this->p->NumberOfMasses(); i++) {
            eval = this->polys[i].D2xi(0, tkc);
            assert(!eval.isError);
            accels[i] = eval.val;
        }

        // Update the loss value, which is the sum of the square error of the
        // accelerations
        for (int i = 0; i < this->p->NumberOfMasses(); i++) {
            (*loss) = (*loss) + pow(accels[i] - trueAccels.val[i], 2);
        }
        return;
    }

    double min;
    double max;
    if (tkcIndex == 0) {
        min = 0;
        max = this->finalT;
    } else {
        if (tkcIndex - 1 < this->p->springs.size()) {
            min = this->p->springs[tkcIndex - 1].kMin;
            max = this->p->springs[tkcIndex - 1].kMax;
        } else {
            min = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMin;
            max = this->p->dampers[tkcIndex - 1 - this->p->springs.size()].cMax;
        }
    }
    for (double val = min; val <= max;
         val += (max - min) / this->discretization) {
        tkc->at(tkcIndex) = val;
        this->LossDfs(tkc, tkcIndex + 1, loss);
    }
}

double Pimodel::Loss() {
    std::vector<double> tkc = std::vector<double>(this->inputSize());
    int tkcI = 1;
    double loss = 0;

    // Add initial displacements loss

    // Add initial velocities loss

    // Add physical loss
    this->LossDfs(&tkc, 0, &loss);
    return loss;
}

std::vector<double> Pimodel::LossGradient() {}