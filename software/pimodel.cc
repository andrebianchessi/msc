#include "pimodel.h"

#include <boost/numeric/ublas/operations.hpp>
#include <cassert>
#include <cmath>
#include <vector>

#include "bounded.h"
#include "maybe.h"
#include "problem.h"
#include "problem_description.h"

Pimodel::Pimodel(ProblemDescription* p, double finalT, int discretization,
                 int order) {
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

// void Pimodel::LossDfs(std::vector<double>* X, int xIndex, double* loss) {
//     if (xIndex == X->size()) {
//         // Array containing values of springs and dampers
//         std::vector<Bounded> DNA = std::vector<Bounded>(X->size() - 1);
//         std::copy(X->begin() + 1, X->end(), DNA.begin());
//         Maybe<Problem> problem = this->p->BuildFromDNA(DNA);
//         assert(!problem.isError);
//         boost::numeric::ublas::vector<double> speedsAndAccel =
//             problem.val.GetXDot()
//     }

//     double minVal;
//     double maxVal;

//     // At this case, we want to set values of time
//     if (xIndex == 0) {
//         minVal = 0;
//         maxVal = this->finalT;
//     }
// }

double Pimodel::Loss() {
    // // vector containing normalized time, followed by
    // // values of springs and dampers
    // std::vector<Bounded> X = std::vector<Bounded>(this->inputSize() - 1);
}

std::vector<double> Pimodel::LossGradient() {}