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
    this->models = bst::matrix<Poly>(p->NumberOfMasses(), 1);
    Poly poly;
    Maybe<Void> r;
    for (int i = 0; i < p->NumberOfMasses(); i++) {
        r = poly.Build(nSprings + nDampers + 1, order);
        assert(!r.isError);
        this->models(i, 0) = poly;
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
        position = this->models(i, 0)(*tkc);
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
        return this->models(0, 0).nMonomials() * this->p->NumberOfMasses();
    }
    int rv = 0;
    for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
        rv += this->models(pIndex, 0).nMonomials();
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
    Maybe<double> Dai;
    for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
        for (int m = 0; m < this->models(pIndex, 0).nMonomials(); m++) {
            this->models(pIndex, 0).monomials[m].a =
                parameters->at(parametersI);
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
    Maybe<double> Dai;
    for (int pIndex = 0; pIndex < this->p->NumberOfMasses(); pIndex++) {
        for (int m = 0; m < this->models(pIndex, 0).nMonomials(); m++) {
            target->at(targetI) = this->models(pIndex, 0).monomials[m].a;
            targetI += 1;
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

    Maybe<Void> err;
    Maybe<double> eval;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        // Fill displacements
        eval = this->models(i, 0)(*tkc);
        assert(!eval.isError);
        X[i] = eval.val;
    }
    Poly dp_dt;
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        // Fill velocities
        dp_dt = this->models(i, 0);
        err = dp_dt.Dxi(0);
        assert(!err.isError);
        eval = dp_dt(*tkc);
        assert(!eval.isError);
        X[Problem::GetMassVelIndex(this->p->NumberOfMasses(), i)] = eval.val;
    }
    return X;
}

std::vector<Polys> Pimodel::getAccelsFromModel() {
    std::vector<Polys> A = std::vector<Polys>(this->p->NumberOfMasses());
    for (int i = 0; i < this->p->NumberOfMasses(); i++) {
        A[i] = Polys(this->models(i, 0));
        assert(!A[i].Dxi(0).isError);
        assert(!A[i].Dxi(0).isError);
    }
    return A;
}

bst::matrix<Polys> Pimodel::getAccelsFromDiffEq(Problem* problem) {
    // See SetXDot at problem.h for reference
    int nMasses = problem->masses.size();
    std::vector<double> coefs = std::vector<double>(4);

    // List of displacements and velocities
    bst::matrix<Poly> Disps = bst::matrix<Poly>(nMasses, 1);
    bst::matrix<Poly> Vels = bst::matrix<Poly>(nMasses, 1);
    for (int i = 0; i < nMasses; i++) {
        Disps(i, 0) = this->models(i, 0);
        Vels(i, 0) = this->models(i, 0);
        assert(!Vels(i, 0).Dxi(0).isError);
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
void Pimodel::PhysicsLossDfs(std::vector<double>* tkc, int tkcIndex,
                             double* loss) {
    if (tkcIndex == int(tkc->size())) {
        Problem problem = this->problemFromTkc(tkc);

        bst::matrix<Polys> AccelsFromDiffEq = getAccelsFromDiffEq(&problem);

        std::vector<Polys> AccelsFromModel = this->getAccelsFromModel();

        Maybe<double> residueEval;
        Polys residue;
        for (int m = 0; m < int(p->masses.size()); m++) {
            residue = AccelsFromModel[m] + (-1.0) * AccelsFromDiffEq(m, 0);
            residueEval = residue(*tkc);
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
    for (int i = 0; i <= discretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / discretization * i;
        this->PhysicsLossDfs(tkc, tkcIndex + 1, loss);
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

// Maybe<Void> Da(std::vector<Poly>* polys, std::vector<double>* X,
//                std::vector<double>* target) {
//     Maybe<Void> r;
//     for (int i = 0; i < int(polys->size()) - 1; i++) {
//         if ((*polys)[i].n != (*polys)[i + 1].n) {
//             r.isError = true;
//             r.errMsg = "all polynomials must have the same number of inputs";
//             return r;
//         }
//     }
//     if (int(X->size()) != (*polys)[0].n) {
//         r.isError = true;
//         r.errMsg = "X of invalid length";
//         return r;
//     }
//     int nTermsTotal = 0;
//     for (int i = 0; i < int(polys->size()); i++) {
//         nTermsTotal += (*polys)[i].nMonomials();
//     }
//     if (int(target->size()) != nTermsTotal) {
//         r.isError = true;
//         r.errMsg = "target must have same length as the number of terms";
//         return r;
//     }

//     int gradIndex = 0;

//     std::vector<double> piGrad;
//     for (int pi = 0; pi < int((*polys).size()); pi++) {
//         piGrad = std::vector<double>((*polys)[pi].nMonomials());
//         (*polys)[pi].Da(X, &piGrad);
//         for (int i = 0; i < int(piGrad.size()); i++) {
//             target->at(gradIndex) += piGrad[i];
//             gradIndex += 1;
//         }
//     }
//     return r;
// }

// Maybe<Void> Da(boost::numeric::ublas::matrix<Poly>* polys,
//                std::vector<double>* X, std::vector<double>* target) {
//     Maybe<Void> r;
//     if ((*polys).size2() != 1) {
//         r.isError = true;
//         r.errMsg = "(*polys) must be column matrix";
//         return r;
//     }
//     for (int i = 0; i < int((*polys).size1()) - 1; i++) {
//         if ((*polys)(i, 0).n != (*polys)(i + 1, 0).n) {
//             r.isError = true;
//             r.errMsg = "all polynomials must have the same number of inputs";
//             return r;
//         }
//     }
//     if (int(X->size()) != (*polys)(0, 0).n) {
//         r.isError = true;
//         r.errMsg = "X of invalid length";
//         return r;
//     }
//     int nTermsTotal = 0;
//     for (int i = 0; i < int((*polys).size1()); i++) {
//         nTermsTotal += (*polys)(i, 0).nMonomials();
//     }
//     if (int(target->size()) != nTermsTotal) {
//         r.isError = true;
//         r.errMsg = "target must have same length as the number of terms";
//         return r;
//     }

//     int gradIndex = 0;

//     std::vector<double> piGrad;
//     for (int pi = 0; pi < int((*polys).size1()); pi++) {
//         piGrad = std::vector<double>((*polys)(pi, 0).nMonomials());
//         (*polys)(pi, 0).Da(X, &piGrad);
//         for (int i = 0; i < int(piGrad.size()); i++) {
//             target->at(gradIndex) += piGrad[i];
//             gradIndex += 1;
//         }
//     }
//     return r;
// }

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
        double residue;
        int gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            residue = (XModel[Problem::GetMassDispIndex(nMasses, m)] -
                       initialX[Problem::GetMassDispIndex(nMasses, m)]);

            std::vector<double> piGrad =
                std::vector<double>(this->models(m, 0).nMonomials());
            Maybe<double> maybeDouble;
            for (int a = 0; a < this->models(m, 0).nMonomials(); a++) {
                maybeDouble = this->models(m, 0).Dai(a, *tkc);
                assert(!maybeDouble.isError);
                piGrad[a] = maybeDouble.val;
            }

            for (int i = 0; i < int(piGrad.size()); i++) {
                grad->at(gradIndex) += 2 * residue * piGrad[i];
                gradIndex += 1;
            }
        }
        Poly pVel;
        gradIndex = 0;
        for (int m = 0; m < nMasses; m++) {
            residue = (XModel[Problem::GetMassVelIndex(nMasses, m)] -
                       initialX[Problem::GetMassVelIndex(nMasses, m)]);

            std::vector<double> piGrad =
                std::vector<double>(this->models(m, 0).nMonomials());
            pVel = this->models(m, 0);
            assert(!pVel.Dxi(0).isError);
            Maybe<double> maybeDouble;
            for (int a = 0; a < pVel.nMonomials(); a++) {
                maybeDouble = pVel.Dai(a, *tkc);
                assert(!maybeDouble.isError);
                piGrad[a] = maybeDouble.val;
            }

            for (int i = 0; i < int(piGrad.size()); i++) {
                grad->at(gradIndex) += 2 * residue * piGrad[i];
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
    for (int i = 0; i <= this->kcDiscretization; i++) {
        tkc->at(tkcIndex) = min + (max - min) / this->kcDiscretization * i;
        this->InitialConditionsLossGradientDfs(tkc, tkcIndex + 1, grad);
    }
}

// // THIS IS WRONG AND MUST BE FIXED
// // XDotDotFromDiffEq also depends on the parameters of the model,
// // thus it can't be treated like a constant as in the initial conditions loss
// // gradient
// void Pimodel::PhysicsLossGradientDfs(std::vector<double>* tkc, int tkcIndex,
//                                      std::vector<double>* grad) {
//     if (tkcIndex == int(tkc->size())) {
//         Problem problem = this->problemFromTkc(tkc);

//         bst::matrix<Poly> AccelsFromDiffEq = getAccelsFromDiffEq(&problem);

//         std::vector<Poly> AccelsFromModel = this->getAccelsFromModel();

//         int nMasses = problem.masses.size();

//         auto residueGrad = std::vector<double>(this->models(0,
//         0).nMonomials()); Poly residue; Maybe<double> residueEval; int
//         gradIndex = 0; for (int m = 0; m < nMasses; m++) {
//             residue = AccelsFromModel[m] + (-1) * AccelsFromDiffEq(m, 0);
//             residueEval = residue(tkc);
//             assert(!residueEval.isError);

//             residue.Da(tkc, &residueGrad);
//             for (int i = 0; i < int(residueGrad.size()); i++) {
//                 grad->at(gradIndex) += 2 * (residueEval.val) *
//                 residueGrad[i]; gradIndex += 1;
//             }
//         }
//         return;
//     }

//     double min;
//     double max;
//     int discretization;
//     if (tkcIndex == 0) {
//         min = 0;
//         max = this->finalT;
//         discretization = this->timeDiscretization;
//     } else {
//         if (tkcIndex - 1 < int(this->p->springs.size())) {
//             min = this->p->springs[tkcIndex - 1].kMin;
//             max = this->p->springs[tkcIndex - 1].kMax;
//         } else {
//             min = this->p->dampers[tkcIndex - 1 -
//             this->p->springs.size()].cMin; max = this->p->dampers[tkcIndex -
//             1 - this->p->springs.size()].cMax;
//         }
//         discretization = this->kcDiscretization;
//     }
//     for (int i = 0; i <= discretization; i++) {
//         tkc->at(tkcIndex) = min + (max - min) / discretization * i;
//         this->PhysicsLossGradientDfs(tkc, tkcIndex + 1, grad);
//     }
// }

std::vector<double> Pimodel::LossGradient() {
    // std::vector<double> tkc = std::vector<double>(this->inputSize());
    // std::vector<double> grad = std::vector<double>(this->nParameters());

    // // AInitial displacements and velocities gradient
    // tkc[0] = 0;  // t = 0
    // this->InitialConditionsLossGradientDfs(&tkc, 1, &grad);

    // // Physics loss gradient
    // this->PhysicsLossGradientDfs(&tkc, 0, &grad);

    // return grad;
}