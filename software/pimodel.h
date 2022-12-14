#pragma once
#include <gtest/gtest.h>

#include "maybe.h"
#include "model.h"
#include "polynomial.h"
#include "problem_description.h"

class Pimodel : public Model {
   public:
    // Creates multivariate polynomial models of a given order that is physics
    // informed, i.e. has physics knowledge embedded into its loss function.
    // Each model corresponds to a function that describes the position of one
    // of the masses as a function of time, the values of the springs
    // and the values of the dampers.
    // Parameters:
    // p: Description of this problem
    // finalT: this model should model the position of the mass from t=0 to
    // finalT
    // timeDiscretization: How many intervals we'll discretize time
    // kcDiscretization: How many intervals we'll discretize
    // springs/dampers values in our loss function oder: order of the
    // multivariate polynomial used in the model
    Pimodel(ProblemDescription* p, double finalT, int timeDiscretization,
            int kcDiscretization, int order);

    // Returns the position of each mass.
    // The input should be an array with:
    // The time
    // The values of the springs
    // The values of the dampers
    // U = [t, k1, k2, ... , kn, c1, c2, ... , cn]
    Maybe<std::vector<double>> operator()(std::vector<double>* tkc);

    int nParameters() override;

   private:
    FRIEND_TEST(PimodelTest, ConstructorTest);
    FRIEND_TEST(PimodelTest, GetParametersTest);
    FRIEND_TEST(PimodelTest, SetParametersTest);
    FRIEND_TEST(PimodelTest, getAccelsFromDiffEqTest);
    FRIEND_TEST(PimodelTest, InitialConditionsLossTest);
    FRIEND_TEST(PimodelTest, PhysicsLossTest);
    FRIEND_TEST(PimodelTest, LossTest);
    FRIEND_TEST(PimodelTest, LossTest2);
    FRIEND_TEST(PimodelTest, LossGradientTest);
    FRIEND_TEST(PimodelTest, LossGradientTest2);

    ProblemDescription* p;

    // column matrix that contains the polynomials that represent the
    // displacement of each mass. Ex: polys[0] = Polynomial that represents the
    // displacement (x) of mass 0, as a function of time and the values of the
    // springs and the masses
    boost::numeric::ublas::matrix<Poly> polys;

    // This model describes the system from t = 0 to t = finalT
    double finalT;

    // Parameter that describes in how many intervals we'll discretize time and
    // the springs/dampers values in our loss function
    int timeDiscretization;
    int kcDiscretization;

    // Returns the size of the input of this model.
    // This will always be the number of springs of the problem
    // plus the number of dampers plus 1 (time).
    int inputSize();

    Maybe<Void> GetParameters(std::vector<double>* target) override;

    Maybe<Void> SetParameters(std::vector<double>* parameters) override;

    double Loss() override;
    // Auxiliary functions used inside Loss() that calculates the losses
    // using "Depth-First-Search"
    void InitialConditionsLossDfs(std::vector<double>* tkc, int tkcIndex,
                                  double* loss);
    void PhysicsLossDfs(std::vector<double>* tkc, int tkcIndex, double* loss);

    std::vector<double> LossGradient() override;
    // Auxiliary functions used inside LossGradient() that calculates the
    // gradients using "Depth-First-Search"
    void InitialConditionsLossGradientDfs(std::vector<double>* tkc,
                                          int tkcIndex,
                                          std::vector<double>* grad);
    void PhysicsLossGradientDfs(std::vector<double>* tkc, int tkcIndex,
                                std::vector<double>* grad);

    // Auxiliary functions:

    // Create problem for given tkc (time, springs and dampers)
    Problem problemFromTkc(std::vector<double>* tkc);

    // Return the state vector with the initial conditions set
    std::vector<double> getInitialX();

    // Calculate the state vector for given tkc (time, springs and dampers)
    // using the polynomials
    std::vector<double> getXModel(std::vector<double>* tkc);

    // Get polynomials that model the accelerations by  taking second time
    // derivative of the polynomials
    std::vector<Poly> getAccelsFromModel();

    // Gets the polynomials that represent the acceleration of each mass
    // by using the differential equation from the discrete element method
    boost::numeric::ublas::matrix<Poly> getAccelsFromDiffEq(Problem* problem);
};

// Sets the gradient of each polynomial at target, one after the other
// Ex:
// p0 = 1*x + 2*1
// p1 = 3*x + 4*1
// polys = [p0, p1]
// grad(p0) = [x, 2]
// grad(p1) = [x, 4]
// for X = [9.0] -> target = [9.0, 2, 9.0, 4]
Maybe<Void> Da(std::vector<Poly>* polys, std::vector<double>* X,
               std::vector<double>* target);
// Same as the above, but input is column matrix
Maybe<Void> Da(boost::numeric::ublas::matrix<Poly>* polys,
               std::vector<double>* X, std::vector<double>* target);