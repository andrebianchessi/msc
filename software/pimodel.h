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
    FRIEND_TEST(PimodelTest, ProblemFromTkcTest);
    FRIEND_TEST(PimodelTest, getXModelTest);
    FRIEND_TEST(PimodelTest, getAccelsFromDiffEqTest);
    FRIEND_TEST(PimodelTest, getInitialXTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTkcTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesTkcTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesTest);
    FRIEND_TEST(PimodelTest, LossTest);
    FRIEND_TEST(PimodelTest, LossGradientTest);

    ProblemDescription* p;

    // column matrix that contains the polynomials that represent the
    // displacement of each mass. Ex: models[0] = Polynomial that represents the
    // displacement (x) of mass 0, as a function of time and the values of the
    // springs and the masses
    boost::numeric::ublas::matrix<Poly> models;
    std::vector<std::vector<double>> modelsCoefficients;

    // Residues of initial displacement
    std::vector<Polys> initialDispResidues;
    // Residues of initial velocity
    std::vector<Polys> initialVelResidues;
    // Values of time, spring coefficients and damper coefficients in which
    // the initial displacement and velocity residues must be calculated to make
    // up the loss function
    std::vector<std::vector<double>> initialConditionsResiduesTkc;

    // Physical residue
    std::vector<Polys> physicsResidues;
    std::vector<std::vector<double>> physicsResiduesTkc;

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

    void AddInitialConditionsResiduesTkc(std::vector<double>* tkc,
                                         int tkcIndex);
    void AddPhysicsResiduesTkc(std::vector<double>* tkc, int tkcIndex);
    void AddResiduesTkc();

    void AddInitialConditionsResidues();
    void AddPhysicsResidues();
    void AddResidues();

    std::vector<double> LossGradient() override;

    // Auxiliary functions:

    // Create problem for given tkc (time, springs and dampers)
    Problem problemFromTkc(std::vector<double>* tkc);

    // Return the state vector with the initial conditions set
    std::vector<double> getInitialX();

    // Calculate the state vector for given tkc (time, springs and dampers)
    // using the polynomials
    std::vector<double> getXModel(std::vector<double>* tkc);

    // Gets the polynomials that represent the acceleration of each mass
    // by using the differential equation from the discrete element method
    boost::numeric::ublas::matrix<Polys> getAccelsFromDiffEq(
        Problem* problem, std::vector<double>& tkc);
};
