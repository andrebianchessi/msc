#pragma once
#include <gtest/gtest.h>

#include "maybe.h"
#include "model.h"
#include "polynomial.h"
#include "problem_description.h"

class Pimodel : public Model {
   public:
    // Parameters:
    // p: Description of this problem
    // finalT: this model models the position of the mass for t in [0,finalT]
    // nModels: the interval [0, finalT] will be divided into nTimeBuckets
    // time buckets. For every mass, there'll be 1 model for each interval.
    // timeDiscretization: How many intervals we'll discretize time in
    // every bucket.
    // kcDiscretization: How many intervals we'll discretize
    // springs/dampers values in our loss function
    // order: order of the multivariate polynomial used in the models
    Pimodel(ProblemDescription* p, double finalT, int nTimeBuckets,
            int timeDiscretization, int kcDiscretization, int order);

    // Returns the position of each mass.
    // The input should be an array with:
    // The time
    // The values of the springs
    // The values of the dampers
    // U = [t, k1, k2, ... , kn, c1, c2, ... , cn]
    Maybe<std::vector<double>> operator()(std::vector<double>* tkc);

    int nParameters() override;

   private:
    friend class PimodelTest;
    FRIEND_TEST(PimodelTest, ConstructorTest);
    FRIEND_TEST(PimodelTest, TimeBucketTest);
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
    FRIEND_TEST(PimodelTrainingTest, TrainTest);

    ProblemDescription* p;

    // Start and end times of all time buckets.
    // Ex: for 2 time buckets -> [0, tFinal/2, tFinal]
    std::vector<double> timeBuckets;
    // Returns index at timeBuckets
    // For the example above:
    // timeBucket(0) -> 0
    // timeBucket(tFinal/4) -> 0
    // timeBucket(tFinal/2) -> 1
    int timeBucket(double t);
    int timeBucket(std::vector<double>* tkc);

    // vector of column matrices that contain the polynomials that represent the
    // displacement of each mass for each time bucket. Ex: models[0][0] =
    // Polynomial that represents the displacement (x) of mass 0, as a function
    // of time and the values of the springs and the masses for the first time
    // bucket.
    std::vector<boost::numeric::ublas::matrix<Poly>> models;
    std::vector<std::vector<std::vector<double>>> modelsCoefficients;

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

    void AddInitialConditionsResiduesTkc(std::vector<double>* tkc,
                                         int tkcIndex);
    void AddPhysicsResiduesTkc(std::vector<double>* tkc, int tkcIndex);
    void AddResiduesTkc();

    void AddInitialConditionsResidues();
    void AddPhysicsResidues();

    // Auxiliary functions that return the weight of each type of loss.
    // Must be called after Add*Residue methods.
    double InitialConditionsWeight();
    double PhysicsWeight();

    double Loss() override;
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
