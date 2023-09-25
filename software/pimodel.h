#pragma once
#include <gtest/gtest.h>

#include "bounded.h"
#include "maybe.h"
#include "model.h"
#include "polynomial.h"
#include "problem_description.h"

class Pimodel : public Model {
   public:
    Pimodel() = default;
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
    Pimodel(ProblemDescription p, double initialT, double finalT,
            int timeDiscretization, int kcDiscretization, int order);

    // Must be called before Train()
    // Builds the residues based on problem description.
    // The arguments define which residues will be considered.
    void SetResidues(bool initialConditions, bool physics);

    // Returns the position of each mass.
    // The input should be an array with:
    // The global time
    // The values of the springs
    // The values of the dampers
    // U = [T, K1, K2, ... , Kn, C1, C2, ... , Cn]
    // All inputs are NOT normalized (i.e. they don't have to be in [0, 1]).
    // E.g. of valid input: U = [100.0, 99.0, ...]
    Maybe<std::vector<double>> operator()(std::vector<double>* TKC);
    // Similar to operator(), but returns the velocity of each mass
    Maybe<std::vector<double>> GetVelocities(std::vector<double>* TKC);

    int nParameters() override;

   private:
    friend class PimodelTest;
    friend class Pimodels;
    friend class PimodelTrainTest;
    FRIEND_TEST(PimodelTest, ModelIdTest);
    FRIEND_TEST(PimodelTest, ConstructorTest);
    FRIEND_TEST(PimodelTest, OperatorTest);
    FRIEND_TEST(PimodelTest, TimeBucketTest);
    FRIEND_TEST(PimodelTest, GetParametersTest);
    FRIEND_TEST(PimodelTest, SetParametersTest);
    FRIEND_TEST(PimodelTest, ProblemFromTkcTest);
    FRIEND_TEST(PimodelTest, getXModelTest);
    FRIEND_TEST(PimodelTest, getAccelsFromDiffEqTest);
    FRIEND_TEST(PimodelTest, getInitialXTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTkcTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTkcZeroDiscTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesTkcTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesTkcZeroDiscTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesTest);
    FRIEND_TEST(PimodelTest, InitialConditionsResiduesNumericalTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesTest);
    FRIEND_TEST(PimodelTest, PhysicsResiduesNumericalTest);
    FRIEND_TEST(PimodelTest, nResiduesTest);
    FRIEND_TEST(PimodelTest, LossGradientTest);
    FRIEND_TEST(PimodelsTest, ConstructorTest);
    FRIEND_TEST(PimodelsTest, setContinuityTest);
    FRIEND_TEST(PimodelsTest, OperatorTest);

    ProblemDescription p;

    int nMasses;

    std::vector<double> normalizeTkc(std::vector<double>* TKC);

    // Vector of column matrices that contain the polynomials that represent the
    // displacement of each mass for each time bucket. Ex, for a system of two
    // masses: models[0] = Polynomial that represents the displacement (x) of
    // mass 0, as a function of time and the values of the springs and the
    // masses.
    boost::numeric::ublas::matrix<Poly> models;
    std::vector<std::vector<double>> modelsCoefficients;
    // Derivatives of each model with respect to t (normalized time)
    std::vector<Poly> modelsD;
    // Second derivatives of each model with respect to t (normalized time)
    std::vector<Poly> modelsDD;

    // Residues of initial displacement
    std::vector<Polys> initialDispResidues;
    // Residues of initial velocity
    std::vector<Polys> initialVelResidues;
    // Values of time, spring coefficients and damper coefficients in which
    // the initial displacement and velocity residues must be calculated to make
    // up the loss function
    std::vector<std::vector<Bounded>> initialConditionsResiduesTkc;

    // Physical residue
    std::vector<Polys> physicsResidues;
    std::vector<std::vector<Bounded>> physicsResiduesTkc;

    // Caches of residues evaluations such as
    // initialDispResidues[i](modelsCoefficients)
    Polys* residueById(int i);
    double residueWeight(int i);
    std::vector<bool> residueIsCached;
    std::vector<double> residueCache;
    void setResidueCache(int i);

    // Gradients of each residue
    // i.e. what we get from initialDispResidues[i].Da()
    std::vector<std::map<std::tuple<int, int>, double>> initialDispResiduesDa;
    std::vector<std::map<std::tuple<int, int>, double>> initialVelResiduesDa;
    std::vector<std::map<std::tuple<int, int>, double>> physicsResiduesDa;
    void setResiduesDa();
    std::map<std::tuple<int, int>, double>* residueDaById(int i);

    // This model describes the system from t = t0 to t = t1
    double t0;
    double t1;
    // Derivative of t (normalized time from [0,1]) w.r.t. T (global time)
    double dtdT;
    // Derivative of T (global time) w.r.t. t (normalized time from [0,1])
    double dTdt;

    // Parameter that describes in how many intervals we'll discretize time and
    // the springs/dampers values in our loss function
    int timeDiscretization;
    int kcDiscretization;

    // Order of the polynomial models
    int order;

    // Returns the size of the input of this model.
    // This will always be the number of springs of the problem
    // plus the number of dampers plus 1 (time).
    int inputSize();

    Maybe<Void> GetParameters(std::vector<double>* target) override;

    Maybe<Void> SetParameters(std::vector<double>* parameters) override;

    void AddInitialConditionsResiduesTkc(std::vector<Bounded>* tkc,
                                         int tkcIndex);
    void AddPhysicsResiduesTkc(std::vector<Bounded>* tkc, int tkcIndex);
    void AddResiduesTkc();

    void AddInitialConditionsResidues();
    void AddPhysicsResidues();

    // Auxiliary functions that return the weight of each type of loss.
    // Must be called after Add*Residue methods.
    double InitialConditionsWeight();
    double PhysicsWeight();

    double Residue(int i) override;
    // AddResidues must be called first for this to work properly
    int nResidues() override;
    std::vector<double> LossGradient(int i) override;

    // Auxiliary functions:

    // Create problem for given tkc (time, springs and dampers)
    Problem problemFromTkc(std::vector<Bounded>* tkc);

    // Return the state vector with the initial conditions set
    std::vector<double> getInitialX();

    // Gets the polynomials that represent the acceleration of each mass
    // by using the differential equation from the discrete element method
    boost::numeric::ublas::matrix<Polys> getAccelsFromDiffEq(
        Problem* problem, std::vector<Bounded>& tkc);
};

class Pimodels {
   public:
    Pimodels() = default;
    Pimodels(ProblemDescription p, double finalT, int nModels,
             int timeDiscretization, int kcDiscretization, int order);

    // initialConditionsLearningRate is the learning rate used to train only the
    // initial conditions, whereas physicsLearningRate is the learning rate used
    // when initial Conditions and physics is considered. In my experience,
    // initialConditionsLR can be much larger (100x) than physicsLR.
    Maybe<double> Train(double initialConditionsLearningRate,
                        double physicsLearningRate, int maxSteps,
                        bool logComplexity, bool logTraining);

    Maybe<std::vector<double>> operator()(std::vector<double>* TKC);
    Maybe<std::vector<double>> GetVelocities(std::vector<double>* TKC);

   private:
    int getTimeBucket(double t) const;
    std::vector<double> timeBuckets;
    std::vector<Pimodel> pimodels;

    // We train the models successively in time: first the model that is valid
    // for the first time interval, then the model valid for the second time
    // interval and etc.
    // We impose continuity between the models at the transition times. However,
    // the models are not only a function of time, but also of the values of the
    // springs and dampers. The question is: what values do we user for the
    // springs and dampers when imposing continuity? This function returns the
    // tkc vector (time followed by spring values followed by damper values)
    // which we use to impose continuity. The k's and c's are set to the average
    // value they can have (based on the problem description). t is set to zero,
    // and must be updated for each point in which we want to impose continuity.
    std::vector<double> continuityTkc() const;

    void setContinuity(int timeBucket, std::vector<double>& TKC);

    void logComplexity();

    FRIEND_TEST(PimodelsTest, ConstructorTest);
    FRIEND_TEST(PimodelsTest, GetTimeBucketTest);
    FRIEND_TEST(PimodelsTest, getContinuityTkcTest);
    FRIEND_TEST(PimodelsTest, setContinuityTest);
    FRIEND_TEST(PimodelsTest, OperatorTest);
};