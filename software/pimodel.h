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
    // finalT discretization: How many intervals we'll discretize time and the
    // springs/dampers values in our loss function oder: order of the
    // multivariate polynomial used in the model
    Pimodel(ProblemDescription* p, double finalT, int discretization,
            int order);

    // Returns the position of each mass.
    // The input should be an array with:
    // The time
    // The values of the springs
    // The values of the dampers
    // U = [t, k1, k2, ... , kn, c1, c2, ... , cn]
    Maybe<std::vector<double>> operator()(std::vector<double>* tkc);

    int nParameters() override;

    // Auxiliary function used inside Loss()
    void LossDfs(std::vector<double>* tkc, int tkcIndex, double* loss);

   private:
    FRIEND_TEST(PimodelTest, ConstructorTest);
    FRIEND_TEST(PimodelTest, GetParametersTest);
    FRIEND_TEST(PimodelTest, SetParametersTest);
    int massIndex;

    ProblemDescription* p;
    std::vector<Poly> polys;

    // This model describes the system from t = 0 to t = finalT
    double finalT;

    // Parameter that describes in how many intervals we'll discretize time and
    // the springs/dampers values in our loss function
    int discretization;

    // Returns the size of the input of this model.
    // This will always be the number of springs of the problem
    // plus the number of dampers plus 1 (time).
    int inputSize();

    Maybe<Void> GetParameters(std::vector<double>* target) override;

    Maybe<Void> SetParameters(std::vector<double>* parameters) override;

    double Loss() override;

    std::vector<double> LossGradient() override;
};