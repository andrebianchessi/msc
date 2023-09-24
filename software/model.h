#include <gtest/gtest.h>

#include <vector>

#include "maybe.h"

class Model {
    // Base class that represent an ML model
   public:
    // Returns the number of parameters this model has
    virtual int nParameters() = 0;

    // Train the model, i.e. set the parameters to the values that minimize
    // the Loss function. This is performed using stochastic gradient descent.
    // The learning rate is automatically reduced by half when needed, until
    // it reaches a value lower than minLearningRate. Set log to true to log
    // the learning. Function returns the latest loss function value.
    // Note: its ok if batchSize>nResidues (the method then sets it to
    // nResidues)
    Maybe<double> Train(double learningRate, int batchSize, int maxSteps,
                        bool log);

    // Returns the total loss, i.e. the sum of the squared residues
    double Loss();

   private:
    // After calling this method, the model's parameters will be set in the
    // target vector
    virtual Maybe<Void> GetParameters(std::vector<double>* target) = 0;

    // After calling this method, the model's parameters will be set with the
    // values from the input vector
    virtual Maybe<Void> SetParameters(std::vector<double>* parameters) = 0;

    // Returns the value of the i-th residue given the model's current
    // parameters.
    // Ex: loss = ((model(x0)-y0)^2 + (model(x1)-y1)^2))
    // Residue(0) -> model(x0)-y0
    // Residue(1) -> model(x1)-y1
    virtual double Residue(int i) = 0;

    // Returns the number of residues in the loss
    // Ex: loss = (model(x0)-y0)^2 + (model(x1) - y1)^2 -> 2 residues
    virtual int nResidues() = 0;

    // Returns the gradient of the i-th loss summand with respect to the model's
    // parameters i.e. if the model has 2 parameters, the returned vector will
    // be the derivative of the i-th squared residue with respect to the first
    // parameter followed by the derivative of the i-th
    // squared residue with respect to the second parameter. It's fine if
    // multiplicative constants are removed. I.e. returning [1,5] if the formal
    // mathematical gradient is [2,10] is absolutely fine, because we scale
    // the gradient anyways during training.
    virtual std::vector<double> LossGradient(int i) = 0;

    // Performs a step in gradient descent considering the i-th residue. If the
    // value of the residue decreases, returns true and its value. Else, resets
    // the parameters to what they previously were, and returns false and the
    // residue's value before the step. The step parameter is what
    // multiplies the gradient: Parameters_new = Parameters_old - step*grad
    void StochasticGradientDescentStep(int i, double stepSize);

    FRIEND_TEST(ModelTest, LossTest);
};