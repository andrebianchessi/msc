#include <vector>

#include "maybe.h"

struct StatusAndValue {
    bool success;
    double value;
};

class Model {
    // Base class that represent an ML model
   public:
    // Returns the number of parameters this model has
    virtual int nParameters() = 0;

    // Train the model, i.e. set the parameters to the values that minimize
    // the Loss function. This is performed using gradient descent.
    // The learning rate is automatically reduced by half when needed, until
    // it reaches a value lower than minLearningRate. Set log to true to log
    // the learning. Function returns the latest loss function value.
    Maybe<double> Train(double learningRate, double minLearningRate, bool log);

    Maybe<double> Train(double learningRate, int maxSteps, bool log);

   private:
    // After calling this method, the model's parameters will be set in the
    // target vector
    virtual Maybe<Void> GetParameters(std::vector<double>* target) = 0;

    // After calling this method, the model's parameters will be set with the
    // values from the input vector
    virtual Maybe<Void> SetParameters(std::vector<double>* parameters) = 0;

    // Returns the loss of the model given it's current parameters
    virtual double Loss() = 0;

    // Returns the gradient of the loss function with respect to it's parameters
    // i.e. if the model has 2 parameters, the returned vector will be the
    // derivative of the loss function with respect to the first parameter
    // followed by the derivative of the loss function with respect to the
    // second parameter
    virtual std::vector<double> LossGradient() = 0;

    // Performs a step in gradient descent. If the Loss function decreases,
    // returns true and its value. Else, resets the parameters to what they
    // previously were, and returns false and the loss function value before the
    // step. The step parameter is what multiplies the gradient: Parameters_new
    // = Parameters_old - step*grad
    StatusAndValue GradientDescentStep(double currentLoss, double step);
};