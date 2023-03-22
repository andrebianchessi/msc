#include "model.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

#include "maybe.h"

const double MIN_STEP_IMPROVEMENT =
    0.0;  // min relative improvement to early stop

StatusAndValue Model::GradientDescentStep(double currentLoss, double step) {
    StatusAndValue status;

    std::vector<double> parameters_0 = std::vector<double>(this->nParameters());
    std::vector<double> parameters_1 = std::vector<double>(this->nParameters());
    this->GetParameters(&parameters_0);

    std::vector<double> grad = this->LossGradient();
    assert(grad.size() == parameters_0.size());

    for (int i = 0; i < int(parameters_0.size()); i++) {
        parameters_1[i] = parameters_0[i] - grad[i] * step;
    }

    this->SetParameters(&parameters_1);
    double newLoss = this->Loss();

    // If the loss function decreased more than MIN_STEP_IMPROVEMENT, return
    // true. Else, reset the parameters and return false
    if (newLoss < currentLoss &&
        currentLoss - newLoss > MIN_STEP_IMPROVEMENT * currentLoss) {
        status.success = true;
        status.value = newLoss;
        return status;
    }
    this->SetParameters(&parameters_0);
    status.success = false;
    status.value = currentLoss;
    return status;
}

Maybe<double> Model::Train(double learningRate, double minLearningRate,
                           bool log) {
    Maybe<double> r;
    if (learningRate <= 0) {
        r.isError = true;
        r.errMsg = "LearningRate must be >0";
        return r;
    }
    StatusAndValue stepStatus;
    auto precision = std::cout.precision();
    double loss = this->Loss();
    while (learningRate >= minLearningRate && learningRate != 0) {
        stepStatus = this->GradientDescentStep(loss, learningRate);
        if (log) {
            std::cout.precision(3);
            std::cout << "Learning Rate: " << learningRate;
            std::cout.precision(20);
            std::cout << " Loss: " << stepStatus.value;
            stepStatus.success ? std::cout << " -> Step Ok"
                               : std::cout << " -> Step Not Ok ";
            std::cout << std::endl;
        }
        if (!stepStatus.success) {
            learningRate = learningRate / 2;
        } else {
            loss = stepStatus.value;
        }
    }
    std::cout.precision(precision);
    r.val = stepStatus.value;
    return r;
}

Maybe<double> Model::Train(double learningRate, int maxSteps, bool log) {
    Maybe<double> r;
    if (learningRate <= 0) {
        r.isError = true;
        r.errMsg = "LearningRate must be >0";
        return r;
    }
    int step = 0;
    StatusAndValue stepStatus;
    double loss = this->Loss();
    auto precision = std::cout.precision();
    while (step < maxSteps && learningRate != 0) {
        stepStatus = this->GradientDescentStep(loss, learningRate);
        if (log) {
            std::cout.precision(3);
            std::cout << "Learning Rate: " << learningRate;
            std::cout.precision(20);
            std::cout << " Loss: " << stepStatus.value;
            stepStatus.success ? std::cout << " -> Step Ok"
                               : std::cout << " -> Step Not Ok ";
            std::cout << std::endl;
        }
        if (!stepStatus.success) {
            learningRate = learningRate / 2;
        } else {
            loss = stepStatus.value;
        }
        step++;
    }
    std::cout.precision(precision);
    r.val = stepStatus.value;
    return r;
}
