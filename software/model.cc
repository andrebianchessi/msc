#include "model.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <vector>

#include "maybe.h"

const double MIN_STEP_IMPROVEMENT =
    0.0;  // min relative improvement to early stop

StatusAndValue Model::GradientDescentStep(double step) {
    StatusAndValue status;

    double loss_0 = this->Loss();

    std::vector<double> parameters_0 = std::vector<double>(this->nParameters());
    std::vector<double> parameters_1 = std::vector<double>(this->nParameters());
    this->GetParameters(&parameters_0);

    std::vector<double> grad = this->LossGradient();
    assert(grad.size() == parameters_0.size());

    for (int i = 0; i < int(parameters_0.size()); i++) {
        parameters_1[i] = parameters_0[i] - grad[i] * step;
    }

    this->SetParameters(&parameters_1);
    double loss_1 = this->Loss();

    // If the loss function decreased more than MIN_STEP_IMPROVEMENT, return
    // true. Else, reset the parameters and return false
    if (loss_1 < loss_0 && loss_0 - loss_1 > MIN_STEP_IMPROVEMENT * loss_0) {
        status.success = true;
        status.value = loss_1;
        return status;
    }
    this->SetParameters(&parameters_0);
    status.success = false;
    status.value = loss_0;
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
    while (learningRate >= minLearningRate && learningRate != 0) {
        stepStatus = this->GradientDescentStep(learningRate);
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
        }
    }
    std::cout.precision(precision);
    r.val = stepStatus.value;
    return r;
}