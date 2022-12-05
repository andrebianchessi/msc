#include "model.h"

#include <cassert>
#include <iostream>
#include <vector>

#include "maybe.h"

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

    // If the loss function decreased, return true
    // else, reset the parameters and return false
    if (loss_1 < loss_0) {
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
    while (learningRate >= minLearningRate && learningRate != 0) {
        stepStatus = this->GradientDescentStep(learningRate);
        if (stepStatus.success) {
            if (log) {
                std::cout << "Learning Rate: " << learningRate
                          << " Loss: " << stepStatus.value << std::endl;
            }
        } else {
            learningRate = learningRate / 2;
        }
    }
    r.val = stepStatus.value;
    return r;
}