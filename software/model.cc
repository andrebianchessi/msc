#include "model.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <iostream>
#include <random>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "maybe.h"
#include "utils.h"

auto rng = std::default_random_engine{};

const double MIN_STEP_IMPROVEMENT =
    0.0;  // min relative improvement to early stop

void Model::StochasticGradientDescentStep(int i, double stepSize) {
    std::vector<double> oldParameters =
        std::vector<double>(this->nParameters());
    this->GetParameters(&oldParameters);

    std::vector<double> newParameters =
        std::vector<double>(this->nParameters());

    std::vector<double> grad = this->LossGradient(i);
    assert(grad.size() == oldParameters.size());
    for (int i = 0; i < int(oldParameters.size()); i++) {
        newParameters[i] = oldParameters[i] - grad[i] * stepSize;
    }

    this->SetParameters(&newParameters);
}

double Model::Loss() {
    double l = 0;
    for (int i = 0; i < this->nResidues(); i++) {
        l += pow(this->Residue(i), 2);
    }
    return l;
}

Maybe<double> Model::Train(double learningRate, int batchSize, int maxSteps,
                           bool log) {
    Maybe<double> r;
    if (learningRate <= 0) {
        r.isError = true;
        r.errMsg = "LearningRate must be >0";
        return r;
    }

    std::vector<double> parametersBeforeStep =
        std::vector<double>(this->nParameters());
    int step = 0;
    while (step < maxSteps) {
        this->GetParameters(&parametersBeforeStep);

        this->StochasticGradientDescentStep(RandomInt(0, this->nResidues() - 1),
                                            learningRate);
        step++;
        // Compute loss only every now and then to improve efficiency because
        // computing the whole loss can be expensive.
        if (step % 50 == 0) {
            if (log) {
                std::cout << "Loss: " << this->Loss() << std::endl;
            }
            if (this->Loss() <= 0.000001) {
                break;
            }
        }
    }
    r.val = this->Loss();
    return r;
}
