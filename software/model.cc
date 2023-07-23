#include "model.h"

#include <cassert>
#include <iostream>
#include <sstream>
#include <unordered_set>
#include <vector>

#include "maybe.h"
#include "utils.h"

const double MIN_STEP_IMPROVEMENT =
    0.0;  // min relative improvement to early stop

void Model::GradientDescentStep(int i, double stepSize) {
    std::vector<double> oldParameters =
        std::vector<double>(this->nParameters());
    this->GetParameters(&oldParameters);

    std::vector<double> newParameters =
        std::vector<double>(this->nParameters());

    std::vector<double> grad = this->ResidueGradient(i);
    assert(grad.size() == oldParameters.size());
    for (int i = 0; i < int(oldParameters.size()); i++) {
        newParameters[i] = oldParameters[i] - grad[i] * stepSize;
    }

    this->SetParameters(&newParameters);
}

double Model::Loss() {
    double l = 0;
    for (int i = 0; i < this->nResidues(); i++) {
        l += this->Residue(i);
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
    auto precision = std::cout.precision();

    batchSize = std::min(batchSize, this->nResidues());

    double lossBeforeStep = this->Loss();
    if (log) {
        std::cout.precision(20);
        std::cout << "Initial Loss: " << lossBeforeStep << std::endl;
        ;
    }
    double lossAfterStep;
    bool stepOk;

    std::vector<double> parametersBeforeStep =
        std::vector<double>(this->nParameters());
    int step = 0;
    while (step < maxSteps && learningRate != 0 && lossBeforeStep > 0) {
        this->GetParameters(&parametersBeforeStep);

        std::unordered_set<int> batch;
        int tries = 0;
        const int maxTriesToCreateRandomBatch = 10;
        while (batch.size() < batchSize &&
               tries < maxTriesToCreateRandomBatch) {
            batch.insert(RandomInt(0, this->nResidues() - 1));
            tries += 1;
        }
        int i = 0;
        while (batch.size() < batchSize) {
            batch.insert(i);
            i++;
        }
        for (int residueIndex : batch) {
            this->GradientDescentStep(residueIndex, learningRate);
        }

        lossAfterStep = this->Loss();
        stepOk = lossAfterStep < lossBeforeStep;
        if (!stepOk) {
            // If the step did not reduce the loss, we reset to the initial
            // conditions and reduce the learning rate
            lossAfterStep = lossBeforeStep;
            this->SetParameters(&parametersBeforeStep);
            learningRate /= 2;
        } else {
            // If it was ok, we update the "before" values to prepare for the
            // next step
            lossBeforeStep = lossAfterStep;
        }
        if (log) {
            std::cout.precision(20);
            std::cout << "Loss: " << lossAfterStep;
            std::cout.precision(3);
            std::cout << " Learning Rate: " << learningRate;
            stepOk ? std::cout << " -> Step Ok"
                   : std::cout << " -> Step Not Ok ";
            std::cout << std::endl;
        }
        step++;
    }
    std::cout.precision(precision);
    r.val = lossAfterStep;
    return r;
}
