#include <memory>
#include <vector>

#include "bounded.h"
#include "evolution.h"
#include "problem_creature.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    int massId = 5;
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 2.0, 2.0);  // m2
    pd.AddMass(150, 3.0, 3.0);  // m3
    pd.AddMass(700, 4.0, 4.0);  // m4
    pd.AddMass(90, 5.0, 5.0);   // m5

    double min = 100000.0;
    double max = 200000.0;
    for (int i = 0; i < massId; i++) {
        for (int j = i + 1; j <= massId; j++) {
            pd.AddSpring(i, j, min, max);
            pd.AddDamper(i, j, min / 2, max / 2);
        }
    }

    pd.SetFixedMass(0);
    pd.AddInitialVel(200.0);
    assert(pd.IsOk());

    // Common parameters
    double finalT = 0.1;
    int popSize = 200;
    double geneticAlgoErrorStop = 1.0 / 100.0;

    // Pimodel based optimization
    int nModels = 5;
    int icPoints = 10;
    int physPoints = 10;
    int order = 3;
    double learningRate = 0.005;
    double earlyStopLoss = 0.1;
    int maxSteps = 20000;
    bool logComplexity = false;
    bool logTraining = true;
    int timeDiscretization = 20;  // used to look for max accel

    // Train models
    Pimodels models =
        Pimodels(pd, finalT / 10, nModels, icPoints, physPoints, order);
    auto start = Now();
    assert(!models
                .Train(learningRate, learningRate / 1000, earlyStopLoss,
                       maxSteps, logComplexity, logTraining)
                .isError);
    std::cout << "Time to train models: " << TimeSince(start) << std::endl;

    // Set all the initial populations to have the max value for springs and
    // dampers. This is a bad solution. We're doing this so that "being lucky"
    // in our initial guess doesn't affect the quality of the optimization
    // algorithm.
    std::vector<ProblemCreature> piPopulation = std::vector<ProblemCreature>();
    std::vector<ProblemCreature> integrationPopulation =
        std::vector<ProblemCreature>();
    for (int i = 0; i < popSize; i++) {
        piPopulation.push_back(
            ProblemCreature(&pd, massId, &models, timeDiscretization));
        integrationPopulation.push_back(ProblemCreature(&pd, massId, finalT));
    }
    for (int i = 0; i < popSize; i++) {
        for (int j = 0; j < integrationPopulation[i].dna.size(); j++) {
            piPopulation[i].dna[j].Set(1.0);
            integrationPopulation[i].dna[j].Set(1.0);
        }
    }

    // Sort once using explicit integration to find out what the initial guess
    // best solution is, so that we can compare;
    Evolution<ProblemCreature> evolution =
        Evolution<ProblemCreature>(&integrationPopulation);
    evolution.SortPopulation();
    auto initialBadGuessDna = evolution.GetCreature(0)->dna;
    Problem initialBadGuess = pd.BuildFromDNA(initialBadGuessDna).val;

    // Pimodel-based optimization
    evolution = Evolution<ProblemCreature>(&piPopulation);
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, true);
    std::cout << "Time to G.A. optimization using Pimodels: "
              << TimeSince(start) << std::endl;
    auto pimodelBestDna = evolution.GetCreature(0)->dna;
    Problem pimodelBest = pd.BuildFromDNA(pimodelBestDna).val;

    // Explicit-integration based optimization
    evolution = Evolution<ProblemCreature>(&integrationPopulation);
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, true);
    std::cout << "Time to G.A. optimization using Explicit integration: "
              << TimeSince(start) << std::endl;
    auto explicitBestDna = evolution.GetCreature(0)->dna;
    Problem explicitBest = pd.BuildFromDNA(explicitBestDna).val;

    std::cout << "Initial bad guess solution:" << std::endl;
    for (int i = 0; i < initialBadGuessDna.size(); i++) {
        std::cout << initialBadGuessDna[i].Get() << std::endl;
    }
    std::cout << "Pimodel-based best solution:" << std::endl;
    for (int i = 0; i < pimodelBestDna.size(); i++) {
        std::cout << pimodelBestDna[i].Get() << std::endl;
    }
    std::cout << "Explicit-based best solution:" << std::endl;
    for (int i = 0; i < explicitBestDna.size(); i++) {
        std::cout << explicitBestDna[i].Get() << std::endl;
    }

    assert(!initialBadGuess.Integrate(finalT).isError);
    assert(!pimodelBest.Integrate(finalT).isError);
    assert(!explicitBest.Integrate(finalT).isError);
    std::cout << "Initial bad guess: " << std::endl;
    initialBadGuess.PrintMassTimeHistory(massId);
    std::cout << "Pimodel-based best: " << std::endl;
    pimodelBest.PrintMassTimeHistory(massId);
    std::cout << "Explicit-based best: " << std::endl;
    explicitBest.PrintMassTimeHistory(massId);
}