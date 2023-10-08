#include <memory>
#include <vector>

#include "bounded.h"
#include "evolution.h"
#include "problem_creature.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    int massId = 4;
    ProblemDescription pd = ProblemDescription();

    pd.AddMass(1.0, 0.0, 0.0);  // m0
    for (int i = 1; i <= massId; i++) {
        pd.AddMass(Random(100, 300), i, 0);
    }

    double min = 100000.0;
    double max = 200000.0;
    for (int i = 0; i < massId; i++) {
        for (int j = i + 1; j <= massId; j++) {
            pd.AddSpring(i, j, min, max);
            pd.AddDamper(i, j, min, max);
        }
    }

    pd.SetFixedMass(0);
    for (int i = 1; i <= massId; i++) {
        pd.AddInitialVel(i, Random(0.0, 200.0));
        pd.AddInitialDisp(i, Random(0.0, 200.0));
    }
    assert(pd.IsOk());

    // Common parameters
    double finalT = 0.01;
    int popSize = 50;
    double geneticAlgoErrorStop = 0.001 / 100.0;  // 0.001%

    // Pimodel based optimization
    int nModels = 18;
    int icPoints = 10;
    int physPoints = 10;
    int order = 4;
    double learningRate = 0.01;
    double earlyStopLoss = 0.1;
    int maxSteps = 10000;
    bool logComplexity = false;
    bool logTraining = true;
    int timeDiscretization = 5;  // used to look for max accel

    // Train models
    Pimodels models =
        Pimodels(pd, finalT, nModels, icPoints, physPoints, order);
    auto start = Now();
    assert(!models
                .Train(learningRate, learningRate / 100, earlyStopLoss,
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
    std::vector<ProblemCreature> randomPopulation =
        std::vector<ProblemCreature>();
    for (int i = 0; i < popSize; i++) {
        piPopulation.push_back(
            ProblemCreature(&pd, massId, &models, timeDiscretization));
        integrationPopulation.push_back(ProblemCreature(&pd, massId, finalT));
        randomPopulation.push_back(ProblemCreature(&pd, massId, finalT));
    }
    for (int i = 0; i < popSize; i++) {
        for (int j = 0; j < integrationPopulation[i].dna.size(); j++) {
            integrationPopulation[i].dna[j].Set(piPopulation[i].dna[j].Get());
            randomPopulation[i].dna[j].Set(piPopulation[i].dna[j].Get());
        }
    }

    // Sort once using explicit integration to find out what the initial guess
    // best solution is, so that we can compare;
    Evolution<ProblemCreature> evolution =
        Evolution<ProblemCreature>(&randomPopulation);
    evolution.SortPopulation();
    auto randomGuessDna = evolution.GetCreature(0)->dna;
    Problem randomGuess = pd.BuildFromDNA(randomGuessDna).val;

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

    std::cout << "Random guess best solution:" << std::endl;
    for (int i = 0; i < randomGuessDna.size(); i++) {
        std::cout << randomGuessDna[i].Get() << std::endl;
    }
    std::cout << "Pimodel-based best solution:" << std::endl;
    for (int i = 0; i < pimodelBestDna.size(); i++) {
        std::cout << pimodelBestDna[i].Get() << std::endl;
    }
    std::cout << "Explicit-based best solution:" << std::endl;
    for (int i = 0; i < explicitBestDna.size(); i++) {
        std::cout << explicitBestDna[i].Get() << std::endl;
    }

    assert(!randomGuess.Integrate(finalT).isError);
    assert(!pimodelBest.Integrate(finalT).isError);
    assert(!explicitBest.Integrate(finalT).isError);
    std::cout << "Random guess best: " << std::endl;
    randomGuess.PrintMassTimeHistory(massId);
    std::cout << "Pimodel-based best: " << std::endl;
    pimodelBest.PrintMassTimeHistory(massId);
    std::cout << "Explicit-based best: " << std::endl;
    explicitBest.PrintMassTimeHistory(massId);
}