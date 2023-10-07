#include <memory>
#include <vector>

#include "bounded.h"
#include "evolution.h"
#include "problem_creature.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 2.0, 2.0);  // m2
    pd.AddMass(150, 3.0, 3.0);  // m3
    pd.AddMass(700, 4.0, 4.0);  // m4
    pd.AddMass(90, 5.0, 5.0);   // m5

    double min = 10000.0;
    double max = 100000.0;
    for (int i = 0; i < 5; i++) {
        for (int j = i + 1; j <= 5; j++) {
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
    int massId = 5;

    // Pimodel based optimization
    int nModels = 5;
    int icPoints = 5;
    int physPoints = 5;
    int order = 4;
    double learningRate = 0.01;
    int maxSteps = 2000;
    bool logComplexity = false;
    bool logTraining = true;
    int timeDiscretization = 20;  // used to look for max accel

    // Train models
    Pimodels models =
        Pimodels(pd, finalT, nModels, icPoints, physPoints, order);
    auto start = Now();
    assert(!models
                .Train(learningRate, learningRate / 50000, maxSteps,
                       logComplexity, logTraining)
                .isError);
    std::cout << "Time to train models: " << TimeSince(start) << std::endl;

    // Create initial populations with the same dna values, so that we can
    // compare the results of each optimization.
    std::vector<ProblemCreature> piPopulation = std::vector<ProblemCreature>();
    std::vector<ProblemCreature> integrationPopulation =
        std::vector<ProblemCreature>();
    for (int i = 0; i < popSize; i++) {
        piPopulation.push_back(
            ProblemCreature(&pd, massId, &models, timeDiscretization));
        integrationPopulation.push_back(ProblemCreature(&pd, massId, finalT));
        for (int j = 0; j < integrationPopulation[i].dna.size(); j++) {
            integrationPopulation[i].dna[j].Set(piPopulation[i].dna[j].Get());
        }
    }

    // Pimodel-based optimization
    Evolution<ProblemCreature> evolution =
        Evolution<ProblemCreature>(&piPopulation);
    Problem nonOptimized = pd.BuildFromDNA(evolution.GetCreature(0)->dna).val;
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, true);
    std::cout << "Time to G.A. optimization using Pimodels: "
              << TimeSince(start) << std::endl;
    Problem pimodelBest = pd.BuildFromDNA(evolution.GetCreature(0)->dna).val;

    // Explicit-integration based optimization
    evolution = Evolution<ProblemCreature>(&integrationPopulation);
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, true);
    std::cout << "Time to G.A. optimization using Explicit integration: "
              << TimeSince(start) << std::endl;
    Problem explicitBest = pd.BuildFromDNA(evolution.GetCreature(0)->dna).val;

    assert(!nonOptimized.Integrate(finalT).isError);
    assert(!pimodelBest.Integrate(finalT).isError);
    assert(!explicitBest.Integrate(finalT).isError);
    std::cout << "Random solution: " << std::endl;
    nonOptimized.PrintMassTimeHistory(massId);
    std::cout << "Pimodel-based best: " << std::endl;
    pimodelBest.PrintMassTimeHistory(massId);
    std::cout << "Explicit-based best: " << std::endl;
    explicitBest.PrintMassTimeHistory(massId);
}