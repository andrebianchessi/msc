// Optimizes problems with increasing complexity using the standard Genetic Algo
// (which evaluates solutions by explicitly integrating them numerically) and
// the Pimodel-Genetic-Algorithm (which evaluates solutions by using Physics
// Informed models which are trained before the genetic optimization). Prints
// out scores to compare both methods.

#include <memory>
#include <vector>

#include "bounded.h"
#include "evolution.h"
#include "problem_creature.h"
#include "utils.h"

struct Results {
    double bestRandomGuessMaxAccel;

    int timeUsecToTrainPimodels;
    int timeUsecToPimodelGa;
    double pimodelGaBestSolutionMaxAccel;

    int timeUsecToExplicitIntegrationGa;
    double explicitIntegrationBestSolutionMaxAccel;
};

void printResults(Results r, int nMasses, int id) {
    double lowestCost = 0;
    double highestCost =
        std::max(r.timeUsecToExplicitIntegrationGa,
                 r.timeUsecToTrainPimodels + r.timeUsecToPimodelGa);
    double lowestAccel = std::min(
        std::min(r.bestRandomGuessMaxAccel, r.pimodelGaBestSolutionMaxAccel),
        r.explicitIntegrationBestSolutionMaxAccel);
    double highestAccel = std::max(
        std::max(r.bestRandomGuessMaxAccel, r.pimodelGaBestSolutionMaxAccel),
        r.explicitIntegrationBestSolutionMaxAccel);

    // lowest accel -> accelScore = 100
    // highest accel -> accelScore = 0
    auto accelScore = [=](double x) {
        double a = -1 / (highestAccel - lowestAccel);
        double b = -a * highestAccel;
        return 100 * (a * x + b);
    };
    // lowest cost -> costScore = 100
    // highest cost -> costScore = 0
    auto costScore = [=](double x) {
        double a = -1 / (highestCost - lowestCost);
        double b = -a * highestCost;
        return 100 * (a * x + b);
    };

    double randomGuessScore = accelScore(r.bestRandomGuessMaxAccel);
    double randomGuessCost = 100;

    double pimodelScore = accelScore(r.pimodelGaBestSolutionMaxAccel);
    double pimodelCost =
        costScore(r.timeUsecToTrainPimodels + r.timeUsecToPimodelGa);
    double explicitScore =
        accelScore(r.explicitIntegrationBestSolutionMaxAccel);
    double explicitCost = costScore(r.timeUsecToExplicitIntegrationGa);

    std::cout << "nMasses,id,method,efficiency score(0 is bad and 100 is "
                 "best),quality "
                 "score(0 is bad "
                 "and 100 is best),mean score\n ";
    std::cout << nMasses << "," << id << ",Pimodel-based GA," << pimodelCost
              << "," << pimodelScore << "," << (pimodelCost + pimodelScore) / 2
              << "\n";
    std::cout << nMasses << "," << id << ",ExplicitIntegration-based GA,"
              << explicitCost << "," << explicitScore << ","
              << (explicitCost + explicitScore) / 2 << "\n";
    std::cout << nMasses << "," << id << ",Random Guess," << randomGuessCost
              << "," << randomGuessScore << ","
              << (randomGuessCost + randomGuessScore) / 2 << "\n";
}

Results optimize(int nMasses) {
    bool verbose = true;

    Results res;

    // Setup problem
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    for (int i = 1; i <= nMasses; i++) {
        pd.AddMass(Random(100, 300), i, 0);
    }
    double min = 100000.0;
    double max = 300000.0;
    for (int i = 0; i < nMasses; i++) {
        for (int j = i + 1; j <= nMasses; j++) {
            pd.AddSpring(i, j, min, max);
            pd.AddDamper(i, j, min, max);
        }
    }
    pd.SetFixedMass(0);
    for (int i = 1; i <= nMasses; i++) {
        pd.AddInitialVel(i, Random(0.0, 200.0));
        pd.AddInitialDisp(i, Random(0.0, 200.0));
    }
    assert(pd.IsOk());

    // Optimization parameters
    double finalT = 0.005;
    int popSize = 70;
    double geneticAlgoErrorStop = 0.001 / 100.0;  // 0.001%
    // Pimodel based optimization
    int nModels = 12;
    int icPoints = 8;
    int physPoints = 8;
    int order = 3;
    double learningRate = 0.045;
    double minImprovementToEarlyStop = 0.05;  // 5%
    int maxSteps = 10000;
    bool logComplexity = false;
    int timeDiscretization = 5;  // used to look for max accel

    // Train models
    Pimodels models =
        Pimodels(pd, finalT, nModels, icPoints, physPoints, order);
    auto start = Now();
    assert(!models
                .Train(learningRate, learningRate / 500,
                       minImprovementToEarlyStop, maxSteps, logComplexity,
                       verbose)
                .isError);
    res.timeUsecToTrainPimodels = TimeUsecSince(start);

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
            ProblemCreature(&pd, nMasses, &models, timeDiscretization));
        integrationPopulation.push_back(ProblemCreature(&pd, nMasses, finalT));
        randomPopulation.push_back(ProblemCreature(&pd, nMasses, finalT));
    }
    for (int i = 0; i < popSize; i++) {
        for (int j = 0; j < int(integrationPopulation[i].dna.size()); j++) {
            integrationPopulation[i].dna[j].Set(piPopulation[i].dna[j].Get());
            randomPopulation[i].dna[j].Set(piPopulation[i].dna[j].Get());
        }
    }

    // Sort once using explicit integration to find out what the initial guess
    // best solution is, so that we can compare;
    Evolution<ProblemCreature> evolution =
        Evolution<ProblemCreature>(&randomPopulation);
    evolution.SortPopulation();
    auto bestRandomGuessDna = evolution.GetCreature(0)->dna;
    Problem bestRandomGuess = pd.BuildFromDNA(bestRandomGuessDna).val;
    assert(!bestRandomGuess.Integrate(finalT).isError);
    res.bestRandomGuessMaxAccel =
        bestRandomGuess.GetMassMaxAbsAccel(nMasses).val;

    // Pimodel-based optimization
    evolution = Evolution<ProblemCreature>(&piPopulation);
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, verbose);
    res.timeUsecToPimodelGa = TimeUsecSince(start);
    auto pimodelBestDna = evolution.GetCreature(0)->dna;
    Problem pimodelBest = pd.BuildFromDNA(pimodelBestDna).val;
    assert(!pimodelBest.Integrate(finalT).isError);
    res.pimodelGaBestSolutionMaxAccel =
        pimodelBest.GetMassMaxAbsAccel(nMasses).val;

    // Explicit-integration based optimization
    evolution = Evolution<ProblemCreature>(&integrationPopulation);
    start = Now();
    evolution.Evolve(geneticAlgoErrorStop, verbose);
    res.timeUsecToExplicitIntegrationGa = TimeUsecSince(start);
    auto explicitBestDna = evolution.GetCreature(0)->dna;
    Problem explicitBest = pd.BuildFromDNA(explicitBestDna).val;
    assert(!explicitBest.Integrate(finalT).isError);
    res.explicitIntegrationBestSolutionMaxAccel =
        explicitBest.GetMassMaxAbsAccel(nMasses).val;

    if (verbose) {
        std::cout << "Random guess best solution:" << std::endl;
        for (int i = 0; i < int(bestRandomGuessDna.size()); i++) {
            std::cout << bestRandomGuessDna[i].Get() << std::endl;
        }
        std::cout << "Pimodel-based best solution:" << std::endl;
        for (int i = 0; i < int(pimodelBestDna.size()); i++) {
            std::cout << pimodelBestDna[i].Get() << std::endl;
        }
        std::cout << "Explicit-based best solution:" << std::endl;
        for (int i = 0; i < int(explicitBestDna.size()); i++) {
            std::cout << explicitBestDna[i].Get() << std::endl;
        }

        std::cout << "Random guess best: " << std::endl;
        bestRandomGuess.PrintMassTimeHistory(nMasses);
        std::cout << "Pimodel-based best: " << std::endl;
        pimodelBest.PrintMassTimeHistory(nMasses);
        std::cout << "Explicit-based best: " << std::endl;
        explicitBest.PrintMassTimeHistory(nMasses);
    }
    return res;
}

int main(int argc, char *argv[]) {
    // run the whole analysis 3 times
    for (int p = 0; p < 3; p++) {
        for (int i = 1; i <= 5; i++) {
            std::cout << "### Optimize problem with " << i << " mass ###\n";
            printResults(optimize(i), i, p);
            std::cout << "#####################\n";
        }
    }
}