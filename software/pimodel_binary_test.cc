#include "pimodel.h"
#include "utils.h"

int main(int argc, char *argv[]) {
    ProblemDescription pd = ProblemDescription();
    pd.AddMass(1.0, 0.0, 0.0);  // m0
    pd.AddMass(300, 1.0, 1.0);  // m1
    pd.AddMass(120, 1.0, 0.0);  // m2
    pd.AddMass(150, 1.0, 3.0);  // m3
    pd.AddMass(700, 2.0, 0.0);  // m4
    pd.AddMass(80, 3.0, 0.0);   // m5

    double min = 100.0;
    double max = 100000;
    pd.AddSpring(0, 1, min, max);  // k01
    pd.AddSpring(1, 2, min, max);  // k12
    pd.AddSpring(1, 3, min, max);  // k13
    pd.AddSpring(1, 4, min, max);  // k14
    pd.AddDamper(1, 4, min, max);  // c14
    pd.AddSpring(0, 2, min, max);  // k02
    pd.AddDamper(0, 2, min, max);  // c02
    pd.AddSpring(2, 4, min, max);  // k24
    pd.AddDamper(2, 4, min, max);  // c24
    pd.AddSpring(0, 3, min, max);  // k03
    pd.AddDamper(0, 3, min, max);  // c03
    pd.AddSpring(3, 4, min, max);  // k34
    pd.AddDamper(3, 4, min, max);  // c34
    pd.AddSpring(4, 5, min, max);  // k45
    pd.AddDamper(4, 5, min, max);  // c45
    pd.SetFixedMass(0);
    pd.AddInitialVel(200.0);  // initial vel
    assert(pd.IsOk());

    // Learning parameters
    double finalT = 0.05;
    int nModels = 20;
    int icPoints = 5;
    int physPoints = 5;
    int order = 3;
    double learningRate = 0.01;
    // Will stop training only if improvement <1%
    double minImprovementToEarlyStop = 0.01;
    int maxSteps = 10000;
    bool logComplexity = true;
    bool logTraining = true;

    // Train all models
    Pimodels models =
        Pimodels(pd, finalT, nModels, icPoints, physPoints, order);
    auto start = Now();
    assert(!models
                .Train(learningRate, learningRate / 100,
                       minImprovementToEarlyStop, maxSteps, logComplexity,
                       logTraining)
                .isError);
    std::cout << "Time to train models: " << TimeSince(start) << std::endl;

    // Get problem using intermediate value for k and c, and integrate it.
    double mean = (min + max) / 2;
    Maybe<Problem> mP = pd.BuildFromVector(
        std::vector<double>{mean, mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean});
    assert(!mP.isError);
    Problem p = mP.val;
    start = Now();
    assert(!p.Integrate(finalT).isError);
    std::cout << "Time to integrate: " << TimeSince(start) << std::endl;

    std::vector<double> tkc =
        std::vector<double>{0.0,  mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean, mean};
    start = Now();
    models(&tkc);
    std::cout << "Time to make inference with model: " << TimeSince(start)
              << std::endl;
    Maybe<std::vector<double>> X;
    Maybe<std::vector<double>> XDot;
    std::cout << "t,x5,x5Dot,modelX5,modelX5Dot" << std::endl;
    for (int i = 0; i < int(p.t.size()); i += 1) {
        tkc[0] = p.t[i];
        X = models(&tkc);
        XDot = models.GetVelocities(&tkc);

        std::cout << p.t[i] << ",";                                // t,
        std::cout << p.XHistory[i][p.GetMassDispIndex(5)] << ",";  // x5,
        std::cout << p.XHistory[i][p.GetMassVelIndex(5)] << ",";   // x5Dot,
        std::cout << X.val[5] << ",";                              // modelX5,
        std::cout << XDot.val[5] << std::endl;                     // modelX5Dot
    }
}