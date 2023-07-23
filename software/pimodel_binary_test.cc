#include "pimodel.h"

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
    pd.AddInitialVel(200.0);  // initial speed
    assert(pd.IsOk());

    // Learning parameters
    double finalT = 0.05;
    int nModels = 3;
    int timeDiscretization = 1;
    int kcDiscretization = 0;
    int order = 3;
    double learningRate = 0.1;
    int batchSize = 1000;
    int maxSteps = 2000;
    bool log = false;

    // Train all models
    Pimodels models = Pimodels(pd, finalT, nModels, timeDiscretization,
                               kcDiscretization, order);
    assert(!models.Train(learningRate, batchSize, maxSteps, log).isError);

    // Get problem using intermediate value for k and c, and integrate it.
    double mean = (min + max) / 2;
    Maybe<Problem> mP = pd.BuildFromVector(
        std::vector<double>{mean, mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean});
    assert(!mP.isError);
    Problem p = mP.val;
    assert(!p.Integrate(finalT).isError);

    std::vector<double> tkc =
        std::vector<double>{0.0,  mean, mean, mean, mean, mean, mean, mean,
                            mean, mean, mean, mean, mean, mean, mean, mean};
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