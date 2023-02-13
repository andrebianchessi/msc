#include "pimodel.h"

int main(int argc, char *argv[]) {
    double mass = 1;
    double kMin = 0.8;
    double kMax = 1.0;
    double cMin = 0.0;
    double cMax = 0.1;
    double tMax = 5.0;
    double initialDisp = 1.0;

    auto pd = ProblemDescription();
    pd.AddMass(mass, 0.0, 0.0);
    pd.AddMass(mass, 1.0, 0.0);
    pd.AddSpring(0, 1, kMin, kMax);
    pd.AddDamper(0, 1, cMin, cMax);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, initialDisp);

    int timeBuckets = 10;
    int timeDiscretization = 2;
    int kcDiscretization = 1;
    int order = 2;
    double learningRate = 0.001;
    int maxSteps = 100;
    bool log = true;
    // Train model
    Pimodel model = Pimodel(&pd, tMax, timeBuckets, timeDiscretization,
                            kcDiscretization, order);
    model.Train(learningRate, maxSteps, log);

    // Get problem using intermediate value for k and c, and integrate it.
    // Then, compare the model's prediction with the problem's result.
    double k = (kMin + kMax) / 2;
    double c = (cMin + cMax) / 2;
    Maybe<Problem> mP = pd.BuildFromVector(std::vector<double>{k, c});
    Problem p = mP.val;
    p.Integrate(tMax);

    std::vector<double> tkc = std::vector<double>{0.0, k, c};
    Maybe<std::vector<double>> X;
    std::cout << "t,x0,modelX0,x1,modelX1" << std::endl;
    for (int i = 0; i < int(p.t.size());
         i += std::max(1, int(p.t.size()) / 20)) {
        tkc[0] = p.t[i];
        X = model(&tkc);

        std::cout << p.t[i] << ",";
        std::cout << p.XHistory[i][0] << ",";
        std::cout << X.val[0] << ",";
        std::cout << p.XHistory[i][1] << ",";
        std::cout << X.val[1] << std::endl;
    }
}