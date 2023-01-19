#include "pimodel.h"

int main(int argc, char *argv[]) {
    auto pd = ProblemDescription();
    pd.AddMass(1, 0.0, 0.0);
    pd.AddMass(1, 1.0, 0.0);
    pd.AddSpring(0, 1, 0.0, 1.0);
    pd.AddDamper(0, 1, 0.0, 1.0);
    pd.SetFixedMass(0);
    pd.AddInitialDisp(1, 1.0);

    double T = 1.0;
    int timeDiscretization = 20;
    int kcDiscretization = 1;
    int order = 4;

    // Train model
    Pimodel model =
        Pimodel(&pd, T, timeDiscretization, kcDiscretization, order);
    double lr = 0.00005;
    model.Train(lr, lr / 2, true);

    // Get problem using intermediate value for k and c, and integrate it
    double k = 0.5;
    double c = 0.5;
    Maybe<Problem> mP = pd.BuildFromVector(std::vector<double>{k, c});
    Problem p = mP.val;
    p.Integrate(T);

    // Compare the model's prediction with the problem's result
    std::vector<double> tkc = std::vector<double>{0.0, k, c};
    Maybe<std::vector<double>> X;
    std::cout << "t,x0,modelX0,x1,modelX1" << std::endl;
    for (int i = 0; i < int(p.t.size()); i++) {
        tkc[0] = p.t[i];
        X = model(&tkc);

        std::cout << p.t[i] << ",";
        std::cout << p.XHistory[i][0] << ",";
        std::cout << X.val[0] << ",";
        std::cout << p.XHistory[i][1] << ",";
        std::cout << X.val[1] << std::endl;
    }
}