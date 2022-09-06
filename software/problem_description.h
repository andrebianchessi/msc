// Represents the description of a problem that we want to optimize.
// This is mainly used so that we can create multiple Problem instances
// (with different values of spring and damper constants) for the same
// set of masses and geometries
#include <vector>

#include "maybe.h"
#include "problem.h"

struct MassDescription {
    // Mass
    double m;
    // Initial location x
    double px;
    // Initial location y
    double py;
};

struct SpringDescription {
    // Ids of the masses the spring is connected to
    int m0, m1;
    // Min and max possible values for the spring's elastic constant
    double kMin, kMax;
};

struct DamperDescription {
    // Ids of the masses the damper is connected to
    int m0, m1;
    // Min and max possible values for the damper's damping coefficient
    double cMin, cMax;
};

// Parameters passed to Problem::SetInitialVel
struct InitialVelDescription {
    double val;
    int massId;  // if set to -1, sets vel of all masses to val (similar to
                 // Problem class)
};
// Parameters passed to Problem::SetInitialDisp
struct InitialDispDescription {
    double val;
    int massId;  // if set to -1, sets vel of all masses to val (similar to
                 // Problem class)
};

class ProblemDescription {
   public:
    void AddMass(double px, double py, double m);

    void AddInitialVel(int massId, double value);
    void AddInitialVel(double value);
    void AddInitialDisp(int massId, double value);
    void AddInitialDisp(double value);

    // Sets mass as fixed in the problem description
    void SetFixedMass(int massId);

    void AddSpring(int m0, int m1, double kMin, double kMax);
    void AddDamper(int m0, int m1, double cMin, double cMax);
    Maybe<std::shared_ptr<Problem>> BuildRandom();

   private:
    std::vector<MassDescription> masses;
    std::vector<SpringDescription> springs;
    std::vector<DamperDescription> dampers;
    std::vector<InitialVelDescription> initialVels;
    std::vector<InitialDispDescription> initialDisps;
    std::vector<int> fixedMasses;
};