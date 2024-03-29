// Represents the description of a problem that we want to optimize.
// This is mainly used so that we can create multiple Problem instances
// (with different values of spring and damper constants) for the same
// set of masses and geometries
#pragma once
#include <vector>

#include "bounded.h"
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
    void AddMass(double m, double px, double py);

    void AddInitialVel(int massId, double value);
    void AddInitialVel(double value);
    void AddInitialDisp(int massId, double value);
    void AddInitialDisp(double value);

    // Sets mass as fixed in the problem description
    void SetFixedMass(int massId);

    void AddSpring(int m0, int m1, double kMin, double kMax);
    void AddDamper(int m0, int m1, double cMin, double cMax);

    // The input is a series of normalized values for springs and dampers.
    // It must have the same length of NumberOfSpringsAndDampers().
    // The first values correspond to the springs, and the last to the dampers.
    // Ex:
    // If the problem has 1 spring and 2 dampers, all with a min of 0 and a max
    // of 10:
    // dna = [0,0.5,1.0] -> k0 = 0, c0 = 5, c1 = 10
    Maybe<Problem> BuildFromDNA(std::vector<Bounded> dna);

    // Similar to BuildFromDNA(std::vector<Bounded> dna),
    // but input is not normalized. Be careful when using this.It's usually
    // better, especially in the context of Genetic Algorithm, to normalize the
    // inputs. This has been implemented to be used in Pimodel (see pimodel.cc)
    Maybe<Problem> BuildFromVector(std::vector<double> springsAndDampers);

    // Returns number of springs + number of dampers
    int NumberOfSpringsAndDampers() const;

    int NumberOfMasses();

    std::vector<MassDescription> masses;
    std::vector<SpringDescription> springs;
    std::vector<DamperDescription> dampers;
    std::vector<InitialVelDescription> initialVels;
    std::vector<InitialDispDescription> initialDisps;
    std::vector<int> fixedMasses;

    // Check if the object is valid
    bool IsOk();
};