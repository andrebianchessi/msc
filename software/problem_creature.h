#pragma once

#include <vector>

#include "bounded.h"
#include "creature.h"
#include "problem_description.h"
#include "utils.h"

class ProblemCreature : public Creature {
   public:
    ProblemDescription* problemDescription;

    // Parameters used to calculate cost
    double t;      // integration duration
    double tStep;  // integration time step

    // Id of the mass we want to minimize the acceleration
    int massId;

    // Random constructor. Receives problem description, id of mass
    // we want to minimize the acceleration and the time parameters
    // used when integrating the problem.
    // WARNING: Make sure to first call IsOk() method to check the problem
    // description.
    ProblemCreature(ProblemDescription* pd, int massId, double t, double tStep);

    double GetCost();

    FRIEND_TEST(ProblemCreatureTest, SimpleTest);
};