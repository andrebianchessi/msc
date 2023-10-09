#pragma once

#include <vector>

#include "bounded.h"
#include "creature.h"
#include "pimodel.h"
#include "problem_description.h"
#include "utils.h"

// Class used in Genetic Algorithm to find values of springs and dampers
// that will minimize the acceleration experienced by a mass. The cost can
// be calculated in two ways: either by explicit time integration, or by
// training a physics informed ML models that describe the displacement of
// each mass of the system as a function of time and of the masses and dampers
// that make up the system.
class ProblemCreature : public Creature {
   public:
    ProblemDescription* problemDescription;

    // Random constructor. Receives problem description, id of mass
    // we want to minimize the acceleration and the time duration
    // used when integrating the problem. The accelerations are found by
    // implicit time integration of the system's differential equation.
    // WARNING: Make sure to first call IsOk() method to check the problem
    // description.
    ProblemCreature(ProblemDescription* pd, int massId, double T);

    // Random constructor that receives a problem description and a pimodels
    // used to calculate the cost. timesToCheck specifies how many time instants
    // will be checked to find the max acceleration inside GetCost.
    ProblemCreature(ProblemDescription* pd, int massId, Pimodels* pimodels,
                    int timesToCheck);

    double GetCost();

    FRIEND_TEST(ProblemCreatureTest, SimpleTest);

   private:
    // Parameters used to calculate cost
    double T;  // integration duration

    // Id of the mass we want to minimize the acceleration.
    int massId;

    Pimodels* pimodels;
    int timesToCheck;
};