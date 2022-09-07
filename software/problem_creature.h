#pragma once

#include <vector>

#include "bounded.h"
#include "creature.h"
#include "problem_description.h"

using namespace std;

class ProblemCreature : public Creature {
   public:
    ProblemDescription* problemDescription;

    // Random constructor
    // Make sure pd->BuildRandom() doesn't return an error
    ProblemCreature(ProblemDescription* pd) {
        this->problemDescription = pd;
        this->dna = vector<Bounded>(pd->NumberOfSpringsAndDampers());
        auto p = pd->BuildRandom().val;
        int i = 0;
        for (auto s : p->springs) {
            this->dna[i] = s.k;
        }
    }

    double GetCost() {
        double xVal = this->GetX();
        double yVal = this->GetY();
        return abs(xVal * xVal + yVal * yVal + 2 * xVal + yVal);
    }
};