#pragma once
#include <vector>
#include <set>
#include "mass.h"
#include "maybe.h"

// Represents a 1D system of masses, springs and dampers
class Problem{
    public:
        // Masses in the problem
        std::vector<Mass> masses;
        // x_0 of all masses in the problem.
        // Note that two masses can't be created with the same x_0
        std::set<float> massesX_0;

        // Returns true if system has mass with x_0 = x
        bool HasMassAtX(float x) const;

        // Add a Mass object to the problem
        Maybe<Void> AddMass(double m, double x_0);

        // Returns degrees of freedom
        int GetDof() const;
};