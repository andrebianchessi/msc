#pragma once
#include <vector>
#include <set>
#include <tuple>
#include "mass.h"
#include "spring.h"
#include "maybe.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;

// Represents a 1D system of masses, springs and dampers
class Problem{
    public:
        // Masses in the problem
        std::vector<Mass> masses;

        // Springs in the problem
        std::vector<Spring> springs;

        // Inverse of Mass Matrix
        mapped_matrix<double> MInv;

        // Stiffness Matrix
        mapped_matrix<double> K;

        // Vector of positions
        vector<double> X;

        // Vector of speeds
        vector<double> XDot;

        // Returns true if system has a mass with initial position (x,y)
        bool HasMassAt(double x, double y) const;

        // Add a Mass object to the problem and returns its id
        Maybe<int> AddMass(double m, double x, double y);

        // Returns mass by id
        Maybe<Mass*> GetMass(int id);

        // Returns degrees of freedom
        int GetDof() const;

        // Add a Spring object to the problem and returns its id
        Maybe<int> AddSpring(int m0, int m1, double k);

        // Returns mass by id
        Maybe<Spring*> GetSpring(int id);

        // Set initial displacement and velocities
        Maybe<Void> SetInitialX(int massId, double value);
        Maybe<Void> SetInitialXDot(int massId, double value);

        // Builds MInv, K
        // Creates X and XDot with zero values
        void Build();
    private:
        // Initial position (x,y) of all masses in the problem.
        // Note that two masses can't be created at the same initial position.
        std::set<std::tuple<double,double>> initialPositions;
};