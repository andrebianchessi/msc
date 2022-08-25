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

        // State Vector
        // Contains all displacements, followed by all velocities
        // X = [x0, x1,... xN, xDot0, xDot1, ..., xDotN]
        vector<double> X;

        // Returns time derivative of state vector
        // XDot = [xDot0, xDot1, ..., xDotN, xDotDot0, ..., xDotDotN]
        vector<double> XDot();

        // Returns true if system has a mass with initial position (px,py)
        bool HasMassAt(double px, double py) const;

        // Add a Mass object to the problem and returns its id
        Maybe<int> AddMass(double m, double px, double py);

        // Returns mass by id
        Maybe<Mass*> GetMass(int id);

        // Returns degrees of freedom
        int GetDof() const;

        // Add a Spring object to the problem and returns its id
        Maybe<int> AddSpring(int m0, int m1, double k);

        // Returns mass by id
        Maybe<Spring*> GetSpring(int id);

        // Set initial displacement and velocities by massId
        Maybe<Void> SetInitialX(int massId, double value);
        Maybe<Void> SetInitialXDot(int massId, double value);

        // Set initial displacement and velocities for all masses
        void SetInitialX(double value);
        void SetInitialXDot(double value);

        // Builds MInv, K
        // Creates X and XDot with zero values
        void Build();
    private:
        // Initial position (px,py) of all masses in the problem.
        // Note that two masses can't be created at the same initial position.
        std::set<std::tuple<double,double>> initialPositions;
        
        bool isBuilt;
        
        // Returns the index of the speed of the mass in the state vector (X)
        int xDotIndex(Mass m);
};