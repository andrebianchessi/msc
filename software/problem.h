#pragma once
#include <vector>
#include <set>
#include <tuple>
#include "mass.h"
#include "spring.h"
#include "maybe.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <gtest/gtest.h>

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
        Maybe<Void> SetInitialDisp(int massId, double value);
        Maybe<Void> SetInitialVel(int massId, double value);

        // Set initial displacement and velocities for all masses
        void SetInitialDisp(double value);
        void SetInitialVel(double value);

        // Set mass with id massId as fixed, i.e. always zero speed
        // Not that this takes precedence over setting initial velocity
        // explicitly with SetInitialVel.
        Maybe<Void> FixMass(int massId);

        // Builds MInv, K
        // Creates X and XDot with zero values
        void Build();
    private:
        // Initial position (px,py) of all masses in the problem.
        // Note that two masses can't be created at the same initial position.
        std::set<std::tuple<double,double>> initialPositions;
        
        // Masses that are fixed, i.e. have always zero velocity
        std::set<Mass*> fixedMasses;
        FRIEND_TEST(ProblemTest, FixMassTest);

        bool massIsFixed(int massId);

        bool isBuilt;

        // Returns the index of the speed of the mass in the state vector (X)
        int xDotIndex(Mass m);
        int xDotIndex(int xIndex);
        
        // Returns row matrix of displacements
        //  [[x0],
        //   [x1],
        //   ... ]
        matrix<double> getDisp();
        // Returns row matrix of Velocities
        //  [[x0Dot],
        //   [x1Dot],
        //      ... ]
        matrix<double> getVel();
        FRIEND_TEST(ProblemTest, GetDispAndVelTest);
};