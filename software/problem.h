#pragma once
#include <vector>
#include <set>
#include <tuple>
#include "mass.h"
#include "spring.h"
#include "damper.h"
#include "maybe.h"
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <gtest/gtest.h>

using namespace boost::numeric::ublas;

// Represents a 1D system of masses, springs and dampers
class Problem {
public:
    // Standard constructor
    Problem();

    // Masses in the problem
    std::vector<Mass> masses;

    // Springs in the problem
    std::vector<Spring> springs;

    // Dampers in the problem
    std::vector<Damper> dampers;

    // Inverse of Mass Matrix
    mapped_matrix<double> MInv;

    // Stiffness Matrix
    mapped_matrix<double> K;

    // Damping Matrix
    mapped_matrix<double> C;

    // State Vector
    // Contains all displacements, followed by all velocities
    // X = [x0, x1,... xN, xDot0, xDot1, ..., xDotN]
    vector<double> X;
    
    // Vector to store the time instants after the system is integrated
    std::vector<double> t;
    // Contains the state vector after Integrate is called and
    // saveWholeStateVector =  true
    std::vector<vector<double>> XHistory;
    // Contains the position and velocity of mass with id = massToSave after
    // Integrate is called and saveWholeStateVector =  false.
    // XiHistory = [[xi, xiDot], [xi, xiDot], ...]
    std::vector<vector<double>> XiHistory;

    // After calling this function, the position and velocity of only the
    // specified mass will be tracked after calling Integrate. In that case,
    // XHistory will remain empty. This is used to avoid saving all the state
    // vector if you're only interested in one of the masses.
    Maybe<Void> TrackOnlyMass(int massId);

    // Sets in XDot the values of the derivatives of the state vector,
    // i.e after calling this function:
    // XDot = [x0Dot, x1Dot, ... , x0DotDot, x1DotDot, ...]
    //
    // odeint integrate expects a function with this signature.
    // That's why X is passed as an argument, and not accessed with
    // this->X.
    void SetXDot(const vector<double> &X, vector<double> &XDot, double t);
    // void operator () (const vector<double> &X, vector<double> &XDot, double t);

    // Returns true if system has a mass with initial position (px,py)
    bool HasMassAt(double px, double py) const;

    // Add a Mass object to the problem and returns its id
    Maybe<int> AddMass(double m, double px, double py);

    // Returns mass by id
    Maybe<Mass *> GetMass(int id);

    // Returns degrees of freedom
    int GetDof() const;

    // Add a Spring object to the problem and returns its id
    Maybe<int> AddSpring(int m0, int m1, double k);

    // Add a Damper object to the problem and returns its id
    Maybe<int> AddDamper(int m0, int m1, double c);

    // Returns mass by id
    Maybe<Spring *> GetSpring(int id);

    // Returns damper by id
    Maybe<Damper *> GetDamper(int id);

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

    // Integrates the system and saves the state vectors at t and XHistory
    void Integrate(double t0, double t1, double timestep);

    // Prints the time history of mass by id
    void PrintMassTimeHistory(int massId);

    // Returns the index of the displacement of the mass in the state vector (X)
    int GetMassDispIndex(Mass m);
    int GetMassDispIndex(int xIndex);

    // Returns the index of the velocity of the mass in the state vector (X)
    int GetMassVelIndex(Mass m);
    int GetMassVelIndex(int xIndex);

private:
    // Initial position (px,py) of all masses in the problem.
    // Note that two masses can't be created at the same initial position.
    std::set<std::tuple<double, double>> initialPositions;

    // Masses that are fixed, i.e. have always zero velocity
    std::set<Mass *> fixedMasses;
    FRIEND_TEST(ProblemTest, FixMassTest);

    bool massIsFixed(int massId);

    bool isBuilt; // default = false
    bool saveWholeStateVector; // default = true
    int massToSave;

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

    // Save current state vector
    void save(const vector<double> &X, double t);
};