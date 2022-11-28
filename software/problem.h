#pragma once
#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <set>
#include <tuple>
#include <vector>

#include "damper.h"
#include "mass.h"
#include "maybe.h"
#include "spring.h"

using namespace boost::numeric::ublas;

// Represents a 1D system of masses, springs and dampers
class Problem {
    friend class ProblemDescription;
    friend class Pimodel;

   public:
    // Standard constructor
    Problem();

    // Assignment
    Problem &operator=(const Problem &rhs);

    // Copy constructor
    Problem(const Problem &rhs);

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
    // Contains the state vector after Integrate is called.
    std::vector<vector<double>> XHistory;
    // Contains the accelerations of each mass after Integrate is called:
    // AccelHistory =
    //     [[x0DotDot at t0, x1DotDot at t0, ...], [x0DotDot at t1, ...], ...]
    std::vector<vector<double>> AccelHistory;

    // Sets in XDot the values of the derivatives of the state vector,
    // i.e after calling this function:
    // XDot = [x0Dot, x1Dot, ... , x0DotDot, x1DotDot, ...]
    //
    // odeint integrate expects a function with this signature.
    // That's why X is passed as an argument, and not accessed with
    // this->X.
    void SetXDot(const vector<double> &X, vector<double> &XDot, double t);

    // This is just a wrapper that instantiates a vector, calls
    // SetXDot on it, and returns it.
    vector<double> GetXDot(const vector<double> &X, double t);

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
    Maybe<Void> SetInitialDisp(double value);
    Maybe<Void> SetInitialVel(double value);

    // Set mass with id massId as fixed, i.e. always zero speed
    // Should be called after Build() and SetInitialVel(value);
    Maybe<Void> FixMass(int massId);

    // Builds MInv, K
    // Creates X and XDot with zero values
    void Build();

    // Integrates the system up to given time; saves the state vectors at
    // this->XHistory and the time instants at this->t
    Maybe<Void> Integrate(double t);

    // Prints the time history of mass by id
    void PrintMassTimeHistory(int massId);

    // Returns the index of the displacement of the mass in the state vector (X)
    int GetMassDispIndex(Mass m);
    int GetMassDispIndex(int xIndex);

    // Returns the index of the velocity of the mass in the state vector (X)
    int GetMassVelIndex(Mass m);
    int GetMassVelIndex(int xIndex);

    // Returns maximum absolute value of acceleration of a mass. Must be called
    // after Integrate.
    Maybe<double> GetMassMaxAbsAccel(int massId);
    // Returns maximum value of acceleration of a mass. Must be called
    // after Integrate.
    Maybe<double> GetMassMaxAccel(int massId);
    // Returns minimum value of acceleration of a mass. Must be called
    // after Integrate.
    Maybe<double> GetMassMinAccel(int massId);

    // Clears the XHistory and AccelHistory arrays. X is set to it's initial
    // value.
    // Should be used to save memory when the past values are no longer needed.
    void ClearHistory();

   private:
    // Initial position (px,py) of all masses in the problem.
    // Note that two masses can't be created at the same initial position.
    std::set<std::tuple<double, double>> initialPositions;

    // Masses that are fixed, i.e. have always zero velocity
    std::set<Mass *> fixedMasses;
    FRIEND_TEST(ProblemTest, FixMassTest);
    FRIEND_TEST(ProblemTest, AssignmentTest);
    FRIEND_TEST(ProblemTest, CopyConstructorTest);

    bool massIsFixed(int massId);

    bool isBuilt;       // default = false
    bool isIntegrated;  // default = false
    FRIEND_TEST(ProblemTest, ClearHistoryTest);

    // Returns row matrix of displacements
    //  [[x0],
    //   [x1],
    //   ... ]
    // from the state vector provided
    static matrix<double> getDisp(vector<double> X, int dof);
    // Returns row matrix of Velocities
    //  [[x0Dot],
    //   [x1Dot],
    //      ... ]
    // from the state vector provided
    static matrix<double> getVel(vector<double> X, int dof);
    FRIEND_TEST(ProblemTest, GetDispAndVelTest);

    // Save current state vector
    void save(const vector<double> &X, double t);
};