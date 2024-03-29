#pragma once
#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix.hpp>

#include "mass.h"

using namespace boost::numeric::ublas;

class Damper {
    friend class Problem;

   public:
    int id;
    const Mass* GetM0();
    const Mass* GetM1();
    // Returns Mass Matrix
    const bounded_matrix<double, 2, 2> GetM();
    // Returns Damping Matrix
    const bounded_matrix<double, 2, 2> GetC();

    double Get_c() { return this->c; };

   private:
    // Damper instances can only be constructed by Problem instances
    Damper(int id, const Mass* m0, const Mass* m1, double c);
    const Mass* m0;
    const Mass* m1;
    double c;
    FRIEND_TEST(ProblemTest, AssignmentTest);
    FRIEND_TEST(ProblemTest, CopyConstructorTest);
};