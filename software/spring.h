#pragma once
#include <gtest/gtest.h>

#include <boost/numeric/ublas/matrix.hpp>

#include "mass.h"

using namespace boost::numeric::ublas;

class Spring {
    friend class Problem;

   public:
    int id;
    const Mass* GetM0();
    const Mass* GetM1();
    // Returns Mass Matrix
    const bounded_matrix<float, 2, 2> GetM();
    // Returns Stiffness Matrix
    const bounded_matrix<float, 2, 2> GetK();

    double Get_k() { return this->k; };

   private:
    // Spring instances can only be constructed by Problem instances
    Spring(int id, const Mass* m0, const Mass* m1, double k);
    const Mass* m0;
    const Mass* m1;
    double k;
    FRIEND_TEST(ProblemTest, AssignmentTest);
    FRIEND_TEST(ProblemTest, CopyConstructorTest);
};