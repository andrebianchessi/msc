#pragma once
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
    const bounded_matrix<float, 2, 2> GetM();
    // Returns Damping Matrix
    const bounded_matrix<float, 2, 2> GetC();

   private:
    // Damper instances can only be constructed by Problem instances
    Damper(int id, const Mass* m0, const Mass* m1, double c);
    const Mass* m0;
    const Mass* m1;
    double c;
};