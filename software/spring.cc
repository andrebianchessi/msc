#include "spring.h"

#include <boost/numeric/ublas/matrix.hpp>

#include "mass.h"

using namespace boost::numeric::ublas;

Spring::Spring(int id, const Mass* m0, const Mass* m1, double k) {
    this->id = id;
    this->m0 = m0;
    this->m1 = m1;
    this->k = k;
};
const Mass* Spring::GetM0() { return this->m0; };
const Mass* Spring::GetM1() { return this->m1; };
const bounded_matrix<double, 2, 2> Spring::GetM() {
    bounded_matrix<double, 2, 2> m;
    m(0, 0) = (this->m0)->m;
    m(0, 1) = 0.0;
    m(1, 0) = 0.0;
    m(1, 1) = (this->m1)->m;
    return m;
};
const bounded_matrix<double, 2, 2> Spring::GetK() {
    bounded_matrix<double, 2, 2> m;
    m(0, 0) = -this->k;
    m(0, 1) = this->k;
    m(1, 0) = this->k;
    m(1, 1) = -this->k;
    return m;
};