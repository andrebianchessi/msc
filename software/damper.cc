#include "damper.h"
#include "mass.h"
#include <boost/numeric/ublas/matrix.hpp>

using namespace boost::numeric::ublas;

Damper::Damper(int id, const Mass* m0, const Mass* m1, double c){
    this->id = id;
    this->m0 = m0;
    this->m1 = m1;
    this->c = c;
};
const Mass* Damper::GetM0(){
    return this->m0;
};
const Mass* Damper::GetM1(){
    return this->m1;
};
const bounded_matrix<float,2,2> Damper::GetM(){
    bounded_matrix<float,2,2> m;
    m(0,0) = (this->m0)->m;
    m(0,1) = 0.0;
    m(1,0) = 0.0;
    m(1,1) = (this->m1)->m;
    return m;
};
const bounded_matrix<float,2,2> Damper::GetC(){
    bounded_matrix<float,2,2> m;
    m(0,0) = -this->c;
    m(0,1) = this->c;
    m(1,0) = this->c;
    m(1,1) = -this->c;
    return m;
};