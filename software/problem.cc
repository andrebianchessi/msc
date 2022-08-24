#include "maybe.h"
#include "problem.h"
#include "mass.h"
#include <tuple>
#include <string>
#include <boost/numeric/ublas/matrix_sparse.hpp>

using namespace boost::numeric::ublas;

int Problem::GetDof() const{
    return this->masses.size();
}

bool Problem::HasMassAt(double x, double y) const{
    auto position = std::tuple<double, double>(x,y);
    return (this->initialPositions.find(position) != this->initialPositions.end());
}

Maybe<int> Problem::AddMass(double m, double x, double y){
    Maybe<int> r;
    if (this->HasMassAt(x, y)){
        r.isError = true;
        r.errMsg = "Mass already present at (" +
            std::to_string(x)+","+std::to_string(y)+")";
        return r;
    }
    int dof = this->GetDof();
    r.val = dof;
    this->masses.push_back(Mass(m, x, y, dof));
    this->initialPositions.insert(std::tuple<double, double>(x,y));
    return r;
}

Maybe<Mass*> Problem::GetMass(int id){
    Maybe<Mass*> r;
    if (id >= static_cast<int>(this->masses.size()) || id < 0){
        r.isError = true;
        r.errMsg = "Tried to get mass with invalid id";
        return r;
    }
    r.val = &(this->masses[id]);
    return r;
}

Maybe<int> Problem::AddSpring(int m0, int m1, double k){
    Maybe<int> r;

    if (m0 == m1){
        r.isError = true;
        r.errMsg = "m0 and m1 must be different";
        return r;  
    }

    auto e = this->GetMass(m0);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid m0";
        return r;
    }
    auto M0 = e.val;

    e = this->GetMass(m1);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid m1";
        return r;
    }
    auto M1 = e.val;
    
    int springId = this->springs.size();
    this->springs.push_back(Spring(springId,M0,M1,k));
    r.val = springId;

    return r;
}

Maybe<Spring*> Problem::GetSpring(int id){
    Maybe<Spring*> r;
    if (id >= static_cast<int>(this->springs.size()) || id < 0){
        r.isError = true;
        r.errMsg = "Tried to get spring with invalid id";
        return r;
    }
    r.val = &(this->springs[id]);
    return r;
}

void Problem::Build(){
    int dof = this->GetDof();

    this->MInv.resize(dof,dof,false);
    this->K.resize(dof,dof,false);

    for (Mass mass : this->masses){
        MInv(mass.xIndex,mass.xIndex) = 1/mass.m;
    }

    for (Spring s : this->springs){
        auto localK = s.GetK();
        this->K(s.m0->xIndex, s.m0->xIndex) = localK(0,0);
        this->K(s.m0->xIndex, s.m1->xIndex) = localK(0,1);
        this->K(s.m1->xIndex, s.m0->xIndex) = localK(1,0);
        this->K(s.m1->xIndex, s.m1->xIndex) = localK(1,1);
    }

    this->X = zero_vector<double>(dof);
    this->XDot = zero_vector<double>(dof);
}

Maybe<Void> Problem::SetInitialX(int massId, double value){
    Maybe<Void> r;
    auto e = this->GetMass(massId);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid massId";
        return r;
    }
    int xIndex = e.val->xIndex;
    if(xIndex >= this->X.size()){
        r.isError = true;
        r.errMsg = "Invalid xIndex. Build() has probably not been called yet.";
        return r;
    }
    this->X[xIndex] = value;
    return r;
}

Maybe<Void> Problem::SetInitialXDot(int massId, double value){
    Maybe<Void> r;
    auto e = this->GetMass(massId);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid massId";
        return r;
    }
    int xIndex = e.val->xIndex;
    if(xIndex >= this->XDot.size()){
        r.isError = true;
        r.errMsg = "Invalid xIndex. Build() has probably not been called yet.";
        return r;
    }
    this->XDot[xIndex] = value;
    return r;
}