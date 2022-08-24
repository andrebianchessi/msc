#include "maybe.h"
#include "problem.h"
#include "mass.h"
#include <tuple>
#include <string>

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