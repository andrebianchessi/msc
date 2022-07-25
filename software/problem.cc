#include "maybe.h"
#include "problem.h"
#include "mass.h"
#include <tuple>
#include <string>

int Problem::GetDof() const{
    return this->masses.size()*2;
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
    int massId = static_cast<int>(this->masses.size());
    r.val = massId;
    int dof = this->GetDof();
    this->masses.push_back(Mass(m, x, y, massId, dof, dof+1));
    this->initialPositions.insert(std::tuple<double, double>(x,y));
    return r;
}

Maybe<Mass*> Problem::GetMass(int id){
    Maybe<Mass*> r;
    if (id>=static_cast<int>(this->masses.size())){
        r.isError = true;
        r.errMsg = "Tried to get mass with invalid id";
        return r;
    }
    r.val = &(this->masses[id]);
    return r;
}
