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

Maybe<Void> Problem::AddMass(double m, double x, double y){
    Maybe<Void> r;
    if (this->HasMassAt(x, y)){
        r.isError = true;
        r.errMsg = "Mass already present at (" +
            std::to_string(x)+","+std::to_string(y)+")";
        return r;
    }

    int dof = this->GetDof();
    this->masses.insert(this->masses.begin(), Mass(m, x, y, dof, dof+1));
    this->initialPositions.insert(std::tuple<double, double>(x,y));
    return r;
}
