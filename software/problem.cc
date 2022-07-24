#include "maybe.h"
#include "problem.h"
#include "mass.h"

int Problem::GetDof() const{
    return this->masses.size()*2;
}

bool Problem::HasMassAtX(float x) const{
    return (this->massesX_0.find(x) != this->massesX_0.end());
}

Maybe<Void> Problem::AddMass(double m, double x_0){
    Maybe<Void> r;
    if (this->HasMassAtX(x_0)){
        r.isError = true;
        r.errMsg = "Tried to add masses with duplicated x_0";
        return r;
    }

    int dof = this->GetDof();
    this->masses.insert(this->masses.begin(), Mass(m, x_0,dof,dof+1));
    this->massesX_0.insert(x_0);
    return r;
}
