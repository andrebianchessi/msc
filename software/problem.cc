#include "maybe.h"
#include "problem.h"
#include "mass.h"
#include <tuple>
#include <string>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/operations.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::ublas;
using namespace boost::numeric::odeint;

int Problem::GetDof() const{
    return this->masses.size();
}

bool Problem::HasMassAt(double px, double py) const{
    auto position = std::tuple<double, double>(px,py);
    return (this->initialPositions.find(position) != this->initialPositions.end());
}

bool Problem::massIsFixed(int massId){
    auto e = this->GetMass(massId);
    if (e.isError){
        return false;
    }
    return (this->fixedMasses.find(e.val) != this->fixedMasses.end());
}

Maybe<int> Problem::AddMass(double m, double px, double py){
    Maybe<int> r;
    if (this->HasMassAt(px, py)){
        r.isError = true;
        r.errMsg = "Mass already present at (" +
            std::to_string(px)+","+std::to_string(py)+")";
        return r;
    }
    int dof = this->GetDof();
    r.val = dof;
    this->masses.push_back(Mass(m, px, py, dof));
    this->initialPositions.insert(std::tuple<double, double>(px,py));
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

    this->X = zero_vector<double>(dof*2);
    this->isBuilt = true;
}

int Problem::xDotIndex(Mass m){
    return this->GetDof()+m.xIndex;
}
int Problem::xDotIndex(int xIndex){
    return this->GetDof()+xIndex;
}

Maybe<Void> Problem::SetInitialDisp(int massId, double value){
    Maybe<Void> r;
    auto e = this->GetMass(massId);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid massId";
        return r;
    }
    if (!this->isBuilt){
        r.isError = true;
        r.errMsg = "Problem not built. First call Build().";
        return r;
    }
    int xIndex = e.val->xIndex;
    if(xIndex >= int(this->X.size())){
        r.isError = true;
        r.errMsg = "Invalid xIndex";
        return r;
    }
    this->X[xIndex] = value;
    return r;
}

Maybe<Void> Problem::SetInitialVel(int massId, double value){
    Maybe<Void> r;
    if (this->massIsFixed(massId)){
        r.isError = true;
        r.errMsg = "Mass is fixed";
        return r;
    }
    auto e = this->GetMass(massId);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid massId";
        return r;
    }
    if (!this->isBuilt){
        r.isError = true;
        r.errMsg = "Problem not built. First call Build().";
        return r;
    }
    int xDotIndex = this->xDotIndex(*e.val);
    if(xDotIndex >= int(this->X.size())){
        r.isError = true;
        r.errMsg = "Invalid xDotIndex.";
        return r;
    }
    this->X[xDotIndex] = value;
    return r;
}

void Problem::SetInitialDisp(double value){
    for (int i = 0; i < int(this->X.size()); i++){
        this->X[i] = value;
    }
}

void Problem::SetInitialVel(double value){
    for (Mass m : this->masses){
        this->X[this->xDotIndex(m)] = value;
    }
}

Maybe<Void> Problem::FixMass(int massId){
    Maybe<Void> r;
    auto e = this->GetMass(massId);
    if (e.isError){
        r.isError = true;
        r.errMsg = "Invalid massId";
        return r;
    }
    this->fixedMasses.insert(e.val);
    return r;
}

matrix<double> Problem::getDisp(){
    matrix<double> s = matrix<double>(this->GetDof(),1);
    for (int i = 0; i < this->GetDof(); i++){
        s(i,0) = this->X[i];
    }
    return s;
}

matrix<double> Problem::getVel(){
    matrix<double> v = matrix<double>(this->GetDof(),1);
    for (int i = 0; i < this->GetDof(); i++){
        v(i,0) = this->X[this->xDotIndex(i)];
    }
    return v;
}

void Problem::SetXDot(const vector<double> &X, vector<double> &XDot, double t){
    // Set first half of XDot
    for (int i = 0; i < this->GetDof(); i++){
        XDot[i] = X[this->xDotIndex(i)];
    }

    // 1D matrix with the second derivatives
    // calculating with discrete element method
    // [x0DotDot, x1DotDot, ...]
    matrix<double> KDisp = prod(this->K, this->getDisp());
    matrix<double> xDotDot = prod(this->MInv,KDisp);

    // Set second half of XDot
    for (int i = 0; i < this->GetDof(); i++){
        XDot[this->xDotIndex(i)] = xDotDot(i,0);
    }

    // Apply fixed conditions
    for(auto m : this->fixedMasses) {
        XDot[m->xIndex] = 0;
        XDot[this->xDotIndex(m->xIndex)] = 0;
    }

}

void Problem::write(const vector<double> &X, double t){
    using namespace std;
    cout.precision(3);
    cout << "t: "<< t << '\t';
    for (int i = 0; i < this->GetDof(); i++){
        cout << "x"<< i << ": " << X[i] << '\t';
    }
    for (int i = this->GetDof(); i < 2*this->GetDof(); i++){
        cout << "x"<< i - this->GetDof() << "Dot:" << " " << X[i] << '\t';
    }
    cout << endl;
}

void Problem::Integrate(double t0, double t1, double timestep){
    auto setXDot = [this](vector<double> const& X, vector<double> &XDot , double t ) {
        this->SetXDot(X, XDot, t);
    };
    auto write = [this](const vector<double> &X, double t) {
        this->write(X,t);
    };
    integrate( setXDot, this->X, t0, t1, timestep, write);
}