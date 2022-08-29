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

Problem::Problem(){
    this -> isBuilt = false;
}

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

Maybe<int> Problem::AddDamper(int m0, int m1, double c){
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
    
    int damperId = this->dampers.size();
    this->dampers.push_back(Damper(damperId,M0,M1,c));
    r.val = damperId;

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

Maybe<Damper*> Problem::GetDamper(int id){
    Maybe<Damper*> r;
    if (id >= static_cast<int>(this->dampers.size()) || id < 0){
        r.isError = true;
        r.errMsg = "Tried to get damper with invalid id";
        return r;
    }
    r.val = &(this->dampers[id]);
    return r;
}

void Problem::Build(){
    int dof = this->GetDof();

    this->MInv.resize(dof,dof,false);
    this->K.resize(dof,dof,false);
    this->C.resize(dof,dof,false);

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

    for (Damper d : this->dampers){
        auto localC = d.GetC();
        this->C(d.m0->xIndex, d.m0->xIndex) = localC(0,0);
        this->C(d.m0->xIndex, d.m1->xIndex) = localC(0,1);
        this->C(d.m1->xIndex, d.m0->xIndex) = localC(1,0);
        this->C(d.m1->xIndex, d.m1->xIndex) = localC(1,1);
    }

    this->X = zero_vector<double>(dof*2);
    this->isBuilt = true;
}

int Problem::GetMassDispIndex(Mass m){
    return m.xIndex;
}
int Problem::GetMassDispIndex(int xIndex){
    return xIndex;
}
int Problem::GetMassVelIndex(Mass m){
    return this->GetDof()+m.xIndex;
}
int Problem::GetMassVelIndex(int xIndex){
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
    int xDotIndex = this->GetMassVelIndex(*e.val);
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
        this->X[this->GetMassVelIndex(m)] = value;
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

matrix<double> Problem::getDisp(vector<double> X, int dof){
    matrix<double> s = matrix<double>(dof,1);
    for (int i = 0; i < dof; i++){
        s(i,0) = X[i];
    }
    return s;
}

matrix<double> Problem::getVel(vector<double> X, int dof){
    matrix<double> v = matrix<double>(dof,1);
    for (int i = 0; i < dof; i++){
        v(i,0) = X[dof + i];
    }
    return v;
}

void Problem::SetXDot(const vector<double> &X, vector<double> &XDot, double t){
    // Set first half of XDot
    for (int i = 0; i < this->GetDof(); i++){
        XDot[i] = X[this->GetMassVelIndex(i)];
    }

    // Calculate xDotDot, which contains the accelerations
    // [x0DotDot, x1DotDot, ...]
    // using the discrete element method
    // 
    // Example for single spring and damper:
    // xDot = [[xDot0],
    //         [xDot1]]
    // xDotDot = [[xDotDot0],
    //            [xDotDot1]]
    // k = [[-k,k],
    //      [k,-k]]
    // c = [[-c,c],
    //      [c,-c]]
    // mInv = [[1/m0, 0],
    //         [0,1/m1]]
    //
    // xDotDot = mInv * ( k*x + c*xDot )
    matrix<double> kx = prod(this->K, Problem::getDisp(X, this->GetDof()));
    matrix<double> cxDot = prod(this->C, Problem::getVel(X, this->GetDof()));
    matrix<double> xDotDot = prod(this->MInv,kx+cxDot);

    // Set second half of XDot
    for (int i = 0; i < this->GetDof(); i++){
        XDot[this->GetMassVelIndex(i)] = xDotDot(i,0);
    }

    // Apply fixed conditions
    for(auto m : this->fixedMasses) {
        XDot[m->xIndex] = 0;
        XDot[this->GetMassVelIndex(m->xIndex)] = 0;
    }
    
}

vector<double> Problem::GetXDot(const vector<double> &X, double t){
    vector<double> XDot = zero_vector<double>(X.size());
    this -> SetXDot(X, XDot, t);
    return XDot;
}

void Problem::save(const vector<double> &X, double t){
    this->t.push_back(t);
    vector<double> xNow  = vector<double>(int(X.size()));
    std::copy(X.begin() , X.end(), xNow.begin());
    this->XHistory.push_back(xNow);
}

void Problem::Integrate(double t0, double t1, double timestep){
    auto setXDot = [this](vector<double> const& X, vector<double> &XDot , double t ) {
        this->SetXDot(X, XDot, t);
    };
    auto save = [this](const vector<double> &X, double t) {
        this->save(X,t);
    };
    integrate( setXDot, this->X, t0, t1, timestep, save);
}

void Problem::PrintMassTimeHistory(int massId){
    using namespace std;
    auto e = this->GetMass(massId);
    if (e.isError){
        cout << "Invalid massId" << endl;
        return;
    }
    auto m = *e.val;
    cout << "t,x,xDot" << endl;

    for (int i = 0; i < int(this->t.size()); i++){
        cout << this->t[i] << ",";
        cout << this->XHistory[i][m.xIndex] << ",";
        cout << this->XHistory[i][this->GetMassVelIndex(m)] << endl;
    }
}