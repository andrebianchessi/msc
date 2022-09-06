#pragma once
// Represents an ideal point mass with 1 DOF (horizontal motion only)
class Mass {
    friend class Problem;

   public:
    double m;    // mass value
    double px;   // horizontal position
    double py;   // vertical position
    int xIndex;  // global index of horizontal degree of freedom

   private:
    // Mass instances can only be constructed by Problem instances
    Mass(double m, double px, double py, int xIndex) {
        this->m = m;
        this->px = px;
        this->py = py;
        this->xIndex = xIndex;
    }
};