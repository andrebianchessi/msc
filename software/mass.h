#pragma once
// Represents an ideal point mass with 1 DOF (horizontal motion only)
class Mass {
    friend class Problem;

    public:
        double m; // mass value
        double x; // horizontal position
        double y; // vertical position
        int xIndex; // global index of horizontal degree of freedom

    private:
        // Mass instances can only be constructed by Problem instances
        Mass(double m, double x, double y, int xIndex){
            this->m = m;
            this->x = x;
            this->y = y;
            this->xIndex = xIndex;
        }
};