#pragma once
// Represents an ideal point mass
class Mass {
    friend class Problem;
    // Mass instances can only be constructed by Problem instances
    private:
        Mass(double m, double x_0, int xIndex, int xDotIndex){
            this->m = m;
            this->x_0 = x_0;
            this->xIndex = xIndex;
            this->xDotIndex = xDotIndex;
        }

    public:
        double m; // mass
        double x_0; // initial location
        int xIndex; // global index of position degree of freedom
        int xDotIndex; // global index of speed degree of freedom
};