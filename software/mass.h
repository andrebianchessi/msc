#pragma once
// Represents an ideal point mass
class Mass {
    friend class Problem;

    public:
        double m; // mass value
        double x; // horizontal position
        double y; // vertical position
        int id; // identifies mass in the system
        int xIndex; // global index of position degree of freedom
        int xDotIndex; // global index of speed degree of freedom

    private:
        // Mass instances can only be constructed by Problem instances
        Mass(double m, double x, double y, int id, int xIndex, int xDotIndex){
            this->m = m;
            this->x = x;
            this->id = id;
            this->xIndex = xIndex;
            this->xDotIndex = xDotIndex;
        }
};