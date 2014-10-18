#ifndef RK4_H
#define RK4_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class RK4 {
    /*
     * Takes in matrices for position (velocity), velocity (acceleration) and acceleration
     * of bodies in a N-body simulation. Rows are the dimensions, columns are the bodies.
     * dt is the step size (constant).
     */

    public:

        RK4(mat, mat, double);

        mat currentPos;
        mat currentVel;
        mat currentAcc;
        double dt;

        mat advancePos();
        mat advanceVel();

        mat newPos;
        mat newVel;
        mat vk1, vk2, vk3, vk4, pk1, pk2, pk3, pk4;


};

#endif // RK4_H
