#ifndef VERLET_H
#define VERLET_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Verlet {

    public:

        Verlet(mat, mat, mat, mat, double);

        mat currentPos;
        mat prevPos;
        mat currentVel;
        mat currentAcc;
        double dt;

        mat advancePos();
        mat advanceVel();

        mat newPos;
        mat newVel;



};

#endif // VERLET_H
