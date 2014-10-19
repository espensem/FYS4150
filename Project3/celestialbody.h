#ifndef CELESTIALBODY_H
#define CELESTIALBODY_H
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class CelestialBody {

    public:

        CelestialBody(string, vec, vec, double);

        string nameOfBody;
        vec initialPos;
        vec currentPos;
        vec currentVel;
        double mass;
        bool hasLockedPos;

        string getName();
        double getMass();
        vec getInitialPos();
        vec getCurrentPos();
        vec getCurrentVel();
        vec setPos(vec);
        vec setVel(vec);
        vec getDistanceBetweenBodies(CelestialBody);

        bool getHasLockedPos();
        void setHasLockedPos(bool);
        //void resetForce();
};

#endif // CELESTIALBODY_H
