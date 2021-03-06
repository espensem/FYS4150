#ifndef SOLARSYSTEM_H
#define SOLARSYSTEM_H

#include <cmath>
#include <time.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <armadillo>
#include <vector>
#include "celestialbody.h"

using std::vector;
using namespace std;
using namespace arma;

class SolarSystem {

    public:

        SolarSystem(int, string,double);

        int dimension;
        string nameOfFile;
        int getNuOfDimensions();
        string getNameOfFile();
        int getNuOfBodies();

        void createCelestialBody(string, bool=false);
        void addCelestialBody(CelestialBody);
        vector<CelestialBody> bodies;
        string initialParametersOfBodies;

        mat getAllCurrentPos();
        mat getAllCurrentVel();
        mat getAllCurrentAcc();
        void setAllPos(mat);
        void setAllVel(mat);
        double totalMass;
        double getTotalMass();
        vec centerOfMass;

        vec getForces(CelestialBody);
        double kineticEnergy;
        double potentialEnergy;
        vec angularMomentum;
        void calculateEnergy();
        void calculateTotalAngularMomentum();
        double totalEnergy();

        double dt;
        mat prevPos;
        mat pPrev;

        void rk4Solver();
        void verletFirstStep();
        void verletSolver(double);

};

#endif // SOLARSYSTEM_H
