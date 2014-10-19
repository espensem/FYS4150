#include <iostream>
#include "solarsystem.h"
#include <string>

using namespace std;

int main(int argc, char* argv[]) {

    if (argc < 5) {
        cout << "Wrong number of inputs" << endl;
        cout << "Usage: /exe dimensions simulationTime timeSteps nameOfFile verlet(optional)" << endl;
        return 1;
    }

    if (argc != 6) {
        cout << "Default solver rk4 chosen. " << endl;
        cout << "To use verlet solver last argument needs to be verlet." << endl;

    }
    int dim = atof(argv[1]);
    double simulationTime = atof(argv[2]);  // in years
    double timeSteps = atof(argv[3]);
    string nameOfFile = argv[4];
    string solver;
    if (argc == 6) {
        solver = argv[5];
    }

    double dt = simulationTime / timeSteps;

    SolarSystem mySolarSystem(dim, nameOfFile, dt);
    mySolarSystem.createCelestialBody("sun");
    mySolarSystem.createCelestialBody("mercury");
    mySolarSystem.createCelestialBody("venus");
    mySolarSystem.createCelestialBody("earth");
    mySolarSystem.createCelestialBody("mars");
    mySolarSystem.createCelestialBody("jupiter");
    mySolarSystem.createCelestialBody("saturn");
    mySolarSystem.createCelestialBody("uranus");
    mySolarSystem.createCelestialBody("neptune");
    mySolarSystem.createCelestialBody("pluto");

    // opening the inputfile, and clearing its content
    stringstream s;
    s << "../Project3/data/" << nameOfFile << ".txt";
    string filename = s.str();

    ofstream outfile(filename.c_str());

    // writing body names in first line of the file
    string name;
    for (int i = 0; i < mySolarSystem.getNuOfBodies(); i++) {
        CelestialBody &body = mySolarSystem.bodies[i];
        name = body.getName();
        outfile << name << " ";
    }
    outfile << endl;


    // integrating

    // first step if solver is verlet
    //mat previousPos;
    if (argc == 6) {
        cout << "verlet" << endl;
        mySolarSystem.verletFirstStep();
    }
    double time = 0;
    while (time <= simulationTime) {

        // calculate energy
        mySolarSystem.calculateEnergy();

        mat posNow = mySolarSystem.getAllCurrentPos();                

        // write to file
        outfile << time << " " << mySolarSystem.kineticEnergy << " " << mySolarSystem.potentialEnergy <<
                   " " << mySolarSystem.totalEnergy() << " ";
        for ( int k = 0; k < mySolarSystem.getNuOfBodies(); k++) {
            outfile << posNow(0,k) << " " << posNow(1,k) <<   " ";
        }
        outfile << endl;

        // advancing
        if (argc == 6) {
            mySolarSystem.verletSolver();
        }
        else {
            // using RK4 as the default
            mySolarSystem.rk4Solver();
        }

        time += dt;
    }

    // close file
    outfile.close();

    return 0;
}

