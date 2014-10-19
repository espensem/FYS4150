#include "solarsystem.h"

SolarSystem::SolarSystem(int dimension, string nameOfFile, double dt) {
    /* OBS: This code is written for ONLY two dimensions, assuming coplanar motion.
     * The argument thus have to be 2.
     */
    this->dt = dt;
    this->dimension = dimension;
    this->nameOfFile = nameOfFile;

    this->initialParametersOfBodies = "../Project3/data/bodies.txt";
}

int SolarSystem::getNuOfDimensions() {
    return this->dimension;
}

string SolarSystem::getNameOfFile() {
    return this->nameOfFile;
}

void SolarSystem::addCelestialBody(CelestialBody newBody) {
    bodies.push_back(newBody);
    return;
}


void SolarSystem::createCelestialBody(string nameOfBody, bool locked) {
    /* Collect initial values of the body from "bodies.txt" and write
     * them to the file created above
     */

    ifstream infile;
    try {
        infile.open(this->initialParametersOfBodies.c_str());
    }
    catch(...) {
        cout << "Failed to open file: " << this->initialParametersOfBodies << "\n" <<
                "No such file found" << endl;
    }

    double rx, ry, vx, vy, mass;
    string line;
    bool bodyFound = false;

    while (getline(infile, line)) {
        vector<string> tokens;
        boost::split(tokens, line, boost::is_any_of(" "));
        //if (strcmp(tokens[0].c_str(), nameOfBody.c_str()) == 0) {
        if (tokens[0] == nameOfBody) {
            rx = atof(tokens[1].c_str());
            ry = atof(tokens[2].c_str());
            vx = atof(tokens[3].c_str());
            vy = atof(tokens[4].c_str());
            mass = atof(tokens[5].c_str());
            bodyFound = true;
            break;
        }
    }

    if (!bodyFound) {
        cout << "Couldn't find " << nameOfBody << " in the file bodies.txt." <<
                "Make shure the file is there and that the name given to the " <<
                "body corresponds to that in the file" << endl;
        return;
    }

    infile.close();

    vec position;
    position.zeros(2);
    position(0) = rx;
    position(1) = ry;
    vec velocity;
    velocity.zeros(2);
    velocity(0) = vx;
    velocity(1) = vy;


    cout << "Celestial body " << nameOfBody << " created." << endl;
    cout << "Mass: " << mass << " " << "Initial pos: " << position(0) << " " << position(1) <<
            " " << "Initial vel: " << velocity(0) << " " << velocity(1) << endl;

    CelestialBody body(nameOfBody, position, velocity, mass);

    if (locked) {
        body.setHasLockedPos(locked);
    }

    addCelestialBody(body);

    return;

}

vec SolarSystem::getForces(CelestialBody body) {
    vec force;
    //double forceX;
    //double forceY;
    // resetting force for each iteration of each body
    force.zeros(2);

    for (int i = 0; i < this->getNuOfBodies(); i++) {
        // Avoid calculating force on itself
        if (!(strcmp(bodies[i].getName().c_str(), body.getName().c_str()) == 0)) {

            double dotProduct = sqrt(dot(body.getDistanceBetweenBodies(bodies[i]),
                                   body.getDistanceBetweenBodies(bodies[i])));

            force -= 4*M_PI*M_PI*(bodies[i].getMass() * body.getMass()) *
                     body.getDistanceBetweenBodies(bodies[i]) / pow(dotProduct,3);
            /*forceX -= 4*M_PI*M_PI*(bodies[i].getMass() * body.getMass() /
                     pow(body.getDistanceBetweenBodies(bodies[i])[0],2));
            forceY -= 4*M_PI*M_PI*(bodies[i].getMass() * body.getMass() /
                     pow(body.getDistanceBetweenBodies(bodies[i])[1],2));
            force[0] = forceX;
            force[1] = forceY;
            */
        }

    }
    return force;

}

mat SolarSystem::getAllCurrentPos() {
    mat p = zeros<mat>(getNuOfDimensions(), getNuOfBodies());

    for (int i = 0; i < getNuOfBodies(); i++) {
        p.col(i) = bodies[i].getCurrentPos();
    }
    return p;
}

mat SolarSystem::getAllCurrentVel() {
    mat v = zeros<mat>(getNuOfDimensions(), getNuOfBodies());

    for (int i = 0; i < getNuOfBodies(); i++) {
        v.col(i) = bodies[i].getCurrentVel();
    }
    return v;
}

mat SolarSystem::getAllCurrentAcc() {
    mat a = zeros<mat>(getNuOfDimensions(), getNuOfBodies());

    for (int i = 0; i < getNuOfBodies(); i++) {
        a.col(i) = getForces(bodies[i]) / bodies[i].getMass();
    }
    return a;
}

void SolarSystem::setAllPos(mat positions) {
    for (int i = 0; i < getNuOfBodies(); i++) {
        bodies[i].setPos(positions.col(i));
    }
}

void SolarSystem::setAllVel(mat velocities) {
    for (int i = 0; i < getNuOfBodies(); i++) {
        bodies[i].setVel(velocities.col(i));
    }
}

void SolarSystem::calculateEnergy() {
    // resetting the values for each iteration
    kineticEnergy = 0;
    potentialEnergy = 0;
    angularMomentum.zeros();

    for (int i = 0; i < getNuOfBodies(); i++) {
        CelestialBody &body1 = bodies[i];
        vec vVector = body1.getCurrentVel();
        double v = norm(vVector);
        kineticEnergy += 0.5 * body1.getMass() * pow(v,2);
        for (int j = i+1; j < getNuOfBodies(); j++) {
            CelestialBody &body2 = bodies[j];
            vec rVector = body1.currentPos - body2.currentPos;
            double r = norm(rVector);
            potentialEnergy -= (body1.getMass() * body2.getMass()) / r;
        }
    }
}

double SolarSystem::totalEnergy() {
    return kineticEnergy + potentialEnergy;
}

void SolarSystem::rk4Solver() {
    int n = getNuOfBodies();

    mat pCurrRk4 = zeros<mat>(2,n);

    mat ak1 = zeros<mat>(2,n);
    mat ak2 = zeros<mat>(2,n);
    mat ak3 = zeros<mat>(2,n);
    mat ak4 = zeros<mat>(2,n);

    mat vk1 = zeros<mat>(2,n);
    mat vk2 = zeros<mat>(2,n);
    mat vk3 = zeros<mat>(2,n);
    mat vk4 = zeros<mat>(2,n);

    pCurrRk4 = getAllCurrentPos();

    ak1 = getAllCurrentAcc();   //k1
    vk1 = getAllCurrentVel();

    setAllVel(vk1 + 0.5 * dt * ak1);    // update
    setAllPos(pCurrRk4 + 0.5 * dt * vk1);

    ak2 = getAllCurrentAcc();   // k2
    vk2 = getAllCurrentVel();

    setAllVel(vk1 + 0.5 * dt * ak2);    //update
    setAllPos(pCurrRk4 + 0.5 * dt * vk2);

    ak3 = getAllCurrentAcc();   //k3
    vk3 = getAllCurrentVel();

    setAllVel(vk1 + dt * ak3);  //update
    setAllPos(pCurrRk4 + dt * vk3);

    ak4 = getAllCurrentAcc();   //k4
    vk4 = getAllCurrentVel();

    setAllVel(vk1 + (1./6) * dt * (ak1 + 2*(ak2 + ak3) + ak4));     //final step
    setAllPos(pCurrRk4 + (1./6) * dt * (vk1 + 2*(vk2 + vk3) + vk4));

//    setAllVel(vk1 + dt * ak1);
//    setAllPos(pCurrRk4 + dt * vk1);

}

void SolarSystem::verletFirstStep() {
    mat pPrev = zeros<mat>(2,getNuOfBodies());
    mat pCurr = getAllCurrentPos();
    mat vCurr = getAllCurrentVel();
    //mat aCurr = getAllCurrentAcc();

    pCurr = getAllCurrentPos();
    vCurr = getAllCurrentVel();
    //aCurr = getAllCurrentAcc();

    //pPrev = pCurr - dt *vCurr + 0.5 * pow(dt,2) * pow(aCurr,2);
    pPrev = pCurr - vCurr * dt;
}

void SolarSystem::verletSolver() {
    //mat pPrevVerlet = prevStep;
    mat pCurrVerlet = getAllCurrentPos();
    mat aCurrVerlet = getAllCurrentAcc();
    mat pNew, vNew;

    pNew = 2 * pCurrVerlet - pPrev + pow(dt,2) * (aCurrVerlet % aCurrVerlet);
    vNew = (pNew - pPrev) / (2*dt);

    setAllPos(pNew);
    setAllVel(vNew);

}

int SolarSystem::getNuOfBodies() {
    return this->bodies.size();
}

