#include "celestialbody.h"

using namespace std;
using namespace arma;

CelestialBody::CelestialBody(string nameOfBody, vec pos, vec vel, double mass) {
    this->nameOfBody = nameOfBody;
    this->initialPos = pos;
    this->currentPos = pos;
    this->currentVel = vel;
    this->mass = mass;

    this->hasLockedPos = false;
}

string CelestialBody::getName() {
    return this->nameOfBody;
}

double CelestialBody::getMass() {
    return this->mass;
}

vec CelestialBody::getInitialPos() {
    return this->initialPos;
}

vec CelestialBody::getCurrentPos() {
    return this->currentPos;
}

vec CelestialBody::getCurrentVel() {
    return this->currentVel;
}

vec CelestialBody::setPos(vec newPos) {
    if (this->hasLockedPos) {
        this->currentPos = this->initialPos;
    }
    else {
        this->currentPos = newPos;
    }
    return newPos;
}

vec CelestialBody::setVel(vec newVel) {
    this->currentVel = newVel;
    return newVel;
}

vec CelestialBody::getDistanceBetweenBodies(CelestialBody otherBody) {
    return (this->currentPos - otherBody.getCurrentPos());
}

bool CelestialBody::getHasLockedPos() {
    return this->hasLockedPos;
}

void CelestialBody::setHasLockedPos(bool locked) {
    this->hasLockedPos = locked;
    return;
}





