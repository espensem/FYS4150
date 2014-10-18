#include "verlet.h"

Verlet::Verlet(mat ppos, mat pos, mat vel, mat acc, double dt) {

    this->prevPos = ppos;
    this->currentPos = pos;
    this->currentVel = vel;
    this->currentAcc = acc;
    this->dt = dt;
}

mat Verlet::advancePos() {
    this->newPos = this->newPos = 2 * this->currentPos - this->prevPos + pow(this->dt,2)*pow(this->currentAcc,2);
    return this->newPos;
}

mat Verlet::advanceVel() {
    this->newVel = this->newVel = (this->newPos - this->prevPos) / (2 * dt);
    return this->newVel;
}
