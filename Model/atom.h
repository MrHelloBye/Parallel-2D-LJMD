#ifndef ATOM_H
#define ATOM_H
#include "vec2.h"

class Atom
{
private:
    double m_mass;
    vec2 m_initial_position;
public:
    vec2 position;
    vec2 velocity;
    vec2 force;
    vec2 num_bndry_crossings;  //this is to keep track of boundary crossings for calculating diffusion coeff.

    Atom(double mass);
    void setInitialPosition(double x, double y);
    void resetForce();
    void resetVelocityMaxwellian(double temperature);
    double mass() { return m_mass; }
    void setMass(double mass) { m_mass = mass; }
    double initial_position(int j){return m_initial_position[j];}
};
#endif
