#include "atom.h"
#include "random.h"
//#include <cmath>
#include <math.h>

Atom::Atom(double mass) :
    m_mass(mass)
{

}


void Atom::setInitialPosition(double x, double y)
{
    m_initial_position.set(x,y);
    position.set(x,y);
}

void Atom::resetForce()
{
    force.zeros();
}

void Atom::resetVelocityMaxwellian(double temperature)
{
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
    double standardDeviation =  sqrt(boltzmannConstant*temperature/m_mass);
    velocity.randomGaussian(0, standardDeviation);  //note: randomGaussian is defined in vec2.cpp. arguments(mean, stdev), so this is Gaussian centered about 0
}
