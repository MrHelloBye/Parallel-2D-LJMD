#include <iostream>
#include "extpotential.h"
#include "vec2.h"
#include <cmath>
#include "system.h"



void ExtPotential::setPosition(double x, double y){
   position.set(x,y);
}

void ExtPotential::setStdev(double input_stdev){
    stdev = input_stdev;
    twostdevSqrd = 2*stdev*stdev;
}

void ExtPotential::setMax(double inputmax){
    max = inputmax;
}


vec2 ExtPotential::getForcefromPotential(vec2 atomPosition, double skin_cutoff_sqrd){
    vec2 displacement = position - atomPosition;
    vec2 force;

    double radiusSqrd = displacement.lengthSquared();

    if(radiusSqrd > 2*skin_cutoff_sqrd ) force.set(0,0);       //  NOTE: if use only skin_cutoffsqrd = 3 sigma for cutoff, then blows up after a bit (i.e. 30k steps)
    else force = (max*exp(-radiusSqrd/twostdevSqrd)/radiusSqrd)*displacement;  //is positive b/c assume repulsive force

    return force;

}

