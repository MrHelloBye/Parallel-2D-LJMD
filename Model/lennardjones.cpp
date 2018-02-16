#include "lennardjones.h"
#include "system.h"
#include <cmath>

double LennardJones::potentialEnergy() const
{
    return m_potentialEnergy;
}

double LennardJones::sigma() const
{
    return m_sigma;
}

void LennardJones::setSigma(double sigma)
{
    m_sigma = sigma;
    m_sigma_sqrd = sigma*sigma;
}

double LennardJones::epsilon() const
{
    return m_epsilon;
}


void LennardJones::setEpsilon(double epsilon)
{
    m_epsilon = epsilon;
    m_four_epsilon = 4.0*m_epsilon;
    m_twntyfour_epsilon = 24.0*m_epsilon;  //must reset this too
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
 m_potentialEnergy = 0;  //reset potential energy

 for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
    Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index

    for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
        Atom *other_atom = system.atoms(other_index);

        //distance and vector btw objects dx should obey min. image convention: only closest distance to particle or its image is considered
        vec2 displacement(0.,0.);
        for(int j=0;j<2;j++){
             displacement[j] = current_atom->position[j] - other_atom->position[j];
             //for cases where the folded back particle will be closer than its image to a given particle
             if (displacement[j] >  system.halfsystemSize(j)) displacement[j] -= system.systemSize(j);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[j] <= -system.halfsystemSize(j)) displacement[j] += system.systemSize(j);
         }

        double radiusSqrd = displacement.lengthSquared();


        //apply interaction range cutoff
        //std::cout<<radius<<std::endl;
       // double radius = sqrt(radiusSqrd);
        //if(radius > 5.0*m_sigma) continue;
        //if(radiusSqrd > 16.0*m_sigma_sqrd) continue;

        double radius = sqrt(radiusSqrd);
        double sigma_over_radius = m_sigma/radius;


       //double total_force_over_r = (m_twntyfour_epsilon/pow(m_sigma,2.))*(2.0*pow(sigma_over_radius,14.)-pow(sigma_over_radius,8.));

        //below version is if make in unit converter L0 = sigma*1e-10;
        //note: epsilon was set to 1, so coeff. in front is just 24
        double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
        //double total_force_over_r = m_twntyfour_epsilon*(2.0*pow(radius,-14.)-pow(radius,-8.));
        //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!


        //find and set force components
        //double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
        for(int j=0;j<2;j++) {
            current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
            other_atom->force[j] -= current_atom->force[j]; //using Newton's 3rd law
        }

        if(system.steps() % system.m_sample_freq ==0){
            //calculate potential energy every m_sample_freq steps
           //m_potentialEnergy += m_four_epsilon*(pow(sigma_over_radius,12.)-pow(sigma_over_radius,6));
            //below version is if make in unit converter L0 = sigma*1e-10;
            //since epsilon is set to 1
            m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
           //m_potentialEnergy += m_four_epsilon*(pow(radius,-12.)-pow(radius,-6));
        }
    }
 }


}



