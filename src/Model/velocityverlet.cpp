#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "send_atoms.h"
#include <mpi.h>
#include "global.h"
#include "unitconverter.h"


void VelocityVerlet::integrate(System &system, double dt) //passing by reference &system
{


    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    double half_dt = 0.5*dt;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    //std::cout <<"in VV" <<std::endl;

    for(Atom *atom : system.atoms()) {
        //this operates on the vectors directly using vec2 class
        atom->velocity += atom->force*half_dt/mass;
        atom->position += atom->velocity*dt;  //NOTE: since v is computed 1st, this is actually x(t+dt) = x(t) + vt+0.5at^2 since v = v+0.5at
    }

    system.applyMirrorBCs(dt);
    //system.applyPeriodicBoundaryConditions();

    if(nprocs >1){
       send_atoms(&system);  //send and recieve atoms which have left their processor's domain
    }

    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*half_dt/mass;  //calculate new velocities
    }
}
