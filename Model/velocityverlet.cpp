#include "velocityverlet.h"
#include "system.h"
#include "atom.h"
#include "send_atoms.h"
#include <mpi.h>


void VelocityVerlet::integrate(System &system, double dt) //passing by reference &system, passes using the address but doesn't create a variable (pointer) whose value = that address
{

    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors


    double half_dt = 0.5*dt;
    if(m_firstStep) {
        system.calculateForces();
        m_firstStep = false;
    }

    for(Atom *atom : system.atoms()) {
        //this operates on the vectors directly using vec3 class
        atom->velocity += atom->force*half_dt/atom->mass();
        atom->position += atom->velocity*dt;  //NOTE: since v is computed 1st, this is actually x(t+dt) = x(t) + vt+0.5at^2 since v = v+0.5at
        //std::cout<< atom->m_initial_position[0] << "atom INITIAL position" <<std::endl;  //initial positions seem fine...
        //std::cout << atom->position[0] <<"atomposition in VV" <<std::endl;  //these positions are blown up...
         //std::cout << atom->force[0] << "atom force" <<std::endl; //force is blown up

    }

    system.applyPeriodicBoundaryConditions();
   // system.applyMirrorBCs_inX(dt);


    //send_atoms accepts a pointer, so declare a pointer(variable whose value = the memory address)
    //System  * pt_system = &system; //assign it the address of system

    if(nprocs >1){
       send_atoms(&system);  //send and recieve atoms which have left their processor's domain
    }


    system.calculateForces(); // New positions, recompute forces

    for(Atom *atom : system.atoms()) {
        atom->velocity += atom->force*half_dt/atom->mass();  //calculate new velocities
    }
}
