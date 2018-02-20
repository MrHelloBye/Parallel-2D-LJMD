#include "lennardjones.h"
#include "system.h"
#include <cmath>
#include <mpi.h>

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

    const double skin_cutoff = 8.*m_sigma;
    const double skin_cutoff_sqrd = skin_cutoff*skin_cutoff;
    const double too_close_sqrd = 0;//(0.5*m_sigma)*(0.5*m_sigma);  //%BE CAREFUL HERE -> if make too large, it will ignore important forces

    vec2 sys_size = system.subsystemSize(); //returns size of LOCAL processors system box
    int decomp_dim = 0;  // 0 or 1, x or y direction of decomposition

    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    // Store number of atoms to send to and receive from the processor on the left and on the right
    int num_to_left = 0;
    int num_to_right = 0;
    int num_from_left, num_from_right;

    // Store atoms to send and receive
    std::vector<Atom> to_left;
    std::vector<Atom> to_right;
    std::vector<Atom> from_left;
    std::vector<Atom> from_right;

    MPI_Request req[4], req2[4];
    MPI_Status stat[4], stat2[4];
    
    
   //First calculate interactions btw. all atoms in the system



    for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
        Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index

        //if(system.steps() > 2200){
        //    std::cout << "current atom force" << current_atom->force[0] <<" "<<  current_atom->force[1] << "atom index" << current_index <<"num_atoms" << system.num_atoms() << "proc" <<rank << std::endl;
       // }

        for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
            Atom *other_atom = system.atoms(other_index);


            //distance and vector btw objects dx
              vec2 displacement(0.,0.);
             for(int j=0;j<2;j++){
                 displacement[j] = current_atom->position[j] - other_atom->position[j];
             }
             //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!
             if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
             if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

             //for case of 1 processor, need to implement x mirror image convention also
             if(nprocs ==1){
                 //for cases where the folded back particle will be closer than its image to a given particle
                 if (displacement[0] >  system.halfsystemSize(0)) displacement[0] -= system.systemSize(0);   //systemSize(j) returns m_systemSize[j] from system class
                 if (displacement[0] <= -system.halfsystemSize(0)) displacement[0] += system.systemSize(0);
             }


            // std::cout << displacement[0] <<"displacement x" <<displacement[1] <<"displacement y" <<std::endl;
             //std::cout<<"numatoms" <<system.num_atoms() <<std::endl;

             //displacements are fine

            double radiusSqrd = displacement.lengthSquared();
            //if(radiusSqrd > skin_cutoff_sqrd ) continue;

            if(radiusSqrd > skin_cutoff_sqrd || radiusSqrd < too_close_sqrd) continue;  //cutoff radius and if 2 particles are too close, don't compute the force to prevent blowup

            double radius = sqrt(radiusSqrd);
            double sigma_over_radius = m_sigma/radius;

            //std::cout <<sigma_over_radius <<"sigmaoverradius" << std::endl;


            double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
            //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!!
            //std::cout << total_force_over_r <<"force" <<std::endl;

            //find and set force components
            //double total_force_over_r = total_force/radius; //precalculate to save 2 FLOPS
            for(int j=0;j<2;j++) {
                current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                other_atom->force[j] -= current_atom->force[j]; //using Newton's 3rd law
            }
            
              //force is blown up here
            //std::cout << current_atom->force[0] <<"force in lj LINE 99" <<std::endl;

            

             if(system.steps() % system.m_sample_freq ==0){
                 //calculate potential energy every m_sample_freq steps
             //m_potentialEnergy += m_four_epsilon*(pow(sigma_over_radius,12.)-pow(sigma_over_radius,6));
                 //below version is if make in unit converter L0 = sigma*1e-10;
                 //since epsilon is set to 1
                 m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
             //m_potentialEnergy += m_four_epsilon*(pow(radius,-12.)-pow(radius,-6));
             }


        }//end of inner loop
        
        //remember if atom needs to be considered as ghost atom for neighboring processors
        // special case for nprocs==2; don't want to send same atom twise
        if (nprocs == 2) {
            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                to_left.push_back(*(system.atoms(current_index)));
                num_to_left++;
            }
            else if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
                to_right.push_back(*(system.atoms(current_index)));
                num_to_right++;
            }
        } else {
            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                to_left.push_back(*(system.atoms(current_index)));
                num_to_left++;
            }
            if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
                to_right.push_back(*(system.atoms(current_index)));
                num_to_right++;
            }
        }
            
    }//end of outer loop
    


    //send ghost atoms
     if (nprocs > 1) {
                // Send ghost atoms to neighboring processor
                MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req); //note: these formulas will pass i.e. from proc 0 to proc (max at right)...--> satisfy PBCs...
                MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);
                MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
                MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
                MPI_Waitall (4, req, stat);

                from_left.resize(num_from_left);
                from_right.resize(num_from_right);
                //verified that these arrays do get filled

                //std::cout <<"num_from_left for ghosts" << num_from_left <<std::endl;

                MPI_Isend(&to_left[0], num_to_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
                MPI_Irecv(&from_left[0], num_from_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
                MPI_Isend(&to_right[0], num_to_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
                MPI_Irecv(&from_right[0], num_from_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
                MPI_Waitall (4, req2, stat2);

                // Calculate interactions between new (ghost) atoms and atoms on processor
                for (int current_index=0; current_index!=system.num_atoms(); ++current_index) {
                    Atom *current_atom = system.atoms(current_index); //*curr... means pointer to current_atom = ... RHS is a pointer
                    //forces with atoms which came from left
                    for (int other_index=0; other_index < num_from_left; ++other_index) {
                        Atom *other_atom = &from_left[other_index];  //*other_atom = &.. means pointer to other_atom = address of from_left[]. other_atom itself is an object.
                        vec2 displacement(0.,0.);
                        for(int j=0;j<2;j++){
                            displacement[j] = current_atom->position[j] - other_atom->position[j];
                        }
                        //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!--> .i.e. think at the corners along skin
                        if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                        if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

                        double radiusSqrd = displacement.lengthSquared();
                        // if(radiusSqrd > skin_cutoff_sqrd ) continue;
                        if(radiusSqrd > skin_cutoff_sqrd || radiusSqrd < too_close_sqrd) continue;  //cutoff radius and if 2 particles are too close, don't compute the force to prevent blowup
                        double radius = sqrt(radiusSqrd);
                        double sigma_over_radius = m_sigma/radius;
                        double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));

                        //find and set force components
                        for(int j=0;j<2;j++) {
                            current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                            //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                        }


                    }
                    //forces with atoms which came from right
                    for (int other_index=0; other_index < num_from_right; ++other_index) {
                        Atom *other_atom = &from_right[other_index]; //other atom is from vector of atoms from right
                        vec2 displacement(0.,0.);
                        for(int j=0;j<2;j++){
                            displacement[j] = current_atom->position[j] - other_atom->position[j];
                        }
                        //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!
                        if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                        if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

                        double radiusSqrd = displacement.lengthSquared();
                       //  if(radiusSqrd > skin_cutoff_sqrd ) continue;
                         if(radiusSqrd > skin_cutoff_sqrd || radiusSqrd < too_close_sqrd) continue;  //cutoff radius and if 2 particles are too close, don't compute the force to prevent blowup
                        double radius = sqrt(radiusSqrd);
                        double sigma_over_radius = m_sigma/radius;
                        double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));

                        //find and set force components
                        for(int j=0;j<2;j++) {
                            current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                            //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                        }
                    }
                }

                //clear the vectors
                from_left.clear();
                from_right.clear();



        }//end of if(nprocs>1)

}



