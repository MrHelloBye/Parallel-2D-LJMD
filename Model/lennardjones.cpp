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
    m_twntyfour_epsilon = 24.0*m_epsilon;
}

void LennardJones::calculateForces(System &system)  //object system is passed by reference (allows changing)
{
    m_potentialEnergy = 0;  //reset potential energy

    const double skin_cutoff = 3.*m_sigma;
    const double skin_cutoff_sqrd = skin_cutoff*skin_cutoff;

    vec2 sys_size = system.systemSize(); //returns size of GLOBAL system box
    int decomp_dim = 0;  // 0 or 1, x or y direction of decomposition

    int nprocs, rank;
    MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    // Store number of atoms to send to and receive from the processor on the left and on the right
    int num_to_left = 0;
    int num_to_right = 0;
    int num_from_left, num_from_right;

    // Store atoms to send and receive
    std::vector<double> to_left;
    std::vector<double> to_right;
    std::vector<double> from_left;
    std::vector<double> from_right;

    MPI_Request req[4], req2[4];
    MPI_Status stat[4], stat2[4];
    
    
   //First calculate interactions btw. all atoms in the processor's domain
    for(int current_index=0; current_index<system.num_atoms()-1; current_index++){  //-1 b/c don't need to calculate pairs when get to last atom
        Atom *current_atom = system.atoms(current_index);  //system.atoms(index) returns the pointer to the atom corresponding to index

        for(int other_index=current_index+1;other_index<system.num_atoms();other_index++){   //to avoid double counting
            Atom *other_atom = system.atoms(other_index);

            //distance and vector btw objects dx
            vec2 displacement(0.,0.);
            for(int j=0;j<2;j++){
                displacement[j] = current_atom->position[j] - other_atom->position[j];
            }
            //NOTE: For 1D decomposition, I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!
            if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
            if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

            //for case of 1 processor, need to implement x mirror image convention also
            if(nprocs ==1){
                //for cases where the folded back particle will be closer than its image to a given particle
                if (displacement[0] >  system.halfsystemSize(0)) displacement[0] -= system.systemSize(0);   //systemSize(j) returns m_systemSize[j] from system class
                if (displacement[0] <= -system.halfsystemSize(0)) displacement[0] += system.systemSize(0);
            }

            double radiusSqrd = displacement.lengthSquared();
            //if(radiusSqrd > skin_cutoff_sqrd ) continue;

            if(radiusSqrd > skin_cutoff_sqrd ) continue;  //cutoff radius and if 2 particles are too close, don't compute the force to prevent blowup

            double radius = sqrt(radiusSqrd);
            double sigma_over_radius = m_sigma/radius;

            //find and set force components
            double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
            //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!
            for(int j=0;j<2;j++) {
                current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                other_atom->force[j] -= total_force_over_r*displacement[j]; //using Newton's 3rd law  NOTE: can't use current atom force b/c that is additive!
            }


            if(system.steps() % system.m_sample_freq ==0){
                //calculate potential energy every m_sample_freq steps
                m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
            }
        }//end of inner loop
        
        //remember if atom needs to be considered as ghost atom for neighboring processors
        // special case for nprocs==2; don't want to send same atom twise
        if (nprocs == 2) {
            //std::cout << "check positions" << system.atoms(current_index)->position[0] << std::endl;

            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                //std::cout << "test if" << rank * sys_size[decomp_dim] / nprocs + skin_cutoff <<std::endl;
                to_left.push_back(system.atoms(current_index)->position[0]);
                to_left.push_back(system.atoms(current_index)->position[1]);
                //note: don't need velocities for ghost atoms
                num_to_left++;
            }
            else if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
               // std::cout << "test if" << (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff <<std::endl;
                to_right.push_back(system.atoms(current_index)->position[0]);
                to_right.push_back(system.atoms(current_index)->position[1]);
                num_to_right++;
            }
        } else {
            if (system.atoms(current_index)->position[decomp_dim] < rank * sys_size[decomp_dim] / nprocs + skin_cutoff) {
                to_left.push_back(system.atoms(current_index)->position[0]);
                to_left.push_back(system.atoms(current_index)->position[1]);
                num_to_left++;
            }
            if  (system.atoms(current_index)->position[decomp_dim] > (rank + 1) * sys_size[decomp_dim] / nprocs - skin_cutoff) {
                to_right.push_back(system.atoms(current_index)->position[0]);
                to_right.push_back(system.atoms(current_index)->position[1]);
                num_to_right++;
            }
        }

        if(system.steps() % 100 ==0){
            //std::cout <<"forces before add ghosts contribution" << current_atom->force[0] << " " << current_atom->force[1] <<std::endl;
        }
            
    }//end of outer loop

    // std::cout <<"ghost position BEFORE send to right" << to_right[0] <<" " << to_right[1] <<"proc" << rank <<std::endl;
    // std::cout <<"ghost position BEFORE send to right" << to_right[2] <<" " <<to_right[3] <<std::endl;
    //std::cout <<"ghost position BEFORE send to left" <<to_left[0] <<" " <<to_left[1]  <<"proc" << rank <<std::endl;
    // std::cout <<"ghost position BEFORE send to left" <<to_left[2]<<" " <<to_left[3] <<std::endl;
    // std::cout << "num to right" <<num_to_right <<"proc" << rank <<std::endl;
    // std::cout << "num to left" <<num_to_left <<"proc" << rank <<std::endl;


    //send ghost atoms
     if (nprocs > 1) {
         // Send ghost atoms to neighboring processor
         MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req); //note: these formulas will pass i.e. from proc 0 to proc (max at right)...--> satisfy PBCs...
         MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);
         MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
         MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
         MPI_Waitall (4, req, stat);

         //resize to fit x and y positions
         from_left.resize(2*num_from_left);
         from_right.resize(2*num_from_right);
         //verified that these arrays do get filled

         //std::cout <<"num_from_left for ghosts" << num_from_left <<std::endl;

         MPI_Isend(&to_left[0], to_left.size(), MPI_DOUBLE, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
         MPI_Irecv(&from_left[0], from_left.size(), MPI_DOUBLE, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
         MPI_Isend(&to_right[0], to_right.size(), MPI_DOUBLE, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
         MPI_Irecv(&from_right[0], from_right.size(), MPI_DOUBLE, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
         MPI_Waitall (4, req2, stat2);

         // Calculate interactions between new (ghost) atoms and atoms on processor
         for (int current_index=0; current_index!=system.num_atoms(); ++current_index) {
             Atom *current_atom = system.atoms(current_index); //*curr... means pointer to current_atom = ... RHS is a pointer

             //forces with atoms which came from left
             for (int ghost_index=0; ghost_index < num_from_left; ++ghost_index) {
                 int index = 2*ghost_index;  //for unpacking recieved data
                 //creating new atoms is VERY SLOW and not necessary--> we just need positions...so just use them directly in displacement

                 if(ghost_index <2 && current_index ==0){
                     //std::cout <<"ghost position AFTER send" << from_left[0] << " " << from_left[0+1] << "proc" << rank <<std::endl;
                     // std::cout <<"ghost position AFTER send" << from_left[0] << " " << from_left[1+1] <<std::endl;
                 }


                 vec2 displacement(0.,0.);
                 for(int j=0;j<2;j++){
                     //displacement[j] = current_atom->position[j] - ghost_atom->position[j];
                     displacement[j] = current_atom->position[j] - from_left[index + j];  //use ghost atoms positions data directly
                 }

                 //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!--> .i.e. think at the corners along skin
                 if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                 if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1); //verified that halfsystemsize and systemsize is correct

                 double radiusSqrd = displacement.lengthSquared();

                 if(radiusSqrd > skin_cutoff_sqrd ) continue;

                 double radius = sqrt(radiusSqrd);
                 double sigma_over_radius = m_sigma/radius;
                 double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));

                 //find and set force components
                 for(int j=0;j<2;j++) {
                     current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                     //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                 }

                 if(system.steps() % system.m_sample_freq ==0) m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
             }


             //forces with atoms which came from right
             for (int ghost_index=0; ghost_index < num_from_right; ++ghost_index) {
                 int index = 2*ghost_index;  //for unpacking recieved data

                 if(ghost_index <2 && current_index ==0){
                     //std::cout <<"ghost position AFTER send" << from_right[0] << " " << from_right[0+1] << "proc" << rank <<std::endl;
                     // std::cout <<"ghost position AFTER send" << from_right[1] << " " << from_right[1+1] <<std::endl;
                 }
                 vec2 displacement(0.,0.);
                 for(int j=0;j<2;j++){
                     displacement[j] = current_atom->position[j] - from_right[index + j];  //use ghost atoms positions data directly
                 }
                 //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!
                 if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                 if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

                 double radiusSqrd = displacement.lengthSquared();
                 if(radiusSqrd > skin_cutoff_sqrd ) continue;

                 double radius = sqrt(radiusSqrd);
                 double sigma_over_radius = m_sigma/radius;
                 double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));

                 //find and set force components
                 for(int j=0;j<2;j++) {
                     current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                     //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                 }
                 if(system.steps() % 100 ==0){
                     //std::cout << "force after ghosts contribution" << current_atom->force[0] << " " <<current_atom->force[1] << std::endl;
                 }
                 if(system.steps() % system.m_sample_freq ==0) m_potentialEnergy += 4.*(pow(radius,-12.)-pow(radius,-6));
             }
         }

                //clear the vectors
                from_left.clear();
                from_right.clear();
                to_left.clear();
                to_right.clear();

                MPI_Barrier(MPI_COMM_WORLD);
        }//end of if(nprocs>1)
}

