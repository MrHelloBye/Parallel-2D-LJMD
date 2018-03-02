#include "lennardjones.h"
#include "system.h"
#include <cmath>
#include <mpi.h>
#include "atom.h"
#include "extpotential.h"

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

    //vec2 system.systemSize() = system.systemSize(); //returns size of GLOBAL system box
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

    //store atom pointers for use in current processor since only these atoms within skin region will interact with neighboring ghosts
    std::vector<Atom*> Atomsto_left;
    std::vector<Atom*> Atomsto_right;

    MPI_Request req[4], req2[4];
    MPI_Status stat[4], stat2[4];

    double l_skin_bndry =(rank-1) * system.systemSize(decomp_dim)/ (nprocs-1) + skin_cutoff;  //boundary of the skin at left of subdomain
    double r_skin_bndry = (rank ) * system.systemSize(decomp_dim) / (nprocs-1) - skin_cutoff;
    
   //First calculate interactions btw. all atoms in the processor's domain
    for(int current_index=0; current_index<system.num_atoms(); current_index++){
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

            if(radiusSqrd > skin_cutoff_sqrd ) continue;  //cutoff radius and if 2 particles are too close, don't compute the force to prevent blowup

            //double radius = sqrt(radiusSqrd);
            //double sigma_over_radius = m_sigma/radius;

            //find and set force components
            double total_force_over_r = 24.*(2.0*pow(radiusSqrd,-7.)-pow(radiusSqrd,-4.));
            //double total_force_over_r = 24.*(2.0*pow(radius,-14.)-pow(radius,-8.));
            //ATTRACTIVE FORCE SHOULD POINT TOWARDS OTHER ATOM. REPULSIVE AWAY FROM OTHER ATOM!!
            for(int j=0;j<2;j++) {
                current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                other_atom->force[j] -= total_force_over_r*displacement[j]; //using Newton's 3rd law  NOTE: can't use current atom force b/c that is additive!
            }

            if(system.steps() % system.m_sample_freq ==0){
                //calculate potential energy every m_sample_freq steps
                m_potentialEnergy += 4.*(pow(radiusSqrd,-6.)-pow(radiusSqrd,-3));
            }
        }//end of inner loop

        //add force due to external potential
        current_atom->force += system.extPotential.getForcefromPotential(current_atom->position, skin_cutoff_sqrd);

        //std::cout << "force form ext potential" << system.extPotential.getForcefromPotential(current_atom->position, skin_cutoff_sqrd) <<std::endl;

        
        //remember if atom needs to be considered as ghost atom for neighboring processors
        if (current_atom->position[decomp_dim] < l_skin_bndry) {
            //std::cout << "test if" << rank * system.systemSize()[decomp_dim] / nprocs + skin_cutoff <<std::endl;
            //special case for ghosts from leftmost proc, come from right side--> must add systemsize from x positions
            if(rank == 1){
                to_left.push_back(current_atom->position[0] + system.systemSize(0));  //verified that this does give the correct numbers
                //std::cout <<system.atoms(current_index)->position[0] + system.systemSize()[0] <<std::endl;
            }else{
            to_left.push_back(current_atom->position[0]);
            }

            to_left.push_back(current_atom->position[1]);
            //note: don't need velocities for ghost atoms

            //also save vectors of pointers to atom for using in current processor for more efficient ghost interaction calculation
            Atomsto_left.push_back(current_atom);  //this positions don't need to be changed b/c are for current proc
            num_to_left++;
        }
        else if  (current_atom->position[decomp_dim] >r_skin_bndry) {
            //special case for ghosts from rightmost proc, must substract systemsize from x positions
            if(rank == nprocs-1){//rightmost proc
                to_right.push_back(current_atom->position[0] - system.systemSize(0));  //THIS DOESN'T GET TRIGGERED
                //std::cout <<system.atoms(current_index)->position[0] - system.systemSize()[0] <<std::endl;
            }else{
                 to_right.push_back(current_atom->position[0]);
            }
            to_right.push_back(current_atom->position[1]);
            Atomsto_right.push_back(current_atom);
            num_to_right++;

        }
    }//end of outer loop

    // std::cout <<"ghost position BEFORE send to right" << to_right[0] <<" " << to_right[1] <<"proc" << rank <<std::endl;
    // std::cout <<"ghost position BEFORE send to right" << to_right[2] <<" " <<to_right[3] <<std::endl;
    //std::cout <<"ghost position BEFORE send to left" <<to_left[0] <<" " <<to_left[1]  <<"proc" << rank <<std::endl;
    // std::cout <<"ghost position BEFORE send to left" <<to_left[2]<<" " <<to_left[3] <<std::endl;
    //std::cout << "num to right" <<num_to_right <<"proc" << rank <<std::endl;
     //std::cout << "num to left" <<num_to_left <<"proc" << rank <<std::endl;

   // std::cout <<"system size" <<system.systemSize()[0] <<std::endl;


    //send ghost atoms
     if (nprocs > 1) {

         int ln =  (rank -1- 1 + nprocs-1) % (nprocs-1)+1;  //left neighbor
         int rn = (rank ) %( nprocs-1)+1;                   //right neighbor

         // Send ghost atoms to neighboring processor
         MPI_Isend(&num_to_left, 1, MPI_INT, ln, 10*rank + ln, MPI_COMM_WORLD, req); //note: these formulas will pass i.e. from proc 0 to proc (max at right)...--> satisfy PBCs...
         MPI_Irecv(&num_from_left, 1, MPI_INT, ln, 10*ln +rank, MPI_COMM_WORLD, req+1);
         MPI_Isend(&num_to_right, 1, MPI_INT, rn, 10*rank + rn, MPI_COMM_WORLD, req+2);
         MPI_Irecv(&num_from_right, 1, MPI_INT, rn, 10*rn +rank, MPI_COMM_WORLD, req+3);
         MPI_Waitall (4, req, stat);

         //resize to fit x and y positions
         from_left.resize(2*num_from_left);
         from_right.resize(2*num_from_right);
         //verified that these arrays do get filled


         MPI_Isend(&to_left[0], to_left.size(), MPI_DOUBLE, ln, 10*rank + ln, MPI_COMM_WORLD, req2);
         MPI_Irecv(&from_left[0], from_left.size(), MPI_DOUBLE, ln, 10*ln +rank, MPI_COMM_WORLD, req2+1);
         MPI_Isend(&to_right[0], to_right.size(), MPI_DOUBLE, rn, 10*rank + rn, MPI_COMM_WORLD, req2+2);
         MPI_Irecv(&from_right[0], from_right.size(), MPI_DOUBLE, rn, 10*rn +rank, MPI_COMM_WORLD, req2+3);
         MPI_Waitall (4, req2, stat2);


         // Calculate interactions between new (ghost) atoms and atoms on processor--> NOTE: need to only interact with atoms within the skin b/c of cutoff radius
         //So i.e. atoms which were sent to left interact with atoms which came from left

         //forces with atoms which came from left
         for (int current_index=0; current_index<system.num_atoms(); current_index++){
             Atom *current_atom = system.atoms(current_index);
             if(current_atom->position[0] < l_skin_bndry){
                 //for (int current_index=0; current_index<num_to_left; current_index++) { //need real current_atom objects to be able to update their forces
                 //Atom *current_atom = system.atoms(current_index);////*curr... means pointer to current_atom = ... RHS is a pointer //these are addresses...

                 for (int ghost_index=0; ghost_index < 2*num_from_left; ghost_index+=2) {
                     //creating new atoms is VERY SLOW and not necessary--> we just need positions...so just use them directly in displacement

                     if(ghost_index <2 && current_index ==0){
                         //std::cout <<"ghost position AFTER send" << from_left[0] << " " << from_left[0+1] << "proc" << rank <<std::endl;
                         // std::cout <<"ghost position AFTER send" << from_left[0] << " " << from_left[1+1] <<std::endl;
                     }

                     vec2 displacement(0.,0.);
                     for(int j=0;j<2;j++){
                         //displacement[j] = current_atom->position[j] - ghost_atom->position[j];
                         displacement[j] = current_atom->position[j] - from_left[ghost_index + j];  //use ghost atoms positions data directly
                         // displacement[j] = Atomsto_left[current_index]->position[j] - from_left[ghost_index + j];
                     }

                     //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!--> .i.e. think at the corners along skin
                     if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                     if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1); //verified that halfsystemsize and systemsize is correct

                     double radiusSqrd = displacement.lengthSquared();

                     if(radiusSqrd > skin_cutoff_sqrd ) continue;

                     //double radius = sqrt(radiusSqrd);
                     //double sigma_over_radius = m_sigma/radius;
                     double total_force_over_r = 24.*(2.0*pow(radiusSqrd,-7.)-pow(radiusSqrd,-4.));

                     //find and set force components
                     for(int j=0;j<2;j++) {
                         current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                         //Atomsto_left[current_index]->force[j] += total_force_over_r*displacement[j];
                         //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                     }

                     if(system.steps() % system.m_sample_freq ==0) m_potentialEnergy += 4.*(pow(radiusSqrd,-6.)-pow(radiusSqrd,-3));
                 }
             }else if(current_atom->position[decomp_dim] > r_skin_bndry){

                 //forces with atoms which came from right
                 // for (int current_index=0; current_index<system.atoms(current_index); current_index++) { //need real current_atom objects to be able to update their forces
                 for (int ghost_index=0; ghost_index < 2*num_from_right; ghost_index+=2) {

                     if(ghost_index <2 && current_index ==0){
                         //std::cout <<"ghost position AFTER send" << from_right[0] << " " << from_right[0+1] << "proc" << rank <<std::endl;
                         // std::cout <<"ghost position AFTER send" << from_right[1] << " " << from_right[1+1] <<std::endl;
                     }
                     vec2 displacement(0.,0.);
                     for(int j=0;j<2;j++){
                         displacement[j] = current_atom->position[j] - from_right[ghost_index + j];  //use ghost atoms positions data directly
                         //displacement[j] = Atomsto_right[current_index]->position[j] - from_right[ghost_index + j];  //use ghost atoms positions data directly
                     }
                     //NOTE: I need minimum image convention along y-direction!--> b/c no passing of atoms/ghosts etc there!
                     if (displacement[1] >  system.halfsystemSize(1)) displacement[1] -= system.systemSize(1);   //systemSize(j) returns m_systemSize[j] from system class
                     if (displacement[1] <= -system.halfsystemSize(1)) displacement[1] += system.systemSize(1);

                     double radiusSqrd = displacement.lengthSquared();
                     if(radiusSqrd > skin_cutoff_sqrd ) continue;

                     //double radius = sqrt(radiusSqrd);
                     //double sigma_over_radius = m_sigma/radius;
                     double total_force_over_r = 24.*(2.0*pow(radiusSqrd,-7.)-pow(radiusSqrd,-4.));

                     //find and set force components
                     for(int j=0;j<2;j++) {
                         current_atom->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                         //Atomsto_right[current_index]->force[j] += total_force_over_r*displacement[j]; //i.e. Fx = (F/r)*x
                         //WE DO NOT want to update the force on other_atom here, b/c it is a ghost atom...--> outside of the system
                     }

                     if(system.steps() % system.m_sample_freq ==0) m_potentialEnergy += 4.*(pow(radiusSqrd,-6.)-pow(radiusSqrd,-3));
                 }
             }
         }//end of outer loop

         //clear the vectors, not clearing these doesn't seem to improve timing
         from_left.clear();
         from_right.clear();
         to_left.clear();
         to_right.clear();

         Atomsto_right.clear();
         Atomsto_left.clear();


         //MPI_Barrier(MPI_COMM_WORLD);
     }//end of if(nprocs>1)
}

