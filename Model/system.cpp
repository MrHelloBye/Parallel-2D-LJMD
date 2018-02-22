#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "random.h"
#include <map>
#include <mpi.h>
#include <algorithm> //need for sort etc.
#include "global.h"

System::System()   //constructor
{
 m_num_atoms = 0;
}

System::~System()  //desctructor: remove all the atoms
{
    for(Atom *atom : m_atoms) {
        delete atom;
    }
    m_atoms.clear();    //vector containing all atoms objects
}

void System::applyPeriodicBoundaryConditions() {
    // Read here: http://en.wikipedia.org/wiki/Periodic_boundary_conditions#Practical_implementation:_continuity_and_the_minimum_image_convention
    //using version A: where fold the atoms back into the simulation cell. The cell has its lower left corner at origin of coord. system.

    for(Atom *atom : atoms()) {

        for(int j=0;j<2;j++){
            //fold atoms back into the box if they escape the box
            //Note: if atom in 1 step moves further than the neighboring image of the simulation cell, it is not fully brought back.
            if (atom->position[j] <  0.) {
                atom->position[j] += m_systemSize[j];
                atom->num_bndry_crossings[j] -= 1;  //crossing left or bottom boundary is counted as -1 crossing.
            }
            if (atom->position[j] >  m_systemSize[j]){
                atom->position[j] -= m_systemSize[j];
                atom->num_bndry_crossings[j] += 1;    //crossing right or top boundary is counted as +1 crossing
            }
        }

    }
}

//
void System::applyMirrorBCs_inX(double dt){
    for(Atom *atom : atoms()) {

             if (atom->position[0] <  0. ) {
                 //we go back to the previous x position, and keep the y position
                 //also we reverse sign on X component of velocity
                 // std::cout << "test mirror bc's BEFORE mirroring" << atom->position[0] <<"time step" << steps() << std::endl;
                 //atom->position[0] -=dt*atom->velocity[0];
                 atom->position[0] = 0.01;  //POOR MAN'S mirror BCs
                // std::cout << "test mirror bc's" << atom->position[0] <<std::endl; // not working properly it seems--> after mirroring, stll have negative positoin...
                 atom->velocity[0] = - atom->velocity[0];

             }

             if(atom->position[0] >  m_systemSize[0]){
                 atom->position[0] = m_systemSize[0] -0.01;  //POOR MAN'S mirror BCs
                // std::cout << "test mirror bc's" << atom->position[0] <<std::endl; // not working properly it seems--> after mirroring, stll have negative positoin...
                 atom->velocity[0] = - atom->velocity[0];
             }
        }

}


// remove atoms that escape the SUB system (corresponding to current processor)
/*
void System::removeEscapedAtoms() {
    for(Atom *atom : atoms()) {
        if(atom->position[0]>m_subsystemSize[0])
            //m_atoms.erase(atom);
            delete atom;
    }
}
*/



void System::rescaleVelocities(StatisticsSampler &statisticsSampler, double currentTemperature, double desiredTemperature, int N_steps){
    //rescale velocities using equipartition theorem: v_desired = sqrt(T_desired/T_actual)*v_actual
    //double rescaling_factor = sqrt(desiredTemperature/statisticsSampler.temperature()); //sqrt(T_desired/T_actual)

    double rescaling_factor = sqrt(1+ (1/N_steps)*(desiredTemperature/currentTemperature-1)); //sqrt(1+(1/N_steps)(T_desired/T_current - 1))
    //Note: if N_steps is set to 1, then just gives the old scaling: sqrt(T_desired/T_actual)
    for(Atom *atom : atoms()) {
        atom->velocity *= rescaling_factor;  //a*=b means a = a*b
    }
    removeTotalMomentum();  //If don't do this, eventually will have drifting issue!
}

void System::removeTotalMomentum() {
    // Find the total momentum and remove momentum equally on each atom so the total momentum becomes zero.
    vec2 total_momentum;
    //double mass;
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
        total_momentum += mass*atom->velocity;
    }
    vec2 Mom_per_atom;   //2D momentum components
    Mom_per_atom = total_momentum/num_atoms(); //this is amount of abs(momentum) per atom that need to remove to get total to 0
    for(Atom *atom : atoms()) {
        for(int j=0;j<2;j++){
            //evenly modify components of velocity of each atom to yield total system momentum of 0
            atom->velocity[j] -= Mom_per_atom[j]/mass;
        }
    }

    //test if total_momentum was rezeroed 0: (a unit test)
    /*
    total_momentum.print("Total Momentum before removeMomentum");
    total_momentum.set(0,0,0);  //reset total-momentum to 0
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
          total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
     }
    total_momentum.print("Total Momentum after removeMomentum");   //print() is fnc in vec2 class, takes in a string input
    */
}


//THIS IS NOT READY TO BE USED FOR MULTIPLE PROCESSORS--> need to adjust in terms of nproc and subsystem size...
/*
void System::createFCCLattice(vec2 numberOfUnitCellsEachDimension, double latticeConstant, double temperature,  double variance, bool input_variance, double mass) {


    vec2 LatticeVector;  //vector which points to the origin of each unit cell
    //Note: 1st unit cell starts at 0,0

    double x,y;
    double halfLatticeConstant=0.5*latticeConstant;

    //each processor finds out what it's ID (AKA rank) is and how many processors there are
    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    for(int i=0;i<numberOfUnitCellsEachDimension[0];i++){
        //i.e. i = 0,1...N_x-1
        for(int j=0;j<numberOfUnitCellsEachDimension[1];j++){

                LatticeVector.set(latticeConstant*i,latticeConstant*j);

                //Place the 4 atoms of each fcc cell into coordinates. Use setInitialPosition(): this will both set position and
                //save the atom's initial position for use later.
                //NOTE: The PBCs will prevent from adding atoms which are beyond the system dimensions when approach the boundaries.
                Atom *atom1 = new Atom(UnitConverter::massFromSI(mass)); //uses mass in kg: mass is correct
                x = LatticeVector[0];
                y = LatticeVector[1];

                atom1->setInitialPosition(x,y);
                atom1->num_bndry_crossings.set(0.,0.);   //make sure initial # of bndry crossings is 0
                atom1->resetVelocityMaxwellian(temperature, variance, input_variance);
                m_atoms.push_back(atom1);     //add element to vector m_atoms 1 element (atom object)

                Atom *atom2 = new Atom(UnitConverter::massFromSI(mass));
                x = halfLatticeConstant + LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                atom2->setInitialPosition(x,y);
                atom2->num_bndry_crossings.set(0.,0.);
                atom2->resetVelocityMaxwellian(temperature, variance, input_variance);
                m_atoms.push_back(atom2);

                /*
                Atom *atom3 = new Atom(UnitConverter::massFromSI(mass));
                x = LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                atom3->setInitialPosition(x,y);
                atom3->num_bndry_crossings.set(0.,0.);
                atom3->resetVelocityMaxwellian(temperature, variance, input_variance);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(mass));
                x = halfLatticeConstant + LatticeVector[0];
                y = LatticeVector[1];
                atom4->setInitialPosition(x,y);
                atom4->num_bndry_crossings.set(0.,0.);
                atom4->resetVelocityMaxwellian(temperature, variance, input_variance);
                m_atoms.push_back(atom4);

                //std::cout << "atom mass = " <<atom1->mass() <<std::endl;

            }
        }


    setSystemSize(latticeConstant*numberOfUnitCellsEachDimension); //system size set by multiply vec2 # of unit cells by latticeConstant
    std::cout<<"system size = " << m_systemSize <<std::endl;
    std::cout<<"num_atoms = " << num_atoms() <<std::endl;
}
*/




void System::createSCLattice(vec2 Total_systemSize, vec2 subsystemSize, double latticeConstant, double temperature, double mass, vec2 subsystemOrigin) {

    double x;
    //each processor finds out what it's ID (AKA rank) is and how many processors there are
    int nprocs, rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

    x =rank*subsystemSize[0]*latticeConstant + 0.5*latticeConstant;  //initializes x at the right position for the processor's domain

    //std::cout <<"x" <<x <<std::endl;  //works correctly


    for(int i=0;i<subsystemSize[0];i++){  // i = 1 and < b/c can't have atoms on both boundaries--> will blow up!
        //i.e. i = 0,1...N_x-1

        double y = 0.5*latticeConstant;

        for(int j=0;j<subsystemSize[1];j++){

                Atom *atom = new Atom(mass); //uses mass already converted to LJ units in main //* atom means create a pointer --> atom is a pointer
                //equivalent to: Atom *atom;  atom = new Atom(..);  //declare a pointer and have it point to new Atom object

                atom->setInitialPosition(x,y);
                atom->num_bndry_crossings.set(0.,0.);   //make sure initial # of bndry crossings is 0
                atom->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (pointer to atom)
                y += latticeConstant;


                //std::cout << "mass line 227 in system" <<atom->m_mass <<std::endl;

        }
        x +=latticeConstant;
    }

    //this sets the TOTAL system size--> is used for PBCs which are at the outer boundaries
    setSystemSize(latticeConstant*Total_systemSize ); //system size set by multiply vec2 # of unit cells by latticeConstant
    //vec2 include_bndry(0.001,0.001);  //will add to SystemSize so include bndry for domain distributions among processors in SimSize
    setSimSize(Total_systemSize*latticeConstant);
    //setSimSize(Total_systemSize*latticeConstant + include_bndry);
    std::cout<<"system size = " << m_systemSize <<std::endl;
    std::cout<<"num_atoms = " << num_atoms() <<std::endl;


     setSubSystemSize(latticeConstant*subsystemSize);

}




void System::createRandomPositions(int num_particles, double side_length, double temperature, double mass){

    //Places particles randomly into a cube
    for(int i=0; i<num_particles; i++) {      //for all atoms (100 right now)
        Atom *atom = new Atom(mass); //uses mass in kg
        //Choose random x,y,z positions
        double x = Random::nextDouble(0, side_length); // random number in the interval [0,10]. nextDouble is defined in random.h
        double y = Random::nextDouble(0, side_length);

        atom->position.set(x,y);    //set atom to position x,y,z

        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)

        vec2 system_size;                                           //define a vector object containing system size
        system_size.set(side_length,side_length);       //set the values of elements of the vector
        setSystemSize(system_size); //system size defines size of system in EACH dimension
        m_num_atoms = num_atoms();
    }
}

void System::increaseTemperature(StatisticsSampler &statisticsSampler, double increment){
    //increases system temperature by factor: meant for use at each MD step to slowly heat system.
    double velocity_rescaling_factor = sqrt((statisticsSampler.temperature()+increment)/statisticsSampler.temperature()); //sqrt(T_new/T_current)=sqrt(factor)
    for(Atom *atom : atoms()) {
        atom->velocity *= velocity_rescaling_factor;  //a*=b means a = a*b
    }
}

void System::calculateForces() {
    for(Atom *atom : m_atoms) {
        atom->resetForce();
    }
    m_potential.calculateForces(*this); // this is a pointer, *this is a reference to this object
}

void System::step(double dt) {
    m_integrator.integrate(*this, dt);  //this calls velocityverlet.cpp
    m_steps++;
    m_time += dt;

}

//---------------------------------------------------------------------
/*!
 Remove atoms from the system.  Needs to sort indices because erase() operation reorders things; also, because of this it is fastest to pop from lowest to highest index.
 Returns the number of atoms deleted.
 \param [in] indices Vector of local indices of atoms to delete from the system
*/
int System::delete_atoms (std::vector <int> indices) {
        std::vector <Atom*>::iterator it = m_atoms.begin();
        //std::map <int, int>::iterator map_it;

        // Sort indices from lowest to highest
        sort (indices.begin(), indices.end());

        // Pop in this order
        int shift = 0;
        int upper;
        for (unsigned int i = 0; i < indices.size(); ++i) {           //finds size = 6
                if (i < indices.size()-1) {
                        upper = indices[i+1]-shift;
                } else {
                        upper = m_atoms.size();
                }

                m_atoms.erase(it - shift + indices[i]); //must pass it a accept a const iterator: note: 'it' is an iterator
          //      glob_to_loc_id_.erase(map_it);  //GETS STUCK AT THIS ERASE!
                ++shift;
        }
        m_num_atoms -= shift;

        //test deletion: loops through  all atoms
        /*
        for(int i=0; i<num_atoms()-1; i++){
            std::cout <<"erase test" << atoms(i)->position[0] << std::endl;
        }
        */



        return shift;
}
/* I think this version is not needed --> is for 1 atom only it seems...
/*
 Attempt to push ONE atom into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms;
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.

std::vector <int> System::add_atoms (const int natoms, Atom &new_atoms) {
        int index = m_atoms.size();
        std::vector <int> update_proc(natoms);  //make vector of size natoms
        for (int i = 0; i < natoms; ++i) {
                m_atoms.push_back(new_atoms[i]);
                glob_to_loc_id_[new_atoms[i].atom_index] = index;
                update_proc[i] = new_atoms[i].atom_index; //filled with indices of the newly added atoms
                index++;
        }
        m_num_atoms += natoms;
        return update_proc;
}
*/


/*

 Attempt to push an atom(s) into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms.
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.

std::vector <int> System::add_atoms (std::vector <Atom> new_atoms) { //accepts POINTER to vector of new_atoms
        int index = m_atoms.size(), natoms = new_atoms.size();  //use . for objects -> for pointers
        //double mass;
        std::vector <int> update_proc(natoms);

        int nprocs, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

        //gets through here

        //std::cout <<"index at line 248"<< index <<std::endl;
        //std::cout <<"natoms at line 248"<< natoms <<std::endl;  //finds that there are 100 new atoms --> ok.. make sense for 2 processors...
        //issue is that it might try to add atoms that already are in the processor AGAIN
        for (int i = 0; i < natoms; ++i) {
               Atom *adding_atom = &new_atoms[i];  //makes a pointer to the individual atom
               m_atoms.push_back(adding_atom); //note: at() is member fnc of c++ vector class: Returns a reference to the element at position n in the vector. is almost same as [] operator but here checks whether i is within bounds of the vector.




               //hard code resseting the mass to correct value
              //adding_atom->setMass(UnitConverter::massFromSI(6.63352088e-26));

               std::cout <<"added atom velocity" << adding_atom->velocity[0] << " " << adding_atom->velocity[1] << "force should be 0 right now" << adding_atom->force[0] << " " <<adding_atom->force[1] << "proc " << rank <<std::endl;
               //glob_to_loc_id_[adding_atom->atom_index] = index;  //BREAKS HERE!  //it has atom_index = 0...
               //std::cout <<"line 270 in addatoms" <<std::endl;
               //works through here
               update_proc[i] = i;  //Later can add atom indices to actually track how the atoms move around to different processors//adding_atom->atom_index; //note: if accessing member property through a pointer, must use -> instead of .
               index++;
\
        }
        m_num_atoms += natoms;
        return update_proc;
}
*/

void System::add_atoms (std::vector <double> new_atoms, double num_recieved) { //accepts POINTER to vector of new_atoms
        //int index = m_atoms.size();//, natoms = new_atoms.lenth()/4; //4 numbers for each atom //use . for objects -> for pointers
        //double mass;
        //std::vector <int> update_proc(natoms);

        int nprocs, rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);   //find ID
        MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

        //gets through here

        //std::cout <<"index at line 248"<< index <<std::endl;
        //std::cout <<"natoms at line 248"<< natoms <<std::endl;  //finds that there are 100 new atoms --> ok.. make sense for 2 processors...
        //issue is that it might try to add atoms that already are in the processor AGAIN
        for (int i = 0; i < num_recieved; ++i) {

            int index =4*i;


              //we will initialize NEW atom object
               Atom *adding_atom = new Atom(mass);  //makes a pointer to the individual atom
               adding_atom->position[0] = new_atoms[index];
               adding_atom->position[1] = new_atoms[index+1];
               //TRY resetting any added atoms velocities according to maxwell boltzmann--> b/c otherwise they blow up...
               //adding_atom->resetVelocityMaxwellian(600.);  //need to change this to variable



               adding_atom->velocity[0] =new_atoms[index+2];
               adding_atom->velocity[1] = new_atoms[index+3];
               adding_atom->force.set(0,0);  //intialize w/ 0 forces since calculated in LJ



               m_atoms.push_back(adding_atom); //note: at() is member fnc of c++ vector class: Returns a reference to the element at position n in the vector. is almost same as [] operator but here checks whether i is within bounds of the vector.

               //hard code resseting the mass to correct value
              //adding_atom->setMass(UnitConverter::massFromSI(6.63352088e-26));

               std::cout <<"added atom position" << adding_atom->position[0] <<" " << adding_atom->position[1] << "added atom velocity" << adding_atom->velocity[0] << " " << adding_atom->velocity[1] << "proc " << rank <<std::endl;
               //glob_to_loc_id_[adding_atom->atom_index] = index;  //BREAKS HERE!  //it has atom_index = 0...
               //std::cout <<"line 270 in addatoms" <<std::endl;
               //works through here
              // update_proc[i] = i;  //Later can add atom indices to actually track how the atoms move around to different processors//adding_atom->atom_index; //note: if accessing member property through a pointer, must use -> instead of .
               //index++;
\
        }
        m_num_atoms += num_recieved;
       // return update_proc;
}




/*!
 Attempt to push ghost atom(s) into the system.  This assigns the map automatically to link the atoms global index to the local storage location.
 This reallocates the internal vector that stores the atoms; if a memory error occurs during such reallocation, an error is given and the system exits.
 Does not change m_num_atoms (the number of atoms a processor is responsible for)
 \param [in] natoms Length of the array of atoms to add to the system.
 \param [in] \*new_atoms Pointer to an array of atoms the user has created elsewhere.
 */

/*
void System::add_ghost_atoms (const int natoms, std::vector <Atom> new_atoms) {
        for (int i = 0; i < natoms; ++i) {
            /* Only add the atom to the system if it is not already contained in the system.
                        Do this by going through the atoms and comparing atom_index of the atoms already in the system and the atoms to be added
            bool found_atom = false;
            for (unsigned int j = 0; j < m_atoms.size(); ++j) {
                if (m_atoms[j]->atom_index == new_atoms.at(i).atom_index) {
                    found_atom = true;
                    break;
                }
            }
            if (!found_atom) {
                Atom *adding_atom = &new_atoms[i];
                m_atoms.push_back(adding_atom);
            }

        }
        return;
}
*/

/*!
 Clears the atoms communicated from neighbouring domains from the list of atoms stored in the system leaving only the atoms the system is responsible for.
 */
void System::clear_ghost_atoms () {
        m_atoms.erase(m_atoms.begin()+m_num_atoms, m_atoms.end());
}
