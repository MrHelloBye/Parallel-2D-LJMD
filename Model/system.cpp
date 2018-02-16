#include "system.h"
#include "velocityverlet.h"
#include "lennardjones.h"
#include "statisticssampler.h"
#include "unitconverter.h"
#include "random.h"

System::System()   //constructor
{

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

/*
// remove atoms that escape the SUB system (corresponding to current processor)
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
    for(Atom *atom : atoms()) { //c++11 way of iterating through  entire vector or array
        total_momentum += atom->mass()*atom->velocity;  //mass() returns value of m_mass (atom's mass)
    }
    vec2 Mom_per_atom;   //2D momentum components
    Mom_per_atom = total_momentum/num_atoms(); //this is amount of abs(momentum) per atom that need to remove to get total to 0
    for(Atom *atom : atoms()) {
        for(int j=0;j<2;j++){
            //evenly modify components of velocity of each atom to yield total system momentum of 0
            atom->velocity[j] -= Mom_per_atom[j]/atom->mass();
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

void System::createSCLattice(vec2 Total_systemSize, vec2 subsystemSize, double latticeConstant, double temperature, double mass, vec2 subsystemOrigin) {

    double x = subsystemOrigin[0];


    for(int i=0;i<subsystemSize[0];i++){
        //i.e. i = 0,1...N_x-1
        x +=latticeConstant;
        double y = subsystemOrigin[1];

        for(int j=0;j<subsystemSize[1];j++){
                y += latticeConstant;

                Atom *atom = new Atom(UnitConverter::massFromSI(mass)); //uses mass in kg: mass is correct

                atom->setInitialPosition(x,y);
                atom->num_bndry_crossings.set(0.,0.);   //make sure initial # of bndry crossings is 0
                atom->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)
        }
    }

    //this sets the TOTAL system size--> is used for PBCs which are at the outer boundaries
    setSystemSize(latticeConstant*Total_systemSize); //system size set by multiply vec2 # of unit cells by latticeConstant
    std::cout<<"system size = " << m_systemSize <<std::endl;
    std::cout<<"num_atoms = " << num_atoms() <<std::endl;

    //setSubSystemSize(latticeConstant*subsystemSize);

}



/*
void System::createFCCLattice(vec2 numberOfUnitCellsEachDimension, double latticeConstant, double temperature, double mass) {


    vec2 LatticeVector;  //vector which points to the origin of each unit cell
    //Note: 1st unit cell starts at 0,0

    double x,y;
    double halfLatticeConstant=0.5*latticeConstant;

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
                atom1->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom1);     //add element to vector m_atoms 1 element (atom object)

                Atom *atom2 = new Atom(UnitConverter::massFromSI(mass));
                x = halfLatticeConstant + LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                atom2->setInitialPosition(x,y);
                atom2->num_bndry_crossings.set(0.,0.);
                atom2->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom2);

                Atom *atom3 = new Atom(UnitConverter::massFromSI(mass));
                x = LatticeVector[0];
                y = halfLatticeConstant + LatticeVector[1];
                atom3->setInitialPosition(x,y);
                atom3->num_bndry_crossings.set(0.,0.);
                atom3->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom3);

                Atom *atom4 = new Atom(UnitConverter::massFromSI(mass));
                x = halfLatticeConstant + LatticeVector[0];
                y = LatticeVector[1];
                atom4->setInitialPosition(x,y);
                atom4->num_bndry_crossings.set(0.,0.);
                atom4->resetVelocityMaxwellian(temperature);
                m_atoms.push_back(atom4);

                //std::cout << "atom mass = " <<atom1->mass() <<std::endl;
            }
        }



    setSystemSize(latticeConstant*numberOfUnitCellsEachDimension); //system size set by multiply vec2 # of unit cells by latticeConstant
    std::cout<<"system size = " << m_systemSize <<std::endl;
    std::cout<<"num_atoms = " << num_atoms() <<std::endl;
}
*/

void System::createRandomPositions(int num_particles, double side_length, double temperature, double mass){

    //Places particles randomly into a cube
    for(int i=0; i<num_particles; i++) {      //for all atoms (100 right now)
        Atom *atom = new Atom(UnitConverter::massFromSI(mass)); //uses mass in kg
        //Choose random x,y,z positions
        double x = Random::nextDouble(0, side_length); // random number in the interval [0,10]. nextDouble is defined in random.h
        double y = Random::nextDouble(0, side_length);

        atom->position.set(x,y);    //set atom to position x,y,z

        atom->resetVelocityMaxwellian(temperature);
        m_atoms.push_back(atom);     //add element to vector m_atoms 1 element (atom object)

        vec2 system_size;                                           //define a vector object containing system size
        system_size.set(side_length,side_length);       //set the values of elements of the vector
        setSystemSize(system_size); //system size defines size of system in EACH dimension
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
