#include "atom.h"
#include "random.h"
//#include <cmath>
#include <math.h>
#include "global.h"
#include "mpiatom.h"


Atom::Atom(double mass) :
    m_mass(mass)
{

}

Atom::Atom()
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

    //extern double mass;
    // Resetting the velocity according to a Maxwell-Boltzmann distribution (see http://en.wikipedia.org/wiki/Maxwell%E2%80%93Boltzmann_distribution )
    double boltzmannConstant = 1.0; // In these units, the boltzmann constant equals 1
    double standardDeviation =  sqrt(boltzmannConstant*temperature/mass);
    velocity.randomGaussian(0, standardDeviation);  //note: randomGaussian is defined in vec2.cpp. arguments(mean, stdev), so this is Gaussian centered about 0
}

/*
void create_MPI_ATOM () {


        //with sending mass and atom index
        /*
         const int num_fields = 4;
        MPI_Datatype type[num_fields] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE, MPI_INT}; //Type of elements in each block (array of handles to data-type objects).
        //MPI_UB is an upper-bound indicator, tells where the datatype needs to end
        int blocklen[num_fields] = {2, 2, 1, 1};


        //w/o sending mass nor atom index
        const int num_fields = 4;
        MPI_Datatype type[num_fields] = {MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE };   //try doing each vector component explicitely since my vectors are not C++ arrays, but are a custom class...
        int blocklen[num_fields] = {1, 1,1,1};



        //no need to send the forces b/c need to be recalculated

        /* To ensure memory padding (done by compiler) in an array (or vector) is
         measured correctly, define a dummy array here to check addresses
         explicitly.  This is much better practice and safer for more compilers.


        //disp  are displacements in bytes

        Atom atom[2];  //dummy array of pointers to atoms
        MPI_Aint disp[num_fields];  //MPI_Aint is C datatype that holds any valid address. holds array of displacements in bytes
        MPI_Aint start_address, address;

        MPI_Get_address(&(atom[0]), &start_address); //get address of start of atom array, store in start_address

        MPI_Get_address(&(atom[0].position[0]), &address);
        disp[0] = address - start_address;     //distance in bytes between data found just by subtracting

        MPI_Get_address(&(atom[0].position[1]), &address);
        disp[1] = address - start_address;


        MPI_Get_address(&(atom[0].velocity[0]), &address);
        disp[2] = address - start_address;

        MPI_Get_address(&(atom[0].velocity[1]), &address);
        disp[3] = address - start_address;



        //MPI_Get_address(&(atom[0].m_mass), &address);
        //disp[2] = address - start_address;

        //MPI_Get_address(&(atom[0].atom_index), &address);
        //disp[3] = address - start_address;

        //MPI_Get_address(&(atom[1]), &address);  //the next atom--> used for the MPI_UB upper bound indicator
        //disp[4] = address - start_address;

        MPI_Type_create_struct(num_fields, blocklen, disp, type, &MPI_ATOM);  //make a structure MPI datatype
        MPI_Type_commit(&MPI_ATOM);
        /*The commit operation commits the data type. A data type is the formal description of a communication buffer,
          NOT the content of that buffer. After a data type has been committed, it can be repeatedly reused to communicate the changing content of a buffer
}
 */

/*
/*!
 * This function utilizes the complementary MPI_Type_free routine
 * to mark the MPI_ATOM type for deallocation at the end of the program.
 * \sa finalize

void delete_MPI_atom() {
        MPI_Type_free (&MPI_ATOM);
}
*/
