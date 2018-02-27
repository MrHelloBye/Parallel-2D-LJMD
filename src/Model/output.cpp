#include "output.h"
#include <mpi.h>



//in c++ can only output a pointer to a C style array from functions
double* output_positions(System system){
    double all_positions [2*(nprocs - 1)*system.num_atoms()];
    int i = 0;

    for(Atom *atom : system.atoms()) {
        all_positions[i] = atom->position.x();
        all_positions[i+1] = atom->position.y();
        i+=2;
    }
    return all_positions;
}

