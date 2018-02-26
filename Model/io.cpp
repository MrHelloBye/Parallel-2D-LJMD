#include "io.h"
#include "system.h"
#include "atom.h"
#include "unitconverter.h"
#include "statisticssampler.h"
#include <cstdlib>
using std::endl; using std::cout;


IO::IO(const char *filename)
{
    open(filename);
}

IO::~IO() {
    close();
}

void IO::open(const char *filename) {
    if(file.is_open()) {
        std::cout << "<IO.cpp> Error, tried to open file " << filename << ", but some file is already open." << endl;
        exit(1);
    }

    file.open(filename);
}

void IO::close() {
    if(file.is_open()) {
        file.close();
    }
}

// This saves the current state to a file following the xyz-standard (see http://en.wikipedia.org/wiki/XYZ_file_format )
// It can easily be opened in Ovito. Note that you can also output more properties than just the position. You can print the
// velocities per particle (or kinetic energy etc), and color the atoms in Ovito based on these properties.

void IO::saveState(System &system, StatisticsSampler &statisticsSampler) //statisticsSampler object of the class StatisticsSampler is passed by reference
{
    if(file.is_open()) {
        file << system.atoms().size() << endl;
        file << "The is an optional comment line that can be empty. The reason we use H is so particles get smaller in Ovito" << endl;
        for(Atom *atom : system.atoms()) {
            file << "H " <<         
                    UnitConverter::lengthToAngstroms(atom->position.x()) << " " <<
                    UnitConverter::lengthToAngstroms(atom->position.y()) << " " << endl;

                    //statisticsSampler.density() << "\n"; //density  is just a function that accesses the member of the class m_density and returns it
                    //UnitConverter::velocityToSI(atom->velocity.x()) <<" "<<
                    //UnitConverter::velocityToSI(atom->velocity.y()) <<" "<<

        }
    }
}

//======================
//the below functions are not members of the class

double * getPositions(System &system){
    double positions[2*system.num_atoms()];
    int i = 0;
     for(Atom *atom : system.atoms()) {
         positions[i] = system.atoms(i)->position[0];
         positions[i+1] = system.atoms(i+1)->position[1];
         i+=2;
     }

     return positions;

}


