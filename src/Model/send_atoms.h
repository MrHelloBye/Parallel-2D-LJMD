#ifndef SEND_ATOMS_H
#define SEND_ATOMS_H

#include "system.h"
#include "atom.h"
#include <math.h>

// Move atoms between processors (domains)
void send_atoms(System *system);

#endif


