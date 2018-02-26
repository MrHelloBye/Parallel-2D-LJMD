#include "send_atoms.h"
#include <mpi.h>
#include <math.h>  //<> for include for files in other directories, i.e. STL
#include "global.h" //tells it that MPI_ATOM is extern variable--> defined elsewhere


using namespace std;

/*
 This function sends the atoms that have left the domain of the processor to the relevant neighboring processor.
*/
void send_atoms(System *system) {

        vec2 sim_size = system->simSize(); //returns total simulation size, needed for finding which proc atoms belong to.
        //NOTE: sim_size is just slightly larger than systemSize in order to include the atoms on the system edges
        int decomp_dim = 0;  // 0 or 1, x or y direction of decomposition

	int nprocs, rank;
	MPI_Comm_size (MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank (MPI_COMM_WORLD, &rank);
	
	// store number of atoms to send to and receive from the processor on the left and on the right
	int num_to_left = 0;
        int num_to_right = 0;
	int num_from_left, num_from_right;

        //try to store as just vectors
        vector<double> to_left;
        vector<double> to_right;
        vector<double> from_left;
        vector<double> from_right;

        int array_index;  //for keeping track of storing values in arrays
	
	// store indices of atoms that have been sent (so we can delete them)
	vector<int> to_delete;

        //get info about the requests to send/recieve --> later in WaitAll can make sure that all send and recieves are complete before moving forward
	MPI_Request req[4], req2[4];  //req[4] are the 4 requests to send/recieve # of atoms, req2[4] are the 4 requests sending/recieving the atoms
	MPI_Status stat[4], stat2[4];  //get status of the requests

	int proc_to;
        for (int i=0; i!=system->num_atoms(); ++i) {   //-> are used b/c
		// calculate the processor for each atom
                proc_to = floor(system->atoms(i)->position[decomp_dim]/ sim_size[decomp_dim] * (nprocs-1));

                if (proc_to == (rank -1- 1 + nprocs-1) % (nprocs-1)) {
                        num_to_left++;
                        to_left.push_back(system->atoms(i)->position[0]);
                        to_left.push_back(system->atoms(i)->position[1]);
                        to_left.push_back(system->atoms(i)->velocity[0]);
                        to_left.push_back(system->atoms(i)->velocity[1]);

                      // std::cout <<"BEFORE SEND position" <<system->atoms(i)->position[0] << " " <<system->atoms(i)->position[1] << "vel" <<system->atoms(i)->velocity[0] << " " <<system->atoms(i)->velocity[1] <<std::endl;

                        to_delete.push_back(i);  //using erase, requires an iterator //never gets sent anywhere-< just for cucrrent proc
		}
                else if (proc_to == (rank ) % (nprocs-1)) {
			num_to_right++;
                        to_right.push_back(system->atoms(i)->position[0]);
                        to_right.push_back(system->atoms(i)->position[1]);
                        to_right.push_back(system->atoms(i)->velocity[0]);
                        to_right.push_back(system->atoms(i)->velocity[1]);
                        to_delete.push_back(i);
		}
                else if (proc_to != rank) {
                        //std::cout <<"Atom moved too many boxes" << proc_to << std::endl;
                }
	}


	// send number of atoms
        int ln =  (rank -1- 1 + nprocs-1) % (nprocs-1)+1;
        int rn = (rank ) %( nprocs-1)+1;

        //synthax: MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm, MPI_Request *request)
        //(starting address of data that sending (called buffer), # of elements in buffer, MPI data type, destination processor, message tag, communicator, pointer to the request)

        MPI_Isend(&num_to_left, 1, MPI_INT, ln, 10*rank + ln, MPI_COMM_WORLD, req); //note: these formulas will pass i.e. from proc 0 to proc (max at right)...--> satisfy PBCs...
        MPI_Irecv(&num_from_left, 1, MPI_INT, ln, 10*ln +rank, MPI_COMM_WORLD, req+1);
        MPI_Isend(&num_to_right, 1, MPI_INT, rn, 10*rank + rn, MPI_COMM_WORLD, req+2);
        MPI_Irecv(&num_from_right, 1, MPI_INT, rn, 10*rn +rank, MPI_COMM_WORLD, req+3);
        MPI_Waitall (4, req, stat);

        //SENDING ATOM INFORMATION AS SIMPLY 4 doubles FOR EACH ATOM!! rx, ry, vx, vy

        // resize vectors of atom data--> 4*# of atoms
        from_left.resize(4*num_from_left);
        from_right.resize(4*num_from_right);

        MPI_Isend(&to_left[0], to_left.size(), MPI_DOUBLE, ln, 10*rank + ln, MPI_COMM_WORLD, req2);
        MPI_Irecv(&from_left[0], from_left.size(), MPI_DOUBLE, ln, 10*ln +rank, MPI_COMM_WORLD, req2+1);
        MPI_Isend(&to_right[0], to_right.size(), MPI_DOUBLE, rn, 10*rank + rn, MPI_COMM_WORLD, req2+2);
        MPI_Irecv(&from_right[0], from_right.size(), MPI_DOUBLE, rn, 10*rn +rank, MPI_COMM_WORLD, req2+3);
        MPI_Waitall (4, req2, stat2);

        // add atoms to system
        system->add_atoms(from_left, num_from_left);  //thes are just arrays of numbers...
        system->add_atoms(from_right, num_from_right);

        // delete atoms that we sent to another system
        system->delete_atoms(to_delete);

        //clear the vectors so can be used fresh again
        to_left.clear();
        to_right.clear();
        from_right.clear();
        from_left.clear();
        to_delete.clear();

        //MPI_Barrier(MPI_COMM_WORLD);
}
