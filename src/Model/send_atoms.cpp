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
	
        // store atoms to send and receive: vectors of atom objects (not pointers)
        /*
        vector<Atom> to_left;
        vector<Atom> to_right;
        vector<Atom> from_left;
        vector<Atom> from_right;
        */

        //try to store as just vectors
        vector<double> to_left;
        vector<double> to_right;
        vector<double> from_left;
        vector<double> from_right;


       int array_index;  //for keeping track of storing values in arrays
	
	// store indices of atoms that have been sent (so we can delete them)
	vector<int> to_delete;

        //NOTE: synthax may need to be changed for openMPI
        //get info about the requests to send/recieve --> later in WaitAll can make sure that all send and recieves are complete before moving forward
	MPI_Request req[4], req2[4];  //req[4] are the 4 requests to send/recieve # of atoms, req2[4] are the 4 requests sending/recieving the atoms
	MPI_Status stat[4], stat2[4];  //get status of the requests


	int proc_to;
        for (int i=0; i!=system->num_atoms(); ++i) {   //-> are used b/c
		// calculate the processor for each atom
                proc_to = floor(system->atoms(i)->position[decomp_dim]/ sim_size[decomp_dim] * nprocs);
                //


               // std::cout << "index  " << i<< "position  " << system->atoms(i)->position[0] << " " <<system->atoms(i)->position[1] << "proc  " <<rank <<std::endl;
                // std::cout << "velocity  " << system->atoms(i)->velocity[0] << "   " <<system->atoms(i)->velocity[1] << "mass  " << system->atoms(i)->m_mass << "proc  " <<rank <<std::endl;

                //print out if atom has wrong mass
                 //if(system->atoms(i)->m_mass < 39.){

                     //std::cout <<"wrong mass atom mass " <<system->atoms(i)->m_mass <<std::endl;
                // }
                //std::cout <<"velocityline 48" <<system->atoms(i)->velocity[0] <<std::endl;
                  //std::cout <<"sim_size" <<sim_size[decomp_dim] <<std::endl;
                 // std::cout <<"proc_to" << proc_to <<std::endl;
                //proc_to works fine

		if (proc_to == (rank - 1 + nprocs) % nprocs) {

                        num_to_left++;

                        //array_index = (num_to_left-1)*4;

                        //SEGMENTATION FAULT IS SOMEWHERE HERE!
                        //FIND THAT IF DON'T USE push_back, GET SEGMENTATON  FAULT! (I.E. IF USE to_left[index] = for deifnitions...

                        to_left.push_back(system->atoms(i)->position[0]);
                        to_left.push_back(system->atoms(i)->position[1]);
                        to_left.push_back(system->atoms(i)->velocity[0]);
                        to_left.push_back(system->atoms(i)->velocity[1]);

                        std::cout <<"BEFORE SEND position" <<system->atoms(i)->position[0] << " " <<system->atoms(i)->position[1] << "vel" <<system->atoms(i)->velocity[0] << " " <<system->atoms(i)->velocity[1] <<std::endl;


                        //to_left[array_index] = system->atoms(i)->position[0];
                       // to_left[array_index+1] = system->atoms(i)->position[1];
                       // to_left[array_index+2] = system->atoms(i)->velocity[0];
                        //to_left[array_index+3] = system->atoms(i)->position[1];

                        //to_left.push_back(*(system->atoms(i)));  //* means pass the value that atoms (which is a pointer) points to
                        //std::cout<<"position to_left send line 53" << system->atoms(i)->position[0] <<endl;

                        to_delete.push_back(i);  //never gets sent anywhere-< jus tfor cucrrent proc
		}
		else if (proc_to == (rank + 1) % nprocs) {
                        //to_right.push_back(*(system->atoms(i)));
			num_to_right++;
                        //array_index = (num_to_right-1)*4;

                        to_right.push_back(system->atoms(i)->position[0]);
                        to_right.push_back(system->atoms(i)->position[1]);
                        to_right.push_back(system->atoms(i)->velocity[0]);
                        to_right.push_back(system->atoms(i)->velocity[1]);


                         std::cout <<"BEFORE SEND position" <<system->atoms(i)->position[0] << " " <<system->atoms(i)->position[1] << "vel" <<system->atoms(i)->velocity[0] << " " <<system->atoms(i)->velocity[1] <<std::endl;


			to_delete.push_back(i);
		}
		else if (proc_to != rank) {
                       // std::cout <<"Atom moved too many boxes" << proc_to << std::endl;
		}			
	}
         //gets through here fine



	// send number of atoms
        //synthax: MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,  MPI_Comm comm, MPI_Request *request)
        //(starting address of data that sending (called buffer), # of elements in buffer, MPI data type, destination processor, message tag, communicator, pointer to the request)
        MPI_Isend(&num_to_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req); // req will input pointer to beggining of req array --> req[0].
        MPI_Irecv(&num_from_left, 1, MPI_INT, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req+1);  //req+1 b/c req is pointer --> pts to req[1]
        MPI_Isend(&num_to_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+2);
        MPI_Irecv(&num_from_right, 1, MPI_INT, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req+3);
        MPI_Waitall (4, req, stat);  //wait for all send and receive requests to be completed


        //TRY SENDING ATOM INFORMATION AS SIMPLY 4 doubles FOR EACH ATOM!! rx, ry, vx, vy

        // resize vectors of atom data--> 4*# of atoms
        from_left.resize(4*num_from_left);
        from_right.resize(4*num_from_right);



         MPI_Isend(&to_left[0], to_left.size(), MPI_DOUBLE, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
         MPI_Irecv(&from_left[0], from_left.size(), MPI_DOUBLE, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
         MPI_Isend(&to_right[0], to_right.size(), MPI_DOUBLE, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
         MPI_Irecv(&from_right[0], from_right.size(), MPI_DOUBLE, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
         MPI_Waitall (4, req2, stat2);


       // std::cout << "got past send" << "to_left[0]" << to_left[0] << std::endl;

	
        // send atoms --> & here is b/c sendsing the vectors by address to MPI_Isend etc..
        /*
        MPI_Isend(&to_left[0], num_to_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2);
        MPI_Irecv(&from_left[0], num_from_left, MPI_ATOM, (rank - 1 + nprocs) % nprocs, 1, MPI_COMM_WORLD, req2+1);
        MPI_Isend(&to_right[0], num_to_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+2);
        MPI_Irecv(&from_right[0], num_from_right, MPI_ATOM, (rank + 1) % nprocs, 1, MPI_COMM_WORLD, req2+3);
        MPI_Waitall (4, req2, stat2);
        */

        // add atoms to system



        // add atoms to system

        /*

        if(num_to_right !=0){
           std::cout <<"before sent to right" << to_right[0].position[0] <<" " << to_right[0].position[1] << "velocity" << to_right[0].velocity[0]
                   <<" " << to_right[0].velocity[1] <<" proc" << rank << std::endl;
        }

        if(num_to_left !=0){
            std::cout <<"before sent to left" << to_left[0].position[0] <<" " << to_left[0].position[1] << "velocity" << to_left[0].velocity[0]
                    <<" " << to_left[0].velocity[1] <<" proc" << rank<< std::endl;
        }

        if(num_from_left !=0){
            //std::cout <<"Send atoms: num_from left " <<num_from_left << "proc" << rank << std::endl;

            std::cout <<"received from left" << from_left[0].position[0] <<" " <<from_left[0].position[1] << "velocity" << from_left[0].velocity[0]
                    <<" " << from_left[0].velocity[1] <<" proc" << rank << std::endl;
        }

         if(num_from_right != 0){
            //std::cout <<"num_from_right" <<num_from_right << "proc" << rank << std::endl;

            std::cout <<"received from right" << from_right[0].position[0] <<" " << from_right[0].position[1]  << "velocity" << from_right[0].velocity[0]
                    <<" " << from_right[0].velocity[1] << " proc" <<rank << std::endl;
         }
         */




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


	
	MPI_Barrier(MPI_COMM_WORLD);
}
