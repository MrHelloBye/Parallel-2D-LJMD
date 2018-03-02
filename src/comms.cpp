#include <mpi.h>

#include "comms.hpp"
#include "Model/system.h"
#include "Controller/controllerState.hpp"

/**********************************************************************
 * Send the model data (positions, maybe velocites)
 * to the view
 **********************************************************************/
int comms_sendModelData(System& system){
  int my_id, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

  int num_local = system.num_atoms();

  //Gather all positions
  float* pos_buf = new float[2*system.num_atoms()];
  int i = 0, j=0;
  for(Atom *atom : system.atoms()) {
    pos_buf[2*i] = atom->position[0];
    pos_buf[2*i+1] = atom->position[1];
    i++;
  }

  //Send atom counts to root
  MPI_Gather(&num_local,1,MPI_INT,NULL,1,MPI_INT,0,MPI_COMM_WORLD);


  //Gatherv allows to recieve varying sizes of data
  MPI_Gatherv(pos_buf, num_local*2, MPI_FLOAT, NULL, NULL, NULL, MPI_FLOAT, 0, MPI_COMM_WORLD);

  //Clean up
  delete pos_buf;

  return 0;
}

/**********************************************************************
 * Get the model data from the model to the view
 *
 * Inputs:
 * pos_buf - the positions that will be gathered from other MPI tasks
 * hue_buf - the hues that will be gathered/set by MPI task
 **********************************************************************/

int comms_gatherModelData(float** pos_buf,float** hue_buf,int* atomCounts,int& numAtoms){



  int my_id, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors



  int num_local = 0;

  int offset=0;
  int displs[nprocs];
  int rcounts[nprocs];

  MPI_Gather(&num_local,1,MPI_INT,atomCounts,1,MPI_INT,0,MPI_COMM_WORLD);

  int num_total = 0;	

  for(int i =0; i < nprocs; i++){
    num_total += atomCounts[i];

    //Get the total count of elements (twice the number of atoms)
    rcounts[i] = atomCounts[i]*2;
    if( i == 0){
      displs[i] = 0;
    } else {
      displs[i] = displs[i-1]+rcounts[i-1];
    }
  }

  //WARNING pos is assumed big enough
  if(numAtoms!=num_total){
    delete *pos_buf;
    delete *hue_buf;
    *pos_buf = new float[num_total*2];
    *hue_buf = new float[num_total];

    numAtoms = num_total;
  }

  int h = 0;
  for(int proc = 1; proc < nprocs; proc++){
    //Skip root: no particles there
    for( int i =0; i < atomCounts[proc]; i++,h++){
      (*hue_buf)[h] = (proc-1)*360./(nprocs-1);
    }
  }

  MPI_Gatherv(NULL,0,MPI_FLOAT, *pos_buf, rcounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);
  return 0;
}

/**********************************************************************
 * Send/Recv the controller state from the view to model
 *
 * For now send and recv call the same bcast, but in the future we may want the
 * view and model to behave differently
 *
 * Also, by the number of static variables, this should go in a class with
 * friend functions
 **********************************************************************/

int comms_bcastControllerState(ControllerState* state){
  int my_id;
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID

  static bool alreadySending = false;
  static ControllerState buf;
  static MPI_Request request;
  static MPI_Datatype stateType =0;


  if(stateType == 0){
    int count = 3;//3 variables
    int blockLengths[] = {2,1,1};
    MPI_Aint displs[] = {0,2*sizeof(float),3*sizeof(float)};
    MPI_Datatype types[] = {MPI_FLOAT,MPI_FLOAT,MPI_UNSIGNED};

    MPI_Type_create_struct(3, blockLengths, displs, types, &stateType);
    MPI_Type_commit(&stateType);
  }

  if(!alreadySending){
    //Start a bcast
    alreadySending = true;

    //Fill the buffer
    if(my_id == 0)
      buf = *state;
    MPI_Ibcast(&buf,1,stateType,0,MPI_COMM_WORLD,&request);
  } else{
    int flag;
    MPI_Status status;
    MPI_Test(&request,&flag,&status);

    if(flag){
      //Bcast successful
      //Copy over the buffer
      if(my_id != 0)
        *state = buf;
      alreadySending = false;

      return 1;
    }
  }
  return 0;

}

int comms_sendControllerState(ControllerState* state){
  return comms_bcastControllerState(state);
}
int comms_recvControllerState(ControllerState* state){
  return comms_bcastControllerState(state);
}
