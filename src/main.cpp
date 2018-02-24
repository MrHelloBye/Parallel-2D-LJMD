
#include <mpi.h>

#include "main.hpp"

int main(int argc,char **argv){
  int my_id, nprocs;

  MPI_Init(&argc, &argv);  
  MPI_Comm_rank(MPI_COMM_WORLD, &my_id);   //find ID
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);  //find # of processors

#ifdef VIEW_ENABLED
  if(my_id == 0){
    //Run the View on root
    MPI_Comm comm_model;
    MPI_Comm_split(MPI_COMM_WORLD,
        MPI_UNDEFINED,my_id,&comm_model);

    viewMain(argc,argv);
  } else{
    MPI_Comm comm_model;
    MPI_Comm_split(MPI_COMM_WORLD,
        1,my_id,&comm_model);

    modelMain(argc,argv,comm_model);
  }
#else
  MPI_Comm comm_model;
  MPI_Comm_split(MPI_COMM_WORLD,
      1,my_id,&comm_model);

  modelMain(argc,argv,comm_model);
#endif
  MPI_Finalize();

}
