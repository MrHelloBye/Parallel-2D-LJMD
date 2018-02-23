
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
    viewMain(argc,argv);
  } else{
    modelMain(argc,argv);
  }
#else
    modelMain(argc,argv);
#endif
  MPI_Finalize();

}
