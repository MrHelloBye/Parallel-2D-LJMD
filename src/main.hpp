#ifndef MAIN_HPP_
#define MAIN_HPP_

#ifdef VIEW_ENABLED
int viewMain(int argc,char** argv);
#endif

int modelMain(int argc,char** argv,MPI_Comm comm_model);

#endif //MAIN_HPP_
