#ifndef COMMS_HPP_
#define COMMS_HPP_


class System;
//Send the model data (positions, maybe velocites)
//To the view
int comms_sendModelData(System& system);

//Get the model data from the model to the view
int comms_gatherModelData(float** pos_buf,float** hue_buf,int* atomCounts,int& numAtoms);

#endif //COMMS_HPP_
