#ifndef COMMS_HPP_
#define COMMS_HPP_


class System;
//Send the model data (positions, maybe velocites)
//To the view
int comms_sendModelData(System& system);

//Get the model data from the model to the view
int comms_gatherModelData(float** pos_buf,float** hue_buf,int* atomCounts,int& numAtoms);

class ControllerState;
//Send/Recv the controller state from the view to model
int comms_sendControllerState(ControllerState* state);
int comms_recvControllerState(ControllerState* state);
//int comms_bcastControllerState(ControllerState* state);
#endif //COMMS_HPP_
