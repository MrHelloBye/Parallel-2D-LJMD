
#include <cstring>
#include <iostream>
#include <GLFW/glfw3.h>

#include "controllerState.hpp"
#include "controller.hpp"
#include "../comms.hpp"


const char XBONE_CONTROLLER_NAME[] = "Xbox Wireless Controller";

Controller::Controller(){
  int i;
  for(i = GLFW_JOYSTICK_1; i <= GLFW_JOYSTICK_LAST;i++){
    int present = glfwJoystickPresent(joystickID);
    if(present){
      joystickID = i;

      const char* name = glfwGetJoystickName(joystickID);

      if(strcmp(name,XBONE_CONTROLLER_NAME) == 0){
        //Found the xbone controller!
        delete name;
        break;
      }

    }
  }

  if(i > GLFW_JOYSTICK_LAST){
    const char* name = glfwGetJoystickName(joystickID);
    std::cerr <<"WARNING: Did not find: '"<<XBONE_CONTROLLER_NAME<< "' using '"
      << name << "'"<<std::endl;
  }
}


//========================================================================
//========================================================================
int Controller::readState(){

  //Check that there is a joystick
  int present = glfwJoystickPresent(joystickID);
  if(present == GLFW_TRUE){

    //**************************************************
    //Check the axes
    //**************************************************
    int j_count;
    const float* j_axes = glfwGetJoystickAxes(joystickID, &j_count);

    float leftX,leftY,rightX,rightY,leftTrigger,rightTrigger,dpadX,dpadY;

    leftX = j_axes[0];
    leftY = j_axes[1];
    rightX = j_axes[2];
    rightY = j_axes[3];

    rightTrigger = j_axes[4];
    leftTrigger = j_axes[5];

    dpadX = j_axes[6];
    dpadY = j_axes[7];

    state.cursorPos[0] = leftX;
    state.cursorPos[1] = leftY;

    state.trigger = rightTrigger;

    //**************************************************
    //Check the Buttons
    //**************************************************

    int count;
    const unsigned char* axes = glfwGetJoystickButtons(joystickID, &count);

    unsigned char aButton, bButton, xButton,yButton,menuButton,
                  leftBumper,rightBumper,leftStick,rightStick;

    aButton = axes[0];
    bButton = axes[1];
    xButton = axes[3];
    yButton = axes[4];

    menuButton = axes[11];

    leftBumper = axes[6];
    rightBumper = axes[7];

    leftStick = axes[13];
    rightStick = axes[14];


    if(menuButton == GLFW_PRESS){
      state.reset = true;
    } else {
      state.reset = false;
    }
  }
  return 0;
}

int Controller::commState(){
  return comms_sendControllerState(&state);
}
