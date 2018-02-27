#ifndef CONTROLLER_HPP_
#define CONTROLLER_HPP_

#include "controllerState.hpp"

class Controller{
  private:
    ControllerState state;
    int joystickID;

  public:
    Controller();

    int readState();
    int commState();

    ControllerState& getState(){return state;}
};

#endif //CONTROLLER_HPP_
