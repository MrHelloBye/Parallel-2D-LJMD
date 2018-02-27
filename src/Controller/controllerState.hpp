#ifndef CONTROLLER_STATE_HPP_
#define CONTROLLER_STATE_HPP_

#include <ostream>

struct ControllerState{
  float cursorPos[2] = {0,0};
  float trigger =0;
  unsigned int reset = false;

  friend std::ostream& operator<< (std::ostream& out, const ControllerState& state); 

};

#endif //CONTROLLER_STATE_HPP_
