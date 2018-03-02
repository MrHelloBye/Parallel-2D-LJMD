#include <ostream>

#include "controllerState.hpp"

std::ostream& operator<<(std::ostream& out, const ControllerState& state) {
    return out << "cursorPos ("<<state.cursorPos[0] << ","<<state.cursorPos[1]<<") "
               << "trigger: "<<state.trigger << "bumper: "<<state.bumper
               <<" reset: "<< ((state.reset)?"true":"false");
}
