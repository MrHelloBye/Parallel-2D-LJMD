
#include <GLFW/glfw3.h>

GLfloat alpha = 210.f, beta = -70.f;
GLfloat zoom = 2.f;

//========================================================================
//========================================================================
const int ID_JOYSTICK = GLFW_JOYSTICK_2;
void joystick_handle(double dt_total){

  //Check that there is a joystick
  int present = glfwJoystickPresent(ID_JOYSTICK);
  if(present == GLFW_TRUE){

    //**************************************************
    //Check the axes
    //**************************************************
    int j_count;
    const float* j_axes = glfwGetJoystickAxes(ID_JOYSTICK, &j_count);

    float leftX,leftY,rightX,rightY,leftTrigger,rightTrigger,dpadX,dpadY;

    leftX = j_axes[0];
    leftY = j_axes[1];
    rightX = j_axes[2];
    rightY = j_axes[3];

    leftTrigger = j_axes[4];
    rightTrigger = j_axes[5];

    dpadX = j_axes[6];
    dpadY = j_axes[7];

    alpha += (GLfloat) dt*leftX * 60.f;
    beta += (GLfloat) dt*leftY * 60.f;

    alpha += (GLfloat) dt*rightX * 60.f;
    beta += (GLfloat) dt*rightY * 60.f;

    zoom -= (float) dt*leftTrigger * 0.5f;
    zoom += (float) dt*rightTrigger * 0.5f;

    alpha += (GLfloat) dt*dpadX * 60.f;
    beta += (GLfloat) dt*dpadY * 60.f;

    //**************************************************
    //Check the Buttons
    //**************************************************

    int count;
    const unsigned char* axes = glfwGetJoystickButtons(ID_JOYSTICK, &count);

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

    if(aButton == GLFW_PRESS)
      beta -= dt*60.f;
    if(bButton == GLFW_PRESS)
      alpha += dt*60.f;
    if(xButton == GLFW_PRESS)
      alpha -= dt*60.f;
    if(yButton == GLFW_PRESS)
      beta += dt*60.f;

    if(menuButton == GLFW_PRESS){
      alpha = 210.f;
      beta = -70.f;
      zoom = 2.f;
    }

    if(leftBumper == GLFW_PRESS){
      zoom += dt*0.5;
    }
    if(rightBumper == GLFW_PRESS){
      zoom -= dt*0.5;
    }

    if(leftStick == GLFW_PRESS){
      zoom += dt*0.5;
    }
    if(rightStick == GLFW_PRESS){
      zoom -= dt*0.5;
    }

    if (zoom < 0)
      zoom = 0;
  }

}


