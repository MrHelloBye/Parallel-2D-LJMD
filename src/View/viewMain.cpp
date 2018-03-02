//
//  main.cpp
//  OpenGL Schrodinger
//
//  Created by Liam Clink on 2/7/18.
//
//

// OpenGL practice
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>
#include <unistd.h>
#include <mpi.h>

#include "circles.hpp"
#include "../comms.hpp"
#include "../Controller/controller.hpp"

void reshape(int, int);
void check_GLSL_compile(GLuint shader);
void check_GLSL_link(GLuint shader_program);
void initShaders(GLuint&);

#define CIRCLE_SPEED 0.1

void MessageCallback( GLenum source,
                      GLenum type,
                      GLuint id,
                      GLenum severity,
                      GLsizei length,
                      const GLchar* message,
                      const void* userParam )
{
  fprintf( stderr, "GL CALLBACK: %s type = 0x%x, severity = 0x%x, message = %s\n",
           ( type == GL_DEBUG_TYPE_ERROR ? "** GL ERROR **" : "" ),
            type, severity, message );
}


int viewMain(int argc,char** argv)
{
  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  /**********************************************************************
   *  Load in some commandline arguments
   **********************************************************************/
  bool drawInstanced = false;

  int c;
  while((c = getopt (argc, argv, "i")) != -1){
    switch(c){
      case 'i':
        drawInstanced=true;
        break;
      case '?':
        if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr,
              "Unknown option character `\\x%x'.\n",
              optopt);
        return 1;
      default:
        abort ();
    }
  }

  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit())
  {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
    return 1;
  }

  // uncomment these lines if on Apple OS X
  if(drawInstanced){
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  }

  GLFWwindow* window = glfwCreateWindow(800, 500, "Particle Demo", NULL, NULL);
  if (!window)
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");

  glfwMakeContextCurrent(window);

  //--------------------------------------------------------//

  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  // Initialize GLEW, must be done after window creation
  GLenum err = glewInit();
  if (GLEW_OK != err)
  {
    /* Problem: glewInit failed, something is seriously wrong. */
    fprintf(stderr, "GLEW Error: %s\n", glewGetErrorString(err));
  }
  fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));




  // get version info
  const GLubyte* renderer = glGetString(GL_RENDERER); // get renderer string
  const GLubyte* version = glGetString(GL_VERSION); // version as a string
  printf("Renderer: %s\n", renderer);
  printf("OpenGL version supported %s\n", version);

  // During init, enable debug output
  //glEnable              ( GL_DEBUG_OUTPUT );
  //glDebugMessageCallback( (GLDEBUGPROC) MessageCallback, 0 );

  //--------------------------------------------------------//
  //Initialize the circles, give them colors based on ID,
  //Position them in a Lattice
  const int INIT_NUM_CIRCLES = 100;

  int numAtoms = INIT_NUM_CIRCLES;
  Circles circles(0.01, 12,numAtoms,drawInstanced);
  circles.setColorsID();
  circles.setPosLattice();
  Circles cursor(0.01, 12, 1,drawInstanced);
  GLfloat clrDum[] ={1,1,1};
  cursor.setColors(clrDum);


  //MPI Buffers
  //--------------------------------------------------------//
  float* pos_buf = new float[numAtoms*2]; 
  float* hue_buf = new float[numAtoms]; 
  int atomCounts[nprocs];

  Controller controller;
  ControllerState& state=controller.getState();


  int width, height;
  //Resize viewport to match window
  glfwGetFramebufferSize(window, &width, &height);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0f, width, height, 0.0f, 0.0f, 1.0f);

  //--------------------------------------------------------//

  // Initialize timer
  double t, dt, t_old;
  t_old = glfwGetTime() - 0.01;

  double t_lastControllerEcho = t_old;

  //Draw in a loop
  while(!glfwWindowShouldClose(window))
  {
    // wipe the drawing surface clear
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    //Update the timer, determine dt
    t = glfwGetTime();
    dt = t-t_old;
    t_old = t;

    if(t - t_lastControllerEcho > 1.0){
      t_lastControllerEcho = t;
      //std::cout <<"View's Controller: "<< state<<std::endl;
    }

    //Move Positions
    //circles.movePos(dt*CIRCLE_SPEED, dt*CIRCLE_SPEED);

    //Get model data
    comms_gatherModelData(&pos_buf,&hue_buf,atomCounts,numAtoms);

    //Send controller data
    controller.readState(dt);
    cursor.setEPos(state.cursorPos);
    cursor.setRadius(0.02 + 0.015*state.trigger);
    float newHSL[] = { 180*(1+state.bumper),.75,.75}; 
    cursor.setHSLs(newHSL);
    controller.commState();


    circles.setPosAndHues(pos_buf,hue_buf,numAtoms);

    //Draw the circles 
    circles.draw();
    cursor.draw();

    // put the stuff we've been drawing onto the display
    glfwSwapBuffers(window);

    // update other events like input handling
    glfwPollEvents();

    if(glfwGetKey(window, GLFW_KEY_ESCAPE)) {
      glfwSetWindowShouldClose(window, 1);
    }
  }

  glfwTerminate();
  return 0;
}

void reshape(int w, int h)
{
  glViewport(0, 0, w, h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0, w, h, 0, 0, 0);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

