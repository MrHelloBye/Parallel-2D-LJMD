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
#include <string>

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>

#include "hsl_to_rgb.hpp"
#include "circles.hpp"
#include "loadShaders.hpp"



/**************************************************
 * Generate a circle as a fan of vertices
 **************************************************/
void Circles::verticesFan()
{
  GLfloat twoPi = 2.0f*M_PI;

  vertices[0] = 0;
  vertices[1] = 0;

  int numSides = numVertices-2;
  for ( int i = 1; i < numVertices; i++ )
  {
    vertices[i*2]     = (radius * cos(i*twoPi/numSides));
    vertices[i*2+1] = ( radius * sin( i * twoPi / numSides ) );
  }

  updateVerticesBuffer();
}

/**************************************************
 * Create the Buffers used by OpenGL, and bind them
 **************************************************/
//When instanced drawing isn't available
int Circles::createBuffers(){

  //Initialize Host memory
  vertices = new GLfloat[numVertices*2];
  pos = new GLfloat[numCircles*2];
  colors = new GLfloat[numCircles*3];

  vertices_copies = new GLfloat[numVertices*numCircles*2];
  pos_copies = new GLfloat[numVertices*numCircles*2];
  colors_copies = new GLfloat[numVertices*numCircles*3];

  //Initialize buffer objects
  glGenBuffers(1, &vertices_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glBufferData(GL_ARRAY_BUFFER, numVertices*numCircles*2*sizeof(GLfloat),
      vertices_copies, GL_STATIC_DRAW);

  glGenBuffers(1, &pos_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glBufferData(GL_ARRAY_BUFFER, numVertices*numCircles*2*sizeof(GLfloat),
      pos_copies, GL_DYNAMIC_DRAW);

  glGenBuffers(1, &colors_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glBufferData(GL_ARRAY_BUFFER, numVertices*numCircles*3*sizeof(GLfloat),
      colors_copies, GL_DYNAMIC_DRAW);

  //Generate the vao
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  //Bind the vbos to the vao
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);

  //Enable them for the vertex shader
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  //Make some CPU side array for specifying where circles start How many
  //verticies each circle has
  first = new GLint[numCircles];
  count = new GLsizei[numCircles];
  for(int i =0; i < numCircles;i++){
    first[i] = i*numVertices;
    count[i] = numVertices;
  }

  return 0;//Success
}
int Circles::createInstancedBuffers(){

  //Initialize Host memory
  vertices = new GLfloat[numVertices*2];
  pos = new GLfloat[numCircles*2];
  colors = new GLfloat[numCircles*3];

  //Initialize buffer objects
  glGenBuffers(1, &vertices_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glBufferData(GL_ARRAY_BUFFER, numVertices*2*sizeof(GLfloat),
      vertices, GL_STATIC_DRAW);

  glGenBuffers(1, &pos_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glBufferData(GL_ARRAY_BUFFER, numCircles*2*sizeof(GLfloat),
      pos, GL_DYNAMIC_DRAW);

  glGenBuffers(1, &colors_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glBufferData(GL_ARRAY_BUFFER, numCircles*3*sizeof(GLfloat),
      colors, GL_DYNAMIC_DRAW);


  //Generate the vao
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  //Bind the vbos to the vao
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glVertexAttribDivisor(1, 1); //Increment the pos with each circle
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribDivisor(2, 1); //Increment the color with each circle
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);

  //Enable them for the vertex shader
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  return 0;//Success
}

int Circles::deleteBuffers(){

  GLuint vbos[] = {vertices_vbo,pos_vbo,colors_vbo};
  glDeleteBuffers(3,vbos);

  delete vertices;
  delete pos;
  delete colors;
  if(!drawInstanced){
    delete vertices_copies;
    delete pos_copies;
    delete colors_copies;
    delete first;
    delete count;
  }
}
int Circles::resizeNumCircles(int newNumCircles){
  //Delete the buffers first
  deleteBuffers();

  numCircles = newNumCircles;

  //Initialize the Buffers
  if(drawInstanced){
    createInstancedBuffers();
  } else{
    createBuffers();
  }

  //Create one circle in a vertex fan, one point in the middle
  verticesFan();
}


/**************************************************
 * Send the CPU side buffers to the device, for drawing
 **************************************************/
void Circles::updateVerticesBuffer(){
  if(drawInstanced){
    glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        2*sizeof(GLfloat)*numVertices, vertices);
  } else{
    //Make copies for each circle
    for(int c =0; c < numCircles;c++){
      for(int v =0; v < numVertices;v++){
        for(int d =0; d < 2;d++){
          vertices_copies[d + 2*(v + numVertices*c)] = vertices[d + 2*v];
        }
      }
    }
    glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        2*sizeof(GLfloat)*numVertices*numCircles, vertices_copies);
  }
}

void Circles::updatePosBuffer(){
  if(drawInstanced){
    glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        2*sizeof(GLfloat)*numCircles, pos);
  } else {
    //Make copies for each circle
    for(int c =0; c < numCircles;c++){
      for(int v =0; v < numVertices;v++){
        for(int d =0; d < 2;d++){
          pos_copies[d + 2*(v + numVertices*c)] = pos[d + 2*c];
        }
      }
    }
    glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        2*sizeof(GLfloat)*numVertices*numCircles, pos_copies);
  }
}

void Circles::updateColorsBuffer(){
  if(drawInstanced){
    glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        3*sizeof(GLfloat)*numCircles, colors);
  } else {
    //Make copies for each circle
    for(int c =0; c < numCircles;c++){
      for(int v =0; v < numVertices;v++){
        for(int d =0; d < 3;d++){
          colors_copies[d + 3*(v + numVertices*c)] = colors[d + 3*c];
        }
      }
    }
    glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 
        3*sizeof(GLfloat)*numVertices*numCircles, colors_copies);
  }
}

/**************************************************
 * Initialize the shader for the circles
 **************************************************/
int Circles::initShaders(){

  shader_program = glCreateProgram();
  //Load the correct shader from file
  std::string shadersPath = SHADERS_PATH;
  if(drawInstanced){
    std::cout<<"Using glsl 4.00 with instanced drawing"<<std::endl;
    loadShaders(shader_program,shadersPath+"/circle4.00.vert",
        shadersPath+"/circle4.00.frag");
  } else{
    std::cout<<"Using glsl 1.10 without instanced drawing"<<std::endl;
    loadShaders(shader_program,shadersPath+"/circle1.10.vert",
        shadersPath+"/circle1.10.frag");
  }

  //Bind the arguments of the shader
  glBindAttribLocation(shader_program, 0, "vertex_position");
  glBindAttribLocation(shader_program, 1, "circle_position");
  glBindAttribLocation(shader_program, 2, "circle_color");

  glLinkProgram(shader_program);

  check_GLSL_link(shader_program);
}

/*********************************************************************
 * Contructors and Destructor
 *********************************************************************/
Circles::Circles(GLfloat radius, GLint numVertices, GLint numCircles,bool drawInstanced):
  radius(radius),numVertices(numVertices),numCircles(numCircles),
  drawInstanced(drawInstanced)
{

  //Initialize the Buffers
  if(drawInstanced){
    createInstancedBuffers();
  } else{
    createBuffers();
  }

  //Create one circle in a vertex fan, one point in the middle
  verticesFan();

  //Initilize the shader
  initShaders();
  glUseProgram(shader_program);
}

Circles::~Circles(){
  deleteBuffers();
}

/*********************************************************************
 * Manipulating vertices
 *********************************************************************/
void Circles::setRadius(float radius){
  this->radius = radius;

  verticesFan();

  updateVerticesBuffer();
}

/*********************************************************************
 * Manipulating positions
 *********************************************************************/
void Circles::setPosLattice(){
  GLint numRow = (GLint) sqrt(numCircles);
  for(int i = 0; i < numCircles; i++){
    pos[i*2]   = 2*((i%numRow)*1.)/numRow -1;
    pos[i*2+1] = 2*((i/numRow)*1.)/numRow -1;
  }

  updatePosBuffer();
}

void Circles::movePos(GLfloat dx, GLfloat dy)
{
  for(int i = 0; i < numCircles; i++){
    pos[i*2]   += dx;
    //pos[i*2+1] += dy;
    //Use this line for more fun
    pos[i*2+1] += (i%2?-1:1)*dy;
  }

  //Periodic Boundaries
  for(int i = 0; i < numCircles*2; i++){
    if(pos[i] >  1) pos[i] -= 2;
    if(pos[i] < -1) pos[i] += 2;
  }

  updatePosBuffer();
}

void Circles::setPos(float* newPos)
{
  for(int i = 0; i < numCircles*2; i++){
//	  std::cout<<"ey: "<<i<< " "<<newPos << " "<< newPos[i]<<std::endl;
    pos[i] = newPos[i]/19-1;
  }

  updatePosBuffer();
}


void Circles::setEPos(float* newPos)
{
  for(int i = 0; i < numCircles*2; i++){
    pos[i] = newPos[i];
  }

  updatePosBuffer();
}

/*********************************************************************
 * Manipulating Colors
 *********************************************************************/


//Set the colors with different Hues, but the same saturation and lightness
void Circles::setColorsID(){
  /*for(int i = 0; i < numCircles; i++){
    colors[i*3]   = .25 + .75*i/(numCircles-1.0); //Mostly red
    colors[i*3+1] = i/(numCircles-1.0); //Some green
    colors[i*3+2] = 1 -i/(numCircles-1.0); //Inverse blue
    }*/
  float hsl[3];
  hsl[1] = 1;hsl[2] = .5;
  for(int i = 0; i < numCircles; i++){
    hsl[0] = i*360./(numCircles-1.);
    hsl_to_rgb(colors+i*3,hsl);
  }

  updateColorsBuffer();
}

void Circles::setColors(GLfloat* newColors)
{
  for(int i = 0; i < numCircles*3; i++){
    colors[i] = newColors[i];
  }

  updateColorsBuffer();
}

void Circles::setHues(float* newHues)
{
  float hsl[3];
  float rgb[3];
  hsl[1] = 1;hsl[2] = .5;
  for(int i = 0; i < numCircles; i++){
    hsl[0] = newHues[i];
    hsl_to_rgb(rgb,hsl);
    for(int j =0; j < 3; j++){
      colors[i*3+j] = rgb[j];
    }
  }

  updateColorsBuffer();
}

void Circles::setPosAndHues(float* newPos,float* newHues, int newSize)
{
  if(newSize != numCircles) {
    resizeNumCircles(newSize);
  }
  setPos(newPos);
  setHues(newHues);
}

void Circles::setHSLs(float* newHSLs)
{
  float rgb[3];
  for(int i = 0; i < numCircles; i++){
    hsl_to_rgb(rgb,newHSLs+i*3);
    for(int j =0; j < 3; j++){
      colors[i*3+j] = rgb[j];
    }
  }

  updateColorsBuffer();
}

/*********************************************************************
 * Draw the Circles
 *********************************************************************/
int Circles::draw(){
  //If we end up using multiple shaders, we'll need this
  //glUseProgram(shader_program);

  glBindVertexArray(vao);
  // draw points from the currently bound VAO with current in-use shader
  if(drawInstanced){
    glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVertices,numCircles);
  }
  else{
    glMultiDrawArrays(GL_TRIANGLE_FAN,first,count,numCircles); 
    /*for(int i =0;i < numCircles;i++){
      glDrawArrays(GL_TRIANGLE_FAN, i*numVertices, numVertices);
      }*/
  }
  //glDrawArrays(GL_TRIANGLE_FAN, 0, numVertices);
  return 0;

}


