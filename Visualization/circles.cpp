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

#include "circles.hpp"
#include "loadShaders.hpp"



/**************************************************
 * Helper Function Generate Circles
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
 * Private OpenGL related
 **************************************************/
int Circles::createBuffers(){
  vertices = new GLfloat[numVertices*2];
  pos = new GLfloat[numCircles*2];
  colors = new GLfloat[numCircles*3];

  vertices_copies = new GLfloat[numVertices*numCircles*2];
  pos_copies = new GLfloat[numVertices*numCircles*2];
  colors_copies = new GLfloat[numVertices*numCircles*3];

  //Initialize vertex buffer object
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


  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  //vertices_vbo
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  //pos_vbo
  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  //colors_vbo
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  first = new GLint[numCircles];
  count = new GLsizei[numCircles];
  for(int i =0; i < numCircles;i++){
    first[i] = i*numVertices;
    count[i] = numVertices;
  }

  return 0;//Success
}
int Circles::createInstancedBuffers(){
  vertices = new GLfloat[numVertices*2];
  pos = new GLfloat[numCircles*2];
  colors = new GLfloat[numCircles*3];

  //Initialize vertex buffer object
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


  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);

  //vertices_vbo
  glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
  glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  //pos_vbo
  glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
  glVertexAttribDivisor(1, 1);
  glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

  //colors_vbo
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribDivisor(2, 1);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  return 0;//Success
}

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


/*********************************************************************
 * Contructors and Destructor
 *********************************************************************/
Circles::Circles(GLfloat &&radius, GLint &&numVertices, GLint &&numCircles,bool drawInstanced):
  radius(radius),numVertices(numVertices),numCircles(numCircles),drawInstanced(drawInstanced)
{

  if(drawInstanced){
    createInstancedBuffers();
  } else{
    createBuffers();
  }

  //Create one circle in a vertex fan, one point in the middle
  verticesFan();

}

Circles::~Circles(){
  delete vertices;
  delete pos;
  delete colors;
  if(!drawInstanced)
    delete vertices_copies;
  delete pos_copies;
  delete colors_copies;
  delete first;
  delete count;
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

void Circles::movePos(GLfloat &&dx, GLfloat &&dy)
{
  for(int i = 0; i < numCircles; i++){
    pos[i*2]   += dx;
    pos[i*2+1] += (i%2?-1:1)*dy;
  }

  //Periodic Boundaries
  for(int i = 0; i < numCircles*2; i++){
    if(pos[i] >  1) pos[i] -= 2;
    if(pos[i] < -1) pos[i] += 2;
  }

  updatePosBuffer();
}

void Circles::setPos(GLfloat* newPos)
{
  for(int i = 0; i < numCircles*2; i++){
    pos[i] = pos[i];
  }

  updatePosBuffer();
}

/*********************************************************************
 * Manipulating Colors
 *********************************************************************/
void hsl_to_rgb( GLfloat* rgb,GLfloat* hsl){
  GLfloat c = (1 - fabs(2*hsl[2] -1))*hsl[1];
  GLfloat hp = hsl[0]/60;
  GLfloat x= c*(1 - fabs(fmod(hp,2) -1));
  if(hsl[0] != hsl[0] || hp < 0){
    rgb[0]=0;rgb[1]=0;rgb[2]=0;
  }
  else if(hp <= 1){
    rgb[0]=c;rgb[1]=x;rgb[2]=0;
  }
  else if(hp <= 2){
    rgb[0]=x;rgb[1]=c;rgb[2]=0;
  }
  else if(hp <= 3){
    rgb[0]=0;rgb[1]=c;rgb[2]=x;
  }
  else if(hp <= 4){
    rgb[0]=0;rgb[1]=x;rgb[2]=c;
  }
  else if(hp <= 5){
    rgb[0]=x;rgb[1]=0;rgb[2]=c;
  }
  else if(hp <= 6){
    rgb[0]=c;rgb[1]=0;rgb[2]=x;
  }
  GLfloat m = hsl[2] - 0.5*c;
  rgb[0] += m; rgb[1] += m; rgb[2] += m;
}


void Circles::setColorsID(){
  /*for(int i = 0; i < numCircles; i++){
    colors[i*3]   = .25 + .75*i/(numCircles-1.0); //Mostly red
    colors[i*3+1] = i/(numCircles-1.0); //Some green
    colors[i*3+2] = 1 -i/(numCircles-1.0); //Inverse blue
    }*/
  GLfloat hsl[3];
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
    colors[i] = colors[i];
  }

  updateColorsBuffer();
}


int Circles::draw(){
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


int Circles::initShaders(GLuint &shader_program)
{
  if(drawInstanced){
    std::cout<<"Using glsl 4.00 with instanced drawing"<<std::endl;
    loadShaders(shader_program,"shaders/vert4.00.glsl","shaders/frag4.00.glsl");
  } else{
    std::cout<<"Using glsl 1.10 without instanced drawing"<<std::endl;
    loadShaders(shader_program,"shaders/vert1.10.glsl","shaders/frag1.10.glsl");
  }

  // insert location binding code here
  glBindAttribLocation(shader_program, 0, "vertex_position");
  glBindAttribLocation(shader_program, 1, "circle_position");
  glBindAttribLocation(shader_program, 2, "circle_color");

  glLinkProgram(shader_program);

  check_GLSL_link(shader_program);
}

