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
	glVertexAttribDivisor(1, numVertices);
	glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 2*sizeof(GLfloat), NULL);

	//colors_vbo
	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
	glVertexAttribDivisor(2, numVertices);
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	glEnableVertexAttribArray(2);

  return 0;//Success
}

void Circles::updateVerticesBuffer(){
	glBindBuffer(GL_ARRAY_BUFFER, vertices_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 
			2*sizeof(GLfloat)*numVertices, vertices);
}

void Circles::updatePosBuffer(){
	glBindBuffer(GL_ARRAY_BUFFER, pos_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 
			2*sizeof(GLfloat)*numCircles, pos);
}

void Circles::updateColorsBuffer(){
	glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
	glBufferSubData(GL_ARRAY_BUFFER, 0, 
			3*sizeof(GLfloat)*numCircles, colors);
}


/*********************************************************************
 * Contructors and Destructor
 *********************************************************************/
Circles::Circles(GLfloat &&radius, GLint &&numVertices, GLint &&numCircles):
	radius(radius),numVertices(numVertices),numCircles(numCircles)
{

  createBuffers();

	//Create one circle in a vertex fan, one point in the middle
	verticesFan();


}

Circles::~Circles(){
  delete vertices;
  delete pos;
  delete colors;
}

/*********************************************************************
 * Manipulating positions
 *********************************************************************/
void Circles::setPosLattice(){
	GLint numRow = (GLint) sqrt(numCircles);
	for(int i = 0; i < numCircles; i++){
		pos[i*2]   = 2*((i%numRow)*1.)/numRow -1;
		pos[i*2+1] = 2*((i/numRow)*1.)/numRow -1;
    std::cout<<pos[i*2]<<" , "<<pos[i*2]<<std::endl;
	}

	updatePosBuffer();
}

void Circles::movePos(GLfloat &&dx, GLfloat &&dy)
{
	for(int i = 0; i < numCircles; i++){
		pos[i*2]   += dx;
		pos[i*2+1] += dy;
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
void Circles::setColorsID(){
	for(int i = 0; i < numCircles; i++){
		colors[i*3]   = 1.0; //Mostly red
		colors[i*3+1] = i/(numCircles-1.0); //Some green
		colors[i*3+2] = 1 -i/(numCircles-1.0); //Inverse blue
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
	glDrawArraysInstanced(GL_TRIANGLE_FAN, 0, numVertices,numCircles);
	//glDrawArrays(GL_TRIANGLE_FAN, 0, numVertices);
  return 0;

}


int Circles::initShaders(GLuint &shader_program)
{
  loadShaders(shader_program,"shaders/vert1.10.glsl","shaders/frag1.10.glsl");
    
    // insert location binding code here
    glBindAttribLocation(shader_program, 0, "vertex_position");
    glBindAttribLocation(shader_program, 1, "circle_position");
    glBindAttribLocation(shader_program, 2, "circle_color");
    
    glLinkProgram(shader_program);
    
    check_GLSL_link(shader_program);
}

