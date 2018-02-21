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


/*********************************************************************
 * Generate a Circles struct of numCircles circles each with numSide sides.
 * Each circle is given a position in a lattice and a color assigned by ID.
 * vertex, position, and color vbos are also created
 *********************************************************************/
class Circles
{
  private:

    GLfloat radius;
    GLint numVertices;
    GLint numCircles;

    bool drawInstanced;

    //GPU side vbo's
    GLuint vertices_vbo = 0;
    GLuint pos_vbo = 0;
    GLuint colors_vbo = 0;

    GLuint vao = 0;


    //CPU side memory
    GLfloat *vertices = NULL;
    GLfloat *pos = NULL;
    GLfloat *colors = NULL;

    //For non-instanced rendering
    GLfloat *vertices_copies = NULL;
    GLfloat *pos_copies = NULL;
    GLfloat *colors_copies = NULL;
    GLint *first = NULL;
    GLsizei *count = NULL;

    //Helper Function to Generate Circles
    void verticesFan();
    
    //Private OpenGL Related
    int createBuffers();
    int createInstancedBuffers();

    void updateVerticesBuffer();
    void updatePosBuffer();
    void updateColorsBuffer();
  public:
    Circles(GLfloat &&radius, GLint &&numVertices, GLint &&numCircles,bool drawInstanced=false);
    ~Circles();

    //Manipulating positions
    void setPosLattice();
    void setPos(GLfloat *newPos);
    void movePos(GLfloat &&dx, GLfloat &&dy);

    //Manipulating colors
    void setColorsID();
    void setColors(GLfloat *newPos);

    //More OpenGL related
    int draw();
    int initShaders(GLuint &shader_program);
};
