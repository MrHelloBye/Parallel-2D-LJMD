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


struct Circle
{
    GLuint points_vbo = 0;
    GLuint vao = 0;
    int numVertices = 2;
    GLfloat radius = 1;
    GLfloat *vertices;
};

void circleFan(GLfloat*, GLfloat, GLfloat, GLfloat, GLfloat, GLint&&);
void generateCircle(Circle&, GLfloat &&x, GLfloat &&y, GLfloat &&z, GLint&&);
void moveCircle(Circle&, GLfloat, GLfloat, GLfloat);
void reshape(int, int);
void check_GLSL_compile(GLuint shader);
void check_GLSL_link(GLuint shader_program);
void initShaders(GLuint&);

int main()
{
    // start GL context and O/S window using the GLFW helper library
    if (!glfwInit())
    {
        fprintf(stderr, "ERROR: could not start GLFW3\n");
        return 1;
    }
    
    // uncomment these lines if on Apple OS X
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    //glfwWindowHint(G#LFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    //glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
    //glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    
    GLFWwindow* window = glfwCreateWindow(500, 500, "Particle Demo", NULL, NULL);
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
        fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
    }
    fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    
    // get version info
    const GLubyte* renderer = glGetString(GL_RENDERER); // get renderer string
    const GLubyte* version = glGetString(GL_VERSION); // version as a string
    printf("Renderer: %s\n", renderer);
    printf("OpenGL version supported %s\n", version);
    
    //--------------------------------------------------------//
    
    Circle circle1;
    circle1.radius = 0.2;
    generateCircle(circle1, 0.25, 0.5, 0, 100);
    Circle circle2;
    circle2.radius = 0.2;
    generateCircle(circle2, 0.25, 0.5, 0, 100);
    
    //--------------------------------------------------------//
    
    GLuint shader_program = glCreateProgram();
    initShaders(shader_program);
    
    
    //Only render the front side of the face
    //glEnable(GL_CULL_FACE); // cull face
    //glCullFace(GL_BACK); // cull back face
    //glFrontFace(GL_CW); // GL_CCW for counter clock-wise
    
    //glClearColor(1.0f, 0.0f, 0.0f, 1.0f);
    
    
    int width, height;
    //Resize viewport to match window
    glfwGetFramebufferSize(window, &width, &height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(0.0f, width, height, 0.0f, 0.0f, 1.0f);
    
    //--------------------------------------------------------//
    
    GLfloat x1 = 0;
    GLfloat x2 = -1;
    
    
    glUseProgram(shader_program);
    //Draw in a loop
    while(!glfwWindowShouldClose(window))
    {
        // wipe the drawing surface clear
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        x1 += 0.01;
        x2 -= 0.03;
        if (x1 > 1) x1 = -1;
        if (x2 < -1) x2 = 1;
        moveCircle(circle1, x1, 0, 0);
        moveCircle(circle2, x2, 0.5, 0);
        
        glBindVertexArray(circle1.vao);
        // draw points from the currently bound VAO with current in-use shader
        glDrawArrays(GL_TRIANGLE_FAN, 0, circle1.numVertices);
        
        glBindVertexArray(circle2.vao);
        glDrawArrays(GL_TRIANGLE_FAN, 0, circle2.numVertices);
        
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

void circleFan(GLfloat *allCircleVertices, GLfloat x, GLfloat y,
               GLfloat z, GLfloat radius, GLint &&numberOfSides)
{
    int numberOfVertices = numberOfSides + 2;
    
    GLfloat twoPi = 2.0f*M_PI;
    
    allCircleVertices[0] = x;
    allCircleVertices[1] = y;
    allCircleVertices[2] = z;
    
    for ( int i = 1; i < numberOfVertices; i++ )
    {
        allCircleVertices[i*3] = x + (radius * cos(i*twoPi/numberOfSides));
        allCircleVertices[(i*3)+1] = y + ( radius * sin( i * twoPi / numberOfSides ) );
        allCircleVertices[(i*3)+2] = z;
    }
}

void generateCircle(Circle &circle, GLfloat &&x, GLfloat &&y, GLfloat &&z,
                    GLint &&numSides)
{
    circle.numVertices = numSides+2;
    circle.vertices = new GLfloat[circle.numVertices*3];
    circleFan(circle.vertices, x, y, z, circle.radius, circle.numVertices-2);
    
    GLfloat *colors = new GLfloat[circle.numVertices*3];
    for (int i = 0; i<circle.numVertices*3; i++)
        colors[i] = 0.5;
    
    //Initialize vertex buffer object
    glGenBuffers(1, &circle.points_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, circle.points_vbo);
    glBufferData(GL_ARRAY_BUFFER, circle.numVertices*3*sizeof(GLfloat),
                 circle.vertices, GL_STATIC_DRAW);
    
    GLuint colors_vbo = 0;
    glGenBuffers(1, &colors_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
    glBufferData(GL_ARRAY_BUFFER, 3*circle.numVertices*sizeof(float), colors, GL_STATIC_DRAW);
    
    glGenVertexArrays(1, &circle.vao);
    glBindVertexArray(circle.vao);
    glBindBuffer(GL_ARRAY_BUFFER, circle.points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3*sizeof(GLfloat), NULL);
    glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
}

void moveCircle(Circle &circle, GLfloat x, GLfloat y, GLfloat z)
{
    circleFan(circle.vertices, x, y, z, circle.radius, circle.numVertices-2);
    glBindBuffer(GL_ARRAY_BUFFER, circle.points_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, 3*sizeof(GLfloat)*circle.numVertices, circle.vertices);
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

void check_GLSL_compile(GLuint shader)
{
    // Check GLSL compile status
    GLint status = 0; glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    std::cout << "Shader Compile status: " << status << std::endl;
    
    if (!status)
    {
        // Print GLSL compile log
        GLint logLen = 0; glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &logLen);
        GLchar* logStr = new GLchar[logLen];
        glGetShaderInfoLog(shader, logLen, 0, logStr); // report logStr delete[] logStr;
        std::cout << "Compile logStr: " << logStr << std::endl;
    }
}

void check_GLSL_link(GLuint shader_program)
{
    // Check GLSL compile status
    GLint status = 0; glGetProgramiv(shader_program, GL_LINK_STATUS, &status);
    std::cout << "Shader Link status: " << status << std::endl;
    
    if (!status)
    {
        // Print GLSL compile log
        GLint logLen = 0; glGetProgramiv(shader_program, GL_INFO_LOG_LENGTH, &logLen);
        GLchar* logStr = new GLchar[logLen];
        glGetShaderInfoLog(shader_program, logLen, 0, logStr); // report logStr delete[] logStr;
        std::cout << "Link logStr: " << logStr << std::endl;
    }
}

void initShaders(GLuint &shader_program)
{
#ifndef RPI
    //Define vertex shader (better to load from text file)
    const char* vertex_shader =
    "#version 400\n"
    "in vec3 vertex_position;"
    "in vec3 vertex_color;"
    "out vec3 color;"
    "void main() {"
    "  color = vertex_color;"
    "  gl_Position = vec4(vertex_position, 1.0);"
    "}";
    
    //Define fragment shader (for drawing surfaces)
    const char* fragment_shader =
    "#version 400\n"
    "in vec3 color;"
    "out vec4 frag_color;"
    "void main() {"
    "  frag_color = vec4(color, 1.0);"
    "}";
#else
    static const char* vertex_shader=
    "#version 110\n"
    "attribute vec3 vertex_color;\n"
    "attribute vec3 vertex_position;\n"
    "varying vec3 color;\n"
    "void main()\n"
    "{\n"
    "    gl_Position = vec4(vertex_position,  1.0);\n"
    "    color = vertex_color;\n"
    "}\n";

    static const char* fragment_shader=
    "#version 110\n"
    "varying vec3 frag_color;\n"
    "void main()\n"
    "{\n"
    "    gl_FragColor = vec4(frag_color, 1.0);\n"
    "}\n";

 #endif
    
    //Compile shaders
    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertex_shader, NULL);
    glCompileShader(vs);
    
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragment_shader, NULL);
    glCompileShader(fs);
    
    check_GLSL_compile(vs);
    check_GLSL_compile(fs);
    
    //Link shaders into program
    glAttachShader(shader_program, vs);
    glAttachShader(shader_program, fs);
    
    // insert location binding code here
    glBindAttribLocation(shader_program, 0, "vertex_position");
    glBindAttribLocation(shader_program, 1, "vertex_color");
    
    glLinkProgram(shader_program);
    
    check_GLSL_link(shader_program);
}
