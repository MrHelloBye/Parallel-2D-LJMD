#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>


using namespace std;

#include "loadShaders.hpp"
string readShaderSource(const string& filename) {
  ifstream file;
  file.open(filename.c_str());

  if (!file) {
    return(0);
  }

  std::stringstream stream;

  stream << file.rdbuf();

  file.close();

  return stream.str();
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

int loadShaders(GLuint &shader_program,
    const string& vertFilename,
    const string& fragFilename)
{
  string vert_source = readShaderSource(vertFilename);
  string frag_source = readShaderSource(fragFilename);
  const char* vert_src = vert_source.c_str();
  const char* frag_src = frag_source.c_str();

  //cout <<vert_src;

  //Compile shaders
  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vs, 1, &vert_src, NULL);
  glCompileShader(vs);

  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fs, 1, &frag_src, NULL);
  glCompileShader(fs);

  cout <<"Checking frag shader compilation"<<endl;
  check_GLSL_compile(vs);
  cout <<"Checking vert shader compilation"<<endl;
  check_GLSL_compile(fs);

  //Link shaders into program
  glAttachShader(shader_program, vs);
  glAttachShader(shader_program, fs);

}



