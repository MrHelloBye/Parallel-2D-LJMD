#ifndef LOAD_SHADERS_HPP_
#define LOAD_SHADERS_HPP_

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <string>

#include <GL/glew.h> // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>

void check_GLSL_compile(GLuint shader);

void check_GLSL_link(GLuint shader_program);

int loadShaders(GLuint &shader_program,
    const std::string& vertFilename,
    const std::string& fragFilename);

#endif //LOAD_SHADERS_HPP_

