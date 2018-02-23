#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(parallel-2D-LJMD, 0.01, forrestglines@gmail.com)
AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
AC_CONFIG_SRCDIR([src/main.cpp])
#AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC([mpicc])
AC_PROG_CXX([mpic++ mpicxx])


# Checks for libraries.
AC_CHECK_LIB([m],[sqrt])

#Determine if the View can be enabled
AC_DEFINE([VIEW_ENABLED],1)
AC_DEFINE_UNQUOTED([SHADERS_PATH],"`pwd`/src/View/shaders/")

AC_CHECK_LIB([glfw], [glfwInit])
AC_CHECK_LIB([GLEW], [glewInit])
AC_CHECK_LIB([GL], [glCreateShader])
AC_CHECK_LIB([GLU], [gluGetString])

# Checks for header files.
AC_CHECK_HEADERS([unistd.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_UINT64_T

# Checks for library functions.
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES([
 Makefile
 src/Makefile
])
AC_OUTPUT
