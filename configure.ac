#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.62])
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

#AC_CHECK_LIB([glfw], [glfwInit])
#AC_CHECK_LIB([glfw3], [glfwCreateWindow])
#AC_CHECK_LIB([GLEW], [glewInit])
#AC_CHECK_LIB([GL], [glCreateShader])


AC_ARG_ENABLE(view,
  [ --disable-view   Disable compilation of view and controller, compile only the model],
  [case "${enableval}" in
     yes | no ) WITH_VIEW="${enableval}" ;;
     *) AC_MSG_ERROR(bad value ${enableval} for --disable-view) ;;
   esac],
  [WITH_VIEW="yes"]
)


AM_CONDITIONAL([WITH_VIEW], [test "x$WITH_VIEW" = "xyes"])

# Define VIEW in config.h if we're going to compile against it
if test "x$WITH_VIEW" = "xyes"; then
    AC_DEFINE([VIEW_ENABLED], 1)
	AC_DEFINE_UNQUOTED([SHADERS_PATH],"`pwd`/src/View/shaders/")

    AC_MSG_NOTICE([View enabled])

	PKG_CHECK_MODULES(GL, [gl glu], ,[AC_MSG_WARN([GL not found: $GL_PKG_ERRORS])])
	PKG_CHECK_MODULES(GLEW, [glew], ,[AC_MSG_WARN([GLEW not found: $GLEW_PKG_ERRORS])])
	PKG_CHECK_MODULES(GLFW, [glfw3], ,[AC_MSG_WARN([GLFW not found: $GLFW_PKG_ERRORS])])
else
    AC_MSG_NOTICE([View disabled])
fi


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
