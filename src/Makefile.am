
modelSources = Model/atom.cpp Model/io.cpp Model/lennardjones.cpp Model/modelMain.cpp Model/send_atoms.cpp Model/statisticssampler.cpp Model/system.cpp Model/unitconverter.cpp Model/vec2.cpp Model/velocityverlet.cpp 

if WITH_VIEW
viewSources = View/circles.cpp View/loadShaders.cpp View/viewMain.cpp View/hsl_to_rgb.cpp Controller/controller.cpp Controller/controllerState.cpp
viewLDADD = $(GL_LIBS) $(GLEW_LIBS) $(GLFW_LIBS) 
viewCFLAGS = $(GL_CFLAGS) $(GLEW_CFLAGS) $(GLFW_CFLAGS) 
endif

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = $(top_builddir)/bin/parallel2DLJMD

__top_builddir__bin_parallel2DLJMD_SOURCES = main.cpp comms.cpp $(modelSources) $(viewSources)
__top_builddir__bin_parallel2DLJMD_LDADD = $(viewLDADD)
__top_builddir__bin_parallel2DLJMD_CFLAGS = $(viewCFLAGS)
