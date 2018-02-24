
modelSources = Model/atom.cpp Model/io.cpp Model/lennardjones.cpp Model/modelMain.cpp Model/send_atoms.cpp Model/statisticssampler.cpp Model/system.cpp Model/unitconverter.cpp Model/vec2.cpp Model/velocityverlet.cpp 

if WITH_VIEW
viewSources = View/circles.cpp View/loadShaders.cpp View/viewMain.cpp 
viewLDADD = -L/usr/local/lib -lXrandr -lXrender -lXi -lGLEW -lGLU -lm -lGL -lm -lpthread -ldl -ldrm -lXdamage -lXfixes -lX11-xcb -lxcb-glx -lxcb-dri2 -lXxf86vm -lXext -lglfw3 -lrt -lm -ldl -lX11 -lpthread -lxcb -lXau -lXdmcp -lglfw3
endif

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile
bin_PROGRAMS = $(top_builddir)/bin/parallel2DLJMD

__top_builddir__bin_parallel2DLJMD_SOURCES = main.cpp $(modelSources) $(viewSources)
__top_builddir__bin_parallel2DLJMD_LDADD = $(viewLDADD)
