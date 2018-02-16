# Commet lines
# Here we define compiler option, libraries and the target
FLAGS = -O2 -Wall

# Here we make the executable file
SRCS = main.cpp atom.cpp io.cpp lennardjones.cpp statisticssampler.cpp system.cpp unitconverter.cpp vec2.cpp velocityverlet.cpp 
OBJS = $(subst .cpp,.o,$(SRCS))
all: MPI-md

# Whereas here we create the object file
MPI-md: $(OBJS) 
	mpic++ ${FLAGS} -o MPI-md $(OBJS) # note: must use mpic++ compiler for MPI with cpp files (mpicc for .c)
#	export OMP_NUM_THREADS=20
	./md

%.o: %.cpp  
	mpic++ $(FLAGS) -c $<

#	g++ ${FLAGS} -c init.cpp
#	g++ ${FLAGS} -c integrator.cpp
#	g++ ${FLAGS} -c MD.cpp

# Clean
clean:
	rm *.o ./MPI-md
