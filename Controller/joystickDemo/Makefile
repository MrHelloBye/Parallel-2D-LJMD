
#Make a copy of this Makefile and change the following three lines as needed
#These should work for any standard linux system, and INCLUDES and LIBRARIES
#probably unnecessary
CC = gcc
INCLUDES = -I/usr/local/include -I/usr/include
LIBRARIES = -L/usr/local/lib -L/usr/lib 

#Make sure all these libraries (and linmath) are installed or otherwise in your
#include/linking path. dl, m, and pthread should come with most systems
LINKING =  -ldl -lm -lpthread -lglfw -lglad -ltinycthread 

all:
	$(CC) joystickDemo.c -o joystickDemo $(INCLUDES) $(LIBRARIES) $(LINKING)

.PHONY : clean
clean:
	rm joystickDemo
