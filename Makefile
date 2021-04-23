LIBS=-larmadillo
CC=g++
COMPILER_OPTIONS=-std=c++11 -g -fno-inline -Wall

OBJECTS=Point_Vortex.o Post_Process.o $(LIBS) 


main: $(OBJECTS)
	$(CC) $(COMPILER_OPTIONS) -o $@ $^


$(OBJECTS): Point_Vortex.hpp Post_Process.hpp 



clean: 
	$(RM) main *.o Vortex*.txt python.txt 
