CXX=g++

FLAGS=-O3

all: main

bfgs.o : bfgs.h bfgs.cpp
	$(CXX) -c $(FLAGS) bfgs.cpp -Wno-write-strings


analysis.o: analysis.h analysis.cpp
	$(CXX) -c $(FLAGS) analysis.cpp -Wno-write-strings


main: analysis.o NGSrelate.cpp analysis.o bfgs.o
	$(CXX)  $(FLAGS) NGSrelate.cpp -lz analysis.o -o NGSrelate bfgs.o

clean: 
	rm -f *.o NGSrelate