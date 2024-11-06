all: main

main: main.o Parameters.o
	g++ main.o Parameters.o -o main

main.o: main.cpp
	g++ -c main.cpp main.o

Parameters.o: Parameters.cpp
	g++ -c Parameters.cpp main.o
