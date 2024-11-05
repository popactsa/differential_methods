all:
	g++ main.o Parameters.o -o main

main.o: main.cpp
	g++ -c main.cpp

Parameters.o: Parameters.cpp
	g++ -c Parameters.cpp
