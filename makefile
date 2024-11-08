all: main

main: main.o Parameters.o auxiliary_functions.o
	g++ -o main main.o Parameters.o auxiliary_functions.o

main.o: main.cpp
	g++ -c -o main.o main.cpp

Parameters.o: Parameters.cpp
	g++ -c -o Parameters.o Parameters.cpp

auxiliary_functions.o: auxiliary_functions.cpp
	g++ -c -o auxiliary_functions.o auxiliary_functions.cpp
