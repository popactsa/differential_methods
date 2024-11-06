all: main

main: Manager.o Solver.o Writer.o
	g++ Manager.o Solver.o Writer.o -o main main.cpp

Manager.o: Manager.cpp
	g++ -c Manager.cpp

Solver.o: Solver.cpp
	g++ -c Solver.cpp

Writer.o: Writer.cpp
	g++ -c Writer.cpp