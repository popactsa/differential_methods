//#include "iSolver.h"
#include "Solver_Lagrange1D.h"
#include "Parameters.h"
//#include "auxiliary_functions.h"
//#include "Wall.h"

int main()
{
	std::string file_name = "Parameters_file.txt";
	Parameters par(file_name); // Parameters initialization by file;
	Solver_Lagrange1D solver(par);
}
