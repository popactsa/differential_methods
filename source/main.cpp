//#include "iSolver.h"
#include "Solver_Godunov1D.h"
#include "Parameters.h"
#include <cstring>
//#include "auxiliary_functions.h"
//#include "Wall.h"

int main(int argv, char** argc)
{
	std::string file_name = "Parameters_file.txt";
	Parameters par(file_name); // Parameters initialization by file;
	Solver_Godunov1D solver(par);
	if (argv > 1) return 0;
	system("python Post.py");
#ifdef WIN32
#else
    system("sxiv graph.png");
#endif // WIN32
}
