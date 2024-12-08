#include "Solver_Godunov1D.h"
#include "Solver_Lagrange1D.h"
#include "Parameters.h"
#include <cstring>

int main(int argv, char** argc)
{
	std::string params_file = "params.txt";
	Parameters par(params_file); // Parameters initialization by file;
	Solver_Lagrange1D solver(par);
	if (argv > 1)
		return 0;
	system("python Post.py");
#ifdef WIN32
#else
    system("sxiv graph.png");
#endif // WIN32
}
