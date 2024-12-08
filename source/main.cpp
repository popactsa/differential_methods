#include "Solver_Godunov1D.h"
#include "Parameters.h"
#include <cstring>

int main(int argv, char** argc)
{
	std::string params_file = "Parameters_file.txt";
	Parameters par(params_file); // Parameters initialization by file;
	Solver_Godunov1D solver(par);
	if (argv > 1)
		return 0;
	system("python Post.py");
#ifdef WIN32
#else
    system("sxiv graph.png");
#endif // WIN32
}
