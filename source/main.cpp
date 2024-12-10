#include "Solver_Godunov1D.h"
#include "Solver_Lagrange1D.h"
#include "Parameters.h"
#include <cstring>

int main(int argv, char** argc)
{
	std::string params_file = "params.txt";
	Parameters par(params_file); // Parameters initialization by file;
	int choose = 0;
	while (!choose)
	{
		std::cout << std::endl << "\tChoose solver:" << std::endl;
		std::cout << "\t\t1.Lagrange1D" << std::endl;
		std::cout << "\t\t2.Godunov1D" << std::endl;
		std::cin >> choose;
		if (choose == 1) Solver_Lagrange1D solver(par);
		else if (choose == 2) Solver_Godunov1D solver(par);
		else
		{
				std::cout << "Choose another value";
				choose = 0;
		}
	}
	system("python Post.py");
#ifdef WIN32
#else
    system("sxiv graph.png");
#endif // WIN32
}
