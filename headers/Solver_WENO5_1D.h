#ifndef Solver_WENO5_1D_H
#define Solver_WENO5_1D_H

#include "Parameters.h"
#include "Solver_Godunov1D.h"
#include <algorithm>
#include <cmath>

class Solver_WENO5_1D: public Solver_Godunov1D
{
public:
	Solver_WENO5_1D(const Parameters& _par);
	void solve_step();
	void get_WENO_reconstruction(double* flux, int n_borders, double* out_arr, int size1);
};

#endif
