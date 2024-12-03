#ifndef SOLVER_GODUNOV1D_H
#define SOLVER_GODUNOV1D_H

#include "Parameters.h"
#include "iSolver.h"

class Solver_Godunov1D: public iSolver
{
	const Parameters par;
	double *p, *rho, *rho_u, *rho_e;
	double *x;
	double t, dt;
	unsigned int step;
	double *F_m, *F_imp, *F_e;
public:
	Solver_Godunov1D(const Parameters& _par);
	void apply_boundary_conditions();
	void solve_step();
	void set_initial_conditions();
	void get_time_step();
	void write_data();
	~Solver_Godunov1D();
};
#endif
