#ifndef SOLVER_LAGRANGE1D_H
#define SOLVER_LAGRANGE1D_H

#include "Parameters.h"
#include "iSolver.h"

class Solver_Lagrange1D: public iSolver {
	public:
		Solver_Lagrange1D(const Parameters& _par);
		void apply_boundary_conditions();
		void solve_step();
		void set_initial_conditions();
		void get_time_step();
		void write_data();
		void check_parameters();
		~Solver_Lagrange1D();
	private:
		const Parameters par;
		double *p, *rho, *U, *m;
		double *v, *x;
		double t, dt;
		int step; // current time step
		double *omega; // viscosity
};
#endif
