#ifndef SOLVER_LAGRANGE1D_H
#define SOLVER_LAGRANGE1D_H

#include "Parameters.h"
#include "iSolver.h"
#include <memory>

class Solver_Lagrange1D: public iSolver {
	public:
		Solver_Lagrange1D(const Parameters& _par);
		void apply_boundary_conditions() override;
		void solve_step() override;
		void set_initial_conditions() override;
		void get_time_step() override;
		void write_data() override;
		void check_parameters();
	private:
		const Parameters par;
		std::shared_ptr<double[]> p, rho, U, m;
		std::shared_ptr<double[]> v, x;
		std::shared_ptr<double[]> omega; //viscosity
		double t, dt;
		int step;
};
#endif
