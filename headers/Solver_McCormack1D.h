#ifndef SOLVER_MCCORMACK1D_H
#define SOLVER_MCCORMACK1D_H

#include "Parameters.h"
#include "Solver_Godunov1D.h"
#include <cmath>

class Solver_McCormack1D: public Solver_Godunov1D
{
	public:
		Solver_McCormack1D(const Parameters& _par);
		void solve_step() override;
		~Solver_McCormack1D();
	private:
		constexpr static double a[2]{0.0, 1.0};
		constexpr static double b[2]{0.5, 0.5};
		double* U[3];
		double* H[3];
		double* h[3];
		double* Flux[3];
		double* u[3];
		double* p_temp;
};

#endif
