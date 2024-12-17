#include "Solver_McCormack1D.h"
#include <algorithm>

Solver_McCormack1D::Solver_McCormack1D(const Parameters& _par):
	Solver_Godunov1D(_par, false)
{
	U[0] = rho;
	U[1] = rho_u;
	U[2] = rho_e;
	for (int i = 0; i < 3; ++i)
	{
		H[i] = new double[par.nx_all]{};
		h[i] = new double[par.nx_all]{};
		Flux[i] = new double[par.nx_all]{};
		u[i] = new double[par.nx_all]{};
	}
	p_temp = new double[par.nx_all]{};

	set_initial_conditions();
	t = 0.0;
	for (step = 1; step <= par.nt; ++step)
	{
		solve_step();
		if (step % par.nt_write == 0) write_data();
	}
	write_exact_solution(*this);
	std::cout << "\nDone" << std::endl;
}

Solver_McCormack1D::~Solver_McCormack1D()
{
	for(int k = 0; k < 3; ++k)
	{
		delete[] H[k];
		delete[] h[k];
		delete[] Flux[k];
		delete[] u[k];
	}
	delete[] p_temp;
}

void Solver_McCormack1D::solve_step()
{
	apply_boundary_conditions();
	get_time_step();
	t += dt;	
	for (int s = 0; s < 2; ++s)
	{
		for (int k = 0; k < 3; ++k)
			for (int j = 0; j < par.nx_all; ++j)
				u[k][j] = U[k][j] + a[s] * h[k][j];
		for (int i = 0; i < par.nx_all; ++i)
		{
			p_temp[i] = (par.gamma - 1) * (u[2][i] - std::pow(u[1][i], 2)/u[0][i]/2.0);
			Flux[0][i] = u[1][i];
			Flux[1][i] = std::pow(u[1][i], 2)/u[0][i] + p_temp[i];
			Flux[2][i] = u[1][i] * u[2][i] / u[0][i] + p_temp[i] * u[1][i] / u[0][i];
		}
		for (int k = 0; k < 3; ++k)
		{
				for (int i = 1; i < par.nx_all - 1; ++i)
				{
					h[k][i] = (s == 0)
						? -dt * (Flux[k][i + 1] - Flux[k][i])/par.dx
						: -dt * (Flux[k][i] - Flux[k][i - 1])/par.dx;
					H[k][i] += b[s] * h[k][i];
				}
			for (int i = 0; i < 2; ++i) // 1D case
			{
				int sign = i == 0 ? 1 : -1;
				switch (par.boundaries[i].b_type)
				{
					case B_FLUX:
						h[k][i * (par.nx_all - 1)] = h[k][i * (par.nx_all - 1) + sign];
					case B_WALL:
						h[k][i * (par.nx_all - 1)] = -h[k][i * (par.nx_all - 1) + sign];
				}
			}
		}
	}
	for (int i = par.nx_fict; i < par.nx+par.nx_fict; ++i)
	{
		rho[i] += H[0][i];
		rho_u[i] += H[1][i];
		rho_e[i] += H[2][i];
		p[i] = (par.gamma - 1) * (rho_e[i] - std::pow(rho_u[i], 2)/rho[i]/2.0);
		H[0][i] = 0.0;
		H[1][i] = 0.0;
		H[2][i] = 0.0;
	}
	// artificial diffusion
	double s_[3][par.nx];
	double um[3][par.nx];
	double u_[3][par.nx];
	double up[3][par.nx];
	double sud[2][3][par.nx];
	double phi[2][3][par.nx];
	double Du[3][par.nx];
	for (int s = 0; s < 2; ++s)
	{
		int is_odd = s%2 ? 1 : 0;
		for (int i = 1; i < par.nx_all - par.nx_fict; ++i)
		{
			for (int k = 0; k < 3; ++k)
			{
				s_[k][i] = U[k][i + is_odd] - U[k][i - 1 + is_odd] > 0 ? 1.0 : -1.0;
				um[k][i] = 0.5 * (U[k][i - 2 + is_odd] - U[k][i - 3 + is_odd]);
				u_[k][i] = 0.5 * (U[k][i - 1 + is_odd] - U[k][i - 2 + is_odd]);
				up[k][i] = 0.5 * (U[k][i + is_odd] - U[k][i - 1 + is_odd]);
				get_diffusion(Du[k], um[k], u_[k], up[k], par.nx);
				sud[0][k][i] = u_[k][i] + Du[k][i];

				um[k][i] = 0.5 * (U[k][i + is_odd] - U[k][i - 1 + is_odd]);
				u_[k][i] = 0.5 * (U[k][i + 1 + is_odd] - U[k][i + is_odd]);
				up[k][i] = 0.5 * (U[k][i + 2 + is_odd] - U[k][i + 1 + is_odd]);
				get_diffusion(Du[k], um[k], u_[k], up[k], par.nx);
				sud[1][k][i] = u_[k][i] + Du[k][i];

				phi[is_odd][k][i] = std::max(
							0.0,
							std::min(
								s_[k][i] * sud[0][k][i],
								std::fabs(par.Q * (u[k][i + is_odd] - u[k][i - 1 + is_odd]))
								)
						);
				phi[is_odd][k][i] = s_[k][i] * std::max(
							phi[is_odd][k][i],
							s_[k][i] * sud[1][k][i]
						);
			}
		}
	}
	double Nu[3][par.nx];
	for (int k = 0; k < 3; ++k)
		for (int i = 0; i < par.nx; ++i)
		{
			Nu[k][i] = phi[0][k][i] - phi[1][k][i];
			U[k][i] = U[k][i] * (1 + Du[k][i])*(1 + Nu[k][i]);
		}
}

void Solver_McCormack1D::get_diffusion(double *Du, double *um, double *u_, double *up, int size)
{
	for (int j = 0; j < size; ++j)
	{
		Du[j] = par.Q * (up[j - 1] - up[j]
				+ up[j + 1] - up[j]
		);
	}
}

