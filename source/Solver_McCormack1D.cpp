#include "Solver_McCormack1D.h"

Solver_McCormack1D::Solver_McCormack1D(const Parameters& _par):
    Solver_Godunov1D(_par, false)
{
	U[0] = rho;
	U[1] = rho_u;
	U[2] = rho_e;
	for (int i = 0; i < 3; ++i)
	{
		H[i] = new double[par.nx_all] {0.0};
		h[i] = new double[par.nx_all] {0.0};
		Flux[i] = new double[par.nx_all] {0.0};
		u[i] = new double[par.nx_all] {0.0};
	}
	p_temp = new double[par.nx_all] {0.0};
	
	for (step = 1; step <= par.nt; ++step)
	{
		solve_step();
		if (step % par.nt_write == 0) write_data();
	}
	write_exact_solution(*this);
	std::cout << "\nDone McCormack" << std::endl;
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
			p_temp[i] = (par.gamma - 1) * (u[2][i] - std::pow(u[1][i], 2) / u[0][i] / 2.0);
			Flux[0][i] = u[1][i];
			Flux[1][i] = std::pow(u[1][i], 2) / u[0][i] + p_temp[i];
			Flux[2][i] = u[1][i] * u[2][i] / u[0][i] + p_temp[i] * u[1][i] / u[0][i];
		}
		if (s == 0)
		{
			for (int k = 0; k < 3; ++k)
			{
				for (int i = 1; i < par.nx_all - 1; ++i)
					h[k][i] = -dt * (Flux[k][i + 1] - Flux[k][i]) / par.dx;
				h[k][0] = h[k][1];
				h[k][par.nx_all - 1] = h[k][par.nx_all - 2];
			}
		}
		else
		{
			for (int k = 0; k < 3; ++k)
			{
				for (int i = 1; i < par.nx_all - 1; ++i)
					h[k][i] = -dt * (Flux[k][i] - Flux[k][i - 1]) / par.dx;
				h[k][0] = h[k][1];
				h[k][par.nx_all - 1] = h[k][par.nx_all - 2];
			}
		}
		for (int k = 0; k < 3; ++k)
			for (int i = 1; i < par.nx_all - 1; ++i)
		      	H[k][i] += b[s] * h[k][i];
	}
	for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
	{
		rho[i] += H[0][i];
		rho_u[i] += H[1][i];
		rho_e[i] += H[2][i];
		p[i] = (par.gamma - 1) * (rho_e[i] - std::pow(rho_u[i], 2) / rho[i] / 2.0);
		H[0][i] = 0.0;
		H[1][i] = 0.0;
		H[2][i] = 0.0;
	}
	//монотонность
	double QD = 0.1;
	double QN = 0.1;
	double Du[par.nx_all];
	double Nu[par.nx_all];
	double NDu[par.nx_all];
	for(int k = 0; k < 3; k++)
	{
		diffusion(QD, U[k], Du);
		anti_diffusion(QN, U[k], Nu);
		anti_diffusion(QN, Du, NDu);
		for(int i = 0; i < par.nx_all; i++)
			U[k][i] += Du[i] + Nu[i] + NDu[i];
	}
	//давление для красоты
	for (int i = par.nx_fict; i < par.nx+par.nx_fict; ++i)
		p[i] = (par.gamma - 1) * (rho_e[i] - std::pow(rho_u[i], 2) / rho[i] / 2.0);
}

void Solver_McCormack1D::diffusion(double Q, double* u, double* Du)
{
	double phi[par.nx_all-1];
	for(int i = 0; i < par.nx_all - 1; i++)
		phi[i] = Q * (u[i + 1] - u[i]);
	for(int i = 1; i < par.nx_all - 1; i++)
		Du[i] = phi[i] - phi[i - 1];
	Du[0] = Du[1];
	Du[par.nx_all - 1] = Du[par.nx_all - 2];
}

void Solver_McCormack1D::anti_diffusion(double Q, double* u, double* Nu)
{
	double phi[par.nx_all - 1];
	double u_b[par.nx_all + 1];
	for(int i = 0; i < par.nx_all - 1; i++)
	{
		phi[i] = Q * (u[i + 1] - u[i]);
		u_b[i + 1] = 0.5 * (u[i + 1] + u[i]);
	}
	u_b[0] = u_b[1];
	u_b[par.nx_all] = u_b[par.nx_all - 1];
	double dud[par.nx_all + 1];
	double Du_b[par.nx_all + 1];
	diffusion(Q, u_b, Du_b);
	for(int i = 0; i < par.nx_all + 1; i++)
		dud[i] = u_b[i] + Du_b[i];

	double phi_c[par.nx_all - 1];
	for(int i = 0; i < par.nx_all - 1; i++)
	{
		double s = phi[i] > 0 ? 1 : -1;
		phi_c[i] = s 
			* std::max(
				0.0, 
				std::min(
					s * dud[i], 
					std::min(
						fabs(phi[i]),
						s * dud[i + 2]
					)
				)
			);
	}
	for(int i = 1; i < par.nx_all - 1; i++)
		Nu[i] = -(phi_c[i] - phi_c[i - 1]);
	Nu[0] = Nu[1];
	Nu[par.nx_all - 1] = Nu[par.nx_all - 2];
}
