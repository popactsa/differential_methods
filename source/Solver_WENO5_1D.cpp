#include "Solver_WENO5_1D.h"

Solver_WENO5_1D::Solver_WENO5_1D(const Parameters& _par):
     Solver_Godunov1D(_par, false)
{
	for (step = 1; step <= par.nt; ++step)
	{
    		solve_step();
		if (step % par.nt_write == 0) write_data();
	}
  	write_exact_solution(*this);
	std::cout << "\nDone WENO5" << std::endl;
}

void Solver_WENO5_1D::get_WENO_reconstruction(double* flux, int n_borders, double* out_arr, int size1)
{
	double h0[n_borders] = {};
	double h1[n_borders] = {};
	double h2[n_borders] = {};
	double h3[n_borders] = {};

	double hl[n_borders] = {};
	double hr[n_borders] = {};

	double epsilon = 1E-40;

	double alpha_tmp[3] = {};
	double omega_tmp[3] = {};
	double gamma_tmp[3] = {3.0/10, 3.0/5, 1.0/10}; // стр. 112
	double beta_tmp[3] = {};

	for (int i = 1; i < n_borders - 1; ++i)
	{
		// Стр. 112
		int i_tmp = i + par.nx_fict - 1;

	  	h0[i] = 1.0/3 * flux[i_tmp] + 5.0/6 * flux[i_tmp + 1] - 1.0/6 * flux[i_tmp + 2];
	  	h1[i] = -1.0/6 * flux[i_tmp - 1] + 5.0/6 * flux[i_tmp] + 1.0/3 * flux[i_tmp + 1];
	  	h2[i] = 1.0/3 * flux[i_tmp - 2] - 7.0/6 * flux[i_tmp - 1] + 11.0/6 * flux[i_tmp];
	  	h3[i] = 11.0/6 * flux[i_tmp + 1] - 7.0/6 * flux[i_tmp + 2] + 1.0/3 * flux[i_tmp + 3];

		// Стр. 113
		beta_tmp[0] = 13.0/12 * pow(flux[i_tmp] - 2 * flux[i_tmp + 1] + flux[i_tmp + 2], 2) + \
			1.0/4 * pow(3 * flux[i_tmp] - 4 * flux[i_tmp + 1] + flux[i_tmp + 2], 2);
		beta_tmp[1] = 13.0/12 * pow(flux[i_tmp - 1] - 2 * flux[i_tmp] + flux[i_tmp + 1], 2) + \
			1.0/4 * pow(flux[i_tmp - 1] - flux[i_tmp + 1], 2);
		beta_tmp[2] = 13.0/12 * pow(flux[i_tmp - 2] - 2 * flux[i_tmp - 1] + flux[i_tmp], 2) + \
			1.0/4 * pow(flux[i_tmp - 2] - 4 * flux[i_tmp - 1] + 3 * flux[i_tmp], 2);

		// hl == u_left
		double sum_alpha = 0.0;
		for (int k = 0; k < 3; ++k)
		{
			alpha_tmp[k] = pow(epsilon + beta_tmp[k], 2);
			alpha_tmp[k] = gamma_tmp[k] / alpha_tmp[k];
			sum_alpha += alpha_tmp[k];
		}
		for (int k = 0; k < 3; ++k)
		{
			omega_tmp[k] = alpha_tmp[k] / sum_alpha;
		}
		hl[i] = omega_tmp[0] * h0[i] + omega_tmp[1] * h1[i] + omega_tmp[2] * h2[i];

		// hr == u_right
		i_tmp = i + par.nx_fict;

		beta_tmp[0] = 13.0/12 * pow(flux[i_tmp] - 2 * flux[i_tmp + 1] + flux[i_tmp + 2], 2) + \
			1.0/4 * pow(3 * flux[i_tmp] - 4 * flux[i_tmp + 1] + flux[i_tmp + 2], 2);
		beta_tmp[1] = 13.0/12 * pow(flux[i_tmp - 1] - 2 * flux[i_tmp] + flux[i_tmp + 1], 2) + \
			1.0/4 * pow(flux[i_tmp - 1] - flux[i_tmp + 1], 2);
		beta_tmp[2] = 13.0/12 * pow(flux[i_tmp - 2] - 2 * flux[i_tmp - 1] + flux[i_tmp], 2) + \
			1.0/4 * pow(flux[i_tmp - 2] - 4 * flux[i_tmp - 1] + 3 * flux[i_tmp], 2);

		sum_alpha = 0.0;
		for (int k = 0; k < 3; ++k)
		{
			alpha_tmp[k] = pow(epsilon + beta_tmp[k], 2);
			alpha_tmp[k] = gamma_tmp[2 - k] / alpha_tmp[k];
			sum_alpha += alpha_tmp[k];
		}
		for (int k = 0; k < 3; ++k)
		{
			omega_tmp[k] = alpha_tmp[k] / sum_alpha;
		}
		hr[i] = omega_tmp[0] * h3[i] + omega_tmp[1] * h0[i] + omega_tmp[2] * h1[i];
	}

	for (int i = 0; i < n_borders; ++i)
	{
		//std::cout << hr[i] << " ";
		out_arr[0*n_borders + i] = hr[i];
		out_arr[1*n_borders + i] = hl[i];
	}
}

void Solver_WENO5_1D::solve_step()
{
	apply_boundary_conditions();
	get_time_step();
	t += dt;

	int num_stencils = 1 + 3;
	double rho_stencil[num_stencils][par.nx_all] = {};
	double rho_u_stencil[num_stencils][par.nx_all] = {};
	double rho_e_stencil[num_stencils][par.nx_all] = {};

	for (int i = 0; i < par.nx_all; ++i)
	{
		rho_stencil[0][i] = rho[i];
		rho_u_stencil[0][i] = rho_u[i];
		rho_e_stencil[0][i] = rho_e[i];
	}

	for (int stencil = 1; stencil < 4; ++stencil)
	{
		// Приращение
		// Потоки в ячейки
		double flux_m[par.nx_all] = {};
		double flux_imp[par.nx_all] = {};
		double flux_e[par.nx_all] = {};

		double u_tmp[par.nx_all] = {};
		double rho_tmp[par.nx_all] = {};
		double p_tmp[par.nx_all] = {};
		// Скорости звука
	  	double A[par.nx_all] = {};

		for (int i = 0; i < par.nx_all; ++i)
		{
			u_tmp[i] = rho_u_stencil[stencil - 1][i] / rho_stencil[stencil - 1][i];
			rho_tmp[i] = rho_stencil[stencil - 1][i];
			p_tmp[i] = (par.gamma - 1) * (rho_e_stencil[stencil - 1][i] - rho_tmp[i] * pow(u_tmp[i], 2) / 2.0);

			flux_m[i] = rho_tmp[i] * u_tmp[i];
			flux_imp[i] = rho_tmp[i] * pow(u_tmp[i], 2) + p_tmp[i];
			flux_e[i] = (rho_e_stencil[stencil - 1][i] + p_tmp[i]) * u_tmp[i];

	  		A[i] = sqrt(par.gamma * p_tmp[i] / rho_tmp[i]); // Какие max тут брать не понятно
		}

		// Расщепление потков по Лаксу-Фридрихсу
		double flux_m_plus[par.nx_all] = {};
		double flux_imp_plus[par.nx_all] = {};
		double flux_e_plus[par.nx_all] = {};

		double flux_m_minus[par.nx_all] = {};
		double flux_imp_minus[par.nx_all] = {};
		double flux_e_minus[par.nx_all] = {};

		for (int i = 0; i < par.nx_all; ++i)
		{
			flux_m_plus[i] = 0.5 * (flux_m[i] + A[i] * rho_stencil[stencil - 1][i]);
			flux_imp_plus[i] = 0.5 * (flux_imp[i] + A[i] * rho_u_stencil[stencil - 1][i]);
			flux_e_plus[i] = 0.5 * (flux_e[i] + A[i] * rho_e_stencil[stencil - 1][i]);

			flux_m_minus[i] = 0.5 * (flux_m[i] - A[i] * rho_stencil[stencil - 1][i]);
			flux_imp_minus[i] = 0.5 * (flux_imp[i] - A[i] * rho_u_stencil[stencil - 1][i]);
			flux_e_minus[i] = 0.5 * (flux_e[i] - A[i] * rho_e_stencil[stencil - 1][i]);
		}

		// Реконструкция WENO5 для потков
		double flowL_m[par.nx + 1] = {};
		double flowL_imp[par.nx + 1] = {};
		double flowL_e[par.nx + 1] = {};

		double flowR_m[par.nx + 1] = {};
		double flowR_imp[par.nx + 1] = {};
		double flowR_e[par.nx + 1] = {};

		double h_arr[2][par.nx + 1] = {};

		get_WENO_reconstruction(flux_m_plus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowL_m[i] = h_arr[1][i];
			//std::cout << flowL_m[i] << " ";
		}
		get_WENO_reconstruction(flux_imp_plus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowL_imp[i] = h_arr[1][i];
			//std::cout << flowL_m[i] << " ";
		}

		get_WENO_reconstruction(flux_e_plus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowL_e[i] = h_arr[1][i];
			//std::cout << flowL_m[i] << " ";
		}

		get_WENO_reconstruction(flux_m_minus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowR_m[i] = h_arr[0][i];
			//std::cout << flowL_m[i] << " ";
		}

		get_WENO_reconstruction(flux_imp_minus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowR_imp[i] = h_arr[0][i];
			//std::cout << flowL_m[i] << " ";
		}

		get_WENO_reconstruction(flux_e_minus, par.nx + 1, &h_arr[0][0], 2);

		for (int i = 0; i < par.nx + 1; ++i)
		{
			// Для flowL по flux_plus
			flowR_e[i] = h_arr[0][i];
			//std::cout << flowL_m[i] << " ";
		}

		//std::cout << std::endl;

		// Граничные условия для потока
		// Потоки для свободной границы
		flowL_m[0] = 0.5 * (flux_m[0] + A[0] * rho_stencil[stencil - 1][0]);
		flowL_imp[0] = 0.5 * (flux_imp[0] + A[0] * rho_u_stencil[stencil - 1][0]);
		flowL_e[0] = 0.5 * (flux_e[0] + A[0] * rho_e_stencil[stencil - 1][0]);

		flowR_m[0] = 0.5 * (flux_m[0] - A[0] * rho_stencil[stencil - 1][0]);
		flowR_imp[0] = 0.5 * (flux_imp[0] - A[0] * rho_u_stencil[stencil - 1][0]);
		flowR_e[0] = 0.5 * (flux_e[0] - A[0] * rho_e_stencil[stencil - 1][0]);


		flowL_m[par.nx] = 0.5 * (flux_m[par.nx_all - 1] + A[par.nx_all - 1] * rho_stencil[stencil - 1][par.nx_all - 1]);
		flowL_imp[par.nx] = 0.5 * (flux_imp[par.nx_all - 1] + A[par.nx_all - 1] * rho_u_stencil[stencil - 1][par.nx_all - 1]);
		flowL_e[par.nx] = 0.5 * (flux_e[par.nx_all - 1] + A[par.nx_all - 1] * rho_e_stencil[stencil - 1][par.nx_all - 1]);

		flowR_m[par.nx] = 0.5 * (flux_m[par.nx_all - 1] - A[par.nx_all - 1] * rho_stencil[stencil - 1][par.nx_all - 1]);
		flowR_imp[par.nx] = 0.5 * (flux_imp[par.nx_all - 1] - A[par.nx_all - 1] * rho_u_stencil[stencil - 1][par.nx_all - 1]);
		flowR_e[par.nx] = 0.5 * (flux_e[par.nx_all - 1] - A[par.nx_all - 1] * rho_e_stencil[stencil - 1][par.nx_all - 1]);

		// Приращение L[s] = k[s] из Википедии
		double L_m[par.nx] = {};
		double L_imp[par.nx] = {};
		double L_e[par.nx] = {};

		for (int i = 0; i < par.nx; ++i)
		{
			L_m[i] = -((flowR_m[i + 1] + flowL_m[i + 1]) - (flowR_m[i] + flowL_m[i])) / par.dx;
			L_imp[i] = -((flowR_imp[i + 1] + flowL_imp[i + 1]) - (flowR_imp[i] + flowL_imp[i])) / par.dx;
			L_e[i] = -((flowR_e[i + 1] + flowL_e[i + 1]) - (flowR_e[i] + flowL_e[i])) / par.dx;
			//std::cout << L_m[i] << " ";
		}
		//std::cout << std::endl;

		// Решение на текущем шаге
		// Рунге-Кутта 3
		double rk[num_stencils][4] = {
			{1.0, 0.0, 0.0, 0.0}, \
			{1.0, 0.0, 0.0, 1.0}, \
			{0.75, 0.25, 0.0, 0.25}, \
			{1.0/3.0, 0.0, 2.0/3.0, 2.0/3.0}
		};


		for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
		{
			rho_stencil[stencil][i] += rk[stencil][0] * rho_stencil[0][i];
			rho_stencil[stencil][i] += rk[stencil][1] * rho_stencil[1][i];
			rho_stencil[stencil][i] += rk[stencil][2] * rho_stencil[2][i];

			rho_u_stencil[stencil][i] += rk[stencil][0] * rho_u_stencil[0][i];
			rho_u_stencil[stencil][i] += rk[stencil][1] * rho_u_stencil[1][i];
			rho_u_stencil[stencil][i] += rk[stencil][2] * rho_u_stencil[2][i];

			rho_e_stencil[stencil][i] += rk[stencil][0] * rho_e_stencil[0][i];
			rho_e_stencil[stencil][i] += rk[stencil][1] * rho_e_stencil[1][i];
			rho_e_stencil[stencil][i] += rk[stencil][2] * rho_e_stencil[2][i];

			rho_stencil[stencil][i] += rk[stencil][3] * L_m[i - par.nx_fict] * dt;
			rho_u_stencil[stencil][i] += rk[stencil][3] * L_imp[i - par.nx_fict] * dt;
			rho_e_stencil[stencil][i] += rk[stencil][3] * L_e[i - par.nx_fict] * dt;

	   		//std::cout << rho_stencil[1][i] << " ";
		}
	  	//std::cout << std::endl;
		// Граничные условия
		for (int i = 0; i < par.nx_fict; ++i)
		{
			rho_stencil[stencil][i] = rho_stencil[stencil][par.nx_fict];
			rho_stencil[stencil][par.nx_all - i - 1] = rho_stencil[stencil][par.nx_all - par.nx_fict - 1];

			rho_u_stencil[stencil][i] = rho_u_stencil[stencil][par.nx_fict];
			rho_u_stencil[stencil][par.nx_all - i - 1] = rho_u_stencil[stencil][par.nx_all - par.nx_fict - 1];

			rho_e_stencil[stencil][i] = rho_e_stencil[stencil][par.nx_fict];
			rho_e_stencil[stencil][par.nx_all - i - 1] = rho_e_stencil[stencil][par.nx_all - par.nx_fict - 1];
		}
	}

	// Новые значения в массивах c stencil=3
	for (int i = 0; i < par.nx_all; ++i)
	{
		rho[i] = rho_stencil[3][i];
		rho_u[i] = rho_u_stencil[3][i];
		rho_e[i] = rho_e_stencil[3][i];
		// Обновим p для красоты
		p[i] = (par.gamma - 1) * (rho_e[i] - pow(rho_u[i], 2) / 2 / rho[i]);
	}
};

