#include "Solver_Godunov1D.h"
#include <cmath>

Solver_Godunov1D::Solver_Godunov1D(const Parameters& _par): par(_par)
{
    p = new double[par.nx_all];
    rho = new double[par.nx_all];
    rho_u = new double[par.nx_all];
    rho_e = new double[par.nx_all];

    x = new double[par.nx_all + 1];
    F_m = new double[par.nx + 1];
    F_imp = new double[par.nx + 1];
    F_e = new double[par.nx + 1];
    // Начальные условия
    set_initial_conditions();
    t = 0.0;
    for (step = 1; step <= par.nt; step++)
    {
        // Решение на текущем шаге
        solve_step();
        if (step % par.nt_write == 0)
            write_data();
        if(step == 49)
            exact_solution(*this);
    }
    // Конец расчета
    std::cout << "Done" << std::endl;
}

Solver_Godunov1D::~Solver_Godunov1D()
{
    delete[] p;
    delete[] rho;
    delete[] rho_u;
    delete[] rho_e;
    delete[] x;
    delete[] F_m;
    delete[] F_e;
    delete[] F_imp;
}

// TODO переписать под выбор типа граничного условия
void Solver_Godunov1D::apply_boundary_conditions()
{
    for (unsigned int i = 0; i < par.nx_fict; ++i)
    {
        rho[i] = rho[par.nx_fict];
        rho[par.nx_all - i - 1] = rho[par.nx_all - par.nx_fict - 1];

        rho_e[i] = rho_e[par.nx_fict];
        rho_e[par.nx_all - i - 1] = rho_e[par.nx_all - par.nx_fict - 1];

        p[i] = p[par.nx_fict];
        p[par.nx_all - i - 1] = p[par.nx_all - par.nx_fict - 1];
    }
    for (unsigned int i = 0; i < 2; ++i)   // 1D -> 2 walls
    {
        if (par.walls[i].type == WALL)
        {
            double v_b = 0.0;
            if(par.ic_preset == IC_TEST2)
            {
                v_b = -2.0 ? (i == 0) : 2.0;
            }
            int i_temp = 0;
            if(i)
                i_temp = par.nx_all-1;
            int sign = 1-2*(i%2);
            for (unsigned int k = 0; k < par.nx_fict; ++k)
                rho_u[i_temp+sign*k] = -rho_u[i_temp+sign*par.nx_fict];
            int f_i_temp = 0;
            if(i%2)
                f_i_temp = par.nx;
            F_m[f_i_temp] = v_b*rho[i_temp];
            F_imp[f_i_temp] = pow(v_b, 2)*rho[i_temp] + p[i_temp];
            F_e[f_i_temp] = (rho_e[i_temp] + p[i_temp]) * v_b * rho[i_temp];
        }
        else if (par.walls[i].type == FLUX)
        {
            int i_temp = 0;
            if(i%2)
                i_temp = par.nx_all-1;
            int sign = 1-2*(i%2);
            for (unsigned int k = 0; k < par.nx_fict; ++k)
                rho_u[i_temp+sign*k] = rho_u[i_temp+sign*par.nx_fict];
            int f_i_temp = 0;
            if(i%2)
                f_i_temp = par.nx;
            F_m[f_i_temp] = rho_u[i_temp];
            F_imp[f_i_temp] = pow(rho_u[i_temp], 2) / rho[i_temp] + p[i_temp];
            F_e[f_i_temp] = (rho_e[i_temp] + p[i_temp]) * rho_u[i_temp] / rho[i_temp];
        }
    }
    /*F_m[0] = rho_u[0];
    F_m[par.nx] = rho_u[par.nx_all - 1];

    F_imp[0] = pow(rho_u[0], 2) / rho[0] + p[0];
    F_imp[par.nx] = pow(rho_u[par.nx_all - 1], 2) / rho[par.nx_all - 1] + p[par.nx_all - 1];

    F_e[0] = (rho_e[0] + p[0]) * rho_u[0] / rho[0];
    F_e[par.nx] = (rho_e[par.nx_all - 1] + p[par.nx_all - 1]) * rho_u[par.nx_all - 1] / rho[par.nx_all - 1];*/
    // density and internal energy
}

void Solver_Godunov1D::solve_step()
{
    get_time_step();
    t += dt;
    apply_boundary_conditions();
    calc_all_flux(par, \
                  F_m, F_imp, F_e, \
                  rho, rho_u, rho_e, p);
    /*// Рунге-Кутта
    double rho_pred[par.nx_all];
    double rho_u_pred[par.nx_all];
    double rho_e_pred[par.nx_all];
    double p_pred[par.nx_all];
    // Пересчет значений предиктор
    for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
    {
        rho_pred[i] = rho[i] - (F_m[i] - F_m[i - 1]) * dt / par.dx;
        rho_u_pred[i] = rho_u[i] - (F_imp[i] - F_imp[i - 1]) * dt / par.dx;
        rho_e_pred[i] = rho_e[i] - (F_e[i] - F_e[i - 1]) * dt / par.dx;
        p_pred[i] = (par.gamma - 1) * (rho_e_pred[i] - pow(rho_u_pred[i], 2) / rho_pred[i] / 2.0);
    }
    calc_all_flux(par, \
                  F_m, F_imp, F_e, \
                  rho_pred, rho_u_pred, rho_e_pred, p_pred);
    // Пересчет значений корректор
    for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
    {
        rho[i] = 0.5*(rho[i] + rho_pred[i]) - 0.5*(F_m[i] - F_m[i - 1]) * dt / par.dx;
        rho_u[i] = 0.5*(rho_u[i] + rho_u_pred[i]) - 0.5*(F_imp[i] - F_imp[i - 1]) * dt / par.dx;
        rho_e[i] = 0.5*(rho_e[i] + rho_e_pred[i]) - 0.5*(F_e[i] - F_e[i - 1]) * dt / par.dx;
        p[i] = (par.gamma - 1) * (rho_e[i] - pow(rho_u[i], 2) / rho[i] / 2.0);
    }*/
    // формула Эйлера
    // Пересчет значений
    for (int i = par.nx_fict; i < par.nx + par.nx_fict; ++i)
    {
        rho[i] = rho[i] - (F_m[i] - F_m[i - 1]) * dt / par.dx;
        rho_u[i] = rho_u[i] - (F_imp[i] - F_imp[i - 1]) * dt / par.dx;
        rho_e[i] = rho_e[i] - (F_e[i] - F_e[i - 1]) * dt / par.dx;
        p[i] = (par.gamma - 1) * (rho_e[i] - pow(rho_u[i], 2) / rho[i] / 2.0);
    }
};

// TODO переписать под передачу массивов в качестве аргументов
void Solver_Godunov1D::set_initial_conditions()
{
    if (par.ic_preset == IC_TEST1)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx < 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                rho[i] = 1.0;
                p[i] = 1.0;
            }
            else
            {
                rho[i] = 0.125;
                p[i] = 0.1;
            }
            rho_e[i] = rho[i] * ((pow(0.0, 2)) / 2.0  + p[i] / (par.gamma - 1.0) / rho[i]);
            rho_u[i] = rho[i] * 0.0;
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        }
    }

    else if (par.ic_preset == IC_TEST2)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            rho[i] = 1.0;
            rho_u[i] = i * par.dx < 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx ? rho[i] * (-2.0) : rho[i] * 2.0;
            p[i] = 0.4;
            rho_e[i] = rho[i] * pow(2.0, 2) / 2.0 + p[i] / (par.gamma - 1.0);
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        }
    }

    else if (par.ic_preset == IC_TEST3)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            p[i] = i * par.dx <= 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx ? 1000.0 : 0.01;
            rho[i] = 1.0;
            rho_e[i] = rho[i] * ((pow(0.0, 2)) / 2.0  + p[i] / (par.gamma - 1.0) / rho[i]);
            rho_u[i] = rho[i] * 0.0;
            x[i] = (i + 0.5 - par.nx_fict) * par.dx;
        }
    }
    else
    {
        std::cerr << "Initial condition is not set: " << par.ic_preset << std::endl;
        exit(1);
    }
}

void Solver_Godunov1D::get_time_step()
{
    double min_dt = 1.0e6;
    for (unsigned int i = 1; i < par.nx_all; ++i)
    {
        double u = rho_u[i] / rho[i];
        double c = (rho[i] != 0.0) ? c = std::sqrt(par.gamma * p[i] / rho[i]) : 0.0;
        double dt_temp = par.CFL * par.dx / (c + fabs(u));
        if (dt_temp < min_dt) min_dt = dt_temp;
    }
    dt = min_dt;
}

void Solver_Godunov1D::write_data()
{
    std::string file_name = "data/" + std::to_string(step) + ".csv";
    std::ofstream file(file_name);
    double rho_out[par.nx];
    double u_out[par.nx];
    double int_e_out[par.nx];
    double p_out[par.nx];
    double x_out[par.nx];

    for(unsigned int i = 0; i < par.nx; i++)
    {
        rho_out[i] = rho[par.nx_fict + i];
        u_out[i] = rho_u[par.nx_fict + i] / rho_out[i];
        p_out[i] = p[par.nx_fict + i];
        int_e_out[i] = rho_e[par.nx_fict + i] / rho_out[i] - pow(u_out[i], 2) / 2.0;
        x_out[i] = x[par.nx_fict + i];

        file << x_out[i] << " " << rho_out[i]  << " " << u_out[i] << " " << p_out[i] << "\n";

    }
}

void choose_solution(Parameters par, \
                     double p_l, double p_r, \
                     double u_l, double u_r, \
                     double rho_l, double rho_r, \
                     double S, double u_s, double p_s, \
                     double& p_temp, double& rho_temp, double& u_temp)
{
    // Отбор решения
    double c_l = (rho_l != 0.0) ? sqrt(par.gamma * p_l / rho_l) : 0.0;
    double c_r = (rho_r != 0.0) ? sqrt(par.gamma * p_r / rho_r) : 0.0;
    double D;
    if (S <= u_s) // Левое
    {
        if (p_s > p_l) // Ударная волна
        {
            D = u_l - c_l * sqrt((par.gamma + 1) / 2 / par.gamma * p_s / p_l + (par.gamma - 1) / 2 / par.gamma); // Toro, 182 (161)
            if (S < D) // Невозмущенное
            {
                rho_temp = rho_l;
                u_temp = u_l;
                p_temp = p_l;
            }
            else // Перед контактным разрывом с ад. Гюгонио
            {
                rho_temp = rho_l * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l) / ((par.gamma - 1) * p_s + (par.gamma + 1) * p_l);
                u_temp = u_s;
                p_temp = p_s;
            }
        }
        else // Волна разрежения
        {
            double c_s = c_l * pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma);
            // Левая граница волны разрежения
            double s_h = u_l - c_l;
            // Правая граница волны разрежения
            double s_t = u_s - c_s;
            if (S < s_h) // Невозмущенное
            {
                rho_temp = rho_l;
                u_temp = u_l;
                p_temp = p_l;
            }
            else
            {
                if (S > s_t) // Снаружи волны разрежения
                {
                    rho_temp = rho_l * pow(p_s / p_l, 1 / par.gamma);
                    u_temp = u_s;
                    p_temp = p_s;
                }
                else // Внутри волны разрежения, Toro, 181 (160)
                {
                    c_s = 2 / (par.gamma + 1) * (c_l + (par.gamma - 1) * (u_l - S) / 2);
                    rho_temp = rho_l * pow(c_s / c_l, 2 / (par.gamma - 1));
                    u_temp = 2 / (par.gamma + 1) * (c_l + (par.gamma - 1) * u_l / 2 + S);
                    p_temp = p_l * pow(c_s / c_l, 2 * par.gamma / (par.gamma - 1));
                }
            }
        }
    }
    else // Правое
    {
        if (p_s > p_r) // Ударная волна
        {
            D = u_r + c_r * sqrt((par.gamma + 1) / 2 / par.gamma * p_s / p_r + (par.gamma - 1) / 2 / par.gamma); // 1-е соотношение Гюгонио, Toro, 182 (161)
            if (S > D) // Невозмущенное
            {
                rho_temp = rho_r;
                u_temp = u_r;
                p_temp = p_r;
            }
            else // Перед контактным разрывом с адиаб. Гюгонио
            {
                rho_temp = rho_r * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r) / ((par.gamma - 1) * p_s + (par.gamma + 1) * p_r);
                u_temp = u_s;
                p_temp = p_s;
            }
        }
        else // Волна разрежения
        {
            double c_s = c_r * pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma);
            double s_h = u_r + c_r; // Правая граница волны разрежения
            double s_t = u_s + c_s; // Левая граница волны разрежения
            if (S > s_h) // Невозмущенное
            {
                rho_temp = rho_r;
                u_temp = u_r;
                p_temp = p_r;
            }
            else
            {
                if (S < s_t) // Снаружи волны разрежения
                {
                    rho_temp = rho_r * pow(p_s / p_r, 1 / par.gamma);
                    u_temp = u_s;
                    p_temp = p_s;
                }
                else // Внутри волны разрежения, Toro, 183 (161)
                {
                    c_s = 2 / (par.gamma + 1) * (c_r - (par.gamma - 1) * (u_r - S) / 2);
                    rho_temp = rho_r * pow(c_s / c_r, 2 / (par.gamma - 1));
                    u_temp = 2 / (par.gamma + 1) * (-c_r + (par.gamma - 1) * u_r / 2 + S);
                    p_temp = p_r * pow(c_s / c_r, 2 * par.gamma / (par.gamma - 1));
                }
            }
        }
    }
}

void godunov_flux(Parameters par, \
                  double p_temp, double rho_temp, double u_temp, \
                  double& F_m, double& F_imp, double& F_e)
{
    F_m = rho_temp * u_temp;
    F_imp = rho_temp * pow(u_temp, 2) + p_temp;
    F_e = (p_temp + rho_temp * pow(u_temp, 2) / 2.0 + p_temp / (par.gamma - 1)) * u_temp;
}

void newton_iteration(Parameters par, \
                      double p_l, double p_r, \
                      double u_l, double u_r, \
                      double rho_l, double rho_r, \
                      double& p_s)
{
    double c_l = (rho_l != 0.0) ? sqrt(par.gamma * p_l / rho_l) : 0.0;
    double c_r = (rho_r != 0.0) ? sqrt(par.gamma * p_r / rho_r) : 0.0;
    double fl, fr, dfl, dfr;
    if (p_s <= p_l) // Волна разрежения
    {
        fl = 2 * c_l / (par.gamma - 1) * (pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma) - 1);
        dfl = c_l / par.gamma * pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma - 1) / p_l;
    }
    else // Ударная волна
    {
        fl = (p_s - p_l) / sqrt(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l));
        int N_iter = 20; // как в Toro, критерия по точности нет
        dfl = rho_l / 4 * ((par.gamma + 1) * p_s + (3 * par.gamma - 1) * p_l) / pow(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l), 1.5);
    }
    if (p_s <= p_r) // Волна разрежения
    {
        fr = 2 * c_r / (par.gamma - 1) * (pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma) - 1);
        dfr = c_r / par.gamma * pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma - 1) / p_r;
    }
    else // Ударная волна
    {
        fr = (p_s - p_r) / sqrt(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r));
        dfr = rho_r / 4 * ((par.gamma + 1) * p_s + (3 * par.gamma + 1) * p_r) / pow(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r), 1.5);
    }
    p_s = p_s - (fl + fr - (u_l - u_r)) / (dfl + dfr);
    if (p_s < 0.0)
    {
        p_s = 1e-6;
    }
}

void exact_solution(Solver_Godunov1D solver)
{
    std::string file_name = "data/exact_solution.csv";
    std::ofstream file(file_name);
    Parameters par = solver.par;
    double p[par.nx_all];
    double rho[par.nx_all];
    double u[par.nx_all];
    double* x = solver.x;

    double p_l, p_r, u_l, u_r,rho_l, rho_r;
    if(par.ic_preset == IC_TEST1)
    {
        p_l = 1.0;
        p_r = 0.1;
        rho_l = 1.0;
        rho_r = 0.125;
        u_l = 0.0;
        u_r = 0.0;
    }
    else if(par.ic_preset == IC_TEST2)
    {
        p_l = 0.4;
        p_r = 0.4;
        rho_l = 1.0;
        rho_r = 1.0;
        u_l = -2.0;
        u_r = 2.0;
    }
    else
    {
        std::cerr << "Initial condition is not set: " << par.ic_preset << std::endl;
        exit(1);
    }
    double c_l = (rho_l != 0.0) ? sqrt(par.gamma * p_l / rho_l) : 0.0;
    double c_r = (rho_r != 0.0) ? sqrt(par.gamma * p_r / rho_r) : 0.0;
    // Условие из Toro, 175 (153). Вывод, 165 (143)
    if (2.0 / (par.gamma - 1) * (c_l + c_r) <= (u_r - u_l))
    {
        std::cerr << "Vaccum!" << std::endl;
        exit(1);
    }
    double p_s = 0.5 * (p_l + p_r);
    // Метод Ньютона - получим давление на контактном разрыве
    double fl, fr;
    double dfl, dfr;
    int N_iter = 20; // как в Toro, критерия по точности нет
    for (int k = 0; k < N_iter; ++k)
    {
        newton_iteration(par, \
                         p_l, p_r, \
                         u_l, u_r, \
                         rho_l, rho_r, \
                         p_s);
    }
    // Волна разрежения
    if (p_s <= p_l) fl = 2 * c_l / (par.gamma - 1) * (pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma) - 1);
    // Ударная волна
    else fl = (p_s - p_l) / sqrt(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l));

    // Волна разрежения
    if (p_s <= p_r) fr = 2 * c_r / (par.gamma - 1) * (pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma) - 1);
    // Ударная волна
    else fr = (p_s - p_r) / sqrt(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r));
    double u_s = 0.5 * (u_l + u_r + fr - fl);
    for(int i = par.nx_fict ; i < par.nx + par.nx_fict; i++)
    {
        double S = (x[i] - 0.5)/solver.t;
        choose_solution(par, \
                        p_l, p_r, \
                        u_l, u_r, \
                        rho_l, rho_r, \
                        S, u_s, p_s, \
                        p[i], rho[i], u[i]);
        file << x[i] << " " << rho[i]  << " " << u[i] << " " << p[i] << "\n";
    }
}

void godunov_reconstruction(Parameters par, \
                            double* rho_left, double* rho_right, \
                            double* rho_u_left, double* rho_u_right, \
                            double* rho_e_left, double* rho_e_right, \
                            double* p_left, double* p_right, \
                            double* rho, double* rho_u, double* rho_e, double* p)
{
    for (int i = 0; i < par.nx; ++i)
    {
        int i_temp = i + par.nx_fict;
        rho_left[i + 1] = rho[i_temp];
        rho_u_left[i + 1] = rho_u[i_temp];
        rho_e_left[i + 1] = rho_e[i_temp];
        p_left[i + 1] = p[i_temp];
        rho_right[i] = rho[i_temp];
        rho_u_right[i] = rho_u[i_temp];
        rho_e_right[i] = rho_e[i_temp];
        p_right[i] = p[i_temp];
    }
}

void calc_all_flux(Parameters par, \
                   double* F_m, double* F_imp, double* F_e, \
                   double* rho, double* rho_u, double* rho_e, double* p)
{
    double rho_left [par.nx_all + 1];
    double rho_u_left[par.nx_all + 1];
    double rho_e_left[par.nx_all + 1];
    double p_left[par.nx_all + 1];
    double rho_right[par.nx_all + 1];
    double rho_u_right[par.nx_all + 1];
    double rho_e_right[par.nx_all + 1];
    double p_right[par.nx_all + 1];
    godunov_reconstruction(par, \
                           rho_left, rho_right, \
                           rho_u_left, rho_u_right, \
                           rho_e_left, rho_e_right, \
                           p_left, p_right, \
                           rho, rho_u, rho_e, p);
    // Решатель Римана - получим параметры на грани ячейки
    for (int i = 1; i < par.nx; ++i)
    {
        double p_l = p_left[i];
        double p_r = p_right[i];
        double u_l = rho_u_left[i] / rho_left[i];
        double u_r = rho_u_right[i] / rho_right[i];
        double rho_l = rho_left[i];
        double rho_r = rho_right[i];
        double c_l = (rho_l != 0.0) ? sqrt(par.gamma * p_l / rho_l) : 0.0;
        double c_r = (rho_r != 0.0) ? sqrt(par.gamma * p_r / rho_r) : 0.0;
        // Условие из Toro, 175 (153). Вывод, 165 (143)
        if (2.0 / (par.gamma - 1) * (c_l + c_r) <= (u_r - u_l))
        {
            std::cerr << "Vaccum!" << std::endl;
            exit(1);
        }
        double p_s = 0.5 * (p_l + p_r);
        // Метод Ньютона - получим давление на контактном разрыве
        double fl, fr;
        double dfl, dfr;
        int N_iter = 20; // как в Toro, критерия по точности нет
        for (int k = 0; k < N_iter; ++k)
        {
            newton_iteration(par, \
                             p_l, p_r, \
                             u_l, u_r, \
                             rho_l, rho_r, \
                             p_s);
        }
        // Волна разрежения
        if (p_s <= p_l) fl = 2 * c_l / (par.gamma - 1) * (pow(p_s / p_l, (par.gamma - 1) / 2 / par.gamma) - 1);
        // Ударная волна
        else fl = (p_s - p_l) / sqrt(rho_l / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_l));

        // Волна разрежения
        if (p_s <= p_r) fr = 2 * c_r / (par.gamma - 1) * (pow(p_s / p_r, (par.gamma - 1) / 2 / par.gamma) - 1);
        // Ударная волна
        else fr = (p_s - p_r) / sqrt(rho_r / 2 * ((par.gamma + 1) * p_s + (par.gamma - 1) * p_r));

        double u_s = 0.5 * (u_l + u_r + fr - fl);
        // Отбор решения
        double u_temp = 0.0;
        double p_temp = 0.0;
        double rho_temp = 0.0;
        double S = 0.0;
        choose_solution(par, \
                        p_l, p_r, \
                        u_l, u_r, \
                        rho_l, rho_r, \
                        S, u_s, p_s, \
                        p_temp, rho_temp, u_temp);
        // Поток Годунова
        godunov_flux(par, \
                     p_temp, rho_temp, u_temp, \
                     F_m[i], F_imp[i], F_e[i]);
    }
}
