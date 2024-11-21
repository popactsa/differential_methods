#include "Solver_Lagrange1D.h"
#include <cmath>
#include <fstream>
Solver_Lagrange1D::Solver_Lagrange1D(const Parameters& _par):
    par(_par)
{
    p = new double[par.nx_all];
    rho = new double[par.nx_all];
    U = new double[par.nx_all];
    m = new double[par.nx_all];
    v = new double[par.nx_all + 1];
    x = new double[par.nx_all + 1];
    omega = new double[par.nx_all];

    char ch;
    set_initial_conditions();
    for (step = 1; step <= par.nt; step++)
    {
        solve_step();
        t += dt;
        if (step % par.nt_write == 0)
        {
            /*std::cout << "====" << step << "====" << std::endl;
            std::cout << "dt\t" << dt << std::endl;
            std::cout << "p\t" << p[par.nx / 2] << std::endl;
            std::cout << "rho\t" << rho[par.nx / 2] << std::endl;
            std::cout << "U\t" << U[par.nx / 2] << std::endl;
            std::cout << "m\t" << m[par.nx / 2] << std::endl;
            std::cout << "v\t" << v[par.nx / 2] << std::endl;
            std::cout << "x\t" << x[par.nx / 2] << std::endl;
            std::cout << std::endl;*/
            write_data();
        };
//		std::cin >> ch;
    };
    std::cout << "done" << std::endl;
};

Solver_Lagrange1D::~Solver_Lagrange1D()
{
    delete[] p;
    delete[] rho;
    delete[] U;
    delete[] m;
    delete[] v;
    delete[] x;
    delete[] omega;
};

void Solver_Lagrange1D::apply_boundary_conditions()
{
    // velocity
    for (unsigned int i = 0; i < 2; ++i)   // 1D -> 2 walls
    {
        if (par.walls[i].type == WALL)
        {
            if(par.ic_preset == IC_TEST1 || par.ic_preset == IC_TEST3  || par.ic_preset == IC_TEST4 )
            {
                v[par.nx_fict] = 0.0;
                v[par.nx_all - par.nx_fict] = 0.0;
            }
            else if(par.ic_preset == IC_TEST2)
            {
                v[par.nx_fict] = -2.0;
                v[par.nx_all - par.nx_fict] = 2.0;
            }
            else
            {
                std::cerr << "For nx_fict = " << par.nx_fict << " bc is not specified" << std::endl;
                exit(1);
            };
            int i_temp = 0;
            if(i)
                i_temp = par.nx_all;
            v[i_temp] = -v[i_temp+(1-2*i)*2] + 2.0 * v[i_temp+(1-2*i)*1];
            // v[0] = -v[2] + 2.0 * v[1];
            // v[par.nx_all] = -v[par.nx_all - 2] + 2.0 * v[par.nx_all - 1];
        }
        else if (par.walls[i].type == FLUX)
        {
            int i_temp = 0;
            if(i)
                i_temp = par.nx_all;
            v[i_temp] = v[i_temp+(1-2*i)*2];
            v[i_temp+(1-2*i)*1] = v[i_temp+(1-2*i)*2];
        }
    };
    // density and internal energy
    for (unsigned int i = 0; i < par.nx_fict; ++i)
    {
        rho[i] = rho[par.nx_fict];
        rho[par.nx_all - par.nx_fict - i] = rho[par.nx_all - par.nx_fict - 1];
        U[i] = U[par.nx_fict];
        U[par.nx_all - par.nx_fict - i] = U[par.nx_all - par.nx_fict - 1];
    };
};

void Solver_Lagrange1D::solve_step()
{
    apply_boundary_conditions();
    get_time_step();
    double* v_last = new double[par.nx_all+1];
    for(unsigned int i = 0; i < par.nx_all + 1; i++)
    {
        v_last[i] = v[i];
    }
    // solver
    // viscosity
    if (par.viscosity == VISC_NONE)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
            omega[i] = 0.0;
    }
    else if (par.viscosity == VISC_NEUMAN)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
            omega[i] = -mu0*rho[i]*fabs(v[i+1] - v[i])*(v[i+1] - v[i]);
    }
    else if (par.viscosity == VISC_LATTER)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
            omega[i] = ((v[i+1] - v[i]) < 0.0) ? mu0*rho[i]*(v[i+1] - v[i])*(v[i+1] - v[i]) : 0.0;
    }
    else
    {
        std::cerr << "Viscosity is not set: " << par.viscosity << std::endl;
        exit(1);
    };
    // velocity
    for (unsigned int i = 1 + par.nx_fict; i < par.nx_all - par.nx_fict; ++i)
    {
        double dx = 0.5 * (m[i] + m[i - 1]);
        v[i] = v[i] - ((p[i] + omega[i]) - (p[i - 1] + omega[i - 1])) * dt / dx;
    };
    // coordinate
    for (unsigned int i = 0; i < par.nx_all + 1; ++i)
    {
        x[i] += v[i] * dt;
    };
    // pressure at the boundaries for conservative law
    double* pc = new double[par.nx_all];
    for(unsigned int i = 1; i < par.nx_all; i++)
    {
        pc[i] = 0.5*(p[i] + omega[i] + p[i-1] + omega[i-1]);
    }
    // density and internal energy
    for (unsigned int i = par.nx_fict; i < par.nx_all - par.nx_fict; ++i)
    {
        rho[i] = rho[i] / (rho[i] * (v[i + 1] - v[i]) * dt / m[i] + 1.0);
        if (par.is_conservative == false)
            U[i] = U[i] / (rho[i] * (v[i + 1] - v[i]) * (par.gamma - 1.0) * dt / m[i] + 1.0);
        else
        {
            double U_temp = U[i];
            U[i] = U[i] - (v[i+1]*pc[i+1] - v[i]*pc[i])*dt/m[i] + (v_last[i+1] + v_last[i])*(v_last[i+1] + v_last[i])/8.0 - (v[i+1] + v[i])*(v[i+1] + v[i])/8.0;
            if (U[i] < 0.0)
            {
                U[i] = U_temp / (rho[i] * (v[i + 1] - v[i]) * (par.gamma - 1.0) * dt / m[i] + 1.0);
            }
        }
    };
    delete[] pc;
    delete[] v_last;
    // pressure
    for (unsigned int i = 0; i < par.nx_all; ++i)
    {
        p[i] = rho[i] * (par.gamma - 1.0) * U[i];
        // std::cout << "p " << p[i] << std::endl;
    };
};

void Solver_Lagrange1D::set_initial_conditions()
{
    if (par.ic_preset == IC_TEST1)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                p[i] = 1.0;
                rho[i] = 1.0;
            }
            else
            {
                p[i] = 0.1;
                rho[i] = 0.125;
            };
            U[i] = p[i] / (par.gamma - 1.0) / rho[i];
        };
        for (unsigned int i = 0; i < par.nx_all + 1; ++i)
        {
            v[i] = 0.0;
            x[i] = par.x_start + (i - par.nx_fict) * par.dx;
        };
    }
    else if (par.ic_preset == IC_TEST2)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            p[i] = 0.4;
            rho[i] = 1.0;
            U[i] = p[i] / (par.gamma - 1.0) / rho[i];
        };
        for (unsigned int i = 0; i < par.nx_all + 1; ++i)
        {
            if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                v[i] = -2.0;
            }
            else
            {
                v[i] = 2.0;
            };
            x[i] = par.x_start - par.nx_fict * par.dx + i * par.dx;
        };
    }
    else if (par.ic_preset == IC_TEST3)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                p[i] = 1000.0;
                rho[i] = 1.0;
            }
            else
            {
                p[i] =  0.01;
                rho[i] = 1.0;
            }
            U[i] = p[i] / (par.gamma - 1.0) / rho[i];
        };
        for (unsigned int i = 0; i < par.nx_all + 1; ++i)
        {
            v[i] = 0.0;
            x[i] = par.x_start - par.nx_fict * par.dx + i * par.dx;
        };
    }
    else if (par.ic_preset == IC_TEST4)
    {
        for (unsigned int i = 0; i < par.nx_all; ++i)
        {
            if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                p[i] = 0.01;
                rho[i] = 1.0;
            }
            else
            {
                p[i] =  100.0;
                rho[i] = 1.0;
            }
            U[i] = p[i] / (par.gamma - 1.0) / rho[i];
        };
        for (unsigned int i = 0; i < par.nx_all + 1; ++i)
        {
            if (i * par.dx <= par.x_start + 0.5 * (par.x_end - par.x_start) + par.nx_fict * par.dx)
            {
                v[i] = 0.0;
            }
            else
            {
                v[i] = 0.0;
            };
            x[i] = par.x_start - par.nx_fict * par.dx + i * par.dx;
        };
    }
    else
    {
        std::cerr << "Initial condition is not set: " << par.ic_preset << std::endl;
        exit(1);
    };
    for (unsigned int i = 0; i < par.nx_all; ++i) m[i] = rho[i] * (x[i + 1] - x[i]);
};

void Solver_Lagrange1D::get_time_step()
{
    double min_dt = 1.0e6;
    for (unsigned int i = 1; i < par.nx_all; ++i)
    {
        double dx = x[i+1] - x[i];
        double V = 0.5*(v[i+1] + v[i]);
        double c = std::sqrt(par.gamma * p[i] / rho[i]);
        double dt_temp = par.CFL * dx / (c + fabs(V));
        if (dt_temp < min_dt) min_dt = dt_temp;
    }
    dt = min_dt; // what's going on????????
    // std::cout << "dt = " << dt << std::endl;
}

void Solver_Lagrange1D::write_data()
{
    std::string file_name = "data/" + std::to_string(step) + ".csv";
    std::ofstream file(file_name);
    double* x_grid = new double[par.nx_all];
    double* v_grid = new double[par.nx_all];
    for(unsigned int i = 0; i < par.nx_all; i++)
    {
        x_grid[i] = 0.5*(x[i+1] + x[i]);
        v_grid[i] = 0.5*(v[i+1] + v[i]);
    }
    double* c = rho;
    for(int i = par.nx_fict; i < par.nx_all - par.nx_fict; i++)
    {
        file << x_grid[i] << " " << rho[i]  << " " << p[i] << " " << v_grid[i] << "\n";
    }
    delete[] x_grid;
    delete[] v_grid;
}
