#include "Solver.h"

Solver::Solver(Manager& man): manager(man)
{
    t = 0;
    dt = 0;
    M = 0;
}

void Solver::apply_boundary_conditions()
{
    int fict = manager.fict;
    double* c = manager.c;
    for(int i = 0; i < fict; i ++)
    {
        c[i] = c[fict];
        c[manager.N_all-i-1] = c[manager.N_all-fict-1];
    }
}

void Solver::get_time_step()
{
    dt = manager.parameters.CFL * manager.dx;
}

void Solver::solve_problem()
{
    int fict = manager.fict;
    double* ctemp = new double[manager.parameters.N+2*fict];
    double* c = manager.c;
    double CFL = manager.parameters.CFL;
    for(int i = fict;i < manager.N_all-fict; i++)
    {
        ctemp[i] = c[i] - CFL*(c[i] - c[i-1]) - 0.5*(1-CFL)*CFL*(c[i+1]-2*c[i]+c[i-1]);
    }
    for(int i = fict;i < manager.N_all-fict; i++)
    {
        c[i] = ctemp[i];
    }

    t += dt;
    M ++;
}
