#ifndef SOLVER_H
#define SOLVER_H

#include "Manager.h"

class Solver
{
public:
    Manager manager;
    double t, dt;
    int M;

    Solver(Manager&);
    void apply_boundary_conditions();
    void get_time_step();
    void solve_problem();
};

#endif
