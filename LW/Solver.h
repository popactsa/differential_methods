#ifndef SOLVER_H
#define SOLVER_H

#include "Manager.h"
#include<math.h>

class Solver
{
public:
    Manager manager;

    Solver(Manager&);
    void apply_boundary_conditions();
    void get_time_step();
    void solve_problem();
};

#endif
