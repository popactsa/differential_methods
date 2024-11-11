#ifndef SOLVER_H
#define SOLVER_H

#include "Manager.h"
#include<math.h>

class Solver
{
public:
    // Manager manager;

    // Solver(Manager& man){manager = man;};
    virtual void apply_boundary_conditions() = 0;
    virtual void get_time_step() = 0;
    virtual void solve_problem() = 0;
};

#endif
