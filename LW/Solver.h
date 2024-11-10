#ifndef SOLVER_H
#define SOLVER_H

#include "Manager.h"
#include<math.h>

class Solver
{
public:
    Manager manager;

    Solver(Manager&);
    virtual void apply_boundary_conditions() const = 0;
    virtual void get_time_step() const  = 0;
    virtual void solve_problem() const = 0;
};

#endif
