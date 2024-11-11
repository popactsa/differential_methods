#ifndef WENO3_H
#define WENO3_H

#include "Manager.h"
#include <math.h>
#include "Solver.h"

class WENO3: public Solver
{
public:
    Manager manager;

    WENO3(Manager&);
    void apply_boundary_conditions();
    void get_time_step();
    void solve_problem();
};

#endif
