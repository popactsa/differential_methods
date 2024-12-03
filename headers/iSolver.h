#ifndef ISOLVER_H
#define ISOLVER_H

class iSolver
{
public:
    virtual void apply_boundary_conditions() = 0;
    virtual void solve_step() = 0;
    virtual void set_initial_conditions() = 0;
    virtual void get_time_step() = 0;
    virtual void write_data() = 0;
};

#endif
