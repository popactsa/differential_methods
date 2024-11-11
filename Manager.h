#ifndef MANAGER_H
#define MANAGER_H

#include "Parameters.h"
#include "Boundary.h"
#include<iostream>
#include<math.h>

class Manager
{
public:
    Parameters parameters;
    int fict;
    int N_all;
    Boundary* boundaries;
    double dt, t;
    int t_done;

    double dx;
    double* xc;
    double* xb;
    double *rho, *u, *v, *w, *e;
    double *f_rho, *f_u, *f_v, *f_w, *f_e;

    Manager(Parameters&);
    void set_boundary_conditions(Boundary*);
    void create_mesh();
    void create_fields();
    void set_initial_conditions();
    ~Manager();
};

#endif
