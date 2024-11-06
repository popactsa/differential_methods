#ifndef MANAGER_H
#define MANAGER_H

#include "Parameters.h"
#include "Boundary.h"
#include<iostream>

class Manager
{
public:
    Parameters parameters;
    int fict;
    int N_all;
    Boundary* boundarys;
    double dx;
    double* xc;
    double* xb;
    double* c;

    Manager(Parameters&);
    void set_boundary_conditions(Boundary*);
    void create_mesh();
    void create_fields();
    void set_initial_conditions();
    ~Manager();
};

#endif
