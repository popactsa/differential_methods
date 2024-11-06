#include "Manager.h"

Manager::Manager(Parameters& param): parameters(param)
{
    if(param.solver_name == "LW")
        fict = 1;
    else
    {
        std::cout << "No such solver" << param.solver_name << std::endl;
        exit(1);
    }
    boundarys = new Boundary[2];
    dx = param.l/double(param.N);
    N_all = param.N + 2*fict;
}

Manager::~Manager()
{
    delete(boundarys);
    delete(xc);
    delete(xb);
    delete(c);
}

void Manager::create_mesh()
{
    xc = new double[parameters.N+2*fict];
    for(int i = 0; i < parameters.N+2*fict; i++)
    {
        xc[i] = dx*(1./2 + i - fict);
    }
    xb = new double[parameters.N+2*fict+1];
    for(int i = 0; i < parameters.N+2*fict+1; i++)
    {
        xb[i] = dx*(i - fict);
    }
}

void Manager::create_fields()
{
    c = new double[parameters.N+2*fict];
}

void Manager::set_initial_conditions()
{
    for(int i = 0; i < parameters.N+2*fict; i++)
    {
        if(xc[i] < 0.4)
            c[i] = 1;
        else
            c[i] = 0;
    }
}

void Manager::set_boundary_conditions(Boundary* bound_array)
{
    for(int i = 0; i < 2; i++)
        boundarys[i] = bound_array[i];
}
