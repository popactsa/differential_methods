#include "Manager.h"

Manager::Manager(Parameters& param): parameters(param)
{
    if(param.solver_name == "WENO3")
        fict = 1;
    else
    {
        std::cout << "No such solver" << param.solver_name << std::endl;
        exit(1);
    }
    boundaries = new Boundary[2];
    dx = (param.x_end - param.x_start)/double(param.N_x);
    N_all = param.N_x + 2*fict;
    t_done = 0;
}

Manager::~Manager()
{
    delete(boundaries);
    delete(xc);
    delete(xb);
    delete(rho, u, v, w, e);
    delete(f_rho, f_u, f_v, f_w, f_e);
}

void Manager::create_mesh()
{
    xc = new double[parameters.N_x+2*fict];
    for(int i = 0; i < parameters.N_x+2*fict; i++)
    {
        xc[i] = dx*(1./2 + i - fict);
    }
    xb = new double[parameters.N_x+2*fict+1];
    for(int i = 0; i < parameters.N_x+2*fict+1; i++)
    {
        xb[i] = dx*(i - fict);
    }
}

void Manager::create_fields()
{
    int N_x = parameters.N_x;
    // поля
    rho = new double[N_all] {0};
    u   = new double[N_all] {0};
    v   = new double[N_all] {0};
    w   = new double[N_all] {0};
    e   = new double[N_all] {0};
    // потоки
    f_rho = new double[N_x+1] {0};
    f_u   = new double[N_x+1] {0};
    f_v   = new double[N_x+1] {0};
    f_w   = new double[N_x+1] {0};
    f_e   = new double[N_x+1] {0};
}

void Manager::set_initial_conditions()
{
    for(int i = 0; i < N_all; i++)
    {
        /*if(xc[i] < 0.4)
        {
            rho[i] = 1.0;
            u[i] = 1*0.1;
            e[i] = 8.31/(parameters.gamma - 1.0) + 1.0*(0.1*0.1)/2.;
        }
        else
        {
            rho[i] = 0.5;
            u[i] = 0.5*0.01;
            e[i] = 8.31/(parameters.gamma - 1.0) + 0.5*(0.1*0.1)/2.;
        }*/
        /*rho[i] = 2.0 + std::sin(xc[i]);
        u[i] = rho[i]*0.1;
        e[i] = 2*8.31/(parameters.gamma - 1.0) + rho[i]*(0.1*0.1)/2.;*/
        rho[i] = 1.0;
        u[i] = 1.0*0.1;
        e[i] = 8.31/(parameters.gamma - 1.0) + 1.0*(0.1*0.1)/2.;
        v[i] = 0;
        w[i] = 0;
    }
}

void Manager::set_boundary_conditions(Boundary* bound_array)
{
    for(int i = 0; i < 2; i++)
    {
        boundaries[i] = bound_array[i];
        if(boundaries[i].name == "IN")
        {
            boundaries[i].e = parameters.Cv*boundaries[i].T \
                        + boundaries[i].u*boundaries[i].u/2.0;
        }
    }
}
