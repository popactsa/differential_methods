#include "WENO3.h"
#include "BoundaryType.h"

WENO3::WENO3(Manager& man):manager(man)
{
}

void WENO3::apply_boundary_conditions()
{
    int fict = manager.fict;
    int N_all = manager.N_all;
    int N_x = manager.parameters.nx;
    double* fields[5] = {manager.rho, manager.u, manager.v, manager.w, manager.e};
    double* flows[5] = {manager.f_rho, manager.f_u, manager.f_v, manager.f_w, manager.f_e};
    double gamma = manager.parameters.gamma;
    for(int k = 0; k < 2; k++)
    {
        Boundary bound = manager.boundaries[k];
        if(bound.name == IN)
        {
            int intern = fict;
            double rho_i = fields[4][intern] - fields[1][intern]*fields[1][intern]/fields[0][intern]/2.0;
            double a = std::sqrt((rho_i*(gamma - 1.0)*gamma/fields[0][intern]));
            bound.p = (gamma - 1.0)*rho_i \
                      + (a*fields[0][intern])*(bound.u - fields[1][intern]/fields[0][intern]);
            bound.rho = bound.p/manager.parameters.R/bound.T;
            /*std::cout << "IN" << std::endl;
            std::cout << bound.p << std::endl;
            std::cout << bound.rho << std::endl;
            std::cout << a << std::endl;*/
            for(int i = 0; i < fict; i ++)
            {
                double a = -1;
                double b = 2;
                fields[0][i] = a*fields[0][intern] + b*bound.rho;
                fields[1][i] = a*fields[1][intern] + b*bound.u*bound.rho;
                fields[2][i] = a*fields[2][intern] + b*bound.v*bound.rho;
                fields[3][i] = a*fields[3][intern] + b*bound.w*bound.rho;
                fields[4][i] = a*fields[4][intern] + b*bound.e*bound.rho;
            }
            flows[0][0] = bound.rho*bound.u;
            flows[1][0] = bound.rho*bound.u*bound.u + bound.p;
            flows[2][0] = bound.rho*bound.v*bound.u;
            flows[3][0] = bound.rho*bound.w*bound.u;
            flows[4][0] = gamma*bound.rho*bound.u*bound.e - (gamma-1)*bound.rho*bound.u*bound.u*bound.u/2.0;
        }
        else if(bound.name == OUT)
        {
            double dp;
            int intern = N_all-fict-1;
            double rho_i = fields[4][intern] - fields[1][intern]*fields[1][intern]/fields[0][intern]/2.0;
            double a = std::sqrt((rho_i*(gamma - 1.0)*gamma/fields[0][intern]));
            dp = bound.p - (gamma - 1.0)*rho_i;
            bound.v = fields[2][intern]/fields[0][intern];
            bound.w = fields[3][intern]/fields[0][intern];
            bound.u = fields[1][intern]/fields[0][intern] - dp/(a*fields[0][intern]);
            bound.rho = fields[0][intern] + dp/(a*a);
            bound.e = bound.p/(gamma - 1.0)/bound.rho + bound.u*bound.u/2.0;
            /*std::cout << "OUT" << std::endl;
            std::cout << bound.u << std::endl;
            std::cout << bound.rho << std::endl;
            std::cout << bound.e << std::endl;
            std::cout << a << std::endl;*/
            for(int i = 0; i < fict; i ++)
            {
                double a = -1;
                double b = 2;
                fields[0][N_all-i-1] = a*fields[0][intern] + b*bound.rho;
                fields[1][N_all-i-1] = a*fields[1][intern] + b*bound.u*bound.rho;
                fields[2][N_all-i-1] = a*fields[2][intern] + b*bound.v*bound.rho;
                fields[3][N_all-i-1] = a*fields[3][intern] + b*bound.w*bound.rho;
                fields[4][N_all-i-1] = a*fields[4][intern] + b*bound.e*bound.rho;
            }
            flows[0][N_x] = bound.rho*bound.u;
            flows[1][N_x] = bound.rho*bound.u*bound.u + bound.p;
            flows[2][N_x] = bound.rho*bound.v*bound.u;
            flows[3][N_x] = bound.rho*bound.w*bound.u;
            flows[4][N_x] = gamma*bound.rho*bound.u*bound.e - (gamma-1)*bound.rho*bound.u*bound.u*bound.u/2.0;
        }
        else if(bound.name == WALL)
        {
            int bc;
            int intern;
            if(k%2 == 0)
            {
                bc = 0;
                intern = fict;
            }
            else
            {
                bc = N_x;
                intern = N_all-fict-1;
            }
            for(int i = 0; i < fict; i ++)
            {
                fields[0][N_all-i-1] = fields[0][intern];
                fields[1][N_all-i-1] = -fields[1][intern];
                fields[2][N_all-i-1] = -fields[2][intern];
                fields[3][N_all-i-1] = -fields[3][intern];
                fields[4][N_all-i-1] = fields[4][intern];
            }
            flows[0][bc] = 0.0;
            flows[1][bc] = bound.p;
            flows[2][bc] = 0.0;
            flows[3][bc] = 0.0;
            flows[4][bc] = 0.0;
        }
        else if(bound.name == FREE_FLUX)
        {
            int bc;
            int intern;
            if(k%2 == 0)
            {
                bc = 0;
                intern = fict;
            }
            else
            {
                bc = N_x;
                intern = N_all-fict-1;
            }
            double dp;
            double rho_i = fields[4][intern] - fields[1][intern]*fields[1][intern]/fields[0][intern]/2.0;
            double a = std::sqrt((rho_i*(gamma - 1.0)*gamma/fields[0][intern]));
            for(int i = 0; i < fict; i ++)
            {
                fields[0][intern] = bound.rho;
                fields[1][intern] = bound.rho*bound.u;
                fields[2][intern] = bound.rho*bound.v;
                fields[3][intern] = bound.rho*bound.w;
                fields[4][intern] = bound.rho*bound.e;
            }
            // потоки надо считать
            double p = 0.5*((gamma - 1.0)*rho_i + bound.p) + \
                        0.5*fields[0][intern]*a*(fields[1][intern]/fields[0][intern] - bound.u);
            double u = 0.5*(bound.u + fields[1][intern]/fields[0][intern]) + \
                        ((gamma - 1.0)*rho_i - bound.p)/2.0/fields[0][intern]/a;
            double v = bound.v;
            double w = bound.w;
            double rho = bound.rho + 1.0/a/a*( \
                        0.5*((gamma - 1.0)*rho_i - bound.p) + \
                        0.5*(fields[0][intern]*a)*(fields[1][intern]/fields[0][intern] - bound.u));
            flows[0][N_x] = rho*u;
            flows[1][N_x] = rho*u*u + p;
            flows[2][N_x] = rho*v*u;
            flows[3][N_x] = rho*w*u;
            flows[4][N_x] = (p*u)/(gamma - 1.0) - rho*u*u*u/2.0 + p*u;
        }
    }
    /*for(int j = 0; j < 5; j++)
    {
        for(int i = 0; i < fict; i++)
        {
            fields[j][i] = fields[j][fict];
            fields[j][N_all-i-1] = fields[j][N_all-fict-1];
        }
        flows[j][0] = 0.0;
        flows[j][N_x] = 0.0;
    }*/
}

void WENO3::get_time_step()
{
    double dt_temp = 1;
    double dt_min = dt_temp;
    for(int i = 0; i < manager.N_all; i++)
    {
        if(manager.u[i] != 0.0)
        {
            dt_temp = std::abs(manager.parameters.CFL * manager.dx * manager.rho[i]/ manager.u[i]);
            if(dt_temp < dt_min)
            {
                dt_min = dt_temp;
            }
        }
    }
    manager.dt = dt_min;
}

void WENO3::solve_problem()
{
    double* fields[5] = {manager.rho, manager.u, manager.v, manager.w, manager.e};
    double* fields_left[5];
    double* fields_right[5];
    int fict = manager.fict;
    int N_x = manager.parameters.nx;
    // подготовка U_left и U_right
    for(int j = 0; j < 5; j++)
    {
        double* c = fields[j];
        fields_left[j]  = new double[N_x+1] {0};
        fields_right[j] = new double[N_x+1] {0};
        // WENO3 реконструкция поехли
        //  1 - левый 0 - правый 2(-1) - суперправый
        double* h0 = new double[N_x+1] {0};
        double* h1 = new double[N_x+1] {0};
        double* h2 = new double[N_x+1] {0};
        double* h_temp[3] = {h2, h0, h1};
        double eps = 1e-40;
        double alpha[2] {0};
        double omega[2] {0};
        double gamma[2] = {1.0/3, 2.0/3}; // стр. 112
        double beta[2] {0};
        for(int i = 1; i < N_x; i++)
        {
            // полиномы в форме Лежандра
            int i_temp = i + fict - 1;    // отступ на фиктивные
            h0[i] =  1./2.*c[i_temp]   + 1./2.*c[i_temp+1];
            h1[i] = -1./2.*c[i_temp-1] + 3./2.*c[i_temp];
            h2[i] =  3./2.*c[i_temp+1] - 1./2.*c[i_temp+2];
            // коэффициенты гладкости
            beta[0] = (c[i_temp+1] - c[i_temp])*(c[i_temp+1] - c[i_temp]);
            beta[1] = (c[i_temp]   - c[i_temp-1])*(c[i_temp]   - c[i_temp-1]);
            double sum_alpha = 0.0;
            // U_left
            for(int k = 0; k < 2; k++)
            {
                alpha[k] = (eps + beta[k])*(eps + beta[k]);
                alpha[k] = gamma[k]/alpha[k];
                sum_alpha += alpha[k];
            }
            for(int k = 0; k < 2; k++)
            {
                omega[k] = alpha[k]/(sum_alpha);
                fields_left[j][i] += omega[k]*h_temp[k+1][i];
            }
            sum_alpha = 0.0;
            // U_right
            beta[0] = (c[i_temp+2] - c[i_temp+1])*(c[i_temp+2] - c[i_temp+1]);
            beta[1] = (c[i_temp+1]   - c[i_temp])*(c[i_temp+1]   - c[i_temp]);
            for(int k = 0; k < 2; k++)
            {
                alpha[k] = (eps + beta[k])*(eps + beta[k]);
                alpha[k] = gamma[1-k]/alpha[k];
                sum_alpha += alpha[k];
            }
            for(int k = 0; k < 2; k++)
            {
                omega[k] = alpha[k]/(sum_alpha);
                fields_right[j][i] += omega[k]*h_temp[k][i];
            }
        }
        delete(h0, h1, h2);
    }

    /*std::cout << "step " <<  manager.t_done + 1 << std::endl;
    for(int i = 0; i < N_x+1; i++)
        std::cout << fields_left[4][i] << " " << fields_right[4][i] << std::endl;
    */

    // вычисление потоков
    double* flows[5] = {manager.f_rho, manager.f_u, manager.f_v, manager.f_w, manager.f_e};
    for(int j = 0; j < 5; j++)
    {
        double* f_c = flows[j];
        for(int i = 1; i < N_x; i++)
        {
            if(j == 0)
            {
                f_c[i] = fields_left[1][i];
            }
            if(j == 1)
            {
                f_c[i] = fields_left[1][i]*fields_left[1][i]/fields_left[0][i] \
                         +(manager.parameters.gamma - 1.0)*(fields_left[4][i] - \
                                 fields_left[1][i]*fields_left[1][i]/fields_left[0][i]/2.0);
            }
            if(j == 2)
            {
                f_c[i] = fields_left[1][i]*fields_left[2][i]/fields_left[0][i];
            }
            if(j == 3)
            {
                f_c[i] = fields_left[1][i]*fields_left[3][i]/fields_left[0][i];
            }
            if(j == 4)
            {
                f_c[i] = manager.parameters.gamma*fields_left[1][i]*fields_left[4][i]/fields_left[0][i] \
                         - (manager.parameters.gamma - 1.0)*fields_left[1][i]*fields_left[1][i]*fields_left[1][i] \
                         /fields_left[0][i]/fields_left[0][i]/2.;
            }
        }
    }

    for(int j = 0; j < 5; j++)
    {
        delete(fields_left[j], fields_right[j]);
    }

    // новое значение U
    for(int j = 0; j < 5; j++)
    {
        for(int i = 0; i < N_x; i++)
        {
            int i_temp = i + fict;
            fields[j][i_temp] += (flows[j][i] - flows[j][i+1])*manager.dt/manager.dx;
        }
    }

    manager.t += manager.dt;
    manager.t_done ++;
}
