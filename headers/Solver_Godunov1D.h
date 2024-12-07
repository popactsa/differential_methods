#ifndef SOLVER_GODUNOV1D_H
#define SOLVER_GODUNOV1D_H

#include "Parameters.h"
#include "iSolver.h"

class Solver_Godunov1D: public iSolver
{
public:
    const Parameters par;
    double *p, *rho, *rho_u, *rho_e;
    double *x;
    double t, dt;
    unsigned int step;
    double *F_m, *F_imp, *F_e;
public:
    Solver_Godunov1D(const Parameters& _par);
    void apply_boundary_conditions();
    void solve_step();
    void set_initial_conditions();
    void get_time_step();
    void write_data();
    ~Solver_Godunov1D();
};

void choose_solution(Parameters par, \
                     double p_l, double p_r, \
                     double u_l, double u_r, \
                     double rho_l, double rho_r, \
                     double S, double u_s, double p_s, \
                     double& p_temp, double& rho_temp, double& u_temp);

void godunov_flux(Parameters par, \
                  double p_temp, double rho_temp, double u_temp, \
                  double& F_m, double& F_imp, double& F_e);

void newton_iteration(Parameters par, \
                      double p_l, double p_r, \
                      double u_l, double u_r, \
                      double rho_l, double rho_r, \
                      double& p_s);

void exact_solution(Solver_Godunov1D solver);

void godunov_reconstruction(Parameters par, \
                            double* rho_left, double* rho_right, \
                            double* rho_u_left, double* rho_u_right, \
                            double* rho_e_left, double* rho_e_right, \
                            double* p_left, double* p_right, \
                            double* rho, double* rho_u, double* rho_e, double* p);
//TODO как нормально назвать функцию
void calc_all_flux(Parameters par, \
                   double* F_m, double* F_imp, double* F_e, \
                   double* rho, double* rho_u, double* rho_e, double* p);
#endif
