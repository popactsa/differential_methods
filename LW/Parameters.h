#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>

class Parameters
{
public:
    double x_start, x_end;
    int N_x;
    double CFL;
    double t_end;
    std::string solver_name;
    int num_boundaries;
    int n_write;

    double gamma = 1.4;
    double R = 8.31;
    double Cv = R/(gamma -1);
    double Cp = gamma*Cv;

};

#endif // PARAMETERS_H
