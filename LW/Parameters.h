#ifndef PARAMETERS_H
#define PARAMETERS_H

#include<string>

class Parameters
{
public:
    double l;
    int N;
    double CFL;
    double t_end;
    std::string solver_name;
    int num_boundarys;
    int n_write;
};

#endif // PARAMETERS_H
