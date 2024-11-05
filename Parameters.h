#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <iostream>
#include <fstream>
enum E_BOUNDARY_TYPE {IN, OUT, FREE_FLUX, WALL, SYM, PERIODIC, NO_REF};
class Parameters {
    double g_x, g_y, g_z;
    double x_start, x_end;
    double y_start, y_end;
    double z_start, z_end;
    int nx, ny, nz;
    double CFL;
    double t_end;
    double dt;
    E_BOUNDARY_TYPE type;

    double Cp, Cv;
    double gamma = Cp / Cv;
    double R;

    double t_write;
    short int nt_write;
    std::string write_file;

public:
    Parameters(std::string);
};

#endif