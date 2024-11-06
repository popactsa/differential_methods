#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <iostream>
#include <fstream>
enum E_BOUNDARY_TYPE {IN, OUT, FREE_FLUX, WALL, SYM, PERIODIC, NO_REF};
const class Parameters {
    float g_x, g_y, g_z;
    float x_start, x_end;
    float y_start, y_end;
    float z_start, z_end;
    unsigned int nx, ny, nz;
    float CFL;
    float t_end;
    float dt;
    E_BOUNDARY_TYPE type;

    unsigned float Cp, Cv;
    unsigned float gamma = Cp / Cv;
    unsigned float R;

    unsigned float t_write;
    unsigned int nt_write;
    std::string write_file;

public:
    Parameters(std::string);
};

#endif
