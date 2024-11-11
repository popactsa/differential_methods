#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "auxiliary_functions.h"
#include "BoundaryType.h"

static std::map<std::string, E_BOUNDARY_TYPE> wall_type_map = {
	{"IN", IN},
	{"OUT", OUT},
	{"FREE_FLUX", FREE_FLUX},
	{"WALL", WALL},
	{"SYM", SYM},
	{"PERIODIC", PERIODIC},
	{"NO_REF", NO_REF}
};


class Parameters {

public:
	double g_x, g_y, g_z;
	double x_start, x_end;
	double y_start, y_end;
	double z_start, z_end;
    unsigned int nx, ny, nz;
	double CFL;
	double t_end;
	double dt;
	unsigned int n_walls;
	std::vector<E_BOUNDARY_TYPE> wall_type;

    double gamma = 1.4;
    double R = 8.31;
    double Cv = R/(gamma-1);
    double Cp = gamma*Cv;

	float t_write;
	unsigned int nt_write;
	std::string write_file;

	Parameters(std::string);
	inline void print_g_x() const {std::cout << g_x << std::endl;}
	inline void print_write_file() const {std::cout << write_file << std::endl;}
};

#endif
