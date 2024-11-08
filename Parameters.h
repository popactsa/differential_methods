#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "auxiliary_functions.h"

enum E_BOUNDARY_TYPE {IN, OUT, FREE_FLUX, WALL, SYM, PERIODIC, NO_REF};

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
	float g_x, g_y, g_z;
	float x_start, x_end;
	float y_start, y_end;
	float z_start, z_end;
	unsigned int nx, ny, nz;
	float CFL;
	float t_end;
	float dt;
	unsigned int n_walls;
	std::vector<E_BOUNDARY_TYPE> wall_type;
	
	float Cp, Cv;
	float gamma = Cp / Cv;
	float R;
	
	float t_write;
	unsigned int nt_write;
	std::string write_file;
	
public:
	Parameters(std::string);
	inline void print_g_x() const {std::cout << g_x << std::endl;}
	inline void print_write_file() const {std::cout << write_file << std::endl;}
};

#endif
