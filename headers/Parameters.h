#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "auxiliary_functions.h"
#include "Wall.h"

const double R = 8.31;

enum test_preset {
	TEST_CUSTOM,
	TEST1,
	TEST2,
	TEST3,
	TEST4
};

enum IC_preset {
	IC_CUSTOM,
	IC_TEST1,
	IC_TEST2,
	IC_TEST3,
	IC_TEST4
};

enum VISC_types {
	VISC_NONE,
	VISC_NEUMAN,
	VISC_LATTER,
	VISC_LINEAR,
	VISC_SUM
};

enum Reconstruction {
	Godunov,
	Kolgan_1972,
	Kolgan_1975,
	Osher_1984
};

struct Parameters {
	double g_x;
	double x_start, x_end;
	unsigned int nx; // amount of x ceils
	double dx;
	unsigned int nx_fict;
	unsigned int nx_all;
	double CFL;
	unsigned int nt;
	test_preset test;
	Wall walls[2]; // only 2 walls for 1D situation
	IC_preset ic_preset;
	VISC_types viscosity;
	Reconstruction reconstruction;
	bool is_conservative;
	double mu0;

	double gamma;
	double Cv;
	double Cp;

	unsigned int nt_write;
	std::string write_file;

	Parameters(std::string);
	Parameters(const Parameters& rhs);
};

#endif
