#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include "auxiliary_functions.h"
#include "Boundary.h"

const double R = 8.31;

enum InitialPreset {
	IC_CUSTOM,
	IC_TEST1,
	IC_TEST2,
	IC_TEST3,
	IC_TEST4,
};

enum ArtificialViscosity {
	V_NONE,
	V_NEUMAN,
	V_LATTER,
	V_LINEAR,
	V_SUM,
};

enum Reconstruction {
	R_GODUNOV,
	R_KOLGAN72,
	R_KOLGAN75,
	R_OSHER84,
};

enum TimeAlgo {
	T_EULER,
	T_RUNGEKUTTA,
};

enum FluxScheme {
	F_GODUNOV,
	F_RLF,
};


struct Parameters {
	double x_start, x_end;
	int nx; // amount of x ceils
	double dx;
	int nx_fict;
	int nx_all;
	double CFL;
	int nt;
	Boundary boundaries[2]; // only 2 walls for 1D situation
	InitialPreset ic_preset;
	ArtificialViscosity viscosity_type;
	Reconstruction reconstruction_type;
	TimeAlgo time_algo;
	FluxScheme flux_scheme;
	bool is_conservative;
	double mu0;

	double gamma;
	double Cv;
	double Cp;

	int nt_write;
	std::string write_file;

	inline void print() {std::cout << viscosity_type << std::endl;}
	Parameters(const Parameters&) = default;
	Parameters(std::string);
//	Parameters(const Parameters& rhs);
};

#endif
