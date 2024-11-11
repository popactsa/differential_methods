#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "BoundaryType.h"

class Boundary
{
public:
	E_BOUNDARY_TYPE name;
	double rho, u, v, w, e, T, p;
};

#endif
