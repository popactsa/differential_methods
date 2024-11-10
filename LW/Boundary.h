#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "BoundaryType.h"

class Boundary
{
public:
	BOUNDARY_TYPE name;
	double rho, u, v, w, e, T, p;
};

#endif
