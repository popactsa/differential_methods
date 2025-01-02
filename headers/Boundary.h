#ifndef BOUNDARY_H
#define BOUNDARY_H

enum BoundaryType {
	B_WALL,
	B_FLUX,
};

struct Boundary {
	BoundaryType b_type;
	double T;
	double v_x;
};
#endif
