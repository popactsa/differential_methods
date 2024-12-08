#ifndef BOUNDARY_H
#define BOUNDARY_H

enum BoundaryType {
	B_WALL,
	B_FLUX,
};

enum BoundaryPreset {
	B_CUSTOM,
	B_TEST1,
	B_TEST2,
	B_TEST3,
	B_TEST4,
};

struct Boundary {
	BoundaryPreset b_preset; 
	BoundaryType b_type;
	double T;
	double v_x;
};
#endif
