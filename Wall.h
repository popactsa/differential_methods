#ifndef WALL_H
#define WALL_H

enum E_BOUNDARY_TYPE {
	WALL,
	FLUX
};

enum BC_preset {
	BC_CUSTOM,
	BC_TEST1,
	BC_TEST2,
	BC_TEST3,
	BC_TEST4
};

struct Wall {
	BC_preset bc_preset; 
	E_BOUNDARY_TYPE type;
	double T;
	double v_x;
};
#endif
