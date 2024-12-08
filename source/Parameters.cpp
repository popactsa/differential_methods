#include "Parameters.h"

Parameters::Parameters(std::string file_name) {
	std::ifstream fin(file_name);

	std::map<std::string, void*> vars = { // maps to identify values from file_name
		{"x_start", &x_start},
		{"x_end", &x_end},
		{"nx", &nx},
		{"nx_fict", &nx_fict},
		{"CFL", &CFL},
		{"nt", &nt},
		{"nt_write", &nt_write},
		{"write_file", &write_file},
		{"gamma", &gamma},
		{"is_conservative", &is_conservative},
		{"mu0", &mu0},
	};
	std::map<std::string, E_BOUNDARY_TYPE> wall_type_map = {
		{"WALL", WALL},
		{"FLUX", FLUX}
	};
	std::map<std::string, unsigned short int> wall_args = {
		{"E_BOUNDARY_TYPE", 0},
		{"T", 1},
		{"v_x", 2}
	};
	std::map<std::string, test_preset> test_preset_type_map {
		{"TEST1", TEST1},
		{"TEST2", TEST2},
		{"TEST3", TEST3},
		{"TEST4", TEST4}
	};
	std::map<std::string, IC_preset> IC_type_map = {
		{"IC_TEST1", IC_TEST1},
		{"IC_TEST2", IC_TEST2},
		{"IC_TEST3", IC_TEST3},
		{"IC_TEST4", IC_TEST4}
	};
	std::map<std::string, VISC_types> VISC_type_map = {
		{"VISC_NONE", VISC_NONE},
		{"VISC_NEUMAN", VISC_NEUMAN},
		{"VISC_LATTER", VISC_LATTER},
		{"VISC_LINEAR", VISC_LINEAR},
		{"VISC_SUM", VISC_SUM}
	};
	std::map<std::string, Reconstruction> Reconstruction_map = {
		{"Godunov", Godunov},
		{"Kolgan_1972", Kolgan_1972},
		{"Kolgan_1975", Kolgan_1975},
		{"Osher_1984", Osher_1984}
	};

	test = TEST_CUSTOM;
	ic_preset = IC_CUSTOM;
	viscosity = VISC_NONE;
	
	for (unsigned int i = 0; i < 2; ++i) walls[i].bc_preset = BC_CUSTOM;
	
	if (fin.is_open()) {
		std::string line;
		const std::string sep = "\t";
		for (std::getline(fin, line); !fin.eof(); std::getline(fin, line)) {
			unsigned int curr_sep_pos = 0;
			std::string var_type, var_name, var_value; // var_number used if index for any
													   // array element is needed
			std::string args_control; // control if there are any additional args will be provided
			get_next_substr_between_sep(var_type, line, sep, curr_sep_pos);
			get_next_substr_between_sep(var_name, line, sep, curr_sep_pos);
			get_next_substr_between_sep(var_value, line, sep, curr_sep_pos);
			
			if (strcmp(var_type.c_str(), "double") == 0)
				*(double*)vars[var_name] = std::stod(var_value);
			
			else if (strcmp(var_type.c_str(), "bool") == 0)
				*(bool*)vars[var_name] = std::stoi(var_value);
            
			else if (strcmp(var_type.c_str(), "unsignedint") == 0)
			    *(unsigned int*)vars[var_name] = std::stoi(var_value);

        	else if (strcmp(var_type.c_str(), "wall") == 0) // specific wall_type case
				get_next_substr_between_sep(args_control, line, sep, curr_sep_pos);

			else if (strcmp(var_type.c_str(), "Initial") == 0) // specific initial_condition case
				ic_preset = IC_type_map[var_value];

			else if (strcmp(var_type.c_str(), "VISC_types") == 0) // specific viscosity case
				viscosity = VISC_type_map[var_value];
            	
			else if (strcmp(var_type.c_str(), "Reconstruction") == 0) // specific viscosity case
				reconstruction = Reconstruction_map[var_value];

            else if (strcmp(var_type.c_str(), "test") == 0) // specific test_preset case
				test = test_preset_type_map[var_value];

            else if (strcmp(var_type.c_str(), "string") == 0)
				*(std::string*)vars[var_name] = var_value;

            else {
				std::cerr << "Variable type is not identified: " << var_name << " " << var_type << std::endl;
				exit(1);
			}

            if (strcmp(args_control.c_str(), "{") == 0) { // assign wall values if provided
				std::string opening, arg_name, arg_value;
				for (std::getline(fin, line);; std::getline(fin, line)) {
					curr_sep_pos = 0;
					get_next_substr_between_sep(opening, line, sep, curr_sep_pos);
					if (strcmp(opening.c_str(), "}") == 0) break;
					get_next_substr_between_sep(arg_name, line, sep, curr_sep_pos);
					get_next_substr_between_sep(arg_value, line, sep, curr_sep_pos);
					if (strcmp(var_type.c_str(), "wall") == 0) {
						if (arg_name == "type")
							walls[std::stoul(var_value)].type = wall_type_map[arg_value];
						
						else if (arg_name == "T")
							walls[std::stoul(var_value)].T = std::stod(arg_value);
						
						else if (arg_name == "v_x")
							walls[std::stoul(var_value)].v_x = std::stod(arg_value);
						
						else {
							std::cerr << "Wall " << var_value << " arg type not found (" << arg_name << ")" << std::endl;
							exit(1);
						}
						std::cout << "Wall " << var_value << " : " << arg_name << "," << arg_value << std::endl;
					}
				}
			}
			std::cout << var_type << " " << var_name << " " << var_value << std::endl;
		}
		dx = (x_end - x_start) / nx;
		nx_all = nx + 2 * nx_fict;
		Cv = R / (gamma - 1.0);
		Cp = Cv + R;
		
		if (test != TEST_CUSTOM) {
			ic_preset = IC_preset(test);
			for (unsigned int i = 0; i < 2; ++i) {
				walls[i].bc_preset = BC_preset(test); // it's easier to specify further boundary values in solver,
										  // than to add many elif cases here //
										  // assigned boundary values don't matter
										  // i'm sorry
			};
		};
	}
	else {
		std::cerr << "Can't open the file: " << file_name << std::endl;
		exit(1);
	}
	fin.close();
};

Parameters::Parameters(const Parameters& rhs):
	x_start(rhs.x_start),
	x_end(rhs.x_end),
	nx(rhs.nx),
	dx(rhs.dx),
	nx_fict(rhs.nx_fict),
	nx_all(rhs.nx_all),
	CFL(rhs.CFL),
	nt(rhs.nt),
	test(rhs.test),
	ic_preset(rhs.ic_preset),
	viscosity(rhs.viscosity),
	mu0(rhs.mu0),	
	reconstruction(rhs.reconstruction),
	is_conservative(rhs.is_conservative),
	gamma(rhs.gamma),
	Cv(rhs.Cv),
	Cp(rhs.Cp),
	nt_write(rhs.nt_write),
	write_file(rhs.write_file)
{
	std::copy(std::begin(rhs.walls), std::end(rhs.walls), std::begin(walls));
};
