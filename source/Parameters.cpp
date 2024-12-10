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
	std::map<std::string, BoundaryType> boundary_type_map = {
		{"WALL", B_WALL},
		{"FLUX", B_FLUX},
	};
	std::map<std::string, int> boundary_args = {
		{"BoundaryType", 0},
		{"T", 1},
		{"v_x", 2},
	};
	std::map<std::string, InitialPreset> ic_type_map = {
		{"TEST1", IC_TEST1},
		{"TEST2", IC_TEST2},
		{"TEST3", IC_TEST3},
		{"TEST4", IC_TEST4},
	};
	std::map<std::string, ArtificialViscosity> visc_type_map = {
		{"NONE", V_NONE},
		{"NEUMAN", V_NEUMAN},
		{"LATTER", V_LATTER},
		{"LINEAR", V_LINEAR},
		{"SUM", V_SUM},
	};
	std::map<std::string, Reconstruction> reconstruction_map = {
		{"GODUNOV", R_GODUNOV},
		{"KOLGAN72", R_KOLGAN72},
		{"KOLGAN75", R_KOLGAN75},
		{"OSHER84", R_OSHER84},
	};
	std::map<std::string, TimeAlgo> time_algo_map = {
		{"EULER", T_EULER},
		{"RUNGEKUTTA", T_RUNGEKUTTA},
	};
	std::map<std::string, FluxScheme> flux_scheme_map = {
		{"GODUNOV", F_GODUNOV},
		{"RLF", F_RLF},
	};
	
	InitialPreset ic_test = IC_CUSTOM;
	
	for (int i = 0; i < 2; ++i) boundaries[i].b_preset = B_CUSTOM;
	
	if (fin.is_open()) {
		std::string line;
		const std::string sep = "\t";
		for (std::getline(fin, line); !fin.eof(); std::getline(fin, line)) {
			int curr_sep_pos = 0;
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

			else if (strcmp(var_type.c_str(), "int") == 0)
			    *(int*)vars[var_name] = std::stoi(var_value);

	        	else if (strcmp(var_type.c_str(), "Boundary") == 0) // specific wall_type case
				get_next_substr_between_sep(args_control, line, sep, curr_sep_pos);

			else if (strcmp(var_type.c_str(), "InitialPreset") == 0) // specific initial_condition case
				ic_preset = ic_type_map[var_value];

			else if (strcmp(var_type.c_str(), "ArtificialViscosity") == 0) // specific viscosity cases
				viscosity_type = visc_type_map[var_value];

			else if (strcmp(var_type.c_str(), "Reconstruction") == 0) // specific reconstruction case
				reconstruction_type = reconstruction_map[var_value];

			else if (strcmp(var_type.c_str(), "TimeAlgo") == 0) // specific time algorithm case
				time_algo = time_algo_map[var_value];

			else if (strcmp(var_type.c_str(), "FluxScheme") == 0) // specific flux calculation case
				flux_scheme = flux_scheme_map[var_value];

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
						if (strcmp(var_type.c_str(), "Boundary") == 0) {
							if (arg_name == "type")
								boundaries[std::stoul(var_value)].b_type = boundary_type_map[arg_value];
							
							else if (arg_name == "T")
								boundaries[std::stoul(var_value)].T = std::stod(arg_value);
							
							else if (arg_name == "v_x")
								boundaries[std::stoul(var_value)].v_x = std::stod(arg_value);
							
							else {
								std::cerr << "Boundary " << var_value << " arg type not found (" << arg_name << ")" << std::endl;
								exit(1);
							}
							std::cout << "Boundary " << var_value << " : " << arg_name << "," << arg_value << std::endl;
						}
					}
				}
				std::cout << var_type << " " << var_name << " " << var_value << std::endl;
			}
			dx = (x_end - x_start) / nx;
			nx_all = nx + 2 * nx_fict;
			Cv = R / (gamma - 1.0);
			Cp = Cv + R;
			
			if (ic_test != IC_CUSTOM) {
				ic_preset = InitialPreset(ic_test);
				for (int i = 0; i < 2; ++i) {
					boundaries[i].b_preset = BoundaryPreset(ic_test); // it's easier to specify further boundary values in solver,
											  // than to add many elif cases here //
											  // assigned boundary values don't matter
											  // i'm sorry
			}
		}
	}
	else {
		std::cerr << "Can't open the file: " << file_name << std::endl;
		exit(1);
	}
	fin.close();
};
