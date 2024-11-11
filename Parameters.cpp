#include "Parameters.h"

Parameters::Parameters(std::string file_name) {
	std::ifstream fin(file_name);
	std::map<std::string, void*> vars = {
		{"g_x", &g_x},
		{"g_y", &g_y},
		{"g_z", &g_z},
		{"x_start", &x_start},
		{"x_end", &x_end},
		{"y_start", &y_start},
		{"y_end", &y_end},
		{"z_start", &z_start},
		{"z_end", &z_end},
		{"nx", &nx},
		{"ny", &ny},
		{"nz", &nz},
		{"n_walls", &n_walls},
		{"CFL", &CFL},
		{"t_end", &t_end},
		{"nt_write", &nt_write},
		{"write_file", &write_file}
	};
	if (fin.is_open()) {
		std::string line;
		const std::string sep = "\t";
		for (std::getline(fin, line); !fin.eof(); std::getline(fin, line)) {
			unsigned int curr_sep_pos = 0;
			std::string var_type, var_name, var_value;
			get_next_substr_between_sep(var_type, line, sep, curr_sep_pos);
			get_next_substr_between_sep(var_name, line, sep, curr_sep_pos);
			get_next_substr_between_sep(var_value, line, sep, curr_sep_pos);
			if (strcmp(var_type.c_str(), "double") == 0) {
				*(double*)vars[var_name] = std::stof(var_value);
		    	}
		    	else if (strcmp(var_type.c_str(), "unsignedint") == 0) {
		    		*(unsigned int*)vars[var_name] = std::stoi(var_value);
				if (var_name == "n_walls") {
					wall_type.resize(n_walls);
				}
		    	}
		    	else if (strcmp(var_type.c_str(), "E_BOUNDARY_TYPE") == 0) { // specific wall_type case
				std::string number;
				get_next_substr_between_sep(number, line, sep, curr_sep_pos);
				if (wall_type.size() > std::stoi(number)) {
					wall_type.at(std::stoi(number)) = wall_type_map[var_value];
				}
				else {
					wall_type.push_back(wall_type_map[var_value]);
					std::cout << "Wall " << number << " type set to " << var_value << std::endl;
				}
		    	}
			else if (strcmp(var_type.c_str(), "string") == 0) {
				*(std::string*)vars[var_name] = var_value;
			}
			else {
				std::cerr << "Variable type is not identified: " << var_name << " " << var_type << std::endl;
				exit(1);
			}
//			std::cout << var_type << " " << var_name << " " << var_value << std::endl;
		}
	}
	else {
		std::cerr << "Can't open a file: " << file_name << std::endl;
		exit(1);
	}
	fin.close();
};
