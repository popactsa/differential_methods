#include "Parameters.h"

Parameters::Parameters(std::string file_name) {
	std::ifstream fin(file_name);
	std::map<const char*, void*> vars = {
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
		{"CFL", &CFL},
		{"t_end", &t_end},
		{"nt_write", &nt_write},
		{"write_file", &write_file}
	};
	if (fin.is_open()) {
		std::string line;
		const std::string sep = "\t";
		std::getline(fin, line);
		for (; !fin.eof(); std::getline(fin, line)) {
			unsigned int curr_sep = 0;
			const char* var_type = get_next_substr_between_sep(line, sep, curr_sep);
			const char* var_name = get_next_substr_between_sep(line, sep, curr_sep);
			const char* var_value = get_next_substr_between_sep(line, sep, curr_sep);
			if (var_type = "float") {
				*(float*)vars[var_name] = std::stof(var_value);
		    	}
		    	else if (var_type == "unsignedint") {
		    		*(unsigned int*)vars[var_name] = std::stoi(var_value);
		    	}
		    	else if (var_type == "E_BOUNDARY_TYPE") {
		    		*(E_BOUNDARY_TYPE*)vars[var_name] = E_BOUNDARY_TYPE(std::stoi(var_value));
		    	}
			else {
				std::cerr << "Variable type is not indentified" << std::endl;
				exit(1);
			}
		}
	}
	else {
		std::cerr << "Can't open a file: " << file_name << std::endl;
		exit(1);
	}
	fin.close();
}
