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
        {"CFL", &CFL},
        {"t_end", &t_end},
        {"nt_write", &nt_write},
        {"write_file", &write_file}
    };
    if (fin.is_open()) {
        std::string line;
        std::string sep = "\t";
        std::getline(fin, line);
        for (; !fin.eof(); std::getline(fin, line)) {
            std::string var_name = line.substr(0, line.find(sep));
            std::string value = line.substr(line.find(sep), '\0');
            const char* val_type = typeid(vars[var_name]).name();
            if (val_type == "double") {
                *(double*)vars[var_name] = std::stof(value);
            }
            else if (val_type == "int") {

            }
        }
    }
    fin.close();
}