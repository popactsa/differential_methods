#include <iostream>
#include "Parameters.h"
#include "auxiliary_functions.h"

int main() {
	std::string file_name = "Parameters_file.txt";
	const Parameters par(file_name);
	par.print_write_file();
	return 0;
}
