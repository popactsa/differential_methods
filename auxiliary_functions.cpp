#include "auxiliary_functions.h"

const char* get_next_substr_between_sep (const std::string str, const std::string sep, unsigned int& curr) {
	try {
		if (curr == std::string::npos) {
			throw "Separator can't be found. Already at the end of the string";
		}
		unsigned int found = str.find(sep, curr);
		std::cout << "curr = " << curr << std::endl;
		std::cout << "found = " << found << std::endl;
		const char* res = (str.substr(curr, found - curr)).c_str();
		std::cout << res << std::endl << std::endl;
		curr = found + 1;
		return res;
	}
	catch (const char* error_message) {
		std::cerr << error_message << std::endl;
		exit(1);
	}
};

