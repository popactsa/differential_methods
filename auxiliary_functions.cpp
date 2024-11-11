#include "auxiliary_functions.h"

void get_next_substr_between_sep (std::string& substr, const std::string& str, const std::string sep, unsigned int& curr) {
	try {
		if (curr == str.length()) {
			throw "Separator can't be found. Already at the end of the string";
		}
		size_t found_temp = str.find(sep, curr);
		unsigned int found = found_temp != std::string::npos ? found_temp : str.length(); // '\n' at the end popped
//		std::cout << "curr = " << curr << std::endl;
//		std::cout << "found = " << found << std::endl;
//		std::cout << str.length() << std::endl;
		substr.assign(str, curr, found - curr);
//		std::cout << substr << std::endl << std::endl;
		curr = found + 1;
	}
	catch (const char* error_message) {
		std::cerr << error_message << std::endl;
		exit(1);
	}
};
