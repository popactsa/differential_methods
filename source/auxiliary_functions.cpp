#include "auxiliary_functions.h"

#ifdef WIN32
#define l 0
#else
#define l 1
#endif // WIN32

void get_next_substr_between_sep (std::string& substr, const std::string& str, const std::string sep, unsigned int& curr) {
	try {
		size_t found_temp = str.find(sep, curr);
		unsigned int found = found_temp != std::string::npos ? found_temp : str.length() - l; // '\n' at the end popped
		substr.assign(str, curr, found - curr);
		curr = found + 1;
	}
	catch (const char* error_message) {
		std::cerr << error_message << std::endl;
		exit(1);
	}
};
