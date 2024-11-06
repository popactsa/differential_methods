#include <iostream>
#include <string>

enum E{test, test2};

int main() {
	E e;
	unsigned int a;
	unsigned short sh;
	float f;
	double dl;
	std::string str;
	std::cout << typeid(a).name() << std::endl;
	std::cout << typeid(sh).name() << std::endl;
	std::cout << typeid(f).name() << std::endl;
	std::cout << typeid(dl).name() << std::endl;
	std::cout << typeid(str).name() << std::endl;
	std::cout << typeid(e).name() << std::endl;
	return 0;
}
