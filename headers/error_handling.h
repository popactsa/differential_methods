#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

#include <exception>
#include <string>

////////Tour of C++ 2022 B.Stroustrup p.49/////////

enum class Error_action
{
	ignore,
	throwing,
	terminating,
	logging
};

constexpr Error_action default_Error_action = Error_action::throwing;

enum class Error_code
{
	negative_value,
	missing,
	bad_syntax,
	incorrect_order
};

std::string Error_code_name[]
{
	"negative value",
	"missing",
	"bad syntax",
	"incorrect order"
};

template<Error_action action = default_Error_action, class C>
constexpr void expect(C cond, Error_code x)
{
	if constexpr (action == Error_action::throwing)
		if (!cond()) throw x;
	if constexpr (action == Error_action::terminating)
		if (!cond()) std::terminate();
	if constexpr (action == Error_action::logging)
		if (!cond()) std::cerr << "expect() failure: " << int(x) << ' ' << Error_code_name[int(x)] << std::endl;
	// ignore -> nothing happens
}

#endif
