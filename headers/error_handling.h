#ifndef ERROR_HANDLING_H
#define ERROR_HANDLING_H

enum class Error_action {ignore, throwing, terminating, logging};
constexpr Error_action default_Error_action = Error_action::throwing;

enum class Error_code {range_error, length_error};

string error_code_name[] {"range error", "length error"};

template<Error_action action = default_Error_action, class C>
constexpr void expect(C _cond, Error_code _ec)
{
	if constexpr (action == Error_action::logging)
		if (!_cond()) std::cerr << "expect() failure: " << int(_ec) << " " << error_code_name[int _ec] << std::endl;
	if constexpr (action == Error_action::throwing)
		if (!_cond()) throw _ec;
	if constexpr (action == Error_action::terminating)
		if (!_cond()) terminate();
	// no action for action == Error_action::ignore
}

#endif
