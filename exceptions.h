#pragma once

#include <exception>
#include <string>

class exception_setofequations : public std::exception
{
	std::string msg;

public:
	exception_setofequations(const std::string& s) :msg(s) {}

	const char* what()const throw()
	{
		return msg.c_str();
	}
};

