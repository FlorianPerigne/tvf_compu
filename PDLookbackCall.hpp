#ifndef PDLOOKBACKCALL_HPP
#define PDLOOKBACKCALL_HPP

#include "PDLookback.hpp"

class PDLookbackCall :public PDLookback
{
	
public:
	PDLookbackCall(int days, double strike = 100);
	
	virtual std::vector<std::vector<double>> operator()(std::vector<std::vector<std::vector<double>>> x);

};

#endif