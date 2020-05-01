#ifndef PDLOOCKBACK_HPP
#define PDLOOCKBACK_HPP

#include "PathDependent.hpp"
class PDLookback : public PathDependent
{
	
public:
	PDLookback(int days);
	
protected:
	int PD_days;
};

#endif