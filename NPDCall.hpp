#ifndef NPDCALL_HPP
#define NPDCALL_HPP

#include "NonPathDependent.hpp"

class NPDCall :public NonPathDependent
{
	
public:
	NPDCall(double strike, std::vector<double> weights, bool weight_bool);
	NPDCall(double strike, std::vector<double> weights);
	NPDCall(double strike);
	
	std::vector<double> get_weights();
	std::vector<double> operator()(const std::vector<std::vector<std::vector<double>>>& x) const;

private:
	std::vector<double> NPD_weights;
	bool NPD_weight_power;
};

#endif