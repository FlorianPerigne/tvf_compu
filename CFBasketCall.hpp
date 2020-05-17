#ifndef CFBASKETCALL_HPP
#define CFBASKETCALL_HPP

#include "ClosedForm.hpp"
#include "Tools.hpp"

class CFBasketCall :public ClosedForm
{
	
public:
	CFBasketCall(double strike, const std::vector<double>& weights, const std::vector<double>& S0, const std::vector<std::vector<double>>& Cov);
	
	std::vector<double> get_weights();
	std::vector<double> operator()(const std::vector<std::vector<std::vector<double>>>& x) const;
	double operator()(const double& S0, const double& mu, const double& TTM, const double& vol) const;
	
private:
	std::vector<double> CF_weights;
	std::vector<double> CF_S0;
	std::vector<std::vector<double>> CF_Cov;
};

#endif