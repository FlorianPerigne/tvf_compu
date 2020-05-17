#include "CFPut.hpp"

//constructor for Closed form solution of a put with a particular leverage 
CFPut::CFPut(double strike, std::vector<double> weights)
	:ClosedForm(strike), CF_weights(weights)
{
	
}
//constructor for closed form of a put
CFPut::CFPut(double strike)
	:ClosedForm(strike)
{
	CF_weights.push_back(1);
}

std::vector<double> CFPut::get_weights()
{
	return CF_weights;
}

//get the closed form value of a put
double CFPut::operator()(const double& S0, const double& mu, const double& TTM, const double& vol) const
{
	double d1 = (std::log(S0/CF_strike) + (mu + 0.5*vol*vol)*TTM)/(vol*std::sqrt(TTM));
	double d2 = d1 - (vol*std::sqrt(TTM));
	
	return CF_weights[0] * (CF_strike*std::exp(-mu*TTM)*normalCDF(-d2) - S0*normalCDF(-d1)) ;
}
	
std::vector<double> CFPut::operator()(const std::vector<std::vector<std::vector<double>>>& x) const
{
	throw std::runtime_error("Error, operator not applicable for closed form formula");

}