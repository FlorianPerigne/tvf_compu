#include "CFBasketCall.hpp"
#include "Tools.hpp"

CFBasketCall::CFBasketCall(double strike, const std::vector<double>& weights, const std::vector<double>& S0, const std::vector<std::vector<double>>& Cov)
	:ClosedForm(strike), CF_weights(weights), CF_Cov(Cov), CF_S0(S0)
{
	
}

std::vector<double> CFBasketCall::get_weights()
{
	return CF_weights;
}

double CFBasketCall::operator()(const double& S0, const double& mu, const double& TTM, const double& vol) const
{
	double X = 1.0;
	double V = 0;
	double var = 0;
	std::vector<double> W;
	
	mult_matrix_vect(CF_Cov, CF_weights,W);
	
	for (int i = 0; i<CF_S0.size(); ++i)
	{
		X *= std::pow(CF_S0[i], CF_weights[i]);
		V += CF_Cov[i][i] * CF_weights[i];
		for (int j = 0; j<=i; ++j)
		{
			V += CF_Cov[i][j]*CF_Cov[i][j] * CF_weights[i];
		}
		var += W[i] * CF_weights[i];
	}
	
	
	double d1 = (std::log(X/CF_strike) + (mu + 0.5*var - 0.5*V)*TTM)/(std::sqrt(var * TTM));
	double d2 = d1 - (std::sqrt(var * TTM));
	
	return (X*normalCDF(d1) - CF_strike*std::exp(-mu*TTM)*normalCDF(d2));
}

std::vector<double> CFBasketCall::operator()(const std::vector<std::vector<std::vector<double>>>& x) const
{
	throw std::runtime_error("Error, operator not applicable for closed form formula");

}