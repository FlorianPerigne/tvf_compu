#include "PDCall.hpp"

//constructor of the path dependent basket call
PDCall::PDCall(double strike, std::vector<double> weights)
	:PathDependent(strike), PD_weights(weights)
{
	
}

//constructo of the vanilla call 
PDCall::PDCall(double strike)
	:PathDependent(strike)
{
	PD_weights.push_back(1);
}

std::vector<double> PDCall::get_weights()
{
	return PD_weights;
}
	
//reutrn the values of the call using the last step of time to evalue the payoff
std::vector<double> PDCall::operator()(const std::vector<std::vector<std::vector<double>>>& x) const
{
	
	std::vector<double> Values;
	Values.resize(x.size());
	
	for (llong i = 0; i < x.size(); i++)
	{

		double m =  0;
		
		for (llong j = 0; j<PD_weights.size(); j++)
		{
			m+= PD_weights[j]*x[i][j][x[0][j].size()-1];
		}
		
		Values[i] = (m > PD_strike)? m - PD_strike : 0;
	}

	return Values;

}

//return the values of the call if we exercise at a specific date of time
std::vector<double> PDCall::operator()(const std::vector<std::vector<std::vector<double>>>& x, llong k) const
{
	
	std::vector<double> Values;
	Values.resize(x.size());
	
	for (llong i = 0; i < x.size(); i++)
	{
		double m =  0;
		
		for (llong j = 0; j<PD_weights.size(); j++)
		{
			m+= PD_weights[j]*x[i][j][k];
		}
		
		Values[i] = (m > PD_strike)? m - PD_strike : 0;
	}

	return Values;

}




