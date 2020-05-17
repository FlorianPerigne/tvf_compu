#include "NPDCall.hpp"


//constructor for log exponential call payoff
NPDCall::NPDCall(double strike, std::vector<double> weights, bool weight_bool)
	:NonPathDependent(strike), NPD_weights(weights)
{
	NPD_weight_power = weight_bool;
}

//constructor for basket call
NPDCall::NPDCall(double strike, std::vector<double> weights)
	:NonPathDependent(strike), NPD_weights(weights)
{
	NPD_weight_power = false;
}

//constructor for call
NPDCall::NPDCall(double strike)
	:NonPathDependent(strike)
{
	NPD_weights.push_back(1);
	NPD_weight_power = false;
}

std::vector<double> NPDCall::get_weights()
{
	return NPD_weights;
}
	
	
std::vector<double> NPDCall::operator()(const std::vector<std::vector<std::vector<double>>>& x) const
{
	
	std::vector<double> Values;
	Values.resize(x.size());
	
	for (llong i = 0; i < x.size(); i++)
	{

		
		double m =  0;
		if (NPD_weight_power){
			m = 1;
		}
		
		for (llong j = 0; j<NPD_weights.size(); j++)
		{
			if (NPD_weight_power){
				m *= std::pow(x[i][j][x[0][j].size()-1],  NPD_weights[j]);
			}
			else{
				m+= NPD_weights[j]*x[i][j][x[0][j].size()-1];
			}
		}
		
		Values[i] = (m > NPD_strike)? m - NPD_strike : 0;
	}

	return Values;

}