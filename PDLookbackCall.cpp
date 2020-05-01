#include "PDLookbackCall.hpp"

PDLookbackCall::PDLookbackCall(int days, double strike)
	:PDLookback(strike, days)
{
}


std::vector<std::vector<double>> PDLookbackCall::operator()(std::vector<std::vector<std::vector<double>>> x)
{
	std::vector<std::vector<double>> Values;
	Values.resize(x.size(), std::vector<double>(x[0].size())); 
	
	for (int i = 0; i < x.size(); i++)
	{
		for (int j = 0; j < x[i].size(); j++)
		{
			double new_strike = 0;
			for (int k = 0; k<std::min(PD_days,(int)x[i][j].size()); k++)
			{
				new_strike += x[i][j][k];
			}
			PD_strike = new_strike/std::min(PD_days,(int)x[i][j].size());
			Values[i][j] = (x[i][j][x[i][j].size()-1]  > PD_strike)? x[i][j][x[i][j].size()-1] - PD_strike : 0;
		}
	}
		
	return Values;
}