#include "StocksStandardFullPath.hpp"
#include <cmath>
#include <iostream>
StocksStandardFullPath::StocksStandardFullPath(ContinuousGenerator* gen, double s0, double mu, double maturity, llong n_steps)
	: StocksFullPath(gen, s0, mu, maturity, n_steps)
{
}

StocksStandardFullPath::StocksStandardFullPath(ContinuousGenerator* gen, double s0, double mu, double maturity, std::vector<double> dividends, std::vector<double> date_dividends, llong n_steps)
	: StocksFullPath(gen, s0, mu, maturity, dividends, date_dividends, n_steps)
{
}

std::vector<std::vector<std::vector<double>>> StocksStandardFullPath::Generate(llong n_sims)
{

	W.resize(n_sims, std::vector<double>(Gen->get_covariance_matrix().size()));
	S.resize(n_sims, std::vector<std::vector<double>>(W[0].size(), std::vector<double>(N_steps)));

	Gen->Generate(W, n_sims);

	double dt = T / (N_steps - 1);
	llong cpt_div = 0;
	double divVal = 0;

	for (llong i = 0; i < n_sims; ++i)
	{
		for (llong j = 0; j < W[0].size(); ++j)
		{
			S[i][j][0] = S0;
		}
	}

	for (llong z = 1; z < N_steps; ++z)
	{
		if (z * dt < Date_Div[cpt_div] && Date_Div[cpt_div] <= (z + 1) * dt)
		{
			divVal = Div[cpt_div];
			cpt_div += 1;

			if (cpt_div == Div.size())
			{
				cpt_div -= 1;
			}
		}
		Gen->Generate(W, n_sims);
		for (llong i = 0; i < n_sims; ++i)
		{
			for (llong j = 0; j < W[0].size(); j++)
			{
				S[i][j][z] = S[i][j][z - 1] + S[i][j][z - 1]*((Mu) * dt + std::sqrt(dt) * W[i][j]) - divVal;
			}
		}
		divVal = 0;
	}
	
	return S;
}
