#include <iostream>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <exception>
#include <chrono>
#include <vector>
#include <string>

#include "LinearCongruential.hpp"
#include "EcuyerCombined.hpp"

#include "Normal.hpp"
#include "NormalMultiVariate.hpp"

#include "Stocks.hpp"
#include "StocksTerminal.hpp"
#include "StocksFullPath.hpp"

#include "Payoff.hpp"
#include "NonPathDependent.hpp"
#include "NPDCall.hpp"
#include "NPDPut.hpp"
#include "NPDBasketCall.hpp"
#include "NPDBasketPut.hpp"

#include "MonteCarloEuropean.hpp"

int main()
{
	try
	{
		UniformGenerator* ecuyer = new EcuyerCombined(1, 1);

		double mu = 3. / 100.;
		std::vector<std::vector<double>> cov(2, std::vector<double>(2));

		cov[0][0] = 0.2*0.2;
		cov[1][1] = 0.3*0.2;
		cov[0][1] = 0.6*std::sqrt(cov[0][0])*std::sqrt(cov[1][1]);
		cov[1][0] = cov[0][1];

		ContinuousGenerator* biv_norm = new NormalMultiVariate(ecuyer, { mu }, cov);

		Stocks* stocks = new StocksTerminal(biv_norm, 100, mu, 1);
		std::vector<std::vector<std::vector<double>>> S(stocks->Generate(5));

		for (int i = 0; i < S.size(); i++)
		{
			for (int j = 0; j < S[i].size(); j++)
			{
				std::cout << S[i][j][0] << ", ";
			}
			std::cout << std::endl;
		}

		std::cout << "----- FULL STOCK PATH -----" << std::endl;

		llong n_steps = 100;
		Stocks* stocksFull = new StocksFullPath(biv_norm, 100, { mu }, 1, n_steps);

		std::vector<std::vector<std::vector<double>>> S_full(stocksFull->Generate(5));

		for (int i = 0; i < S_full.size(); i++)
		{
			for (int j = 0; j < S_full[i].size(); j++)
			{
				std::cout << S_full[i][j][n_steps-1] << ", ";
			}
			std::cout << std::endl;
		}

		std::cout<<"----- Test Payoff ------" << std::endl;

		Payoff* V1 = new NPDPut(100);
		std::vector<std::vector<double>> Call = (*V1)(S);
		for (int i = 0; i < S.size(); i++)
		{
			for (int j = 0; j < S[i].size(); j++)
			{
				std::cout << Call[i][j] << ", ";
			}
			std::cout << std::endl;
		}
		
		std::cout<<"----- Test Basket Payoff ------" << std::endl;

		std::vector<double> weights = {0.5, 0.5};
		Payoff* V2 = new NPDBasketCall(100, weights);

		std::vector<std::vector<double>> CallBkt = (*V2)(S);

		for (int i = 0; i < S.size(); i++)
		{
			std::cout << CallBkt[i][0] << ", ";
			std::cout << std::endl;
		}

		std::cout << "----- MONTECARLO ------" << std::endl;

		NonPathDependent* call_payoff = new NPDBasketCall(100, weights);
		StocksTerminal* stocksT = new StocksTerminal(biv_norm, 100, { mu }, 1);
		MonteCarlo* mc_solver = new MonteCarloEuropean(stocksT, call_payoff, 100000);

		mc_solver->Solve();
		std::cout << "Price = " << mc_solver->get_price() << std::endl;
	}
	catch (std::exception & e)
	{
		std::cout << e.what() << std::endl;
	}
	
	
	return 0;
}
