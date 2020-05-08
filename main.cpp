#include <iostream>
#include <iomanip>
#include <algorithm>
#include <numeric>
#include <fstream>
#include <exception>
#include <chrono>
#include <vector>
#include <string>

#include "LinearCongruential.hpp"
#include "EcuyerCombined.hpp"
#include "VanDerCorput.hpp"

#include "Normal.hpp"
#include "NormalMultiVariate.hpp"

#include "StandardAntithetic.hpp"

#include "StocksAntitheticTerminal.hpp"
#include "StocksStandardTerminal.hpp"
#include "StocksAntitheticFullPath.hpp"
#include "StocksStandardFullPath.hpp"

#include "Payoff.hpp"
#include "NonPathDependent.hpp"
#include "NPDBasketCall.hpp"
#include "NPDBasketPut.hpp"
#include "PathDependent.hpp"
#include "PDBasketCall.hpp"
#include "CFCall.hpp"
#include "CFPut.hpp"
#include "NPDCall.hpp"
#include "PDCall.hpp"

#include "MonteCarloEuropean.hpp"
#include "MonteCarloLSM.hpp"

#include "BasisLaguerre.hpp"

#include "Simulation.hpp"

#include "Tools.hpp"

// void display_mat(const std::vector<std::vector<double>>& A);

int main()
{
	try
	{
		UniformGenerator* ecuyer = new EcuyerCombined(1, 1);
		UniformGenerator* ecuyer2 = new EcuyerCombined(1, 1);

		UniformGenerator* vdc = new VanDerCorput(1, 2);
		UniformGenerator* vdc2 = new VanDerCorput(1, 3);

		double mu = 0./ 100.;
		std::vector<double> weights = { 0.5, 0.5};

		std::vector<std::vector<double>> cov(2, std::vector<double>(2));

		cov[0][0] = 0.2*0.2;
		cov[1][1] = 0.2*0.2;
		cov[0][1] = 0.9999*std::sqrt(cov[0][0])*std::sqrt(cov[1][1]);
		cov[1][0] = cov[0][1];

		ContinuousGenerator* biv_norm = new NormalMultiVariate(ecuyer, 0, cov);
		ContinuousGenerator* biv_norm2 = new NormalMultiVariate(ecuyer2, 0, cov);
		ContinuousGenerator* biv_norm3 = new NormalMultiVariate(vdc, vdc2, 0, cov);

		std::cout << "----- MONTECARLO ------" << std::endl;

		llong n_simu = 5000;

		NonPathDependent* call_payoff = new NPDBasketCall(100, weights);
		NonPathDependent* put_payoff = new NPDBasketPut(100, weights);
		PathDependent* call_payoff_PD = new PDBasketCall(100, weights);

		R3R1Function* antithetic_function = new StandardAntithetic();

		ClosedForm* call_payoff_CF = new CFCall(100);
		double prix_bs = (*call_payoff_CF)(100, 0 ,1,0.2);
		double prix_bs_div = (*call_payoff_CF)(93, 0 ,1,0.2);
		NonPathDependent* Call_clasic = new NPDCall(100);
		PathDependent* Call_clasic_PD = new PDCall(100);
		
		// --- Standard MC Terminal

		StocksTerminal* stocksT = new StocksStandardTerminal(biv_norm, 100, mu, 1);
		MonteCarlo* mc_solver = new MonteCarloEuropean(stocksT, call_payoff, n_simu);

		// --- MC Terminal with Antitethic variance reduction

		StocksTerminal* stocksT2 = new StocksAntitheticTerminal(biv_norm, 100, mu, 1, antithetic_function);
		MonteCarlo* mc_solver_antithetic = new MonteCarloEuropean(stocksT2, call_payoff, llong(n_simu/2));

		// --- MC Terminal with Quasi random numbers

		StocksTerminal* stocksT3 = new StocksStandardTerminal(biv_norm3, 100, mu, 1);
		MonteCarlo* mc_solver_quasi = new MonteCarloEuropean(stocksT3, call_payoff, n_simu);

		// --- MC Terminal with Control Variate variance reduction

		StocksTerminal* stocksT4 = new StocksStandardTerminal(biv_norm, 100, mu, 1);
		MonteCarlo* mc_solver_CV = new MonteCarloEuropean(stocksT4, call_payoff, llong(n_simu/2), Call_clasic, prix_bs);
		
		// --- MC Terminal with antithetic and quasi random numbers variance reduction

		StocksTerminal* stocksT5 = new StocksAntitheticTerminal(biv_norm3, 100, mu, 1, antithetic_function);
		MonteCarlo* mc_solver_quasi_anti = new MonteCarloEuropean(stocksT5, call_payoff, llong(n_simu/2));
		
		// --- MC Terminal with Antithetic and Control Variate variance reduction

		StocksTerminal* stocksT6 = new StocksAntitheticTerminal(biv_norm, 100, mu, 1, antithetic_function);
		MonteCarlo* mc_solver_anti_CV = new MonteCarloEuropean(stocksT6, call_payoff, llong(n_simu/2), Call_clasic, prix_bs);
		
		// --- MC Terminal  with Quasi random numbers and Control Variate variance reduction

		StocksTerminal* stocksT7 = new StocksStandardTerminal(biv_norm3, 100, mu, 1);
		MonteCarlo* mc_solver_quasi_CV = new MonteCarloEuropean(stocksT7, call_payoff, llong(n_simu/2), Call_clasic, prix_bs);

		// --- MC Terminal  with Quasi random numbers anthitetic and Control Variate variance reduction

		StocksTerminal* stocksT8 = new StocksAntitheticTerminal(biv_norm3, 100, mu, 1, antithetic_function);
		MonteCarlo* mc_simul_quasi_anti_CV = new MonteCarloEuropean(stocksT8, call_payoff, llong(n_simu/2), Call_clasic, prix_bs);


		// ------ MC Simulation ------

		Simulation* MC_simul_standard = new Simulation(mc_solver);
		Simulation* MC_simul_quasi = new Simulation(mc_solver_quasi);
		Simulation* MC_simul_standard_anti = new Simulation(mc_solver_antithetic);
		Simulation* MC_simul_CV = new Simulation(mc_solver_CV);
		Simulation* MC_simul_quasi_anti = new Simulation(mc_solver_quasi_anti);
		Simulation* MC_simul_anti_CV = new Simulation(mc_solver_anti_CV);
		Simulation* MC_simul_quasi_CV = new Simulation(mc_solver_quasi_CV);
		Simulation* MC_simul_quasi_anti_CV = new Simulation(mc_simul_quasi_anti_CV);

		llong n_sims = 100;

		std::cout << "EXPECTATION AND VARIANCE IN FUNCTION OF N_SIMULATION of paths" << std::endl;

		//MC_simul_standard->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000}, "var_MC_simul_standard.csv");
		//MC_simul_quasi->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_quasi.csv");
		//MC_simul_standard_anti->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_standard_anti.csv");
		//MC_simul_quasi_anti->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_quasi_anti.csv");
		//MC_simul_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_CV.csv");
		//MC_simul_anti_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_anti_CV.csv");
		//MC_simul_quasi_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_quasi_CV.csv");
		//MC_simul_quasi_anti_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "var_MC_simul_quasi_anti_CV.csv");
		

		clock_t start, end;

		start = clock();
		MC_simul_standard->expectation_by_sims(200, { 10000 }, "exp_MC_simul_standard.csv");
		end = clock();
		
		std::cout << "Standard : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_quasi->expectation_by_sims(200, { 10000 }, "exp_MC_simul_quasi.csv");
		end = clock();

		std::cout << "Quasi : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_standard_anti->expectation_by_sims(200, { 10000 }, "exp_MC_simul_standard_anti.csv");
		end = clock();

		std::cout << "Standard Antithetic : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_quasi_anti->expectation_by_sims(200, { 10000 }, "exp_MC_simul_quasi_anti.csv");
		end = clock();

		std::cout << "Quasi Antithetic : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_CV->expectation_by_sims(200, { 10000 }, "exp_MC_simul_CV.csv");
		end = clock();

		std::cout << "Standard CV : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_anti_CV->expectation_by_sims(200, { 10000 }, "exp_MC_simul_anti_CV.csv");
		end = clock();

		std::cout << "Antithetic CV : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_quasi_CV->expectation_by_sims(200, { 10000 }, "exp_MC_simul_quasi_CV.csv");
		end = clock();

		std::cout << "Quasi CV : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

		start = clock();
		MC_simul_quasi_anti_CV->expectation_by_sims(200, { 10000 }, "exp_MC_simul_quasi_anti_CV.csv");
		end = clock();

		std::cout << "Quasi Antithetic CV : " << double(end - start) / double(CLOCKS_PER_SEC) / 200 << std::setprecision(5) << std::endl;

	}
	catch (std::exception & e)
	{
		std::cout << e.what() << std::endl;
	}
	
	return 0;
}

// void display_mat(const std::vector<std::vector<double>>& A)
// {
	// for (int ii =0; ii<A.size();++ii)
	// {
		// for(int jj = 0; jj<A[0].size(); ++jj)
		// {
			// std::cout<<A[ii][jj]<<", ";
		// }
		// std::cout<<std::endl;
	// }
	
// }

void test_functions()
{
	double mu = 3. / 100.;
	std::vector<std::vector<double>> cov(2, std::vector<double>(2));

	cov[0][0] = 0.2 * 0.2;
	cov[1][1] = 0.3 * 0.2;
	cov[0][1] = 0.6 * std::sqrt(cov[0][0]) * std::sqrt(cov[1][1]);
	cov[1][0] = cov[0][1];

	UniformGenerator* ecuyer = new EcuyerCombined(1, 1);
	ContinuousGenerator* biv_norm = new NormalMultiVariate(ecuyer, { mu }, cov);

	std::cout << "----- FULL STOCK PATH -----" << std::endl;

	llong n_steps = 100;
	Stocks* stocksFull = new StocksStandardFullPath(biv_norm, 100, { mu }, 1, n_steps);

	std::vector<std::vector<std::vector<double>>> S_full(stocksFull->Generate(5));

	for (int i = 0; i < S_full.size(); i++)
	{
		for (int j = 0; j < S_full[i].size(); j++)
		{
			std::cout << S_full[i][j][n_steps - 1] << ", ";
		}
		std::cout << std::endl;
	}

	std::cout << "----- Test Basket Payoff ------" << std::endl;

	std::vector<double> weights = { 0.5, 0.5 };
	Payoff* V2 = new NPDBasketCall(100, weights);

	/*std::vector<double> CallBkt = (*V2)(S);

	for (int i = 0; i < S.size(); i++)
	{
		std::cout << CallBkt[i] << ", ";
		std::cout << std::endl;
	}*/
	
	
	
	std::cout<<"----- Test Cholesky --------"<<std::endl;
	std::vector<std::vector<double>> MAtrix;
	MAtrix.resize(4, std::vector<double>(4));
	for(int i = 0; i<4; i++)
	{
		for(int j = 0; j<4; j++)
		{
			MAtrix[i][j] = -1.0;
			if (i==j){
				MAtrix[i][i] = 5.0;
			}	
		}
	}

	std::vector<std::vector<double>> L(MAtrix);
	Cholesky(L);
	std::vector<std::vector<double>> Matrix2(MAtrix);
	inv_sym_defpos(MAtrix, Matrix2);
	
	std::vector<std::vector<double>> multi;
	
	for(int i = 0; i<4; i++)
	{
		std::cout<< Matrix2[i][0] <<", "<<  Matrix2[i][1] <<", "<<  Matrix2[i][2] <<", "<<  Matrix2[i][3] <<std::endl;
	}
	
	mult_matrix(MAtrix, Matrix2, multi);
	for(int i = 0; i<4; i++)
	{
		std::cout<< multi[i][0] <<", "<<  multi[i][1] <<", "<<  multi[i][2] <<", "<<  multi[i][3] <<std::endl;
	}
	
		// std::vector<std::vector<double>> MAtrix;
		// MAtrix.resize(4, std::vector<double>(4));
		// std::vector<std::vector<double>> Matrixinv;
		// for(int i = 0; i<4; i++)
		// {
			// for(int j = 0; j<4; j++)
			// {
				// MAtrix[i][j] = -1.0;
				// if (i==j){
					// MAtrix[i][i] = 5.0;
				// }	
			// }
		// }
		
		// std::cout<< inverse(MAtrix, Matrixinv)<<std::endl;
		// display_mat(Matrixinv);
		
		// std::ofstream f;
		// f.open("XQnorm2.csv");

		// for (llong i = 0; i < S.size(); ++i)
		// {
			// f <<S[i][0][S[0][0].size()-2]<< "\n";
		// }

		// f.close();
		
		// std::ofstream f1;
		// f1.open("XQnorm1.csv");

		// for (llong i = 0; i < S.size(); ++i)
		// {
			// f1 <<S[i][0][S[0][0].size()-1]<< "\n";
		// }

		// f1.close();
	
	
}
	
