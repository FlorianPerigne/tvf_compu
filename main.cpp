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
#include "CFCall.hpp"
#include "CFPut.hpp"
#include "NPDCall.hpp"
#include "NPDPut.hpp"
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

		cov[0][0] = 0.3*0.3;
		cov[1][1] = 0.15*0.15;
		cov[0][1] = 0.6*std::sqrt(cov[0][0])*std::sqrt(cov[1][1]);
		cov[1][0] = cov[0][1];

		ContinuousGenerator* biv_norm = new NormalMultiVariate(ecuyer, 0, cov);
		ContinuousGenerator* biv_norm2 = new NormalMultiVariate(ecuyer2, 0, cov);
		ContinuousGenerator* biv_norm3 = new NormalMultiVariate(vdc, vdc2, 0, cov);

		std::cout << "----- MONTECARLO ------" << std::endl;

		llong n_simu = 5000;

		NonPathDependent* call_payoff = new NPDCall(100, weights); //si weights alors basket
		NonPathDependent* put_payoff = new NPDPut(100, weights);
		PathDependent* call_payoff_PD = new PDCall(100, weights);

		R3R1Function* antithetic_function = new StandardAntithetic();

		ClosedForm* call_payoff_CF = new CFCall(100);
		double prix_bs = (*call_payoff_CF)(100, 0 ,1,std::sqrt(cov[0][0]));
		// double prix_bs_div = (*call_payoff_CF)(93, 0 ,1,0.2);//uniquement utile pour control d'un stock avec 7 de div cumulé et 0 de mu
		NonPathDependent* Call_clasic = new NPDCall(100);//si pas weights alors c'est un call classic sur le premier stock path
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

		llong n_sims = 300;

		std::cout << "EXPECTATION AND VARIANCE IN FUNCTION OF N_SIMULATION of paths" << std::endl;

		if (false) {
			MC_simul_standard->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\var_MC_simul_standard_2A.csv");
			MC_simul_quasi->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\var_MC_simul_quasi_2A.csv");
			MC_simul_standard_anti->variance_by_sims(150, { 5, 25, 50, 125, 250, 500, 1000, 2500, 5000 }, "..\\Graphs\\var_MC_simul_standard_anti_2A.csv");
			MC_simul_quasi_anti->variance_by_sims(150, { 5, 25, 50, 125, 250, 500, 1000, 2500, 5000 }, "..\\Graphs\\var_MC_simul_quasi_anti_2A.csv");
			MC_simul_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\var_MC_simul_CV_2A.csv");
			MC_simul_anti_CV->variance_by_sims(150, { 5, 25, 50, 125, 250, 500, 1000, 2500, 5000 }, "..\\Graphs\\var_MC_simul_anti_CV_2A.csv");
			MC_simul_quasi_CV->variance_by_sims(150, { 10, 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\var_MC_simul_quasi_CV_2A.csv");
			MC_simul_quasi_anti_CV->variance_by_sims(150, { 5, 25, 50, 125, 250, 500, 1000, 2500, 5000 }, "..\\Graphs\\var_MC_simul_quasi_anti_CV_2A.csv");
		}

		std::cout << "CV Var in function of Correlation" << std::endl;

		if (false)
		{
			std::vector<double> corr = { -0.999, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0,
									0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.999 };

			std::ofstream f;
			f.open("..\\Graphs\\var_MC_CV_corr.csv");

			for (int i = 0; i < corr.size(); ++i)
			{
				std::cout << "Corr = " << corr[i] << std::endl;
				cov[0][1] = corr[i] * std::sqrt(cov[0][0]) * std::sqrt(cov[1][1]);
				cov[1][0] = cov[0][1];

				ContinuousGenerator* biv_norm_corr = new NormalMultiVariate(vdc, vdc2, 0, cov);
				StocksTerminal* stocks_corr = new StocksStandardTerminal(biv_norm_corr, 100, mu, 1);
				MonteCarlo* mc_solver_quasi_CV_corr = new MonteCarloEuropean(stocks_corr, call_payoff, llong(n_simu / 2), Call_clasic, prix_bs);

				Simulation* MC_simul_corr = new Simulation(mc_solver_quasi_CV_corr);
				MC_simul_corr->compute_variance(300);

				f << corr[i] << "," << MC_simul_corr->get_V() << "\n";
			}

			f.close();
		}
		
		std::cout << "Variance of MC LSM" << std::endl;

		Basis* base = new BasisLaguerre(3);

		// --- Standard MC Terminal

		StocksFullPath* stocksF = new StocksStandardFullPath(biv_norm, 100, mu, 1, 12);
		MonteCarlo* mc_solver_fp = new MonteCarloLSM(stocksF, call_payoff_PD, n_simu, base);

		// --- MC Terminal with Antitethic variance reduction

		StocksFullPath* stocksF2 = new StocksAntitheticFullPath(biv_norm, 100, mu, 1, 12, antithetic_function);
		MonteCarlo* mc_solver_antithetic_fp = new MonteCarloLSM(stocksF2, call_payoff_PD, llong(n_simu / 2), base);
		MonteCarlo* mc = new MonteCarloEuropean(stocksF2, call_payoff_PD, llong(n_simu / 2));
		
		// --- MC Terminal with Control Variate variance reduction

		StocksFullPath* stocksF4 = new StocksStandardFullPath(biv_norm, 100, mu, 1, 12);
		MonteCarlo* mc_solver_CV_fp = new MonteCarloLSM(stocksF4, call_payoff_PD, llong(n_simu), base, Call_clasic, prix_bs);

		// --- MC Terminal with Antithetic and Control Variate variance reduction

		StocksFullPath* stocksF6 = new StocksAntitheticFullPath(biv_norm, 100, mu, 1, 12, antithetic_function);
		MonteCarlo* mc_solver_anti_CV_fp = new MonteCarloLSM(stocksF6, call_payoff_PD, llong(n_simu / 2), base, Call_clasic, prix_bs);

		Simulation* MC_simul_fp = new Simulation(mc_solver_fp);
		Simulation* MC_simul_antithetic_fp = new Simulation(mc_solver_antithetic_fp);
		Simulation* MC_simul_CV_fp = new Simulation(mc_solver_CV_fp);
		Simulation* MC_simul_anti_CV_fp = new Simulation(mc_solver_anti_CV_fp);

		// mc_solver_fp->Solve();
		// std::cout << mc_solver_fp->get_price() << std::endl;
		// mc_solver_antithetic_fp->Solve();
		// mc_solver_fp->Solve();
		// std::cout << mc_solver_antithetic_fp->get_price() << std::endl;
		// std::cout << mc_solver_fp->get_price() << std::endl;

		// mc_solver_CV_fp->Solve();
		// std::cout << mc_solver_CV_fp->get_price() << std::endl;
		// mc_solver_anti_CV_fp->Solve();
		// std::cout << mc_solver_anti_CV_fp->get_price() << std::endl;


		if (true) {
			//MC_simul_fp->variance_by_sims(150, { 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\Variance-BasketCall-LSM\\var_MC_simul_standard.csv");
			//MC_simul_antithetic_fp->variance_by_sims(150, {25, 50, 125, 250, 500, 1000, 2500, 5000 }, "..\\Graphs\\Variance-BasketCall-LSM\\var_MC_anti.csv");
			//MC_simul_CV_fp->variance_by_sims(150, { 50, 100, 250, 500, 1000, 2000, 5000, 10000 }, "..\\Graphs\\Variance-BasketCall-LSM\\var_MC_simul_standard_cv.csv");
			// MC_simul_anti_CV_fp->variance_by_sims(150, {250, 500, 1000, 2500, 5000 }, "..\\Graphs\\Variance-BasketCall-LSM\\var_MC_simul_anti_cv.csv");
		}

		// --- MC Terminal with Quasi random numbers

		//StocksFullPath* stocksF3 = new StocksStandardFullPath(biv_norm3, 100, mu, 1, 10);
		//MonteCarlo* mc_solver_quasi = new MonteCarloLSM(stocksF3, call_payoff_PD, n_simu, base);

		// --- MC Terminal with antithetic and quasi random numbers variance reduction

		//StocksFullPath* stocksF5 = new StocksAntitheticFullPath(biv_norm3, 100, mu, 1, 10, antithetic_function);
		//MonteCarlo* mc_solver_quasi_anti = new MonteCarloLSM(stocksF5, call_payoff_PD, llong(n_simu / 2), base);


		// --- MC Terminal  with Quasi random numbers and Control Variate variance reduction

		//StocksFullPath* stocksF7 = new StocksStandardFullPath(biv_norm3, 100, mu, 1, 10);
		//MonteCarlo* mc_solver_quasi_CV = new MonteCarloLSM(stocksF7, call_payoff_PD, llong(n_simu / 2), base, Call_clasic, prix_bs);

		// --- MC Terminal  with Quasi random numbers anthitetic and Control Variate variance reduction

		//StocksFullPath* stocksF8 = new StocksAntitheticFullPath(biv_norm3, 100, mu, 1, 10, antithetic_function);
		//MonteCarlo* mc_simul_quasi_anti_CV = new MonteCarloLSM(stocksF8, call_payoff_PD, llong(n_simu / 2), base, Call_clasic, prix_bs);

		if(false)
		{
			
			PathDependent* basket_call_payoff_PD = new PDCall(100, weights);
			
			ClosedForm* call_payoff_CF2 = new CFCall(100);
			double prix_bs = (*call_payoff_CF2)(100, 0 ,1,std::sqrt(cov[0][0]));
			NonPathDependent* Call_clasic_1asset = new NPDCall(100);
			//ici on control avec un call classic sur la premiere path 
			//ainsi plus le payoff qu'on estime en proche du call sur le premier asset alors plus on a un control variate proche de ce que l'on cherche
			
			StocksFullPath* stocksF00 = new StocksStandardFullPath(biv_norm, 100, mu, 1, 12);
			MonteCarlo* mc_solver_CV_fp_test = new MonteCarloLSM(stocksF00, basket_call_payoff_PD, llong(n_simu), base, Call_clasic_1asset, prix_bs);
			
			mc_solver_CV_fp_test->Solve();
			std::cout << mc_solver_CV_fp_test->get_price() << std::endl;
			
		}
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
	Payoff* V2 = new NPDCall(100, weights);

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
	
