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
#include "Sobol.hpp"

#include "Normal.hpp"
#include "NormalMultiVariate.hpp"

#include "StandardAntithetic.hpp"

#include "StocksAntitheticTerminal.hpp"
#include "StocksStandardTerminal.hpp"
#include "StocksAntitheticFullPath.hpp"
#include "StocksStandardFullPath.hpp"

#include "Payoff.hpp"
#include "NonPathDependent.hpp"
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

int main()
{
	try
	{
		UniformGenerator* ecuyer = new EcuyerCombined(1, 1);

		UniformGenerator* vdc = new VanDerCorput(1, 2);
		UniformGenerator* vdc2 = new VanDerCorput(1, 3);
		
		UniformGenerator* sobol = new Sobol(1);
		
		// Parameters of the assets dynamics

			double mu = 0. / 100.; // risk free rate

			double rho = 0.6; // correlation between the two assets
			std::vector<double> weights = { 0.5, 0.5 }; // weights in the basket

			std::vector<std::vector<double>> cov(2, std::vector<double>(2)); // Variance-Covariance matrix

			cov[0][0] = 0.3 * 0.3; // Variance of the asset on the diagonal
			cov[1][1] = 0.15 * 0.15;
			cov[0][1] = rho * std::sqrt(cov[0][0]) * std::sqrt(cov[1][1]); // Covariance between asset i and j --> correlation * sqrt(variance_i * variance_j)
			cov[1][0] = cov[0][1];

		// Random numbers generation

			ContinuousGenerator* norm_standard = new NormalMultiVariate(ecuyer, 0, cov);
			ContinuousGenerator* norm_quasi = new NormalMultiVariate(vdc, vdc2, 0, cov);
			ContinuousGenerator* norm_quasi_sob = new NormalMultiVariate(sobol, 0, cov);

		// Payoff definition

			NonPathDependent* call_payoff = new NPDCall(100, weights); // if weights are provided, then it's behaving as a basket option
			NonPathDependent* Call_classic = new NPDCall(100);

			PathDependent* call_payoff_PD = new PDCall(100, weights);

		// For Monte Carlo algorithm

			llong n_simu = 5000;

			// Specific for the Longstaff-Schwartz

			Basis* base = new BasisLaguerre(3);


		// For the variance reduction methods

			R3R1Function* antithetic_function = new StandardAntithetic(); // Antithetic

			ClosedForm* call_payoff_CF = new CFCall(100); // For CV
			double prix_bs = (*call_payoff_CF)(100, 0, 1, std::sqrt(cov[0][0]));


		// --- Standard MC Terminal

		StocksTerminal* stocksT = new StocksStandardTerminal(norm_standard, 100, mu, 1);
		MonteCarlo* mc_solver = new MonteCarloEuropean(stocksT, call_payoff, n_simu);

		// --- MC Terminal with Antitethic variance reduction

		StocksTerminal* stocksT2 = new StocksAntitheticTerminal(norm_standard, 100, mu, 1, antithetic_function);
		MonteCarlo* mc_solver_antithetic = new MonteCarloEuropean(stocksT2, call_payoff, llong(n_simu/2));

		// --- MC Terminal with Quasi random numbers

		StocksTerminal* stocksT3 = new StocksStandardTerminal(norm_quasi_sob, 100, mu, 1);
		MonteCarlo* mc_solver_quasi = new MonteCarloEuropean(stocksT3, call_payoff, n_simu);

		// --- MC Terminal with Control Variate variance reduction

		StocksTerminal* stocksT4 = new StocksStandardTerminal(norm_standard, 100, mu, 1);
		MonteCarlo* mc_solver_CV = new MonteCarloEuropean(stocksT4, call_payoff, n_simu, Call_classic, prix_bs);

		// Solving the Monte Carlo object...

		mc_solver->Solve();
		std::cout << "Monte Carlo price : " << mc_solver->get_price() << std::endl;

		// Or using the Simulation class to have more possibilites and verify the properties of the solver

		Simulation* MC_simul_standard = new Simulation(mc_solver);

		llong n_sims_var = 150;

		std::vector<llong> n_sims_std = { 100, 250, 500, 1000, 2000, 5000, 7500, 10000 };

		if (false) {
			MC_simul_standard->variance_by_sims(n_sims_var, n_sims_std, "var_MC_simul_standardcsv");
		}

		// --- Standard MC LSM

		StocksFullPath* stocksF = new StocksStandardFullPath(norm_standard, 100, mu, 1, 12);
		MonteCarlo* mc_solver_fp = new MonteCarloLSM(stocksF, call_payoff_PD, n_simu, base);

		// --- MC LSM with Antitethic variance reduction

		StocksFullPath* stocksF2 = new StocksAntitheticFullPath(norm_standard, 100, mu, 1, 12, antithetic_function);
		MonteCarlo* mc_solver_antithetic_fp = new MonteCarloLSM(stocksF2, call_payoff_PD, llong(n_simu / 2), base);
		
		// --- MC LSM with Quasi random numbers

		StocksFullPath* stocksF3 = new StocksStandardFullPath(norm_quasi, 100, mu, 1, 12);
		MonteCarlo* mc_solver_quasi_fp = new MonteCarloLSM(stocksF3, call_payoff_PD, n_simu, base);

		// --- MC LSM with Control Variate variance reduction

		StocksFullPath* stocksF4 = new StocksStandardFullPath(norm_standard, 100, mu, 1, 12);
		MonteCarlo* mc_solver_CV_fp = new MonteCarloLSM(stocksF4, call_payoff_PD, n_simu, base, Call_classic, prix_bs);

		Simulation* MC_simul_fp = new Simulation(mc_solver_fp);

		mc_solver_fp->Solve();
		std::cout << "Monte Carlo LSM price : " << mc_solver_fp->get_price() << std::endl;
		
		if (false) {
			MC_simul_fp->variance_by_sims(n_sims_var, n_sims_std, "var_MC_simul_standard.csv");
		}

	}
	catch (std::exception & e)
	{
		std::cout << e.what() << std::endl;
	}
	
	return 0;
}
