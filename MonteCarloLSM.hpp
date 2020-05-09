#ifndef MONTECARLOLSM_HPP
#define MONTECARLOLSM_HPP

#include "MonteCarlo.hpp"
#include "Basis.hpp"
#include "Tools.hpp"
#include <typeinfo>

class MonteCarloLSM : public MonteCarlo
{
public:
	MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction);
	MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction, PathDependent* payoff_CV, double closedFormValue);
	MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction, NonPathDependent* payoff_CV, double closedFormValue);
	
	void Solve();
	
	
	std::vector<double> Regression_C_Hat(const std::vector<double>& X, const std::vector<double>& Y);

private:
	Basis* LSM_Basis;
	bool PD_or_NPD;

};

#endif