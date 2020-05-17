#include "MonteCarloLSM.hpp"
#include <algorithm>
#include <numeric>
#include <iomanip>


//Constructor for classical LSM
MonteCarloLSM::MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction)
	: MonteCarlo(stocks, payoff, n_sims), LSM_Basis(BasisFunction)
{
}


//for control variate with a PathDependent payoff
MonteCarloLSM::MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction, PathDependent* payoff_CV, double closedFormValue)
	: MonteCarlo(stocks, payoff, n_sims, payoff_CV, closedFormValue), LSM_Basis(BasisFunction)
{
	PD_or_NPD = 1;
}

//for control variate with a NonPathDependent payoff
MonteCarloLSM::MonteCarloLSM(StocksFullPath* stocks, PathDependent* payoff, llong n_sims, Basis* BasisFunction, NonPathDependent* payoff_CV, double closedFormValue)
	: MonteCarlo(stocks, payoff, n_sims, payoff_CV, closedFormValue), LSM_Basis(BasisFunction)
{
	PD_or_NPD = 0;
}



void MonteCarloLSM::Solve()
{
	
	std::vector<std::vector<std::vector<double>>> S(mc_stocks->Generate(N_sims)); //generate stock path
	std::vector<double> P((*mc_payoff)(S, S[0][0].size()-1));  //get the payoff at maturity
	std::vector<double> itm_path(P.size());
	std::vector<double> Weights(mc_payoff->get_weights());
	
	std::vector<double> P2;
	
	//get payoff at maturity for the control variate case
	if (PD_or_NPD)
	{
		P2 = ((*MC_payoff_CV)(S, S[0][0].size()-1));
	}
	else
	{
		P2 = ((*MC_payoff_CV)(S));
	}
	
	std::vector<double> itm_path2(P2.size());
	std::vector<double> Weights2(MC_payoff_CV->get_weights());

	price = 0;
	double price2 = 0;
	llong nb_steps = S[0][0].size();
	double dt = (mc_stocks->get_maturity())/(nb_steps-1);
	double r = mc_stocks->get_mu();

	for (llong k = nb_steps-2; k>0; k--)
	{
		std::vector<double> X;
		std::vector<double> Y;
		
		std::vector<double> X2;
		std::vector<double> Y2;
		
		for(llong i = 0; i<P.size(); i++)
		{
			
			//seperate the path in the money
			if (P[i]>0)
			{
				itm_path[i] = 1;
				double xi = 0;
				for (llong j = 0; j<Weights.size(); j++)
				{
					xi += Weights[j]*S[i][j][k];
				}
				X.push_back(xi);
				Y.push_back(std::exp(-r * dt)*P[i]);
			}
			else
			{
				itm_path[i] = 0;
			}
			
			//for the control variate case
			if ((P2[i]>0) && PD_or_NPD)
			{
				itm_path2[i] = 1;
				double xi = 0;
				for (llong j = 0; j<Weights2.size(); j++)
				{
					xi += Weights2[j]*S[i][j][k];
				}
				X2.push_back(xi);
				Y2.push_back(std::exp(-r * dt)*P2[i]);
			}
			else
			{
				itm_path2[i] = 0;
			}
			
		}

		std::vector<double> C_hat(Regression_C_Hat(X, Y)); //Regression of X on Y
		std::vector<double> C_hat2;
		std::vector<double> G2;
		std::vector<double> G((*mc_payoff)(S, k)); //compute the payoff at this step
		
		
		if (X2.size()>0)
		{
			C_hat2 = Regression_C_Hat(X2, Y2);
			G2 = (*MC_payoff_CV)(S, k);
		}
		
		
		llong c = 0;
		llong c2 = 0;
		
		for (llong i = 0; i<P.size(); i++)
		{
			if (itm_path[i] == 1) //if we are in the money we can look if it is better to continue or exercise
			{
				
				if (G[i]>C_hat[c])
				{
					P[i] = G[i];
				}
				c +=1;
			}
			else
			{
				P[i] = std::exp(-r*dt)*P[i];
			}
			
			//for the control variate case
			if ((itm_path2[i] == 1) && PD_or_NPD)
			{	

				if (G2[i]>C_hat2[c2])
				{
					P2[i] = G2[i];
				}
				c2 +=1;
			}
			else
			{
				P2[i] = std::exp(-r*dt)*P2[i];
			}
		}

	}
	
	price2 = std::accumulate(P2.begin(), P2.end(), 0.0)/(P2.size());
	
	if (PD_or_NPD)
	{
		price2 *= std::exp(-r*dt);
	}
	else{
		price2 *= std::exp(-mc_stocks->get_mu() * mc_stocks->get_maturity());
	}
	
	price = std::exp(-r*dt)*std::accumulate(P.begin(), P.end(), 0.0)/(P.size());

	price = price + (closedFormValue - price2);
}

//Regression of X on Y 
std::vector<double> MonteCarloLSM::Regression_C_Hat(const std::vector<double>& X, const std::vector<double>& Y)
{
	
	std::vector<std::vector<double>> L(LSM_Basis->get_matrix_L(X));
	std::vector<std::vector<double>> LT;
	std::vector<std::vector<double>> A;
	std::vector<std::vector<double>> Ainv;
	std::vector<std::vector<double>> C;
	std::vector<double> Beta;	
	std::vector<double> C_hat;
	

	transpose_matrix(L, LT);	//L'
	mult_matrix(LT,L,A);	//L'L
	inv_sym_defpos(A, Ainv);   // (L'L)^-1

	mult_matrix(Ainv,LT,C);		//(L'L)^-1 L'
	mult_matrix_vect(C, Y, Beta); //beta = (L'L)^-1 L'Y

	mult_matrix_vect(L, Beta, C_hat); //C_hat = beta L 
	
	return C_hat;
	
	
}

