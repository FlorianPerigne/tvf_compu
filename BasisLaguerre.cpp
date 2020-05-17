#include "BasisLaguerre.hpp"

BasisLaguerre::BasisLaguerre(llong number)
	:Basis(number)
{
}

//get the B_nb degree laguerre polynome
double BasisLaguerre::get_Lk(double X)
{

	double L = 0;
	for (llong i = 0; i<=B_nb; ++i)
	{
		L += factorial(B_nb)/(factorial(i)*factorial(B_nb-i)) * std::pow(-1.0, i)/factorial(i) * std::pow(X, i);

	}

	return L;
}

//get the kth degree laguerre polynome
double BasisLaguerre::get_Lk(double X, llong k)
{

	double L = 0;
	for (llong i = 0; i<=k; ++i)
	{
		L += factorial(k)/(factorial(i)*factorial(k-i)) * std::pow(-1.0, i)/factorial(i) * std::pow(X, i);

	}

	return L;
}

//create the matrix (L_j(X_i))
std::vector<std::vector<double>> BasisLaguerre::get_matrix_L(const std::vector<double>& X)
{
	std::vector<std::vector<double>> L;
	L.resize(X.size(), std::vector<double>(B_nb+1));

	for(llong i = 0; i<X.size();++i)
	{
		for(llong j = 0; j<=B_nb;++j)
		{
			L[i][j] = get_Lk(X[i], j);
			
		}
	}
	return L;
	
}