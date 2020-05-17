#include "Sobol.hpp"
#include <cmath>
#include <iostream>
#include <ostream>
#include <algorithm>
#include <vector>
#include <limits>
#include <math.h>


//code of numerical finance


Sobol::Sobol(llong seed)
	:QuasiGenerator(seed)
{
	n = -1;
	std::cout<<"Sobol Sequence generation"<<std::endl;
	sequence(n);
	std::cout<<"Sobol Sequence generation Done"<<std::endl;
	n = 2;
	simul = true;
	
};

double Sobol::Generate()
{
	if(simul == true) 
	{
		// std::cout<<"Sobol Sequence generation"<<std::endl;
		sequence(n);
		// std::cout<<"Sobol Sequence generation Done"<<std::endl;
		simul = false;
		// std::cout<<"0D elemente: "<<x[0]<<std::endl;
		// std::cout<<"1D elemente: "<<x[1]<<std::endl;
		// std::cout<<"2D elemente: "<<x[2]<<std::endl;
		// std::cout<<"size: "<<x.size()<<std::endl;
		return x[0];
	}
	else
	{
		simul = true;
		return x[1];
		
	}
};

void Sobol::sequence(int inputn)
{
	/*
	When n is negative, internally initializes a set of MAXBIT direction numbers for each of MAXDIM
	different Sobol’ sequences. When n is positive (but ≤MAXDIM), returns as the vector x[0..n-1]
	the next values from n of these sequences. (n must not be changed between initializations.)
	*/
// x.resize(n)
	const int MAXBIT=30,MAXDIM=6;
	int j,k,l;
	unsigned long i,im,ipp;
	static double fac;
	// static unsigned long in;
	// std::vector<static unsigned long> ix(MAXDIM);
	// std::vector<static unsigned long*> iu(MAXBIT);
	// (*iu).resize(MAXBIT+1);
	static unsigned long in,ix[MAXDIM+1],*iu[MAXBIT+1];
	static int mdeg[MAXDIM]={1,2,3,3,4,4};
	static unsigned long ip[MAXDIM] ={0,1,1,2,1,4};
	static unsigned long iv[MAXDIM*MAXBIT] = {1,1,1,1,1,1,3,1,3,3,1,1,5,7,7,3,3,5,15,11,5,15,13,9};


	if (n < 0) 
	{

		for (k=0;k<MAXDIM;k++) ix[k]=0;
		in=0;
		if (iv[0] != 1) return;
		fac=1.0/(1L << MAXBIT);
		for (j=0,k=0;j<MAXBIT;j++,k+=MAXDIM) 
		{
			iu[j] = &iv[k];
		}

		// To allowboth 1D and 2D addressing.
		for (k=0;k<MAXDIM;k++) 
		{
			for (j=0;j<mdeg[k];j++) iu[j][k] <<= (MAXBIT-1-j);
			for (j=mdeg[k];j<MAXBIT;j++) 
			{
				ipp=ip[k];
				i=iu[j-mdeg[k]][k];
				i ^= (i >> mdeg[k]);
				for (l=mdeg[k]-1;l>=1;l--) 
				{
				if (ipp & 1) i ^= iu[j-l][k];
				ipp >>= 1;
				}
				iu[j][k]=i;
			}
		}
	}
	else 
	{
		x.resize(n+1);
		im=in++;
		for (j=0;j<MAXBIT;j++) 
		{
			if (!(im & 1)) break;
			im >>= 1;
		}
		if (j >= MAXBIT) throw std::runtime_error("MAXBIT too small in sobseq");
		im=j*MAXDIM;
		for (k=0;k<std::min(n,MAXDIM);k++)
		{
			ix[k] ^= iv[im+k];
			x[k]=ix[k]*fac;
			// std::cout<<fac<<", "<<ix[k]<<std::endl;
		}
	}
};