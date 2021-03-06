#ifndef normal_HPP
#define normal_HPP

#include "ContinuousGenerator.hpp"
#include "uniform_generator.hpp"
#include "discrete_generator.hpp"

enum NormalAlgo {
	BoxMuller,
	CLT,
	N_RejectionSamp
};


class Normal : public ContinuousGenerator
{
public:
	Normal(UniformGenerator* gen, NormalAlgo type, double mu, double sigma);
	Normal(UniformGenerator* gen, UniformGenerator* gen2, double mu, double sigma);

	double Generate();

private:
	UniformGenerator* gen1;
	UniformGenerator* gen2;

	DiscreteGenerator* ht_gen;
	ContinuousGenerator* exp_gen;

	double mu;
	double sigma;
	NormalAlgo type;
	double x, y;
	bool isGenerated = false;
};


#endif