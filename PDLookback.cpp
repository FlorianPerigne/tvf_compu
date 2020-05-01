#include "PDLookback.hpp"

PDLookback::PDLookback(double strike, int days)
	:PathDependent(strike), PD_days(days)
{
}