#include "Utilities.hh"
#include <cmath>
#include "cctk_Parameters.h"

double smoothR(double rr)
{
	DECLARE_CCTK_PARAMETERS;
	return pow(pow(rr, 4.0)+pow(zero_epsilon, 4.0), 0.25);
}

