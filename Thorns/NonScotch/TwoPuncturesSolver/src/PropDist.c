// TwoPunctures:  File  "Equations.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

// U.d11[ivar]  = U[ivar]_xx; 

double
propdist (double m_plus_i, double m_minus_i, derivs v)
{
  DECLARE_CCTK_PARAMETERS;
  double U, x, dx, hor1, hor2, x_grid_moved, rp_tmp, rm_tmp, r_plus, r_minus,
 	 gxx, Stot;
  int i, N;

//fprintf(stderr, "\nEntered Propdist.c");
  N = 1000;
  hor1 = - par_b + m_minus_i/2 + move_origin_x;  // Move origin here
  hor2 = par_b - m_plus_i/2 + move_origin_x;     // Move origin here
  // Corrections for the apparent horizon location
  hor1 += - 0.024*m_minus_i;
  hor2 += 0.024*m_plus_i;
 
  dx = (hor2 - hor1)/(N-1);

fprintf(stderr, "\nProper Separation Calculation: Integrating from %f to %f, dx = %f",hor1, hor2, dx);

  Stot=0.0;

  for ( i=0; i<(N-1); i++ )
  {

     // Find the point along the x axis
     x = hor1 + i*dx;

     // Calculate the metric component, automatically evaluate at each point.

	  U = PunctTaylorExpandAtArbitPosition
	      (0, 1, npoints_A, npoints_B, npoints_phi, v, x, 0, 0);
	  //U = PunctIntPolAtArbitPosition
	  //    (0, 1, npoints_A, npoints_B, npoints_phi, v, x, 0, 0);

        double rp_tmp = sqrt(pow(x - par_b,2));
        double rm_tmp = sqrt(pow(x + par_b,2));

        /* modify r_plus and r_minus to get smooth behaviour on puncture */
        double r_plus = pow(pow(rp_tmp,4)+pow(epsilon,4),0.25);
        double r_minus = pow(pow(rm_tmp,4)+pow(epsilon,4),0.25);

        double psi1 = 1                                                //
          + 0.5 * m_plus_i / r_plus                                          //
          + 0.5 * m_minus_i / r_minus + U;

	gxx = 0;
	// gxx = pow(psi1,4);

     // Integrate
     Stot += pow(psi1,2) * dx;
 
  }

fprintf(stderr, "\nProper Separation Calculation: S = %f", Stot);

  return Stot;

}

