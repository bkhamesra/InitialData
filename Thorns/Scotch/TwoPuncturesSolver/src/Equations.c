// TwoPunctures:  File  "Equations.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <assert.h>
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

#define DEBUG

// utility macros
#define SQR(x) ((x)*(x))
#define QUAD(x) (SQR(x)*SQR(x))
#define CUBE(x) ((x)*SQR(x))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

// U.d0[ivar]   = U[ivar];  (ivar = 0..nvar-1) 
// U.d1[ivar]   = U[ivar]_x;  
// U.d2[ivar]   = U[ivar]_y;  
// U.d3[ivar]   = U[ivar]_z;  
// U.d11[ivar]  = U[ivar]_xx; 
// U.d12[ivar]  = U[ivar]_xy; 
// U.d13[ivar]  = U[ivar]_xz;
// U.d22[ivar]  = U[ivar]_yy;
// U.d23[ivar]  = U[ivar]_yz;
// U.d33[ivar]  = U[ivar]_zz;

//---------------------------------------------------------------
//*******           Nonlinear Equations                *********/
//---------------------------------------------------------------
void
NonLinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
                 double A, double B, double X, double R,
		 double x, double r, double phi,
		 double y, double z, derivs U, double *values, double time)
//BK: Compute the Hamiltonian constraint violation - Equation 62 from the Michael's paper
{
  DECLARE_CCTK_PARAMETERS;

  double Aij[6], psiguess, rhohat, junk, junk3[3];

  getSourceInfo(x-move_origin_x, y, z, &rhohat, junk3,  &junk, &psiguess, Aij);

  double psi = psiguess + U.d0[0];
  double AA = traceSquare(Aij);

  
  double lhs = U.d11[0] + U.d22[0] + U.d33[0] + 0.125 * AA /( QUAD(psi)*CUBE(psi)); 
  double rhs = -2.0*Pi*rhohat/pow(psi, conformal_density_power-5.);
  values[0]=lhs-rhs;
//	CCTK_VInfo(CCTK_THORNSTRING, "Trace-square info at (%g, %g, %g): Axx = %g, Axy = %g, Axz = %g, Ayy = %g, Ayz = %g, Azz = %g, trace A^2= %g", x-move_origin_x, y, z, Aij[0], Aij[1], Aij[2], Aij[3], Aij[4], Aij[5], AA);
	/*
   if(abs(x-par_b) < .5 || abs(x+par_b) < .5)
   {
	   CCTK_VInfo(CCTK_THORNSTRING, "NonLinEquations at (%g, %g, %g): U.d0 = %e, psi-1 = %e, AijAij = %e, Laplacian: %e, rhs = %e, lhs = %e, residual = %e", x-move_origin_x, y, z, U.d0[0], psi-1.0, AA, U.d11[0] + U.d22[0] + U.d33[0], rhs, lhs, values[0]);
   }
   */
}

//---------------------------------------------------------------
//*******               Linear Equations                *********/
//---------------------------------------------------------------
void
LinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
              double A, double B, double X, double R,
	      double x, double r, double phi,
	      double y, double z, derivs dU, derivs U, double *values, double time)
{
  DECLARE_CCTK_PARAMETERS;
  double Aij[6], psiguess, rhohat, junk, junk3[3];

  getSourceInfo(x-move_origin_x, y, z, &rhohat, junk3, &junk, &psiguess, Aij);

  double psi = psiguess + U.d0[0];
  double AA = traceSquare(Aij);
  double lhs = dU.d11[0] + dU.d22[0] + dU.d33[0]-.875*(dU.d0[0]*AA)/( QUAD(psi)*QUAD(psi)); 
  double rhs = 2.*(conformal_density_power-5.)*Pi/pow(psi, conformal_density_power-4.)*rhohat*dU.d0[0];
  values[0] = lhs - rhs;
}

//---------------------------------------------------------------
//*******               makeInitialGuess                *********/
//---------------------------------------------------------------
void makeInitialGuess_Psi (double x, double y, double z, double *u, double time)
{
	DECLARE_CCTK_PARAMETERS;
  double psiguess, psi, junk, junk3[3], junk6[6];

  getSourceInfo(x-move_origin_x, y, z, &junk, junk3, &psi, &psiguess, junk6);
  u[0] = psi-psiguess;
}
