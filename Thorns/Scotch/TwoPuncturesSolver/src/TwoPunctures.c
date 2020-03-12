// TwoPunctures:  File  "TwoPunctures.c"
// $Id: TwoPunctures.c,v 1.6 2004/10/25 13:41:10 herrmann Exp $

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

static inline double pow2 (const double x)
{
  return x*x;
}

static inline double pow4 (const double x)
{
  return x*x*x*x;
}

static int grid_data_unset=1;
static int derivs_unset=1;
static double *F = NULL, *AijAij = NULL; //BK Comments: What is F?       
 //BK: v = U 
static derivs u, v, w, A; /* order of Aij components is 11,12,13,22,23,33 */
static derivs v_coeffs, w_coeffs, Aij_coeffs;

// -------------------------------------------------------------------
//

CCTK_REAL TwoPuncturesSolver_GetU(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
	DECLARE_CCTK_PARAMETERS;
	//BK: Returns value of U at arbitrary point on the grid from spectral grid using taylor expansion or interpolation.
	//BK: Returns value of U at arbitrary point on the grid from spectral grid using taylor expansion or interpolation.
	return ValueAtArbitPosition(0, 1, npoints_A, npoints_B, npoints_phi, v, v_coeffs, x+move_origin_x, y, z);
}

CCTK_REAL TwoPuncturesSolver_GetConformalFactor(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z)
{
	DECLARE_CCTK_PARAMETERS;
	double psiguess, junk, junk3[3], Aij_flattened[6];
	//BK: gets initial guess for conformal factor and traceless conformal curvature from BowenID
	getSourceInfo(x, y, z, &junk, junk3, &junk, &psiguess, Aij_flattened);

	//BK: Returns value of U at arbitrary point on the grid from spectral grid using taylor expansion or interpolation.
	const double U = ValueAtArbitPosition(0, 1, npoints_A, npoints_B, npoints_phi, v, v_coeffs, x+move_origin_x, y, z);
	return U + psiguess;
}

void
TwoPuncturesSolver (CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  //BK: grid setup method - Taylor expansion or interpolation
  gsm_type gsm;
  gsm = getGSM();

  int nvar_psi = 1, nvar_Aij = 4, ncomp_Aij = 6;
  //bk: total points on spectral grid for each spectral coordinate a, b, phi
  int n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

  int ntotal_psi  = n1 * n2 * n3 * nvar_psi; 
  int ntotal_Aij = n1 * n2 * n3 * nvar_Aij;  //BK Comments: why nvar_Aij and not ncomp_Aij
  int ntotal_AijAij = n1 * n2 * n3 * 1;
  int precomputed=0;
  int success;             /* return value of Newton solver */

  /*
  if (grid_data_unset || derivs_unset) {
    *par_b_pcmp=-1.;
	*move_origin_x_pcmp=-1.;
  }

  if(*par_b_pcmp != par_b || *move_origin_x_pcmp != move_origin_x)
  {
	  grid_data_unset = 1;
	  derivs_unset = 1;
  }
  */

  TwoPuncturesSolver_ForceResolve(derivs_unset);

    //F_of_v (nvar_psi, n1, n2, n3, v, F, u, cctk_time);       //

    /* print out ADM mass, eq.: \Delta M_ADM=2*r*u=4*b*V for A=1,B=0,phi=0 */
  if(grid_data_unset)
  {
     *admMass = (- 4*par_b*PunctEvalAtArbitPosition(v.d0, 0, 1, 0, 0, 1, n1, n2, n3));
    CCTK_VInfo (CCTK_THORNSTRING, "non-puncture ADM mass is %g", *admMass);

    /* BK: Evaluate the ADM mass at radius given by parameter admrad */
    CCTK_REAL masstmp,radtmp;
    masstmp =PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, admrad, 0.0, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, 0.0, admrad, 0.0);
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, 0.0, 0.0, admrad);
    radtmp=1/sqrt(3.0)*admrad;
    masstmp+=PunctIntPolAtArbitPosition(0, nvar_psi, n1, n2, n3, v, radtmp, radtmp, radtmp);
    masstmp/=4.0;

    *admMassRad = 2.0*admrad*masstmp;

    /* BK: Evaluate u_plus and u_minues via interpolation */
    CCTK_VInfo (CCTK_THORNSTRING, "non-puncture ADM mass from r= %.19g is %.19g", admrad,*admMassRad);
    switch (gsm)
    {
      case GSM_Taylor_expansion:
        *u_plus = PunctTaylorExpandAtArbitPosition
           (0, nvar_psi, n1, n2, n3, v, par_b, 0.0, 0.0);
        *u_minus = PunctTaylorExpandAtArbitPosition
           (0, nvar_psi, n1, n2, n3, v, -par_b, 0.0, 0.0);
        break;
      case GSM_evaluation:
        *u_plus = PunctIntPolAtArbitPosition
            (0, nvar_psi, n1, n2, n3, v, par_b, 0.0, 0.0);
        *u_minus = PunctIntPolAtArbitPosition
            (0, nvar_psi, n1, n2, n3, v, -par_b, 0.0, 0.0);
        break;
      default:
        fprintf(stderr,"we don't have a known gsm: gsm=%d\n",gsm);
        fprintf(stderr,"  and taylor=%d and eval=%d\n",GSM_Taylor_expansion,GSM_evaluation);
        assert (0);
    }
	CCTK_VInfo(CCTK_THORNSTRING, "u+ = %g, u-=%g", *u_plus, *u_minus);
  }
  else
  {
	  CCTK_VInfo(CCTK_THORNSTRING, "no solve performed; see parameters or precomputed values if this is unexpected behavior");
  }

  // Set standard BY puncture ID with U=0 
  if (bypass_solve) {
       CCTK_INFO ("NOT solving the constraints.");
       *admMass = 0.;
       *admMassRad = 0.;
       CCTK_VInfo (CCTK_THORNSTRING, "ADM mass is %g", *admMass);
	   *u_plus = 0.;
	   *u_minus = 0.;
  }


  CCTK_INFO ("Interpolating result");

  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING,"Moving origin in x direction by %f: "
                 "Punctures to appear at %f and %f\n",
                 move_origin_x,-par_b-move_origin_x,par_b-move_origin_x);
  }

#pragma omp parallel for
  for (int ind = 0; ind < cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; ++ind)
  {
        CCTK_REAL x_grid_moved;
	/* BK: Align the puncture on carpet grid to the one on spectral grid */
	/* make sure we are not exactly on a puncture */
	if ((x[ind]+move_origin_x==par_b || x[ind]+move_origin_x==-par_b) 
		&& y[ind]==0 && z[ind]==0) {
	  x_grid_moved = x[ind] + move_origin_x +1e-15;
	}
	else {
	  x_grid_moved = x[ind] + move_origin_x;
	}
/*        fprintf(stderr,"x[ind] %f move_origin_x %f x_grid_move %f\n",x[ind],move_origin_x,x_grid_moved);*/

	// avoid NaNs
        double U;
	if (bypass_solve) {
		makeInitialGuess_Psi(x_grid_moved, y[ind], z[ind], &U, 0.);  //BK: Compute U from BowenID (psi - psi_guess)
	} else {

	        switch (gsm)
	        {
	        case GSM_Taylor_expansion:
	          U = PunctTaylorExpandAtArbitPosition
	            (0, nvar_psi, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
	          break;
	        case GSM_evaluation:
                  if (use_full_transform) {
                    U = PunctIntPolCoeffsAtArbitPosition
                      (0, nvar_psi, n1, n2, n3, v_coeffs, x_grid_moved, y[ind], z[ind]);
                  } else {
                    U = PunctIntPolAtArbitPosition
                      (0, nvar_psi, n1, n2, n3, v, x_grid_moved, y[ind], z[ind]);
                  }
	          break;
	        default:
	          assert (0);
	        }

			/*
			if(abs(x_grid_moved + par_b) < .5 && abs(y[ind]) < .5 && abs(z[ind]) < .5)
			{
				CCTK_VInfo(CCTK_THORNSTRING, "U(%g, %g, %g) = %g", x[ind], y[ind], z[ind], U);
			}
			*/

	}

	double psiguess, junk, junk3[3], Aij_flattened[6];
	getSourceInfo(x[ind], y[ind], z[ind], &junk, junk3, &junk, &psiguess, Aij_flattened);
        const double psi1 = psiguess + U;

        if(keep_u_around)
        {   
	    puncture_u[ind] = U;
            conf_fac[ind] = psi1;
            conf_fac_bowen[ind] = psiguess;
            lorentz_fac[ind] = getLorentzFactor(x[ind], y[ind], z[ind]);
            polytrope_cons[ind] = getConfPolyConstant(x[ind], y[ind], z[ind])*pow(psi1, 8.);
	}

		double Aij[3][3];
		unflattenAij(Aij_flattened, Aij);

        gxx[ind] = pow4 (psi1);
        gxy[ind] = 0;
        gxz[ind] = 0;
        gyy[ind] = pow4 (psi1);
        gyz[ind] = 0;
        gzz[ind] = pow4 (psi1);
		
        kxx[ind] = Aij[0][0] / pow2(psi1);
        kxy[ind] = Aij[0][1] / pow2(psi1);
        kxz[ind] = Aij[0][2] / pow2(psi1);
        kyy[ind] = Aij[1][1] / pow2(psi1);
        kyz[ind] = Aij[1][2] / pow2(psi1);
        kzz[ind] = Aij[2][2] / pow2(psi1);

        if (CCTK_Equals(initial_lapse,"utb")) {
          alp[ind] = 1./( pow2(psi1) );
        }
        else if (CCTK_Equals(initial_lapse,"2/(1+psi_{BL}^4)")) {
          alp[ind] = 2./( 1. + pow4(psi1) );
        }
        else if (CCTK_Equals(initial_lapse,"psi^2")) {
	  double w = getLorentzFactor(x[ind],y[ind],z[ind]);
          alp[ind] =  pow2(psi1*w) ;
        }
  }
  CCTK_INFO("done with initial data setup for now");
  if (0) {
    /* Keep the result around for the next time */
    free_dvector (F, 0, ntotal_psi - 1);
    free_derivs (&u, ntotal_psi);
    free_derivs (&v, ntotal_psi);
  }

  grid_data_unset = 0;
  derivs_unset = 0;
  *par_b_pcmp = par_b;
  *move_origin_x_pcmp = move_origin_x;

}


//BK: Initializes u, v derivatives and defines other functions (NonLinequations, LinEquations, MakeInitialGuess)
void TwoPuncturesSolver_ForceResolve(CCTK_INT resolve)
{
	DECLARE_CCTK_PARAMETERS;
	int nvar_psi = 1, nvar_Aij = 4, ncomp_Aij = 6;
	int n1 = npoints_A, n2 = npoints_B, n3 = npoints_phi;

	int ntotal_psi  = n1 * n2 * n3 * nvar_psi; 
	int ntotal_Aij = n1 * n2 * n3 * nvar_Aij;
	int ntotal_AijAij = n1 * n2 * n3 * 1;
	int success;             /* return value of Newton solver */
  if (  (F == NULL || derivs_unset==1 || resolve == 1) && bypass_solve==0 ) {
    /* Solve only when called for the first time and this set
     * of parameters has not already been computed */
    // free old values if they exists
    if(F != NULL) // assume u,v,w are allocated at the same time!
    {
	free_dvector (F, 0, ntotal_psi - 1);
        free_derivs (&u, ntotal_psi);
	free_derivs (&v, ntotal_psi);
    }

    F = dvector (0, ntotal_psi - 1);
    //BK: Assign vector arrays for derivatives of u with ntotal_psi elements
    allocate_derivs (&u, ntotal_psi);
    allocate_derivs (&v, ntotal_psi);

    CCTK_INFO ("Solving puncture equation (psi)");
    NonLinEquations = NonLinEquations_Psi;
    LinEquations = LinEquations_Psi;
    makeInitialGuess = makeInitialGuess_Psi;
    success = Newton (nvar_psi, n1, n2, n3, v, Newton_tol, Newton_maxit, 0.);
    if(Newton_abort_if_unsuccessful && !success) {
      CCTK_WARN(0, "Newton solver failed to converge when solving for psi");
    }
	if(use_full_transform) {
		v_coeffs.d0 = dvector(0, ntotal_psi - 1);
		TP_FullTransform(v.d0, v_coeffs.d0, nvar_psi, n1, n2, n3);
	}
	derivs_unset = 0;
  }
}
