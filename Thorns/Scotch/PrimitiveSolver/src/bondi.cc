#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <new>

#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "util_String.h"
#include "util_Table.h"

#define INITVALUE (42)
#define M_PI     3.14159265358979323846  /* pi */
#define PSIMAX 10.
#define MAXV 0.98
#define TINY 0.000001

#define MAX_NPTS 2500
#define MAX_SOLNS 20
#define SOLN_TO_POPULATE 0

#define DIM 1
#define INTERP_NPOINTS 1
#define NUM_INPUT_ARRAYS 2
#define NUM_OUTPUT_ARRAYS 2

/***************
 ** FUNCTIONS **
 ***************/

/* Scheduled routines */
extern "C" void MatterLibs_Bondi_Solve(CCTK_ARGUMENTS);
extern "C" void MatterLibs_IsoBondi_Populate(CCTK_ARGUMENTS);
extern "C" void MatterLibs_Bondi_WriteSolns(CCTK_ARGUMENTS);
extern "C" void MatterLibs_Bondi_Cleanup(CCTK_ARGUMENTS);

/* Externally visible functions */
extern "C"
void calcBondiRho( CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *psi, CCTK_REAL *rho, CCTK_REAL *eps, 
                   CCTK_REAL *press, CCTK_REAL *vel_x, CCTK_REAL *vel_y, CCTK_REAL *vel_z,
                   CCTK_REAL *B_x, CCTK_REAL *B_y, CCTK_REAL *B_z );

/* Inline helper functions */
static inline CCTK_REAL SQR (const CCTK_REAL x)
{
  return ((x)*(x));
}
static inline CCTK_REAL QUAD (const CCTK_REAL x)
{
  const CCTK_REAL y = ((x)*(x));
  return ((y)*(y));
}


/*********************************
 ** GLOBAL BONDI SOLUTION CLASS **
 *********************************/

/* Begin definition of a Bondi solution class (Search End for end ...) */
class BondiSolution {

  private:
   /* Boolean for whether solution already exists */
   bool have_solution;

  public:
   /* Storage for solution-specific parameters */
   CCTK_REAL Mdot, Qdot;                       // Conserved quantities
   CCTK_REAL rSonic, rhoSonic, hSonic, aSonic; // Sonic-pt information 
   CCTK_REAL eos_K, M, bmag;                   // Extra parameters
   CCTK_INT  npts, iSonic;                     // Size of solution vectors, location of sonic pt
   CCTK_REAL rmin, rmax, logrmin, dlogr;       // r-axis info
   CCTK_REAL bh_x, bh_y, bh_z;                 // Black hole location

   /* Information for expansion around some point */
   CCTK_REAL rminus, rho_minus, u_minus, uIso_minus, drhodr, dudr;

   /* The solution itself */
   std::vector<CCTK_REAL> rSch;
   std::vector<CCTK_REAL> log_rSch;
   std::vector<CCTK_REAL> rho;
   std::vector<CCTK_REAL> u_r;
   std::vector<CCTK_REAL> int_energy;

   /* Constructor with parameter index as input */
   BondiSolution( int param_idx = SOLN_TO_POPULATE );

   /* Function for getting the Bondi solution at a given r(Schwarzschild) */
   CCTK_INT find_Bondi_Solution();
   CCTK_INT get_Bondi_at_r ( CCTK_REAL r, CCTK_REAL * rho_at_r, CCTK_REAL * u_at_r);
   void print_to_file( char * filename );

};


/************************************
 * Constructor, General             *
 *   Initiates parameters, vectors, * 
 *   and r-axes according to params *
 ************************************/

BondiSolution::BondiSolution( int param_idx ) {

  DECLARE_CCTK_PARAMETERS;

  have_solution=0;

  /******************************************************** 
   * Process parameters, storing anything useful in class *
   ********************************************************/
  M = bondi_central_mass[param_idx];
  Mdot = bondi_mdot[param_idx];

  rSonic = bondi_sonicpt_radius[param_idx];
  rmin = M * bondi_rmin[param_idx]; 
  rmax = M * bondi_rmax[param_idx];
  npts = bondi_npts[param_idx]; 

  bh_x = bondi_bh_x0[param_idx];
  bh_y = bondi_bh_y0[param_idx];
  bh_z = bondi_bh_z0[param_idx];

  bmag = bondi_bmag[param_idx];

  /* Sound speed squared at sonic point */
  CCTK_REAL aSonic2, uSonic2, uSonic;
  aSonic2 = M / ( 2*rSonic - 3*M );
  if ( aSonic2 > ( eos_gamma - 1. ) ) {
     CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                "Problem with Bondi solution %d. aSonic^2 > gam-1.", param_idx);
  }
  aSonic = sqrt(aSonic2); 
  uSonic2 = M / (2*rSonic);
  uSonic = sqrt(uSonic2);
    
  /* EoS K from solution at sonic point */
  hSonic = 1./( 1. - aSonic2/(eos_gamma-1.));
  rhoSonic = Mdot / ( 4. * M_PI * uSonic * SQR(rSonic) );
  eos_K = hSonic * aSonic2 * pow( rhoSonic, 1.-eos_gamma ) / eos_gamma;
  Qdot = SQR(hSonic) * ( 1. - 3.*uSonic2 );

  CCTK_VInfo(CCTK_THORNSTRING,"Bondi Solution %d:\t eos_K=%g, Mdot=%g, Qdot=%g",
             param_idx, eos_K, Mdot, Qdot );
  CCTK_VInfo(CCTK_THORNSTRING,"      at Sonic Pt: r=%g, rho=%g, a=%g, h=%g, u=%g",
             rSonic, rhoSonic, aSonic, hSonic, uSonic );
  

  /* Actually solve the Bondi solution */
  logrmin = log10(rmin);
  dlogr   = ( log10(rmax) - logrmin) / ( 1.*(npts-1) );
  
  /* Resize vectors to proper size */
  rSch.assign(npts,INITVALUE);
  log_rSch.assign(npts,INITVALUE);
  rho.assign(npts,INITVALUE);
  u_r.assign(npts,INITVALUE);
  int_energy.assign(npts,INITVALUE);

  /* Create r_Schwarzschild axis vectors {rSch, log_rSch} */
  iSonic = 0; 
  CCTK_REAL dist_from_rSonic = 1e30; 
  for ( int rn=0; rn < rSch.size(); rn++ ) {
      log_rSch[rn] = logrmin + dlogr*(rn-1);
      rSch[rn] = pow( 10, log_rSch[rn] );
      if ( fabs(rSch[rn] - rSonic) < dist_from_rSonic ) {
         dist_from_rSonic = fabs(rSch[rn] - rSonic);
         iSonic = rn;
      }
  }

  CCTK_INT ierr = find_Bondi_Solution();
  if ( ierr == 0 ) {
     have_solution = 1;
  } else {
     CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
        "find_Bondi_Solution(%d) was unsuccessful at %d/%d points.",param_idx,ierr,npts);
  }

  /* Store rho and its derivative at rSch ~ 2.25 M */
  CCTK_REAL rplus, rho_plus, u_plus, uIso_plus;
  CCTK_INT jminus, jplus;
  
  rminus = 2.25 * M;
  jminus = floor( ( log10(rminus) - logrmin ) / dlogr );
  rho_minus = rho[jminus];
  u_minus   = u_r[jminus];
  ierr = get_Bondi_at_r ( rminus, &rho_minus, &u_minus );
  uIso_minus = (4./3.)*u_minus;
  
  rplus  = 0.25 * SQR(3.02) * M/1.01; // ~2.2575
  jplus = floor( ( log10(rplus) - logrmin) / dlogr );
  rho_plus = rho[jplus];
  u_plus   = u_r[jplus];
  ierr = get_Bondi_at_r ( rplus, &rho_plus, &u_plus );
  uIso_plus = u_minus / (1. - 1./(2.*M)) / (1. + 1./(2.*M));

  drhodr = 100*( rho_plus - rho_minus )/M;
  dudr   = 100*( uIso_plus - uIso_minus ); 

};

/*************************************************************
 * BondiSolution::find_Bondi_Solution                        *
 *   Fill in actual Bondi solution: rho(r), u(r)             *
 *************************************************************/ 
CCTK_INT BondiSolution::find_Bondi_Solution() {

  DECLARE_CCTK_PARAMETERS;
  CCTK_REAL uSonic = sqrt( M/(2.*rSonic) );
  CCTK_INT  ierr=0;

  /* Start solving from sonic point outwards,                    */
  /*   We want the larger rho solution outside the sonic points */
  CCTK_REAL rhoL = rhoSonic;
  CCTK_REAL uL = uSonic;
  for ( int i = iSonic; i < rSch.size(); i++ ) { 
      rhoL = 1.02*rhoL; /* Approach from above to get higher branch */
      ierr += get_Bondi_at_r( rSch[i], &rhoL, &uL );
      if ( rhoL < rho_abs_min ) {
         rhoL = rho_abs_min;
         uL   = eos_K * pow( rhoL, eos_gamma ) / ( eos_gamma - 1. );
      } 
      rho[i] = rhoL;
      u_r[i]   = uL;
      int_energy[i]   = eos_K * eos_gamma * pow( rhoL, eos_gamma ) / ( eos_gamma - 1.);
  }

  /* Solve from sonic point inwards                              */
  /*   We want the smaller rho solution outside the sonic points */
  rhoL = rhoSonic;
  uL = uSonic;
  for ( int i=iSonic-1; i>=0; i-- ) {
      rhoL = 0.98*rhoL; /* Approach from below to get lower branch */
      ierr += get_Bondi_at_r( rSch[i], &rhoL, &uL );
      if ( rhoL < rho_abs_min ) {
         rhoL = rho_abs_min;
         uL   = eos_K * pow( rhoL, eos_gamma ) / ( eos_gamma - 1. );
      } 
      rho[i] = rhoL;
      u_r[i]   = uL;
      int_energy[i]   = eos_K * eos_gamma * pow( rhoL, eos_gamma ) / ( eos_gamma - 1.);
  }

  return ierr;

};

/*************************************************************
 * BondiSolution::get_Bondi_at_r                             *
 *   Given an initial guess of rho, u at r_Schwarzschild,    *
 *   return the iteratively solved solution Bondi.           *
 *************************************************************/ 
CCTK_INT BondiSolution::get_Bondi_at_r( CCTK_REAL r, CCTK_REAL * rho_at_r, CCTK_REAL * u_at_r) { 

   DECLARE_CCTK_PARAMETERS;

   /* Store initial guesses ... */
   CCTK_REAL rho_0, u_0, rhoL, uL;
   rho_0 = *rho_at_r;
   u_0 = *u_at_r;

   CCTK_INT errout;
   /* Iteratively solve for initial condition and if too close to edges of solution */
   if ( !have_solution || !bondi_interp_after_solution || 
       (have_solution && bondi_interp_after_solution && 
         (log10(r)<log_rSch[5] || log10(r)>log10(rmax)-2*dlogr) )  ) {
   
      /* Newton-Raphson Parameters */
      const CCTK_INT  num_rounds_max = 1;
      const CCTK_REAL reset_rho = rho_abs_min;
   
      /* Case where we had a bad initial rho. */
      if ( rho_0 < 0. ) {
         if ( r < 1.1*rSonic && r > 0.9*rSonic ) {
            r = rSonic;
         } else {
            CCTK_REAL u_r;
            if ( r < rSonic ) {
               u_r = 1./sqrt(r);
            } else {
               u_r = 0.5*pow(r, -1.5);
            }
            rho_0 = Mdot/(4.*M_PI*SQR(r)*u_r);
         }
      }
   
      /* ... using Newton-Raphson, ... */
      /* Note: there are 2 branches of solutions we might end up on. */
      CCTK_REAL rho_old;
      CCTK_REAL rhoL = rho_0;
      CCTK_REAL uL = u_0;
      CCTK_INT  itn_count = 0;
      CCTK_INT  num_rounds = 0;
      CCTK_INT  itn_count_extra = 0;
      CCTK_INT  satisfied = 0;
      CCTK_INT  keep_going = 1;
      CCTK_INT  doing_extra = 0;
      CCTK_REAL err = 1e3;
      CCTK_REAL QdotL, errQdot, MdotL, errMdot;
   
   
      while ( keep_going==1 ) {
   
         /* Travel along Mdot=const */
         CCTK_REAL h, dhdrho, dudrho, uterm, resid;
         dhdrho = eos_K * eos_gamma * pow( rhoL, eos_gamma-2. );
         //h = 1. + eos_K * eos_gamma * pow( rhoL, eos_gamma-1.)/(eos_gamma - 1.);
         h = 1. + dhdrho*rhoL/(eos_gamma - 1.);
         uL = Mdot / ( 4. * M_PI * SQR(r) * rhoL );
         dudrho = -uL/rhoL;
         uterm = 1. - 2.*M/r + SQR(uL);
   
         CCTK_REAL f, df, drho;
         f = SQR(h) * ( 1. - 2.*M/r + SQR(uL) ) - Qdot;
         df = 2.*h*( dhdrho*uterm + h*uL*dudrho );
         drho = -f/df; 
   
   
         /* Save & Update rho */
         rho_old = rhoL;
         rhoL += drho;
   
         /* Error: Absolute error in Qdot(rho,u) */
         QdotL = SQR(h)*(1.-2.*M/r + SQR(uL)); 
         errQdot = (QdotL - Qdot)/Qdot;
         MdotL = 4.*M_PI*SQR(r)*rhoL*uL;
         errMdot = (MdotL-Mdot)/Mdot;
   
         /* Error and convergence: relative change in rho over iteration */
         if ( rhoL <= 0 ) {
            err = fabs(drho);
            rhoL = reset_rho;
            if ( verbose >= 2 ) CCTK_VInfo(CCTK_THORNSTRING, "   Reset rho at iteration %d-%d.",num_rounds, itn_count);
         } else {
            err = fabs(drho/rhoL);
         }
   
         /* Relative error satisfied. Extra iterations & extra rounds. */
         if ( ( fabs(err) < rel_tolerance || fabs(f) < abs_tolerance ) && doing_extra==0 && extra_iterations>0  ) {
            satisfied = 1;
            doing_extra = 1;
         }
         if ( doing_extra == 1 ) itn_count_extra++;
   
         /* Reset/check keep_going if we're to head into (are in) extra iterations */
         if ( (satisfied==1 && doing_extra==0) || ( itn_count_extra > extra_iterations) ) { 
            keep_going=0;
         }
   
         /* Start a new round */
         if ( itn_count+1>=max_iterations && fabs(err) > rel_tolerance && num_rounds<num_rounds_max ) {
            num_rounds++;
            itn_count=-1;
            rhoL=1.05*rhoL; /* tweak rho to hopefully get better convergence */
         }
   
         /* Generic fall out if itn & group maxes are met */
         if ( (fabs(err) > rel_tolerance || fabs(f) > abs_tolerance) && num_rounds+1>=num_rounds_max && itn_count+1>=max_iterations ) {
            keep_going=0;
         }
   
         itn_count++;
   
         /* Very verbose */
         if ( verbose >= 3 ) {
            CCTK_VInfo(CCTK_THORNSTRING, "     Itn (%d-%d): (rho,u)=(%g,%g)->(%g,%g), drho=%g, (f,df)=(%g,%g), err=%g; delta(Mdot,Qdot)=(%g,%g)",
                num_rounds, itn_count, rho_old, uL, rhoL, Mdot/(4.*M_PI*SQR(r)*rhoL), drho, f, df, err, errMdot, errQdot);
         }
      }
   
      if ( satisfied == 0 ) {
   
         CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
             "Failed to find Bondi solution for r=%g, EoS=(%g,%g), Mdot=%g, Qdot=%g, rSonic=%g, (rho,u)=(%g,%g)->(%g,%g) with err=%g; delta(Mdot,Qdot)=(%g,%g).",
             r, eos_K, eos_gamma, Mdot, Qdot, rSonic, rho_0, u_0, rhoL, uL, err, errMdot, errQdot); 
   
      }
   
   
      *rho_at_r = rhoL;
      *u_at_r   = Mdot / ( 4. * M_PI * SQR(r) * rhoL );
   
      if ( verbose >=2 ) {
         CCTK_VInfo(CCTK_THORNSTRING," get_Bondi_at_r(r=%g,%g,%g) returning (rho,u)=(%g,%g) with err %g, delta(Mdot,Qdot)=(%g,%g). Success=%d",
            r, rho_0, u_0, *rho_at_r, *u_at_r, err, errMdot, errQdot, satisfied );
      }
   
      errout = ( satisfied > 0 ? 0 : 1 ); 

  } else { 

    /* Get interpolated values over existing 1D solution, uniform in logr space */
    CCTK_INT operand_indices[2] = { 0, 1 };
    CCTK_INT opcodes[2] = { 0, 0 }; // Direct interpolation. No derivative
    const int operator_handle = CCTK_InterpHandle(interpolator_name);
    if (operator_handle < 0)
       CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
               "couldn't find interpolator \"%s\"!",
               interpolator_name);

    const int param_table_handle = Util_TableCreateFromString(interpolator_pars);
    if (param_table_handle < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "bad interpolator parameter(s) \"%s\"!",
                 interpolator_pars);
    }

    Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                          operand_indices, "operand_indices");
  
    Util_TableSetIntArray(param_table_handle, NUM_OUTPUT_ARRAYS,
                          opcodes, "opcodes");


    const CCTK_REAL logr = log10(r);
    const void *interp_coords[DIM] = { (const void *) &logr };

    const CCTK_REAL coord_origin[DIM] = { logrmin };
    const CCTK_REAL coord_delta[DIM] = { dlogr };
    const CCTK_INT input_array_dims[DIM] = { npts };

    const CCTK_INT input_array_types[NUM_INPUT_ARRAYS] = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };
    const void * input_arrays[NUM_INPUT_ARRAYS] = { (const void *) &rho,
                                                    (const void *) &u_r };

    CCTK_INT output_array_types[NUM_OUTPUT_ARRAYS] = { CCTK_VARIABLE_REAL, CCTK_VARIABLE_REAL };
    void * output_arrays[NUM_OUTPUT_ARRAYS] = { (void *) &rhoL, 
                                                (void *) &uL };

    errout = CCTK_InterpLocalUniform(DIM,
                 operator_handle, param_table_handle,
                 coord_origin, coord_delta,
                 INTERP_NPOINTS, CCTK_VARIABLE_REAL, interp_coords,
                 NUM_INPUT_ARRAYS,
                    input_array_dims, input_array_types, input_arrays,
                 NUM_OUTPUT_ARRAYS,
                    output_array_types, output_arrays);

    Util_TableDestroy(param_table_handle);

    *rho_at_r = rhoL;
    *u_at_r   = uL;

  }

  return errout;

}

void BondiSolution::print_to_file( char * filename )
{

   DECLARE_CCTK_PARAMETERS;

   char out_format[5];
   sprintf(out_format, ".19g");

   size_t length_written; 
   FILE *file;

   file = fopen (filename,"w");
   if (!file) {
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
          "Could not open bondi solution output file '%s'", filename); 
      return;
   }
       
   /* Header information */
   fprintf(file, "# Bondi Solution, Numerical\n");
   fprintf(file, "#   BH of M = %g at (%g,%g,%g)\n", M, bh_x, bh_y, bh_z );
   fprintf(file, "#   Mdot = %g, r_sonic = %g\n", Mdot, rSonic ); 
   fprintf(file, "#   %d points, r = [ %g : %g ]\n", npts, rmin, rmax); 
   fprintf(file, "# Data headers:\n");
   fprintf(file, "# 1:r(Iso) 2:r(Sch) 3:rho 4:u_Sch 5:v_Sch 6:v_Iso 7:delta(Mdot) 8:delta(Qdot)\n");
      

   /* Print data */
   char format_str[2048];
   size_t len_written = snprintf (format_str, sizeof(format_str)/sizeof(format_str[0]),
                           "%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\t%%%s\n", out_format, out_format, out_format,
                           out_format, out_format, out_format, out_format, out_format);
   assert( len_written < sizeof(format_str)/sizeof(format_str[0])); 
   for ( int rIdx=0; rIdx<npts; rIdx++ ) {

       /* Convert some quantities from Schwarzschild to Isotropic coordinates */
       CCTK_REAL rSchL = rSch[rIdx];
       CCTK_REAL uSch = u_r[rIdx];

       /* Derive other quantities for output */
       CCTK_REAL rIso, psi, vSch, vIso, uIso, invwIso;
       if ( rSchL >=2.*M ) {
          rIso = 0.5*( sqrt( SQR(rSchL) - 2.*M*rSchL) + rSchL - M );
          psi = 1. + 0.5*M/rIso;
       } else {
          /* FIX */
          continue;
          //rIso = 1.e-5;
          //psi = ;
       }
       vSch = -uSch / sqrt( 1. - 2.*M/rSchL + SQR(uSch) );
       uIso = uSch / ( 1. - M/(2.*rIso) ) / ( 1. + M/(2.*rIso) );
       //uIso = uSch / (1. - SQR( M/(2.*rIso) ) );
       invwIso = 1./sqrt( 1. + QUAD(psi)*SQR(uSch) );
       vIso = uIso * invwIso; 

       /* Fractional error, Mdot */
       CCTK_REAL errMdot = (4.*M_PI*SQR(rSchL)*rho[rIdx]*uSch - Mdot)/Mdot;

       /* Fractional error, Qdot */
       CCTK_REAL QdotL, errQdot;
       QdotL = SQR( 1. + eos_gamma*eos_K*pow(rho[rIdx],eos_gamma-1.)/(eos_gamma-1.))*( 1. - 2.*M/rSchL + SQR(uSch) );
       errQdot = ( QdotL - Qdot ) / Qdot;

       /* Final output */ 
       fprintf(file, format_str, rIso, rSchL, rho[rIdx], uSch, vSch, vIso, errMdot, errQdot ); 

   }
   fprintf(file,"\n");
   fclose(file);

}
// End BondiSolution class

/****************************************
 ** Begin Externally visible Functions **
 ****************************************/

/* Static, global pointers to bondi solutions, to be allocated later */
static BondiSolution * Bondi_Soln[MAX_SOLNS]; 

/**************************************************************
 * MatterLibs_Bondi_Solve                                     *
 *   Populate static structures containing Bondi solutions    * 
 *************************************************************/ 
extern "C"
void MatterLibs_Bondi_Solve ( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  /* Initialize pointers to class */
  for ( int nsol=0; nsol<MAX_SOLNS; nsol++ ) {
      Bondi_Soln[nsol] = NULL;
  }

  /* Solve all solutions */
  for( int nsol=0; nsol<bondi_num_solutions; nsol++ ) {
     /* Construct class instance & a ref to the storage */
     Bondi_Soln[nsol] = new BondiSolution(nsol);
  }

}

/*******************************************************************
 * calcBondiRho                                                    *
 *   Return primitives for a Bondi solution at a specific location *
 *******************************************************************/ 
extern "C"
void calcBondiRho( CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL *psi, CCTK_REAL *rho, CCTK_REAL *eps, 
                   CCTK_REAL *press, CCTK_REAL *vel_x, CCTK_REAL *vel_y, CCTK_REAL *vel_z,
                   CCTK_REAL *B_x, CCTK_REAL *B_y, CCTK_REAL *B_z ) {

  DECLARE_CCTK_PARAMETERS;

  BondiSolution * BondiN = Bondi_Soln[SOLN_TO_POPULATE]; 
  assert(BondiN);

  /* We're here */
  CCTK_REAL xL = x-BondiN->bh_x;
  CCTK_REAL yL = y-BondiN->bh_y;
  CCTK_REAL zL = z-BondiN->bh_z;
  CCTK_REAL rIso = sqrt( SQR(xL) + SQR(yL) + SQR(zL) + SQR(r_epsilon) );
  CCTK_REAL rbyM = rIso/BondiN->M;

  /* Calculate the equivalent schwarzschild radius from the closest 
     solution point where r < r_i, then find rho,uIso there           */
  CCTK_REAL uIso, rhoL;
  //if ( rSch > BondiN->rminus ) {
  //if ( rIso > 0.7*BondiN->M ) {
  if ( rbyM > 1.-TINY ) {
     /* Far from BH, don't need any tricks, just a decent guess would be nice 
        if we are not interpolating */
     CCTK_INT rIdx;
     CCTK_REAL rSch, uSch;
     rSch = 0.25 * SQR(2.*rIso + BondiN->M)/rIso;
     rIdx = std::min<int>( floor( 0.5 + ( log10(rSch) - BondiN->logrmin ) / BondiN->dlogr ), BondiN->npts-1 );
     if ( rIdx >= BondiN->npts ) rIdx = BondiN->npts-1;
     //rhoL = BondiN->rho[rIdx];
     rhoL = BondiN->rho[rIdx] + (BondiN->rho[rIdx+1]-BondiN->rho[rIdx])*
              (rSch-BondiN->rSch[rIdx])/(BondiN->rSch[rIdx+1]-BondiN->rSch[rIdx]);
     uSch = BondiN->u_r[rIdx];

     /* Populate the (rho,u) guesses with the real solution */
     CCTK_INT ierr = BondiN->get_Bondi_at_r ( rSch, &rhoL, &uSch );
     uIso = uSch / ( 1. - BondiN->M/(2.*rIso) ) / ( 1. + BondiN->M/(2.*rIso) );
     //uIso = uSch / ( 1. - SQR( BondiN->M/(2.*rIso) ) );

  } else {

     CCTK_REAL umin, a_uIso, b_uIso;
     if ( rbyM > 0.5*(1.-TINY) ) {
        /* Just outside horizon, linearize rho */
        rhoL = BondiN->rho_minus + BondiN->drhodr*rbyM*(rIso-BondiN->M);
     } else {
        /* Junk rho within the horizon, to be excised later. Must be smooth */
        rhoL = ( BondiN->rho_minus - BondiN->drhodr*BondiN->M/4.)*0.5*( 1.-cos(2.*M_PI*rIso/BondiN->M) );
     }
      
     a_uIso = 1.5*BondiN->uIso_minus - 0.5*BondiN->dudr;
     b_uIso = 0.5*( BondiN->dudr - BondiN->uIso_minus );
     uIso = rbyM*( a_uIso + b_uIso*SQR(rbyM) ); 

  }

  if ( rhoL < rho_abs_min ) {
     rhoL = rho_abs_min;
  }

  // Calculate and return full psi. Used for w later in populate, or overwritten in 2P.
  //if ( rIso < 0.5*BondiN->M ) {
  //   *psi = 2.875 - 5.*SQR(rIso/BondiN->M) + 6.*QUAD(rIso/BondiN->M);
  //} else {
     *psi = 1.+ BondiN->M/(2.*rIso);
  //}
  CCTK_REAL invw = 1./sqrt( 1. + QUAD(*psi)*SQR(uIso) );
  CCTK_REAL xhat = xL/fmax(rIso,r_epsilon);
  CCTK_REAL yhat = yL/fmax(rIso,r_epsilon);
  CCTK_REAL zhat = zL/fmax(rIso,r_epsilon);
  CCTK_REAL vxL = -uIso * invw * xhat;
  CCTK_REAL vyL = -uIso * invw * yhat;
  CCTK_REAL vzL = -uIso * invw * zhat;

  /* Check physicality with psi */
  CCTK_REAL speed = QUAD(*psi)*( SQR(vxL) + SQR(vyL) + SQR(vzL) );
  if ( speed > MAXV ) {
      CCTK_REAL limit_v = sqrt( MAXV/speed );
      CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                 "At (%g,%g,%g), rescaling v=(%g,%g,%g) to v=(%g,%g,%g) due to psi4=%g and mag(v)=%g",
                 xL, yL, zL, vxL, vyL, vzL, limit_v*vxL, limit_v*vyL, limit_v*vzL, QUAD(*psi), speed );
      vxL = limit_v*vxL;
      vyL = limit_v*vyL;
      vzL = limit_v*vzL;
  }

  /* Populate input pointers */
  *rho = rhoL;
  *press = BondiN->eos_K * pow( rhoL, eos_gamma );
  *eps = (*press)/( rhoL*(eos_gamma - 1.));
  *vel_x = vxL;
  *vel_y = vyL;
  *vel_z = vzL;

  if ( B_x != NULL ) {
     /* FIX: Cap here as well? Bad at center. */
     *B_x = BondiN->bmag*xhat;
     *B_y = BondiN->bmag*yhat;
     *B_z = BondiN->bmag*zhat;
  } 

}; 

/**************************************************************
 * MatterLibs_IsoBondi_Populate                               *
 *   Populate the hydro & spacetime GFs with a Bondi solution *
 *   in isotropic coordinates. The result is constraint       *
 *   violating if the density of the solution is too high.    *
 *************************************************************/ 
extern "C"
void MatterLibs_IsoBondi_Populate( CCTK_ARGUMENTS ) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  const CCTK_INT set_bondi_lapse = CCTK_EQUALS(initial_lapse,"Bondi");

  /* Check for mhd activation */
  bool have_bfield = 0;
  if ( Bvec != NULL ) {
     have_bfield = 1;
  }
     

  for( int k=0; k<cctk_lsh[2]; k++ ) {
     for( int j=0; j<cctk_lsh[1]; j++ ) {
        for(int i=0; i<cctk_lsh[0]; i++ ) {
 
           CCTK_INT idx=CCTK_GFINDEX3D(cctkGH, i, j, k);

           /* Get Bondi solution */
           BondiSolution * BondiN = Bondi_Soln[SOLN_TO_POPULATE];
           assert(BondiN);

           /* Assume coordinates are isotropic */
           CCTK_REAL xL, yL, zL, rIso;
           xL = x[idx];
           yL = y[idx];
           zL = z[idx];
           rIso = sqrt( SQR(xL-BondiN->bh_x) + SQR(yL-BondiN->bh_y) + SQR(zL-BondiN->bh_z) + SQR(r_epsilon));

           CCTK_REAL psi, psi4, rhoL, epsL, pressL, vxL, vyL, vzL, bxL, byL, bzL;
           calcBondiRho( xL, yL, zL, &psi, &rhoL, &epsL, &pressL, &vxL, &vyL, &vzL, &bxL, &byL, &bzL );

           /* Fill in geometry, smoothing over the puncture as in Faber et al. */
           /* Note: this is not constraint-satisfying                          */
           //if ( rIso < 0.5*BondiN->M ) {
              //psi = 2.875 - 5.*SQR(rIso/BondiN->M) + 6.*QUAD(rIso/BondiN->M);
           //   psi = 2.875 - 5.*pow(rIso/BondiN->M,2.) + 6.*pow(rIso/BondiN->M,4.);
           //} else {
           //   psi = 1.+ BondiN->M/(2.*rIso);
           //}

           /* Overide psi with proper puncture psi */
           psi = fmin( 1. + 0.5*BondiN->M/rIso, PSIMAX );
           psi4 = QUAD(psi);
 
           /* Spacetime GFs. */
           gxx[idx] = psi4;
           gyy[idx] = psi4;
           gzz[idx] = psi4;
           gxy[idx] = 0.;
           gxz[idx] = 0.;
           gyz[idx] = 0.;
           if ( set_bondi_lapse > 0 ) {
              alp[idx] = 1./SQR(psi);
           }
           kxx[idx] = 0.;
           kyy[idx] = 0.;
           kzz[idx] = 0.;
           kxy[idx] = 0.;
           kxz[idx] = 0.;
           kyz[idx] = 0.;

           /* Check physicality with new psi */
          CCTK_REAL speed = psi4*SQR(vxL) + psi4*SQR(vyL) + psi4*SQR(vzL);
          if ( speed > MAXV ) {
              CCTK_REAL limit_v = sqrt( MAXV/speed );
              CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                         "At (%g,%g,%g), rescaling v=(%g,%g,%g) to v=(%g,%g,%g) due to psi4=%g and mag(v)=%g",
                         xL, yL, zL, vxL, vyL, vzL, limit_v*vxL, limit_v*vyL, limit_v*vzL, psi4, speed );
              vxL = limit_v*vxL;
              vyL = limit_v*vyL;
              vzL = limit_v*vzL;
           }

           /* Populate hydro variables */
           rho[idx] = rhoL;
           eps[idx] = epsL;
           press[idx] = pressL;
           velx[idx] = vxL;
           vely[idx] = vyL;
           velz[idx] = vzL;

           if ( have_bfield ) {
              CCTK_INT idBx=CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 0);
              CCTK_INT idBy=CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 1);
              CCTK_INT idBz=CCTK_VECTGFINDEX3D(cctkGH, i, j, k, 2);
              Bvec[idBx] = bxL;
              Bvec[idBy] = byL;
              Bvec[idBz] = bzL;
           }

        /* Check for NaNs here. */
        if ( std::isnan(rhoL) || std::isnan(epsL) || std::isnan(pressL) 
           || std::isnan(vxL) || std::isnan(vyL) || std::isnan(vzL) ) {
           CCTK_VWarn(1,__LINE__,__FILE__,CCTK_THORNSTRING,
                      "ERROR: NaN found in populating grids with Bondi at (%g,%g,%g): (rho,eps,press,vx,vy,vz)=(%g,%g,%g,%g,%g,%g)",
                       xL, yL, zL, rhoL, epsL, pressL, vxL, vyL, vzL); 
           //CCTK_WARN(0, "ERROR: NaN found in populating grids with Bondi.");
        }

        }
     }
  }

  /* Clear static memory here! */
  

}

/*************************************************************
 * MatterLibs_Bondi_WriteSolns                               *
 *   Write all Bondi solutions to file                       *
 *************************************************************/ 
extern "C"
void MatterLibs_Bondi_WriteSolns(CCTK_ARGUMENTS)
{

   DECLARE_CCTK_ARGUMENTS;
   DECLARE_CCTK_PARAMETERS;

   for ( int nsol=0; nsol<bondi_num_solutions; nsol++ ) {

       BondiSolution * BondiN = Bondi_Soln[nsol];
       assert (BondiN);

       char * filename;
       Util_asprintf ( &filename, "%s/Bondi_Solution_%d.asc", out_dir, nsol );
       assert(filename);

       BondiN->print_to_file( filename );
       free(filename);

   }

}

extern "C" void MatterLibs_Bondi_Cleanup(CCTK_ARGUMENTS) {

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for( int nsol=0; nsol<bondi_num_solutions; nsol++ ) {
     /* Free memory from the solutions */
     delete Bondi_Soln[nsol];
  }

}

