#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "EOS_Base.h"
#include "primsolver.h"
#include "steffenson.h"

#define DEBUG

static inline double pow2 (const double x)
{
  return x*x;
}
static inline double pow4 (const double x)
{
  const double y = x*x;
  return y*y;
}

// used to go from energy density to matter density
struct param_t {
  double E, S2;
  int eos_handle;
};
static double e_to_rho_RHS_f(double x, void *params);
static double e_to_rho_RHS_df(double x, void *params);
static void e_to_rho_RHS_fdf(double x, void *params, double *f, double *df);

// this is used a couple of times
static inline double rho_times_h(int eos_handle, double rho)
{
  return rho*(1+EOS_SpecificIntEnergy(eos_handle, rho, primsolver_poison))+EOS_Pressure(eos_handle, rho, primsolver_poison);
}


int prims_E_to_rho(int eos_handle, double dE, double dJx, double dJy, 
    double dJz, double psi, double *drho, double *dlorentz)
{
  DECLARE_CCTK_PARAMETERS;
  int ierr;
  struct param_t params;
  gsl_function_fdf RHS = {e_to_rho_RHS_f, e_to_rho_RHS_df, e_to_rho_RHS_fdf,
                           &params};

  params.E = dE;
  params.S2 = pow4(psi) * (pow2(dJx) + pow2(dJy) + pow2(dJz));
  params.eos_handle = eos_handle;
  
  if(dE == 0.) // vacuum
  {
    *drho = 0.;
    *dlorentz = 1.;
  }
  else
  {
    *drho = steffenson_solve(&RHS, *drho, max_iter, abs_err, rel_err, &ierr);
    if(ierr != 0)
	return ierr;

    *dlorentz = sqrt( 0.5 * ( 1 + 
	sqrt( 1 + 4*params.S2/pow2(rho_times_h(eos_handle, *drho)) )));
	  
    // force the lorentz factor to unity if it is unity up to
    // the error tolerance
    if( *dlorentz < 1. && 
	 GSL_SUCCESS == prims_gsl_root_test_delta(1., *dlorentz, 
	     abs_err, rel_err))
    {
      *dlorentz = 1.;
    }
  }

  return 0;
}

// ok these all follow from:
// S^j = rho h W^2 v^j
// E = rho h W^2 - P
// h = 1 + eps + P/rho
// W^2 = 1 / ( 1 - v^2 )
// they only ever work for a polytropic EOS (I believe)
// specifically they assume that denthalpy/drho = dP/drho/rho
// Polytype:
// P = (Gamma-1) * eps * rho
// P = K * rho^Gamma

static double e_to_rho_RHS_f(double rho, void *params)
{
  struct param_t *p = (struct param_t *)params;
  const double P = EOS_Pressure(p->eos_handle, rho, primsolver_poison);
  const double rho_h = rho_times_h(p->eos_handle, rho);

  return pow2(p->E+P)-(p->E+P)*rho_h - p->S2;
}
static double e_to_rho_RHS_df(double rho, void *params)
{
  struct param_t *p = (struct param_t *)params;
  const double P = EOS_Pressure(p->eos_handle, rho, primsolver_poison);
  const double rho_h = rho_times_h(p->eos_handle, rho);
  const double DPbyDRho = EOS_DPressByDRho(p->eos_handle, rho, primsolver_poison);

  return (p->E+P-rho_h)*DPbyDRho - rho_h/rho*(p->E+P);
}
static void e_to_rho_RHS_fdf(double rho, void *params, double *f, double *df)
{
  struct param_t *p = (struct param_t *)params;
  const double P = EOS_Pressure(p->eos_handle, rho, primsolver_poison);
  const double rho_h = rho_times_h(p->eos_handle, rho);
  const double DPbyDRho = EOS_DPressByDRho(p->eos_handle, rho, primsolver_poison);

  *f  =  pow2(p->E+P)-(p->E+P)*rho_h - p->S2;
  *df = (p->E+P-rho_h)*DPbyDRho - rho_h/rho*(p->E+P);

  return;
}
