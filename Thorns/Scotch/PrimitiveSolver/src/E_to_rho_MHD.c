#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "EOS_Base.h"
#include "matterlibs.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#define DEBUG

static inline double SQR (const double x)
{
  return x*x;
}
static inline double QUAD (const double x)
{
  const double y = x*x;
  return y*y;
}

// used to go from energy density to matter density
struct param_t {
  int solver_dim, eos_handle;
  double E, S2, B2, SB;
};

static int ej_to_rhoW_f  (const gsl_vector * x, void * params, gsl_vector * f );
static int ej_to_rhoW_df (const gsl_vector * x, void * params, gsl_matrix * J );
static int ej_to_rhoW_fdf(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J);
static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s );

int matterlibs_E_to_rho_MHD(int eos_handle, double dE, double dJx, double dJy, 
    double dJz, double dBx, double dBy, double dBz, double psi, double *drho, double *dlorentz)
{
  DECLARE_CCTK_PARAMETERS;

  // If vacuum, this is trivial. Don't bother with the rest of this function.
  if(dE == 0.)
  {
    *drho = 0.;
    *dlorentz = 1.;
    return 0;
  }

  // Set up parameters for solve
  int ierr, solver_dim;
  struct param_t params;
  params.B2 = QUAD(psi) * (SQR(dBx) + SQR(dBy) + SQR(dBz));
  params.S2 = QUAD(psi) * (SQR(dJx) + SQR(dJy) + SQR(dJz));
  params.SB = dJx*dBx + dJy*dBy + dJz*dBz;
  params.E = dE;
  params.eos_handle = eos_handle;
  solver_dim = params.S2 > 0 ? 2 : 1;
  params.solver_dim = solver_dim;

  // If S2 != 0, we can't simply eliminate w_lorentz anymore in favor of rho*h.  This is
  // because we now use the S2 equation to eliminate the unknown v.B (though this should 
  // vanish in MHD).  For (v.B != 0), we need another relation. 
  if ( params.S2 > 0 && params.SB != 0 ) {
     CCTK_WARN(0,"Currently MatterLibs can't handle non-stationary MHD initial data without v.B=0 assumption. Aborting.");
  }

  // Set up initial guesses (in passed arguments)
  double rho_init = *drho;
  double W2_init = SQR(*dlorentz); 
  double rhoW_init[2] =  { rho_init, W2_init };
  gsl_vector_view x = gsl_vector_view_array( rhoW_init, solver_dim );

  /* Solver type info -- hardcode Newton */
  const gsl_multiroot_fdfsolver_type *my_solver_choice;
  my_solver_choice = gsl_multiroot_fdfsolver_newton;

  /* Set up solver -- hardcoded derivative type */
  size_t iter = 0;
  int status;
  gsl_multiroot_function_fdf ej_to_rhoW_system = {ej_to_rhoW_f, ej_to_rhoW_df,
      ej_to_rhoW_fdf, solver_dim, &params};
  gsl_multiroot_fdfsolver *s = gsl_multiroot_fdfsolver_alloc(my_solver_choice, solver_dim);
  gsl_multiroot_fdfsolver_set ( s, &ej_to_rhoW_system, &x.vector ); 
  if ( verbose > 1 ) 
     print_state(iter,s);

  /* Iterate manually */
  do {
     iter++;
     gsl_vector_memcpy(&x.vector, s->x);
     status = gsl_multiroot_fdfsolver_iterate(s);

     // Check that range 
     int reset = 0;
     if ( gsl_vector_get(s->x, 0) < 0 ) { // rho < 0 
        gsl_vector_set( s->x, 0, fabs(gsl_vector_get (s->x, 0)) );
        reset = 1;
     }
     if ( solver_dim == 2 ) {
        // First try just requiring W2>0. Imposing physical W during solve may be too 
        // restrictive in finding solution
        if ( gsl_vector_get(s->x, 1) < 0 ) {
           gsl_vector_set( s->x, 1, fabs(1-gsl_vector_get(s->x,1)) );
           reset = 1;
        }
     }
     if(reset>0) {
        gsl_vector_memcpy(&x.vector, s->x);
        gsl_multiroot_fdfsolver_set ( s, &ej_to_rhoW_system, &x.vector );
     }

     if ( verbose > 1 )
        print_state ( iter, s );

     if (status) break; /* Is the solver stuck or failed? */

     /* Test solution */
     status = gsl_multiroot_test_delta ( s->dx, s->x, abs_err, rel_err);
 
  } while ( status == GSL_CONTINUE && iter < (size_t)max_iter );

  if ( status == GSL_SUCCESS ) {
     *drho = gsl_vector_get( s->x, 0);
     if ( solver_dim == 2 ) {
        *dlorentz = sqrt( gsl_vector_get( s->x, 1) );
     } else {
        *dlorentz = 1.;
     }
  }

  gsl_multiroot_fdfsolver_free(s);

  return status != GSL_SUCCESS;

}

static void print_state ( size_t iter, gsl_multiroot_fdfsolver * s )
{
  size_t sdim = s->x->size;
  CCTK_VInfo(CCTK_THORNSTRING,"Solver state in E_to_rho_MHD: iter=%d (rho,W2)=(% e % e ) "
     "f(rho,W2)=(% e % e) J(rho,W2) = ((% e,% e),(% e % e))" , (int)iter,
     gsl_vector_get(s->x,0), sdim < 2 ? -1 : gsl_vector_get(s->x,1),
     gsl_vector_get(s->f,0), sdim < 2 ? -1 : gsl_vector_get(s->f,1),
     gsl_matrix_get(s->J,0,0), sdim < 2 ? -1 : gsl_matrix_get( s->J, 0,1 ),
     sdim < 2 ? -1 : gsl_matrix_get(s->J,1,0), sdim < 2 ? -1 : gsl_matrix_get( s->J, 1,1 ));
}

// The system to solve here is more complicated.
// Given the following definitions
//   S_j = ( rho h W2 + B2 )v_j - (v.B) B_j
//   h   = 1 + eps + P/rho   (hydro specific-enthalpy only)
//   E   = rho h W^2 - P + (1-1/(2 W2)) B2 - (v.B)^2/(2 W2)
//   W2  = 1 / ( 1 - v^2 )   -- Use as independent variable
// And assume a polytropic/gamma-law EOS.
//
// If S2=0, W=1 and solve for rho
//    f = ebar + P - rho h W2 + B2/2
//
// If S2>0 (2D solve) 

static int ej_to_rhoW_f  (const gsl_vector * x, void * params, gsl_vector * f )
{
  // Don't duplicate too much code
  return ej_to_rhoW_fdf(x, params, f, NULL);
}

static int ej_to_rhoW_df (const gsl_vector * x, void * params, gsl_matrix * J )
{
  // Don't duplicate too much code
  return ej_to_rhoW_fdf(x, params, NULL, J);
}

static int ej_to_rhoW_fdf(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J)
{
  DECLARE_CCTK_PARAMETERS;

  // Process incoming information
  const double rho = gsl_vector_get(x,0);
  const double w2  = x->size < 2 ? 1 : gsl_vector_get(x,1);
  struct param_t *pars = (struct param_t *)params;

  // Shorthands
  const double P     = EOS_Pressure(pars->eos_handle, rho, matterlibs_poison);
  const double eps   = EOS_SpecificIntEnergy(pars->eos_handle, rho, matterlibs_poison);
  const double enthalpy = rho*(1+eps)+P; // rho*h
  const double B2 = (pars->B2); 

  const double dPdRho = EOS_DPressByDRho(pars->eos_handle, rho, matterlibs_poison);
  const double dPdEps = EOS_DPressByDEps(pars->eos_handle, rho, matterlibs_poison);
  const double dEpsdRho = dPdRho/dPdEps;
  const double denthalpydRho = 1 + eps + rho*dEpsdRho + dPdRho;
  const double invw2 = 1./w2;

  if ( x->size < 2 ) {

      const double f_ebar = rho*(1+eps) + 0.5*B2 - pars->E;
      const double df_ebar_drho = 1 + eps + rho*dEpsdRho;
      gsl_vector_set ( f, 0, f_ebar );
      gsl_matrix_set ( J, 0,0, df_ebar_drho );

  } else {

     const double f_ebar = enthalpy*w2 + (1 - 0.5*invw2)*B2 - pars->E - P;
     const double df_ebar_drho = w2*denthalpydRho + (w2-1)*dPdRho;
     const double df_ebar_dW2 = enthalpy + 0.5*SQR(invw2)*B2;

     const double f_S2   = (2*enthalpy*B2 - SQR(enthalpy))*SQR(w2) 
                   + (SQR(B2)-2*enthalpy*B2-(pars->S2))*w2 - SQR(B2);
     const double df_S2_drho = 2*denthalpydRho*( B2*SQR(w2)- SQR(w2)*enthalpy - B2*w2 );
     const double df_S2_dW2 = 2*(2*enthalpy*B2-SQR(enthalpy))*w2 + SQR(B2)-2*enthalpy*B2-(pars->S2);

     // Set the respective gsl vectors/matrices
     gsl_vector_set ( f, 0, f_ebar );
     gsl_vector_set ( f, 1, f_S2 );
     gsl_matrix_set ( J, 0,0, df_ebar_drho );
     gsl_matrix_set ( J, 0,1, df_ebar_dW2 );
     gsl_matrix_set ( J, 1,0, df_S2_drho );
     gsl_matrix_set ( J, 1,1, df_S2_dW2 );

  }
     
  return GSL_SUCCESS;

}
