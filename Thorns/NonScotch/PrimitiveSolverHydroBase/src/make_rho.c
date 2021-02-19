
#include "matterlibs.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#define SQR(x) ((x)*(x))
#define REPS 1e-8

void MatterLibs_FeedRho( CCTK_ARGUMENTS );
void MatterLibs_SetRho( CCTK_ARGUMENTS );
void MatterLibs_CalculateRho ( CCTK_ARGUMENTS );
void MatterLibs_SetLapse( CCTK_ARGUMENTS );

void matterlibs_calcMatterRho( double x, double y, double z, 
  double *psi, double *rho, double *eps, double *P, double *vel_x, double *vel_y, double *vel_z,
  double *bvec_x, double *bvec_y, double *bvec_z ) 
{

  DECLARE_CCTK_PARAMETERS;

  // Initialize the B-field to either zero or the uniform B-field. Any other B-field 
  // needs to currently subsequently overwrite the zero B-field.
  if ( bvec_x != NULL ) {
     if ( CCTK_EQUALS(rho_profile,"UniformB") ) {
        *bvec_x = b0_uniform[0];
        *bvec_y = b0_uniform[1];
        *bvec_z = b0_uniform[2];
     } else {
        *bvec_x = 0;
        *bvec_y = 0;
        *bvec_z = 0;
     }
  }

  // Choose the function based on the keyword parameter
  if ( CCTK_EQUALS( rho_profile, "Gaussian" ) ) {
	calcGaussRho( x, y, z, rho, vel_x, vel_y, vel_z );
	// Gaussian assumes you specify \hat{rho} and therefore any choice 
	// of the EOS polytrope kappa should depend on the solved rho, not
	// the parameter rho
	*psi = 1.;						// a fudge
	*P = 0;
	*eps = 0;

  } else if ( CCTK_EQUALS( rho_profile, "Constant" ) ) {
	*psi = 1.;						// a fudge
	*rho = calcCloudRho( x, y, z );
	*P = eos_k*pow( *rho, eos_gamma );
	*eps = (1./(eos_gamma-1.))*(*P/ *rho);
	*vel_x = *vel_y = *vel_z = 0.;

  } else if ( CCTK_EQUALS( rho_profile, "CB Disk" ) ) {
	*psi = 1.;						// a fudge
	calcCBDiskRho( x, y, z, rho, P, vel_x, vel_y, vel_z );
	*eps = *P / ( *rho * (eos_gamma-1.) ); // Assumes ideal gas.

  } else if ( CCTK_EQUALS (rho_profile, "TOV") ) {
	calcTOVRho( x, y, z, psi, rho, eps, P, vel_x, vel_y, vel_z );
  } else if ( CCTK_EQUALS (rho_profile, "Bondi") ) {
        calcBondiRho( x, y, z, psi, rho, eps, P, vel_x, vel_y, vel_z, bvec_x, bvec_y, bvec_z );
        // Override psi by matter-only psi 
        *psi = 1.; 
  } else {
	CCTK_VWarn(0,__LINE__, __FILE__, CCTK_THORNSTRING,"The choice '%s' of rho_profile is not yet implemented.", rho_profile);
  } 

}

void matterlibs_calcConformalRho( double x, double y, double z, double *psi, double *rho, double *Ebar, double *Jbarx, double *Jbary, double *Jbarz, double *bvec_x, double *bvec_y, double *bvec_z, double *vB, double *w_lorentz ) 
{

  DECLARE_CCTK_PARAMETERS;

  // Initialize the B-field to either zero or the uniform B-field. Any other B-field 
  // needs to currently subsequently overwrite the zero B-field.
  if ( bvec_x != NULL ) {
     if ( CCTK_EQUALS(rho_profile,"UniformB") ) {
        *bvec_x = b0_uniform[0];
        *bvec_y = b0_uniform[1];
        *bvec_z = b0_uniform[2];
     } else {
        *bvec_x = 0;
        *bvec_y = 0;
        *bvec_z = 0;
     }
  }

  if ( CCTK_EQUALS (rho_profile, "TOV") ) {
	calcConformalTOVRho( x, y, z, psi, rho, Ebar, Jbarx, Jbary, Jbarz, w_lorentz );
  } else 
  {
	CCTK_VWarn(0,__LINE__, __FILE__, CCTK_THORNSTRING,"The choice '%s' of rho_profile is not yet supported for setting momenta in accordance with Bowen-type data.", rho_profile);
  } 
  *vB = 0; //TODO: find some way to guess what this is

}

// Gaussian-based Rho
void calcGaussRho( double x, double y, double z, double *drho, double *dvel_x, double *dvel_y, double *dvel_z )
{

  DECLARE_CCTK_PARAMETERS;
  double rho=0, crho, jx = 0, jy = 0, jz = 0;
  double xx, yy, zz, chi2;
  double rho_min;

  for ( int i=0; i<num_clouds; i++) {

    xx = x-cloud_x0[i];
    yy = y-cloud_y0[i];
    zz = z-cloud_z0[i];

    if ( cloud_sig_x[i] > 0 || cloud_sig_y[i] > 0 || cloud_sig_z[i] > 0 ) {

	double sigx = (cloud_sig_x[i]>0) ? cloud_sig_x[i] : cloud_sig[i];
	double sigy = (cloud_sig_y[i]>0) ? cloud_sig_y[i] : cloud_sig[i];
	double sigz = (cloud_sig_z[i]>0) ? cloud_sig_z[i] : cloud_sig[i];

	chi2 = SQR(xx/sigx) + SQR(yy/sigy) + SQR(zz/sigz);
	if ( cloud_rho[i] > 0 ) {
		double msig = sqrt( SQR(sigx) + SQR(sigy) + SQR(sigz) );
		chi2 = chi2 + cloud_r0[i]/msig * ( cloud_r0[i]/msig - sqrt(chi2) );
	}



    } else {

	double r2, rr;
	r2 = xx*xx+yy*yy+zz*zz;
	rr = sqrt(r2) - cloud_r0[i];
	chi2 = SQR( rr / cloud_sig[i] );

    }
    crho = cloud_rho[i] * exp(-0.5*chi2);

    if ( old_rho )
	crho = crho / cloud_sig[i];

    // calculate atmospheric density mimicking what whisky does
    if (initial_rho_abs_min > 0.0)
      rho_min = initial_rho_abs_min;
    else if (initial_rho_rel_min > 0.0)
      rho_min = cloud_rho[i] * initial_rho_rel_min;
    else if (rho_abs_min > 0.0)
      rho_min = rho_abs_min;
    else
      rho_min = cloud_rho[i] * rho_rel_min;

    if (initial_atmosphere_factor > 0.0)
      rho_min *= initial_atmosphere_factor;

    rho = rho + crho;
    if ( crho > rho_min ) { // only add velocity if not in atmosphere
    
      // add up the current densities instead of just the velocities
      jx = jx + crho * (cloud_omega_y[i]*zz - cloud_omega_z[i]*yy - Omega * y);
      jy = jy + crho * (cloud_omega_z[i]*xx - cloud_omega_x[i]*zz + Omega * x);
      jz = jz + crho * (cloud_omega_x[i]*yy - cloud_omega_y[i]*xx);
    
    }

  }

  // set return values
  *drho = rho;
  *dvel_x = jx / rho;
  *dvel_y = jy / rho;
  *dvel_z = jz / rho;

  if ( (*dvel_x)*(*dvel_x) + (*dvel_y)*(*dvel_y) + (*dvel_z)*(*dvel_z) >= 1. ) {

    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, 
	    "Faster-than-light matter at (%g,%g,%g): v^2 = %g", 
	    x,y,z, (*dvel_x)*(*dvel_x) + (*dvel_y)*(*dvel_y) + (*dvel_z)*(*dvel_z));

  }

} 

// Constant density Rho
double calcCloudRho( double x, double y, double z )
{

  DECLARE_CCTK_PARAMETERS;
  double rho=0;
  double r2;

  for ( int i=0; i<num_clouds; i++) {

      r2 = SQR(x-cloud_x0[i]) + SQR(y-cloud_y0[i]) + SQR(z-cloud_z0[i]);

      if ( r2 <= SQR(cloud_radius[i]) ) {
	   rho = rho + cloud_rho[i];
      }

  } 

  return rho;

}

// Circumbinary-type Disk aka O'Neil
void calcCBDiskRho( double x, double y, double z, double *drho, double *dpress, double *dvel_x, double *dvel_y, double *dvel_z )
{

  DECLARE_CCTK_PARAMETERS;
  double rho=0, jx = 0, jy = 0, jz = 0, press=0;
  double xx, yy, zz, rr, rcyl, cot2theta;
  double rho_min;

  xx = x - disk_center[0];
  yy = y - disk_center[1];
  zz = z - disk_center[2];

  rcyl = sqrt( SQR(xx) + SQR(yy) ) + REPS;
  rr = sqrt( SQR(xx) + SQR(yy) + SQR(zz) );
  cot2theta = SQR(zz) / SQR(rcyl);

  // calculate atmospheric density mimicking what whisky does
  if (initial_rho_abs_min > 0.0)
    rho_min = initial_rho_abs_min;
  else if (initial_rho_rel_min > 0.0)
    rho_min = disk_rho_inner * initial_rho_rel_min;
  else if (rho_abs_min > 0.0)
    rho_min = rho_abs_min;
  else
    rho_min = disk_rho_inner * rho_rel_min;

  if (initial_atmosphere_factor > 0.0)
    rho_min *= initial_atmosphere_factor;

  // Actual distribution
  double rho_of_r = ( disk_tanh_width > 0) ? 0.5 * disk_rho_inner *
	(tanh( (rr - disk_inner_radius)/disk_tanh_width ) + 1.) : disk_rho_inner;

  if ( (rcyl < disk_inner_radius) && (disk_tanh_width <= 0) ) {
     // Set the gas to atmosphere: rho_of_r does not apply everywhere.
     rho = rho_min;
     press = eos_k * pow( rho, eos_gamma );
  } else {
     rho = rho_of_r * exp( - 0.5*cot2theta / SQR(h_by_r) );
     press = ( rr * SQR(h_by_r) * SQR(rcyl/rr) ) * rho / SQR(rr - 2. + REPS);
  }

  if ( rho > rho_min ) { // only add velocity if not in atmosphere
    int vsign = 1;
    double vmag = sqrt(rcyl)/(rcyl-2.);
    if ( fabs(vmag) > disk_max_vel ) vmag = disk_max_vel;
    if ( disk_retrograde )
       vsign = -1;
    jx = rho * vsign * vmag * (-yy/rcyl);
    jy = rho * vsign * vmag * (xx/rcyl);
    jz = 0;
  }
 
  // set return values
  *drho = rho;
  *dpress = press;
  *dvel_x = jx / rho;
  *dvel_y = jy / rho;
  *dvel_z = jz / rho;

  if ( (*dvel_x)*(*dvel_x) + (*dvel_y)*(*dvel_y) + (*dvel_z)*(*dvel_z) >= 1. ) {

    CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING, 
	    "Faster-than-light matter at (%g,%g,%g): v^2 = %g", 
	    x,y,z, (*dvel_x)*(*dvel_x) + (*dvel_y)*(*dvel_y) + (*dvel_z)*(*dvel_z));

  }

} 

