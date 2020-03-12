
#include "matterlibs.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"

double matterlibs_calcMatterPsi( double x, double y, double z ) 
{

  DECLARE_CCTK_PARAMETERS;
  double psi;

  // Choose the function based on the keyword parameter
  if ( CCTK_EQUALS( rho_profile, "Gaussian" ) ) {
       	psi = 1.; // 1. is "flat"

  } else if ( CCTK_EQUALS( rho_profile, "Constant" ) ) {
       	psi = 1.; // 1. is "flat"

  } else if ( CCTK_EQUALS (rho_profile, "TOV") ) {
	double dummy;

	calcTOVRho( x, y, z, &psi, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy );
  } else if ( CCTK_EQUALS (rho_profile, "CB Disk") ) {
	psi = 1.;
  } else if ( CCTK_EQUALS (rho_profile, "Bondi") ) {
        psi = 1.; // Hackish 
  } else {
	CCTK_VWarn(0,__LINE__, __FILE__, CCTK_THORNSTRING,"The choice '%s' of rho_profile is not yet implemented.", rho_profile);
  } 

  return psi;

}
