#include <assert.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

//*******************************//
//****** External Routines ******//
//*******************************//
void TotalRestMass_Startup(CCTK_ARGUMENTS);

void TotalRestMass_Startup(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;

  for (int vol=0; vol<intvolumes; vol++) {
        int my_calc_mass_every;

        /* ignore any volume that ist not actually computed */
        my_calc_mass_every = calc_mass_every_vol[vol] >= 0 ? calc_mass_every_vol[vol] : calc_mass_every;
        assert(my_calc_mass_every >= 0);
        if (my_calc_mass_every == 0)
        {
                continue;
        }
	/* Attempt to get variable index */
	int ierr=CCTK_VarIndex(weight_var[vol]);
	if (ierr < 0)
	{
		CCTK_PARAMWARN("TotalRestMass: Error in checking weight variable.");
	}
        /* check if spherical surface index is valid */
        if (CCTK_Equals(centre_from[vol], "spherical surface"))
        { 
                int sn, sf_error;

                sf_error = 0;
                sn = surface_index[vol];
                if (sn < 0) {
                  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                             "surface number sn=%d is invalid for volume %d", sn,vol);
                  sf_error = 1;
                } else if (sn >= sphericalsurfaces_nsurfaces) {
                  CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                             "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                             sn,sphericalsurfaces_nsurfaces);
                  sf_error = 1;
                }
                if (sf_error)
		  CCTK_PARAMWARN("TotalRestMass: Error in checking spherical surfaces.");
        }
  }

}
