#define INITVALUE (42)

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#define KRANC_C
#include "GenericFD.h"
#include "Symmetry.h"

#define SQR(x) ((x) * (x))
//#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
//#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
//#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
//*******************************//
//****** External Routines ******//
//*******************************//
void TotalRestMass_Scalars(CCTK_ARGUMENTS);
void TotalRestMass_Init(CCTK_ARGUMENTS);
void ZeroMasses(CCTK_ARGUMENTS);
void SumRestMass(CCTK_ARGUMENTS);
void TotalRestMass_output(CCTK_ARGUMENTS);

//*******************************//
//****** Internal Routines ******//
//*******************************//
static int sumLocalMass(CCTK_ARGUMENTS, CCTK_REAL * localVolMass, CCTK_REAL * localVol_dMdt);
static void totalmass_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc);
static CCTK_INT file_created;

void TotalRestMass_Init(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  CCTK_INT   is_symbnd[6];

  *SymmetryFactor = INITVALUE;
  for (int vol=0; vol<intvolumes; vol++) {
	VolumeMass[vol] = -1e234;
	Volume_dMassdt[vol] = -1e234;
  } 
  file_created=0;

  /* Get symmetry info */
  int ierr = GetSymmetryBoundaries(cctkGH,6,is_symbnd);
  if ( ierr < 0 ) {
	CCTK_WARN(1,"Error in returning the symmetries.");
  }
  assert(!ierr);

  CCTK_INT symfactL=1;
  for (int l = 0; l < 6 ; l++)
  {
	if (is_symbnd[l]) symfactL = symfactL * 2;
  }

  CCTK_VInfo(CCTK_THORNSTRING,"Due to symmetries, we multiply by a factor of %d",symfactL);
  *SymmetryFactor = symfactL;
  /*********************/
}

void ZeroMasses(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  for ( int vol=0; vol<intvolumes; vol++) {
      VolumeMass[vol] = 0;
      Volume_dMassdt[vol] = 0;
  }
}


void SumRestMass(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL sumVolMass[intvolumes], sumVol_dMdt[intvolumes];
  double localVolMass[intvolumes];
  double localVol_dMdt[intvolumes];
  int vol, ierr;
  int do_calc_mass, my_calc_mass_every;
 
  do_calc_mass = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_mass_every = calc_mass_every_vol[vol] >= 0 ? calc_mass_every_vol[vol] : calc_mass_every;
    assert(my_calc_mass_every >= 0);
    if ( (my_calc_mass_every > 0) && (cctk_iteration % my_calc_mass_every == 0) )
    {
        do_calc_mass = 1;
        break;
    }
  }
  if (do_calc_mass) 
  {

     if (verbose)
	CCTK_VInfo(CCTK_THORNSTRING,"Calculating total rest mass at iteration %d",cctk_iteration);

     for ( vol=0; vol<intvolumes; vol++)
     {
	localVolMass[vol] = -1e44;
	localVol_dMdt[vol] = -1e44;
     }

     if (verbose) CCTK_INFO("Ready for local loops");
     const int ierri = sumLocalMass(CCTK_PASS_CTOC, localVolMass, localVol_dMdt);
     assert(!ierri);
     if (ierri) CCTK_WARN (CCTK_WARN_ABORT, "Could not sum total masses");

     /* Get the summation operator */
     const int sum_op = CCTK_ReductionHandle ("sum");
     if (sum_op < 0)  CCTK_WARN (CCTK_WARN_ABORT, "error");

     /* do all reductions at once to reduce comminication overhead */
     ierr = CCTK_ReduceLocalArray1D (cctkGH, 0, sum_op, localVolMass, sumVolMass, intvolumes, CCTK_VARIABLE_REAL);
     assert(!ierr);
     if (ierr) CCTK_WARN (CCTK_WARN_ABORT, "error");
     /* Ok, so we lose to some more overhead with another array to be summed */
     ierr = CCTK_ReduceLocalArray1D (cctkGH, 0, sum_op, localVol_dMdt, sumVol_dMdt, intvolumes, CCTK_VARIABLE_REAL);
     assert(!ierr);
     if (ierr) CCTK_WARN (CCTK_WARN_ABORT, "error");
     for (vol=0; vol<intvolumes; vol++) {
     	     VolumeMass[vol] += *SymmetryFactor*sumVolMass[vol];
     	     Volume_dMassdt[vol] += *SymmetryFactor*sumVol_dMdt[vol];
     }

     if (verbose) {
	   CCTK_VInfo(CCTK_THORNSTRING,"Multiplied level results by %d for Symmetries: ",*SymmetryFactor);
	   if (intvolumes>0)
	      CCTK_VInfo(CCTK_THORNSTRING,"   VolumeMasses=[%g+%g,%g+%g,%g+%g,%g+%g,%g+%g].",
	   	VolumeMass[0],*SymmetryFactor*sumVolMass[0], VolumeMass[1],*SymmetryFactor*sumVolMass[1],
	   	VolumeMass[2],*SymmetryFactor*sumVolMass[2], VolumeMass[3],*SymmetryFactor*sumVolMass[3],
	   	VolumeMass[4],*SymmetryFactor*sumVolMass[4]);
     }

  }

}
     
static void totalmass_get_centre(CCTK_ARGUMENTS, int c, CCTK_REAL *xc, CCTK_REAL *yc, CCTK_REAL *zc)
{
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  assert(c >= 0);
  assert(c < intvolumes);
  assert(xc != 0);
  assert(yc != 0);
  assert(zc != 0);

  if (!CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    if (surface_index[c] != -1)
    {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "Parameter %s::surface_index[%d] has been set, but centre %d "
                  "is not being set from a spherical surface, so this parameter will be "
                  "ignored.", CCTK_THORNSTRING, c, c);
    }
  }

  if (CCTK_EQUALS(centre_from[c], "parameter"))
  {
    *xc = centre_x[c]; *yc = centre_y[c]; *zc = centre_z[c];
  }
  else if (CCTK_EQUALS(centre_from[c], "spherical surface"))
  {
    int sn=surface_index[c];
    if (sn < 0) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d is invalid for volume %d", sn,c);
    } else if (sn>=sphericalsurfaces_nsurfaces) {
      CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "surface number sn=%d too large, increase SphericalSurface::nsurfaces from its current value %d",
                 sn,sphericalsurfaces_nsurfaces);
    }
    if (sf_valid[sn]<=0) {
      CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                 "didn't find valid volume centre for sn=%d, vol=%d",sn,c);
    }

    *xc = sf_centroid_x[surface_index[c]];
    *yc = sf_centroid_y[surface_index[c]];
    *zc = sf_centroid_z[surface_index[c]];
  }
  else
  {
    assert(0); // Should never get here
  }

  if (verbose > 0)
  {
    CCTK_VInfo(CCTK_THORNSTRING,
               "centre %d, at (%g,%g,%g) radius %g", c, *xc, *yc, *zc, int_wi_rad[c]);
  }
}

static int sumLocalMass(CCTK_ARGUMENTS, CCTK_REAL * localVolMass, CCTK_REAL * localVol_dMdt )
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  CCTK_REAL rhoL,vxL,vyL,vzL,wL;
  CCTK_REAL gxxL,gyyL,gzzL,gxyL,gxzL,gyzL,detg;
  CCTK_REAL densL,weightL,MassL,dMassdtL,dDensdtL;
  CCTK_REAL dx,dy,dz,dV;
  CCTK_REAL xL,yL,zL,r2L;
  CCTK_REAL x0L[intvolumes],y0L[intvolumes],z0L[intvolumes];
  CCTK_REAL *VolWeights[intvolumes], *dDensVar;
  int vol, my_calc_mass_every;
  int num_negative_mass_points;

  /* Make sure things are initialized */
  for (vol=0; vol<intvolumes; vol++)
  {
    /* Make sure things are initialized */
    localVolMass[vol] = 0;
    localVol_dMdt[vol] = 0;
    /* Get Weight Pointers */
    VolWeights[vol] = NULL;
    VolWeights[vol] = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,weight_var[vol]));
    assert(VolWeights[vol]);
    /* set up local origins for individual integration volumes */
    totalmass_get_centre(CCTK_PASS_CTOC, vol, &x0L[vol], &y0L[vol], &z0L[vol]);
  }
  dDensVar = NULL;
  if(!CCTK_Equals(dDdt_variable, ""))
  {
    dDensVar = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH,0,dDdt_variable));
    assert(dDensVar);
  }

  /* Summation limits */
  CCTK_INT   is_symbnd[6], is_physbnd[6], is_ipbnd[6];
  CCTK_INT   nboundaryzones[6], is_internal[6], is_staggered[6], shiftout[6];
  CCTK_INT   imin[3], imax[3], npoints;

  GenericFD_GetBoundaryInfo(cctkGH, cctk_ash, cctk_lsh, cctk_bbox,
                            cctk_nghostzones,
                            imin, imax, is_symbnd, is_physbnd, is_ipbnd);

  if (CCTK_IsFunctionAliased ("GetBoundarySpecification")) {
    int ierr = GetBoundarySpecification(6, nboundaryzones, is_internal, is_staggered, shiftout);
    if (ierr != 0)
      CCTK_WARN(0, "Could not obtain boundary specification");
  } else {
    /* assume that nboundaryzones = nghostzones
     * leave is_internal, is_staggered, shiftout undefined
     */
    for (int dir = 0; dir<3; dir++)
    {
      for (int face = 0; face < 2; face++)
      {
         CCTK_INT indx = dir*2 + face;
         nboundaryzones[indx] = cctk_nghostzones[dir];
      }
    }
  }

  for (int dir = 0; dir<3; dir++)
  {
    for (int face = 0; face < 2; face++)
    {

       CCTK_INT indx = dir*2 + face;
       if (is_physbnd[indx]) {
   	   npoints=0;
       } else if (is_symbnd[indx]) {
	   npoints=nboundaryzones[indx];
       } else {
  	   npoints=cctk_nghostzones[dir];
       }

       switch(face)
       {
       case 0: /*lower*/
          imin[dir]=npoints;
	  break;
       case 1: /*upper*/
          imax[dir]=cctk_lsh[dir] - npoints;
	  break;
       default:
          CCTK_WARN(0,"internal error");
       }

    }
  }

  if (verbose) { 
	CCTK_VInfo(CCTK_THORNSTRING,"Extent: (%d,%d,%d) to (%d,%d,%d)",imin[0],imin[1],imin[2],imax[0],imax[1],imax[2]);
  }

  /* Loop over everything except ghost zones */ 
  num_negative_mass_points = 0; /* count number of nonsensical points in here */
  for(int k = imin[2] ; k < imax[2] ; k++)
  {
    for(int j = imin[1] ; j < imax[1] ; j++)
    {
      for(int i = imin[0] ; i < imax[0] ; i++)
      {

        /* Initiate variables */
        densL = INITVALUE;
	dDensdtL = INITVALUE;
        rhoL = INITVALUE;
        MassL = INITVALUE;
        weightL = INITVALUE;
	gxxL = INITVALUE;
	gxyL = INITVALUE;
	gxzL = INITVALUE;
	gyyL = INITVALUE;
	gyzL = INITVALUE;
	gzzL = INITVALUE;
	vxL = INITVALUE;
	vyL = INITVALUE;
	vzL = INITVALUE;

        dx = INITVALUE;
        dy = INITVALUE;
        dz = INITVALUE;
        dV = INITVALUE;

	xL = INITVALUE;
	yL = INITVALUE;
	zL = INITVALUE;
	r2L = INITVALUE;

	/* Information from the grid */
        int idx = CCTK_GFINDEX3D(cctkGH, i,j,k);

	rhoL = rho[idx];
        dDensdtL = dDensVar ? dDensVar[idx] : 0.;

	xL = x[idx];
	yL = y[idx];
	zL = z[idx];

	gxxL = gxx[idx];
	gxyL = gxy[idx];
	gxzL = gxz[idx];
	gyyL = gyy[idx];
	gyzL = gyz[idx];
	gzzL = gzz[idx];

	vxL = velx[idx];
	vyL = vely[idx];
	vzL = velz[idx];

        /* Get volume element */
        dx = CCTK_DELTA_SPACE(0);
        dy = CCTK_DELTA_SPACE(1);
        dz = CCTK_DELTA_SPACE(2);
        dV = dx*dy*dz;

	/* Calculate dens = sqrt(detg) rho w */
	detg = 2*gxyL*gxzL*gyzL + gzzL*(gxxL*gyyL - SQR(gxyL)) - gyyL*SQR(gxzL) - gxxL*SQR(gyzL);
	wL = 1. / sqrt( 1. - (gxxL*SQR(vxL) + gyyL*SQR(vyL) + gzzL*SQR(vzL) + 2*gxyL*vxL*vyL 
		+ 2*gxzL*vxL*vzL + 2*gyzL*vyL*vzL) );
	if ( rhoL <= rho_min ) {
		dDensdtL = 0;
		densL = 0;
	} else {
        	densL = sqrt(detg) * wL * rhoL;
		if (debug) {
			densL = 1;
			dDensdtL = 0.1;
		}
		
	}
        MassL = densL*dV;
	dMassdtL = dDensdtL*dV;

	if ( MassL < 0 ) {
	
            num_negative_mass_points += 1;
	    if ( verbose > 1 ) {
		CCTK_VInfo(CCTK_THORNSTRING,"Negative MassL found in TotalRestMass at (%g,%g,%g).  MassL=%g.", xL, yL, zL, MassL);
	    }
	
	}


	/* Add the integrand where appropriate */ 
	for (vol=0; vol<intvolumes; vol++)
	{
            my_calc_mass_every = calc_mass_every_vol[vol] >= 0 ? calc_mass_every_vol[vol] : calc_mass_every;
            assert(my_calc_mass_every >= 0);
            if ( (my_calc_mass_every == 0) || (cctk_iteration % my_calc_mass_every != 0) )
                continue;

	    weightL = VolWeights[vol][idx];
	    if ( weightL == 0 ) 
		continue;

	    r2L = SQR(xL-x0L[vol]) + SQR(yL-y0L[vol]) + SQR(zL-z0L[vol]) ;
	    if ( r2L < SQR(int_wi_rad[vol]) )
	    {
		localVolMass[vol] += weightL*MassL;
		localVol_dMdt[vol] += weightL*dMassdtL;
	    }
        }

      }
    }
  }

  if ( num_negative_mass_points>0 ) {
	CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "Local mass integrand in TotalRestMass is negative in %d points.  This implies there's something wrong.", 
                   num_negative_mass_points);
  }

 if (verbose) {
	CCTK_VInfo(CCTK_THORNSTRING,"After local integration ....");
        if (intvolumes>0)
	    CCTK_VInfo(CCTK_THORNSTRING,"      Integrated volume masses are (%g,%g,%g,%g,%g)",localVolMass[0],
		localVolMass[1], localVolMass[2], localVolMass[3], localVolMass[4]);
 }

 return 0;

}

/* IO */
void TotalRestMass_output(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;

  const char *fmode = (file_created>0) ? "a" : "w";  // file mode: append if already written
  char *filename;
  char varname[1024]; // XXX fixed size
  char file_extension[5]=".asc";
  char format_str_real[2048]; // XXX fixed size
  char addformat[50];
  int vol, do_calc_mass, my_calc_mass_every;
  FILE *file;

  const CCTK_INT myproc= CCTK_MyProc(cctkGH);

  do_calc_mass = 0;
  for (vol=0 ; vol<intvolumes; vol++)
  {
    my_calc_mass_every = calc_mass_every_vol[vol] >= 0 ? calc_mass_every_vol[vol] : calc_mass_every;
    assert(my_calc_mass_every >= 0);
    if ( (my_calc_mass_every > 0) && (cctk_iteration % my_calc_mass_every == 0) )
    {
        do_calc_mass = 1;
        break;
    }
  }
  if ( (myproc == 0) && do_calc_mass )
  {

    if (verbose) {
      CCTK_VInfo(CCTK_THORNSTRING, "writing output");
    }

    // filename
    sprintf(varname, "TotalRestMass");

    filename = (char *) malloc (strlen (out_dir) + strlen (varname) +
                                strlen (file_extension) +2);
    assert(filename);
    sprintf (filename, "%s/%s%s", out_dir, varname, file_extension);

    // open file
    file = fopen (filename, fmode);
    if (!file) {
      CCTK_VWarn (1, __LINE__, __FILE__, CCTK_THORNSTRING,
                  "write_totalrestmass: Could not open scalar output file '%s'",
                  filename);
    }

    // write header on startup
    if (file_created<=0) {
      fprintf(file,"# TotalRestMass mass totals\n");
      fprintf(file,"# gnuplot column index:\n");
      fprintf(file,"# 1:it 2:t \n");
      for (vol=0;vol<intvolumes;vol++) 
      {
	  fprintf(file,"#\t( r=%g ) %d:Mass %d:Int[dD/dt]\n",int_wi_rad[vol],2*vol+3,2*vol+4);
      }
      fprintf(file,"\n");
    }

    // format string
    sprintf (format_str_real, "%%d\t%%%s",out_format);
    sprintf (addformat,"\t%%%s\t%%%s",out_format,out_format);

    // write data
    fprintf(file, format_str_real, cctk_iteration, cctk_time);
    for (vol=0;vol<intvolumes;vol++)
    {
	fprintf(file,addformat,VolumeMass[vol],Volume_dMassdt[vol]);
    }
    fprintf(file,"\n");

    fclose(file);
    free(filename);

    if (file_created==0) {
      file_created=1;
    }

  }
}
