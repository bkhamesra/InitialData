#include <stdio.h>

#include "assert.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "EOS_Base.h"
#include "primsolver.h"
#include "steffenson.h" // for gsl_test_interval
#include "HydroBase.h"

#define DEBUG
#define MINB2 1e-24
#define velx (&vel[0*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define vely (&vel[1*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])
#define velz (&vel[2*cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]])




static inline double pow2 (const double x)
{
  return x*x;
}
static inline double pow3 (const double x)
{
  return x*x*x;
}
static inline double pow4 (const double x)
{
  const double y = x*x;
  return y*y;
}

CCTK_REAL PrimitiveSolver_getRho0W(CCTK_INT polyhandle, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psiL)
{
    DECLARE_CCTK_PARAMETERS;
    double dEbar, dJbar[3];
    double junk, junk6[6];
    getSourceInfo( x, y, z, &dEbar, dJbar, &junk, &junk, junk6);
    
    double dJbarx, dJbary, dJbarz;
    dJbarx = dJbar[0];
    dJbary = dJbar[1];
    dJbarz = dJbar[2];
    //CCTK_VInfo(CCTK_THORNSTRING, "PRIMSOLVER_GETTING_RH0");
    // double psiL = getConformalFactor(x, y, z);
    // Calculate appropriate values for E and S_j after baring and unbaring 
    // as in gr-qc/0606104.
    
    double dE, dJx, dJy, dJz;
    dE = dEbar/pow(psiL, conformal_density_power);
    dJx = dJbarx/pow(psiL, 10.);
    dJy = dJbary/pow(psiL, 10.);
    dJz = dJbarz/pow(psiL, 10.);
    
    // for a polytropic EOS D is linked to E and S^j
    	
    double drho = dE;
    double dlorentz;
    int Etorho_err;
    assert(polyhandle>=0);
    Etorho_err = prims_E_to_rho(polyhandle, dE, dJx, dJy, dJz, psiL, &drho, &dlorentz);
    
    if( Etorho_err != 0 )
    {
        CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
                "initial recovery of primitive variables failed at "
                "(%g,%g,%g), where (E,J^j)=[%g,%g,%g,%g], "
                "psiL = %g, (rho,w_lorentz) = (%g,%g)",
                x,y,z, dE, dJx,dJy,dJz, psiL, drho,dlorentz);
    }
    
    return drho*dlorentz;
}

void PrimitiveSolver_SetPrims(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
  
    //comented  const char * eos_type_name;
    double dpsi, drho, deps, dpress, dvelx, dvely, dvelz; 	// primitive initial data
    //comented  double dE, dD, dtau, dJx, dJy, dJz;    
    double dE, dJx, dJy, dJz;                   			// conservative initial data
    double dEbar, dJbar[3];					// conformal initial data
    double dlorentz2, dlorentz, dh;				// auxiliary variables
    double xx, yy, zz;
    double psiL, chiL, rootdetg;
    double twelfth, third;
    double g11, g21, g31, g22, g32, g33, detg;
    //comente  double uxx, uxy, uxz, uyy, uyz, uzz;
    static double linmom[3];                                      // total linear momentum of matter
    static double angmom[3];                                      // total angular momentum of matter 
    FILE *testsuitefh=NULL;
    static int firsttime = 1;				        // chain output together 
  
  
  
    CCTK_VInfo(CCTK_THORNSTRING, "PRIMITIVE SOLVER WORKING  testsuite");
    // open file to dump internal state to
    if(testsuite)
    {
      char fn[1024];
  
      snprintf(fn, sizeof(fn)/sizeof(fn[0]), "%s/primsolver_internaldata.asc", out_dir);
      if((testsuitefh = fopen(fn, firsttime ? "w" : "a")) == NULL)
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,  "Cannot open file %s to dump internal data", fn);
  
      if(firsttime) {
            fprintf(testsuitefh, "# 1:x 2:y 3:z 4:psi 5:rho 6:eps 7:press 8:velx 9:vely 10:velz 11:psi_final\n");
         }
  
      firsttime = 0;
    }
  
    // Whisky kindly sets the initial rho_min in Whisky_Rho_Minima_Setup
  
    twelfth = 1.0/12.0;
    third = 1.0/3.0;
  
    CCTK_REAL *CarpetWeights = (CCTK_REAL *)(CCTK_VarDataPtr(cctkGH, 0, "CarpetReduce::weight"));
    assert(CarpetWeights);
    for ( int kk = 0 ; kk < cctk_lsh[2] ; kk++) {
      for ( int j = 0 ; j < cctk_lsh[1] ; j++) {
        for (int i = 0 ; i < cctk_lsh[0] ; i++) {
  
          const int ind = CCTK_GFINDEX3D (cctkGH, i, j, kk);
  
          xx = x[ind];
          yy = y[ind];
          zz = z[ind];
  
          g11 = gxx[ind];
          g22 = gyy[ind];
          g33 = gzz[ind];
          g21 = gxy[ind];
          g31 = gxz[ind];
          g32 = gyz[ind];
          detg  =  2*g21*g31*g32 + g33*(g11*g22 - pow2(g21)) - g22*pow2(g31) - g11*pow2(g32);
  
  
  	  chiL = fmax(1/pow(detg,third),chiEps);
  	  // psiL is used in the rescaling from Ebar to E, it should not be
  	  // capped from above since that makes the initial data potentially
  	  // large. rootdetg on the other hand is only used to densitize the
  	  // conservative varibles at the end, here we want to cap it from above
  	  // to limit the magnitude of the ID
          // psiL = 1/pow(chiL,1.0/4.0); 
  	  rootdetg = 1/pow(chiL,3.0/2.0);
          psiL = pow(detg,twelfth);
  	  /*rootdetg = sqrt(detg);*/
  
          // obtain guesses for the primitive variables
  	  double junk, junk6[6];
          getSourceInfo( xx, yy, zz, &dEbar, dJbar, &junk, &junk, junk6);
  
  
          // calculate auxiliary variables
          // Calculate linear and angular momentum through raw primitives
  	  double dJbarx, dJbary, dJbarz;
  	  dJbarx = dJbar[0];
  	  dJbary = dJbar[1];
  	  dJbarz = dJbar[2];
          if (CarpetWeights[ind] == 1.0 &&
               i >= cctk_nghostzones[0] &&  i < cctk_lsh[0]-cctk_nghostzones[0] &&
               j >= cctk_nghostzones[1] &&  j < cctk_lsh[1]-cctk_nghostzones[1] &&
              kk >= cctk_nghostzones[2] && kk < cctk_lsh[2]-cctk_nghostzones[2]) {
  
             // Calculate appropriate values for E and S_j after baring and unbaring 
             // as in gr-qc/0606104.
  
             linmom[0] += dJbarx*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
             linmom[1] += dJbary*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
             linmom[2] += dJbarz*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2); 
             angmom[0] += (yy*dJbarz-zz*dJbary)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
             angmom[1] += (zz*dJbarx-xx*dJbarz)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
             angmom[2] += (xx*dJbary-yy*dJbarx)*CCTK_DELTA_SPACE(0)*CCTK_DELTA_SPACE(1)*CCTK_DELTA_SPACE(2);
  
          }
          // DEBUG
  
          // Calculate appropriate values for E and S_j after baring and unbaring 
          // as in gr-qc/0606104.
          dE = dEbar/pow(psiL, conformal_density_power);
  	  dJx = dJbarx/pow(psiL, 10.);
  	  dJy = dJbary/pow(psiL, 10.);
  	  dJz = dJbarz/pow(psiL, 10.);
  
  	  // for a polytropic EOS D is linked to E and S^j
  
  	  drho = dE;
  
                 
  
  	  int Etorho_err;
  	  // assert(*PrimitiveSolver_polytrope_handle>=0);
  	  // assert(EOS_Handle("2D_Polytrope">=0);
 	  //comente	Etorho_err = prims_E_to_rho(*PrimitiveSolver_polytrope_handle, dE, dJx, dJy, dJz, psiL, &drho, &dlorentz);Esto quiere decir que solo va a funcionar con 2D_Polytrope agregar un parametro para que funcione con otras EOShandles
  	  Etorho_err = prims_E_to_rho(EOS_Handle("2D_Polytrope"), dE, dJx, dJy, dJz, psiL, &drho, &dlorentz);
  	  if( Etorho_err != 0 )
  	  {
  	      CCTK_VWarn(0,__LINE__,__FILE__,CCTK_THORNSTRING,
  		      "initial recovery of primitive variables failed at "
  		      "(%g,%g,%g), where (E,J^j)=[%g,%g,%g,%g], "
  		      "psiL = %g, dpsi = %g => (rho,w_lorentz) = (%g,%g)",
  		      xx,yy,zz, dE, dJx,dJy,dJz, psiL,dpsi, drho,dlorentz);
  	      return; /* NOTREACHED */
  	  }
  
  	  if( drho < PrimitiveSolver_rho_abs_min*(1.0+PrimitiveSolver_atmo_tolerance) ) // enforce at least atmospheric density
  	  {
  	    drho = PrimitiveSolver_rho_abs_min;
  	    dpress = EOS_Pressure(EOS_Handle("2D_Polytrope"), drho, primsolver_poison); // atmosphere always used a polytrope
  	    deps = EOS_SpecificIntEnergy(EOS_Handle("2D_Polytrope"), drho, primsolver_poison);
  	    dvelx = dvely = dvelz = 0.;
  	    dlorentz2 = dlorentz = 1.;
                     
  
  	    dJx = dJy = dJz = 0.;
  	  }
  	  else
  	  {
  	    dpress = EOS_Pressure(EOS_Handle("2D_Polytrope"), drho, primsolver_poison);
  	    deps = EOS_SpecificIntEnergy(EOS_Handle("2D_Polytrope"), drho, primsolver_poison);
  	    dlorentz2 = pow2(dlorentz);
            dh = 1 + deps + dpress/drho;
  	    dvelx = ( dJx ) / (drho * dh * dlorentz2 );
  	    dvely = ( dJy ) / (drho * dh * dlorentz2 );
  	    dvelz = ( dJz ) / (drho * dh * dlorentz2 );
  
  
  	    // calculate a suitable matter density (this is the reason we go 
  	    // through the whole rigmalore...)
  	    // hand crafted to guarantee positivity
  
            if(verbose >= 1 &&
                GSL_SUCCESS != prims_gsl_root_test_interval(dE, drho*dh*dlorentz2-dpress, abs_err, rel_err)) 
            {
                CCTK_VWarn(1, __LINE__, __FILE__, CCTK_THORNSTRING,  
                        "Back and forth conversion of relativistic energy density E failed: "
                        "E(rho) != E_orig; E(rho) = %g, rho = %g, E_orig = %g, diff = %e", 
                        drho*dh*dlorentz2-dpress, 
                        drho, dE, dE - (drho*dh*dlorentz2-dpress));
            }
  	  }
             
          // factor in the determinant that Whisky uses
          dJx *= rootdetg;
          dJy *= rootdetg;
          dJz *= rootdetg;
  	
  
          rho[ind] = drho;
          velx[ind] = dvelx;
          vely[ind] = dvely;
          velz[ind] = dvelz;
          eps[ind] = deps;
          press[ind] = dpress;
          w_lorentz[ind] = dlorentz;
  
  
  
  	  }
      }
   }
  
   if(verbose >= 2) {
  
      CCTK_VInfo (CCTK_THORNSTRING, 
          "Total matter linear momentum P^i = (%g, %g, %g)",
          linmom[0],linmom[1],linmom[2]);
      CCTK_VInfo (CCTK_THORNSTRING, 
          "Total matter angular momentum J^i = (%g, %g, %g)", 
          angmom[0],angmom[1],angmom[2]);
  
    }
  
  
   if(testsuitefh)
   fclose(testsuitefh);

}



