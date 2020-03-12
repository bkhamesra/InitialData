
#include "matterlibs.h"
#include "math.h"
#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#define SQR(x) ((x)*(x))

/* When to add a poloidal Bfield? */
/*   If using twopunctures ... Add post-solving the tov, but before feeding to 2P (Add to getRho). However we need d_j rho and d_j press! */
/*   Simply apply it to a set of hydro variables after the fact (after SetMetric). Need to modify Tmunu then as well, and call Prim2Con */ 

/* Set both B and A if necessary, everywhere */
/* Requires rho, press GFs to have already been set and synced */
/* Requires an additional sync afterwards, including TmunuBase */
void MatterLibs_addPoloidalB( CCTK_ARGUMENTS )
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  

  /* HACK: Don't have the distance to the nearest star center on me here. Assume origin? */

  /* Does not impose Prim2Con. Please make sure Prim2ConInitial is run. */
  int set_Afield = CCTK_Equals(initial_Avec,"AforPoloidalB");
  int set_Bfield = CCTK_Equals(initial_Bvec,"poloidalBfield");

  CCTK_REAL dx, dy, dz, nx, ny, nz;
  CCTK_REAL dx_rho, dy_rho, dz_rho;
  CCTK_REAL dx_press, dy_press, dz_press;

  dx = CCTK_DELTA_SPACE(0);
  dy = CCTK_DELTA_SPACE(1);
  dz = CCTK_DELTA_SPACE(2);
  
  nx =cctk_lsh[0];
  ny =cctk_lsh[1];
  nz =cctk_lsh[2];

  CCTK_WARN(1,"Introduction of a poloidal magnetic field like this assumes compact sources located by parameters cloud_(x|y|z)0[N].");

  for ( int k = 2; k < nz-1; k++ )
     for ( int j = 2; j < ny-1; j++ )
        for ( int i = 2; i < nx-1; i++ ) {

           CCTK_INT idx = CCTK_GFINDEX3D( cctkGH, i, j, k );

           CCTK_INT x_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 0 );
           CCTK_INT y_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 1 );
           CCTK_INT z_idx = CCTK_VECTGFINDEX3D( cctkGH, i, j, k, 2 );

           CCTK_INT iplus1 = CCTK_GFINDEX3D( cctkGH, i+1, j,   k   );
           CCTK_INT jplus1 = CCTK_GFINDEX3D( cctkGH, i,   j+1, k   );
           CCTK_INT kplus1 = CCTK_GFINDEX3D( cctkGH, i,   j  , k+1 );
           CCTK_INT iminus1 = CCTK_GFINDEX3D( cctkGH, i-1, j,   k   );
           CCTK_INT jminus1 = CCTK_GFINDEX3D( cctkGH, i,   j-1, k   );
           CCTK_INT kminus1 = CCTK_GFINDEX3D( cctkGH, i,   j  , k-1 );

           CCTK_REAL rhoL, pressL, xL, yL, zL;
           CCTK_INT have_Bfield = 1;
           rhoL = rho[idx];
           pressL = press[idx];
           xL = x[idx];
           yL = y[idx];
           zL = z[idx];
           if ( press[idx] < add_Bfield_minP ) {
              have_Bfield = 0;
           }

           CCTK_REAL AxL, AyL, AzL, BxL, ByL, BzL;
           CCTK_REAL rhofact = 1. - rhoL/poloidal_rho_max;
           CCTK_REAL delP = pressL - add_Bfield_minP;
           CCTK_REAL maxP_Pcut = fmax( delP, 0. );
           CCTK_REAL AdotPhihat = poloidal_A_b*pow(rhofact,poloidal_n_p)*maxP_Pcut;
 
           CCTK_REAL xx, yy, zz, inv_rr;
           xx=yy=zz=30000;
           for ( int msrce=0; msrce < num_clouds; msrce++ ) {
               xx = fmin( xx, xL - cloud_x0[msrce] );
               yy = fmin( yy, yL - cloud_y0[msrce] );
               zz = fmin( zz, zL - cloud_z0[msrce] );
           }
           inv_rr = 1./sqrt( SQR(xx) + SQR(yy) + SQR(zz) );

           if ( have_Bfield > 0 ) {
              AxL = - yy * AdotPhihat;
              AyL =   xx * AdotPhihat;
              AzL = 0.;
           } else {
              AxL = 0.;
              AyL = 0.;
              AzL = 0.;
           }

           if ( set_Bfield > 0 && have_Bfield > 0 ) {

              // HACK: Hard-coded 2nd order
              CCTK_REAL dx_rho, dy_rho, dz_rho;
              CCTK_REAL dx_press, dy_press, dz_press;
              dx_rho = (rho[iplus1]-rho[iminus1])/dx;
              dy_rho = (rho[jplus1]-rho[jminus1])/dy;
              dz_rho = (rho[kplus1]-rho[kminus1])/dz;
              dx_press = (press[iplus1]-press[iminus1])/dx;
              dy_press = (press[jplus1]-press[jminus1])/dy;
              dz_press = (press[kplus1]-press[kminus1])/dz;
    
              CCTK_REAL dAPhi_dx, dAPhi_dy, dAPhi_dz;
              dAPhi_dx = poloidal_A_b*pow(rhofact,poloidal_n_p-1)*maxP_Pcut
                       * ( rhofact*dx_press/delP - poloidal_n_p*dx_rho/poloidal_rho_max);
              dAPhi_dy = poloidal_A_b*pow(rhofact,poloidal_n_p-1)*maxP_Pcut
                       * ( rhofact*dy_press/delP - poloidal_n_p*dy_rho/poloidal_rho_max);
              dAPhi_dz = poloidal_A_b*pow(rhofact,poloidal_n_p-1)*maxP_Pcut
                       * ( rhofact*dz_press/delP - poloidal_n_p*dz_rho/poloidal_rho_max);
   
              CCTK_REAL dy_dAx, dz_dAx, dx_dAy, dz_dAy;
              dy_dAx = - AdotPhihat - yy*dAPhi_dy;
              dz_dAx =              - yy*dAPhi_dz;
              dx_dAy =   AdotPhihat + xx*dAPhi_dy;
              dz_dAy =                xx*dAPhi_dz;
   
 
              CCTK_REAL detg  = 2*gxy[idx]*gxz[idx]*gyz[idx] + gzz[idx]*(gxx[idx]*gyy[idx] 
                              - SQR(gxy[idx])) - gyy[idx]*SQR(gxz[idx]) - gxx[idx]*SQR(gyz[idx]); 
              CCTK_REAL inv_sqrtdet = 1./sqrt(detg);

              if ( have_Bfield > 0 ) {
                 BxL =  dz_dAy*inv_sqrtdet;
                 ByL = -dz_dAx*inv_sqrtdet;
                 BzL = ( dy_dAx - dx_dAy )*inv_sqrtdet;
              } else {
                 BxL = 0;
                 ByL = 0;
                 BzL = 0;
              }

              Bvec[x_idx] = BxL;
              Bvec[y_idx] = ByL;
              Bvec[z_idx] = BzL;
    
           }
 
           if ( set_Afield ) {
              Avec[x_idx] = AxL;
              Avec[y_idx] = AyL;
              Avec[z_idx] = AzL;
              if ( Aphi != NULL ) {
                 Aphi[idx] = 0.;
              }
           }

  }

}

