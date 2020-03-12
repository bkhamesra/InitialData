// TwoPunctures:  File  "Equations.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "cctk.h"
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

void addteuk (double x, double y, double z, double time, int w_ind, double Aij[3][3])
{
  DECLARE_CCTK_PARAMETERS;
  int i, j;
  double efftime, xx, yy, zz, rad, kktw[3][3], cofs[5], orderr;

  /* Setting up useful vars to increase efficiency*/
  xx = x - move_origin_x - wave_orig_x[w_ind];  /* Sets the origin of the wave: */
  yy = y - wave_orig_y[w_ind];
  zz = z - wave_orig_z[w_ind];

  efftime = time - rad0[w_ind];

  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {
      kktw[i][j]=0;
    }
  }

  /* Check parameters for Teuk wave and get the wave*/

  if( w_ind >= 25) {
	CCTK_WARN (0, "Trying to add more waves than the current hard-coded limit.  Not a hard fix if it's what you really want.");
  }

  if( CCTK_EQUALS( ffshape[w_ind], "eppley" ) || CCTK_EQUALS( ffshape[w_ind], "eppleykij" ) || CCTK_EQUALS( ffshape[w_ind], "eppleysmod" ) ) {
	if ( even[w_ind] == 0 ) {
	  CCTK_WARN (0, "Sorry, odd Teukolsky wave is not yet implemented ..." );
	  } 
	else {
	  /* r=0 requires special treatment */
  	  rad = sqrt( xx*xx + yy*yy + zz*zz );
          if (rad < 1.0e-10 ) {
	    CCTK_INFO ("Entering origin");
	    CCTK_WARN (1, "The origin region of the Teukolsky wave thorn needs fixing if not at origin or far from origin!  Coming soon ...");
	    Origin_Eppley(xx, yy, zz, efftime, w_ind, kktw);
	    orderr = exp(-pow(rad0[w_ind]/sigma[w_ind],2));
	    if ( (texprad > 0) && (orderr >= 1e-20) && (rad0[w_ind] != 0) ) {
		CCTK_WARN (1, "You're using Taylor expansion for a wave that isn't at or far from the origin.  Be aware of errors near the origin.");
		fprintf(stderr,"\tErrors around the origin are at least of order %.19g with the current expansion scheme.\n", orderr);
	    }
	  } 
  	  else {

	    // Get TN tensor coefficients
	    if ( rad <= texprad ) {
		//CCTK_INFO("Going into the Taylor expansion.");
		Texp_ABC_Eppley(xx, yy, zz, time, w_ind, cofs);
	    }
	    else {
	    	ABC_Eppley(xx, yy, zz, time, w_ind, cofs);
	    }

	    // Evaluate the exact tensor
	    ExtCurv(xx, yy, zz, w_ind, cofs, kktw);
	  }
	}
  }
  else { CCTK_WARN (0, "Invalid function type."); }

  /* Add the wave to the traceless-transverse part of Kij */
  for (i = 0; i < 3; i++)
  {
    for (j = 0; j < 3; j++)
    {	
      Aij[i][j] += kktw[i][j];
    }
  }

}

/* Deal with the origin separately */
void Origin_Eppley(double x, double y, double z, double acttime, int w_ind, double kktw[3][3])
{

  DECLARE_CCTK_PARAMETERS;
  double tos, tos2, sigma2, coa, kp2;

  tos = acttime / sigma[w_ind];
  sigma2 = sigma[w_ind] * sigma[w_ind];
  tos2 = tos * tos;
  kp2 = kwave[w_ind]*kwave[w_ind];

  // A, B, C are all equal: we only use A

  // Hard coded in time-symmetric origin
  if ( CCTK_EQUALS( ffshape[w_ind], "eppley" ) ) {
      coa = 8 * amp[w_ind] * exp(-tos2) / (5*sigma2*sigma2) * (8*tos2*tos2*tos2 - 
	    60 * tos2*tos2 + 90*tos2 - 15); 
  } else if ( CCTK_EQUALS( ffshape[w_ind], "eppleykij" ) ) {
      coa = - 16 * amp[w_ind] * acttime * exp(-tos2) / (5*sigma2*sigma2*sigma2) *
	    (8 * tos2 * tos2 * tos2 - 84 * tos2 * tos2 + 210 * tos2 - 105); 
  } else if ( CCTK_EQUALS( ffshape[w_ind], "eppleysmod" ) ) {
      coa = 2*amp[w_ind]*exp(-tos2)/(5*sigma2*sigma2) * ( kwave[w_ind]*acttime*(80*tos2*tos2 - 400*tos2 - 40*tos2*kp2*sigma2 + 300 + 60*kp2*sigma2 
	    + kp2*kp2*sigma2*sigma2)*sin(kwave[w_ind]*acttime) + 
	    (32*tos2*tos2*tos2 -240*tos2*tos2+360*tos2 -80*kp2*tos2*tos2*sigma2 -60*kp2*sigma2 + 10*kp2*kp2*tos2*sigma2*sigma2 -5*kp2*kp2*sigma2*sigma2
	    + 240*kp2*tos2*sigma2-60)*cos(kwave[w_ind]*acttime) );
  }

// Debugging:
//fprintf(stderr,"coa = %.19g, acttime = %.19g, tos = %.19g, sigma = %.19g, w_ind = %d\n",coa,acttime,tos,sigma[w_ind],w_ind);

// Override coefficient?  Sometimes the origin is a problem ...
//coa = 0;

  if (even[w_ind] == 1) {
     switch(mm[w_ind]) {
       case -2:
	  kktw[0][0] = 0;
          kktw[0][1] = coa;
	  kktw[0][2] = 0;
          kktw[1][1] = 0;
	  kktw[1][2] = 0;
          kktw[2][2] = 0;
	  break;
       case -1:
          kktw[0][0] = 0;
	  kktw[0][1] = 0;
	  kktw[0][2] = 0;
          kktw[1][1] = 0;
          kktw[1][2] = coa;
          kktw[2][2] = 0;
	  break;
       case 0:
	  kktw[0][0] = - coa;
	  kktw[0][1] = 0;
	  kktw[0][2] = 0;
          kktw[1][1] = - coa;
          kktw[1][2] = 0;
          kktw[2][2] = 2 * coa;
	  break;
       case 1:
          kktw[0][0] = 0;
          kktw[0][1] = 0;
          kktw[0][2] = coa;
          kktw[1][1] = 0;
          kktw[1][2] = 0;
          kktw[2][2] = 0;
	  break;
       case 2:
          kktw[0][0] = coa;
          kktw[0][1] = 0;
          kktw[0][2] = 0;
          kktw[1][1] = - coa;
          kktw[1][2] = 0;
          kktw[2][2] = 0;
	  break;
     }}
  else {
     switch(mm[w_ind]) {
       case -2:
	  break;
       case -1:
	  break;
       case 0:
	  break;
       case 1:
	  break;
       case 2:
	  break;
     }
   }

  /* Apply symmetries here */
  kktw[1][0] = kktw[0][1];
  kktw[2][0] = kktw[0][2];
  kktw[2][1] = kktw[1][2];

}

/* Calculate the coeffecients */
void ABC_Eppley(double x, double y, double z, double time, int w_ind, double cofs[5])
{

  DECLARE_CCTK_PARAMETERS;
  double rad, rad2, rad3, rad4, rad5, u, dff[5];
  int kij;

  rad2 = x*x + y*y + z*z;
  rad = sqrt( rad2 );
  rad3 = rad2 * rad;
  rad4 = rad3 * rad;
  rad5 = rad4 * rad;

  kij = 0;
  // Set toggle for kij case
  if ( CCTK_EQUALS( ffshape[w_ind], "eppleykij" ) ) {
     kij = 1;
  }

  if (CCTK_EQUALS(wave_dir[w_ind], "outgoing")) {
	  /* Outgoing Part */
	  CCTK_WARN(0, "Beware: Pure outgoing wave isn't regular at origin (unless rad0 is significantly large and you don't mind extra noise).");

	  u = time - (rad - rad0[w_ind]);
	  if (CCTK_EQUALS(ffshape[w_ind],"eppleysmod")) {
	  	Shape_SMEppley( dff, u, amp[w_ind], sigma[w_ind], kwave[w_ind] );
	  } else if (CCTK_EQUALS(ffshape[w_ind],"eppley") || CCTK_EQUALS(ffshape[w_ind],"eppleykij") ) {
	  	Shape_Eppley( dff, u, amp[w_ind], sigma[w_ind], kij );
	  }	
	  cofs[0] = 3 * (dff[2]/rad3 + 3*dff[1]/rad4 + 3*dff[0]/rad5);
	  cofs[1] = -(dff[3]/rad2 + 3*dff[2]/rad3 + 6*dff[1]/rad4 + 6*dff[0]/rad5);
	  cofs[2] = (dff[4]/rad + 2*dff[3]/rad2 + 9*dff[2]/rad3 + 21*dff[1]/rad4 
		    + 21*dff[0]/rad5)/4;
	  cofs[3] = dff[3]/rad + 2*dff[2]/rad2 + 3*dff[1]/rad3 + 3*dff[0]/rad4;
          cofs[4] = 0; 

   }   else if (CCTK_EQUALS(wave_dir[w_ind], "ingoing"))  {
	  /* Ingoing Part */
	  CCTK_WARN(0, "Beware: Pure ingoing wave isn't regular at origin (unless rad0 is sufficiently large and you don't mind extra noise).");

	  u = time + (rad - rad0[w_ind]);
	  if (CCTK_EQUALS(ffshape[w_ind],"eppleysmod")) {
	  	Shape_SMEppley( dff, u, amp[w_ind], sigma[w_ind], kwave[w_ind] );
	  } else if (CCTK_EQUALS(ffshape[w_ind],"eppley") || CCTK_EQUALS(ffshape[w_ind],"eppleykij") ) {
	  	Shape_Eppley( dff, u, amp[w_ind], sigma[w_ind], kij );
	  }	
	  cofs[0] = 3*(dff[2]/rad3 - 3*dff[1]/rad4 + 3*dff[0]/rad5);
	  cofs[1] = dff[3]/rad2 - 3*dff[2]/rad3 + 6*dff[1]/rad4 - 6*dff[0]/rad5;
	  cofs[2] = (dff[4]/rad - 2*dff[3]/rad2 + 9*dff[2]/rad3 - 21*dff[1]/rad4 
		     + 21*dff[0]/rad5) /4;
	  cofs[3] = (-dff[3]/rad + 2*dff[2]/rad2 - 3*dff[1]/rad3 + 3*dff[0]/rad4);
	  cofs[4] = rad * cofs[0] / 3; 

   }   else if   (CCTK_EQUALS(wave_dir[w_ind], "timesym"))   {
	  /* Outgoing Part */
	  u = time - (rad - rad0[w_ind]);
	  if (CCTK_EQUALS(ffshape[w_ind],"eppleysmod")) {
	  	Shape_SMEppley( dff, u, amp[w_ind], sigma[w_ind], kwave[w_ind] );
	  } else if (CCTK_EQUALS(ffshape[w_ind],"eppley") || CCTK_EQUALS(ffshape[w_ind],"eppleykij") ) {
	  	Shape_Eppley( dff, u, amp[w_ind], sigma[w_ind], kij );
	  }	
	  cofs[0] = 3 * (dff[2]/rad3 + 3*dff[1]/rad4 + 3*dff[0]/rad5);
	  cofs[1] = -(dff[3]/rad2 + 3*dff[2]/rad3 + 6*dff[1]/rad4 + 6*dff[0]/rad5);
	  cofs[2] = (dff[4]/rad + 2*dff[3]/rad2 + 9*dff[2]/rad3 + 21*dff[1]/rad4 
		    + 21*dff[0]/rad5)/4;
	  cofs[3] = dff[3]/rad + 2*dff[2]/rad2 + 3*dff[1]/rad3 + 3*dff[0]/rad4; 

	  
	  /* Subtract Ingoing Part */
	  u = time + (rad - rad0[w_ind]);
	  if (CCTK_EQUALS(ffshape[w_ind],"eppleysmod")) {
	  	Shape_SMEppley( dff, u, amp[w_ind], sigma[w_ind], kwave[w_ind] );
	  } else if (CCTK_EQUALS(ffshape[w_ind],"eppley") || CCTK_EQUALS(ffshape[w_ind],"eppleykij") ) {
	  	Shape_Eppley( dff, u, amp[w_ind], sigma[w_ind], kij );
	  }	
	  cofs[0] += - 3*(dff[2]/rad3 - 3*dff[1]/rad4 + 3*dff[0]/rad5);
	  cofs[1] += -dff[3]/rad2 + 3*dff[2]/rad3 - 6*dff[1]/rad4 + 6*dff[0]/rad5;
	  cofs[2] += - (dff[4]/rad - 2*dff[3]/rad2 + 9*dff[2]/rad3 - 21*dff[1]/rad4 
		     + 21*dff[0]/rad5) /4;
	  cofs[3] += - (-dff[3]/rad + 2*dff[2]/rad2 - 3*dff[1]/rad3 + 3*dff[0]/rad4);
	  cofs[4] = rad * cofs[0] / 3;

  }   else   {
	CCTK_WARN(0, "Incorrect parameter for wave_dir" );
  }

/* Debugging */
//  if (x <= 2 && x >= -2 && y<=2 && z<=2 ) {
//     fprintf(stderr, "\nIn ABC_Eppley: (x,y,z,amp)=(%f,%f,%f,%f)\t(A,B,C)=(%f,%f,%f)",
//	     x,y,z,amp,cofs[0],cofs[1],cofs[2]);
//  }

}

/* Calculate the wave portion */
void ExtCurv( double x, double y, double z, int w_ind, double cofs[5], double kktw[3][3])
{

  DECLARE_CCTK_PARAMETERS;
  double x2, y2, z2, rad, rad2, rad3, rad4, rho2, coa, cob, coc, col, cok;
  int sx;

  coa = cofs[0];
  cob = cofs[1];
  coc = cofs[2];
  col = cofs[3];
  cok = cofs[4];

  x2 = x*x;
  y2 = y*y;
  z2 = z*z;
  rad2 = x2 + y2 + z2;
  rad = sqrt(rad2);
  rad3 = rad2 * rad;
  rad4 = rad3 * rad;
  rho2 = x2 + y2;
  sx = x/sqrt(x2);

  if (even[w_ind] == 1) {
     switch(mm[w_ind]) { 
       case -2:
	  kktw[0][0] = 2 * y * x / rad4 * (x2 * coa + (y2 + z2 - x2) * cob 
			- (y2 + z2) * coc);
          kktw[0][1] = (-z2 / rad2 + 2 * x2 * y2 / rad4) * coa 
			+ (1 - z2 / rad2 - 4 * x2 * y2 / rad4) * cob
			+ 2 * (z2 / rad2 + x2 * y2 / rad4) * coc;
          kktw[0][2] = y * z / rad4 * ((rad2 + 2 * x2) * coa 
			+ (rad2 - 4 * x2) * cob - 2 * (y2 + z2) * coc);
          kktw[1][1] = 2 * x * y / rad4 * (y2 * coa + (x2 + z2 - y2) * cob 
			- (x2 + z2) * coc);
          kktw[1][2] = x * z / rad4 * ((rad2 + 2 * y2) * coa + (rad2 - 4 * y2) 
			* cob - 2 * (x2 + z2) * coc);
          kktw[2][2] = 2 * x * y / rad4 * (-(x2 + y2) * coa - 2 * z2 * cob + 
			(rad2 + z2) * coc);
	  break;
       case -1:
	  kktw[0][0] = 2 * y * z * sx / rad4 * (-(y2 + z2) * coa - 2 * x2 * cob + 
			(rad2 + x2) * coc); 
          kktw[0][1] = z * sqrt(x2) / rad4 * ((rad2 + 2 * y2) * coa + (rad - 4 * y2) 
			* cob - 2 * (x2 + z2) * coc);
          kktw[0][2] = sqrt(x2) * y / rad4 * ((rad2 + 2 * z2) * coa + (rad2 
			- 4 * z2) * cob - 2 * (x2 + y2) * coc);
          kktw[1][1] = 2 * y * z * sx / rad4 * (y2 * coa + (rad2 - 2 * y2)
			* cob - (x2 + z2) * coc);
          kktw[1][2] = sx * (-x2 / rad2 + 2 * y2 * z2 / rad4) * coa
			+ sx * (1 - x2 / rad2 - 4 * y2 * z2 / rad4) * cob
			+ sx * 2 * (x2 / rad2 + y2 * z2 / rad4) * coc;
          kktw[2][2] = 2 * y * z / (rad4*sx) * (z2 * coa + (rad2 - 2 * z2)
			* cob - (x2 + y2) * coc);
	  break;
       case 0:
	  kktw[0][0] = (-1 + 3 * y2 / rad2 + 3 * x2 * z2 / rad4) * coa
			- 6 * z2 * x2 / rad4 * cob
			+ 3 * (-y2 / rad2 + x2 * z2 / rad4) * coc;
          kktw[0][1] = 3 * x * y / rad4 * (-rho2 * coa - 2 * z2 * cob
			+ (rad2 + z2) * coc);
          kktw[0][2] = 3 * x * z / rad4 * (z2 * coa + (rho2 - z2) * cob
			- rho2 * coc);
          kktw[1][1] = (-1 + 3 * x2 / rad2 + 3 * y2 * z2 / rad4) * coa
			- 6 * y2 * z2 / rad4 * cob
			+ 3 * (-x2 / rad2 + y2 * z2 / rad4) * coc;
          kktw[1][2] = 3 * y * z / rad4 * (z2 * coa + (rho2 - z2) * cob
			- rho2 * coc);
          kktw[2][2] = (-1 + 3 * z2*z2 / rad4) * coa + 6 * z2 * rho2
			/ rad4 * cob + 3 * rho2*rho2 / rad4 * coc;
	  break;
       case 1:
	  kktw[0][0] = 2 * sqrt(x2) * z / rad4 * (x2 * coa + (rad2 - 2 * x2)
			* cob - (y2 + z2) * coc);
          kktw[0][1] = y * z / (rad4*sx) * ((rad2 + 2 * x2) * coa + (rad2
			- 4 * x2) * cob - 2 * (y2 + z2) * coc);
          kktw[0][2] = ((-y2 / rad2 + 2 * x2 * z2 / rad4) * coa
			+ (1 - y2 / rad2 - 4 * x2 * z2 / rad4) * cob
			+ 2 * (y2 / rad2 + x2 * z2 / rad4) * coc)/sx;
          kktw[1][1] = 2 * sqrt(x2) * z / rad4 * (-(x2 + z2) * coa - 2 * y2
			* cob + (rad2 + y2) * coc);
          kktw[1][2] = sqrt(x2) * y / rad4 * ((rad2 + 2 * z2) * coa + (rad2
			- 4 * z2) * cob - 2 * (x2 + y2) * coc);
          kktw[2][2] = 2 * sqrt(x2) * z / rad4 * (z2 * coa + (rad2 - 2 * z2)
			* cob - (x2 + y2) * coc);
	  break;
       case 2:
	  kktw[0][0] = ((x2 - z2) / rad2 - x2 * (z2 + 2 * y2) / rad4)
			* coa + 2 * x2 * (z2 + 2 * y2) / rad4 * cob
                        + ((y2 + 2 * z2) / rad2 - x2 * (z2 + 2 * y2)
                        / rad4) * coc;
          kktw[0][1] = y * x * (x2 - y2) / rad4 * (coa - 2 * cob + coc);
          kktw[0][2] = x * z / rad4 * ((2 * x2 + z2) * coa + (-x2
			+ 3 * y2 + z2) * cob - (x2 + 2 * z2 + 3 * y2)
                        * coc);
          kktw[1][1] = ((z2 - y2) / rad2 + y2 * (2 * x2 + z2) / rad4)
			* coa - 2 * y2 * (2 * x2 + z2) / rad4 * cob
                        + (-(x2 + 2 * z2) / rad2 + y2 * (2 * x2 + z2)
                        / rad4) * coc;
          kktw[1][2] = y * z / rad4 * (-(z2 + 2 * y2) * coa - (3 * x2
			- y2 + z2) * cob + (3 * x2 + 2 * z2 + y2) * coc);
          kktw[2][2] = (y2*y2 - x2*x2) / rad4 *coa - 2 * z2 * (x2 - y2)
			/ rad4 * cob + (x2 - y2) * (rad2 + z2) / rad4
                        * coc;
	  break;
     }}
  else {
     switch(mm[w_ind]) {
       case -2:
	  break;
       case -1:
	  break;
       case 0:
	  kktw[0][0] = x * y * z * (8 * cok + 2 * col) / rad3; 
          kktw[0][1] = - z * (x2 - y2) * (4 * cok + col) / rad3;
          kktw[0][2] = - y * (col * rho2 - 4 * cok * z2) / rad3;
          kktw[1][1] = - x * y * z * (8 * cok + 2 * col) / rad;
          kktw[1][2] = x * (col * rho2 - 4 * cok * z2) / rad3;
          kktw[2][2] = 0;
	  break;
       case 1:
	  break;
       case 2:
	  break;
     }
  }

  /* Symmetrize */
  kktw[1][0] = kktw[0][1];
  kktw[2][0] = kktw[0][2];
  kktw[2][1] = kktw[1][2];


}

/* You know ... the wave */
void Shape_SMEppley( double dff[5], double u, double amp, double sigma, double k)
{

  double uos, uos2, ex2s2, dffex, sku, cku, sig2, kp2;  

  uos = u / sigma;
  uos2 = uos * uos;
  ex2s2 = amp * exp(-uos2); /* F(t+-r) / u */
  sku = sin(k*u);
  cku = cos(k*u);
  sig2 = sigma*sigma;
  kp2 = k*k;

  dff[0] = cku * u * ex2s2;
  dff[1] = ((1-2*uos2)*cku - k*u*sku) * ex2s2;
  dff[2] = (-2*k*(1-2*uos2)*sku + u*((4*uos2-6)/sig2-kp2)*cku) * ex2s2;
  dff[3] = (k*u*((18-12*uos2)/sig2 + kp2)*sku + 
	   ( (24*uos2-8*uos2*uos2-6)/sig2 -3*kp2 + 6*uos2*kp2)*cku)* ex2s2;
  dff[4] = (u*( (16*uos2*uos2-80*uos2+60)/(sig2*sig2) + (36*kp2-24*kp2*uos2)/sig2 + kp2*kp2 )*cku + 
	   4*k*( (6-24*uos2+8*uos2*uos2)/sig2 - 2*kp2*uos2 + kp2)*sku) * ex2s2;

}

/* You know ... the wave */
void Shape_Eppley( double dff[5], double u, double amp, double sigma, int kij )
{

  double uos, uos2, ex2s2, dffex;  

  uos = u / sigma;
  uos2 = uos * uos;
  ex2s2 = amp * exp(-uos2); /* F(t+-r) / u */

  dff[0] = u * ex2s2;
  dff[1] = (1 - 2 * uos2) * ex2s2;
  dff[2] = uos / sigma * (-6 + 4 * uos2) * ex2s2;
  dff[3] = 1 / (sigma*sigma) * (-6 + 24 * uos2 - 8 * uos2*uos2) * ex2s2;
  dff[4] = uos / (sigma*sigma*sigma) * (60 - 80*uos2 + 16*uos2*uos2) * ex2s2;

  if (kij) {
  	dffex = 1 / (sigma*sigma*sigma*sigma) * (60 - 360*uos2 + 240*uos2*uos2 - 
	        32*uos2*uos2*uos2)*ex2s2; 
  	dff[0] = - dff[1]/2;
	dff[1] = - dff[2]/2;
	dff[2] = - dff[3]/2;
	dff[3] = - dff[4]/2;
	dff[4] = - dffex/2;
  }

}

/* Calculate the coeffecients via Taylor Expansion */
void Texp_ABC_Eppley(double x, double y, double z, double time, int w_ind, double cofs[5])
{

  DECLARE_CCTK_PARAMETERS;
  double rad, rad2, rad3, rad4, rad5, ex2s2;
  double r0, r0p2, r0p4, r0p6, r0p8;
  double sigp2, sigp4, sigp6, sigp8;
  double k, kp2, kp4, kp6, kp8, skr, ckr; 

  rad2 = x*x + y*y + z*z;
  rad = sqrt( rad2 );
  rad3 = rad2 * rad;
  rad4 = rad3 * rad;
  rad5 = rad4 * rad;

  r0=rad0[w_ind];
  r0p2=rad0[w_ind]*rad0[w_ind];
  r0p4=r0p2*r0p2;
  r0p6=r0p4*r0p2;
  r0p8=r0p6*r0p2;

  sigp2=sigma[w_ind]*sigma[w_ind];
  sigp4=sigp2*sigp2;
  sigp6=sigp4*sigp2;
  sigp8=sigp6*sigp2;
  ex2s2=exp(r0p2/sigp2);

  k=kwave[w_ind];
  kp2=k*k;
  kp4=kp2*kp2;
  kp6=kp4*kp2;
  kp8=kp6*kp2;
  skr=sin(k*rad0[w_ind]);
  ckr=cos(k*rad0[w_ind]);

  // Parameter checking and warnings
  if ( CCTK_EQUALS( ffshape[w_ind], "eppleykij" ) ) {
	CCTK_WARN(0,"Taylor expansion has not been implemented for the eppleykij packet.  Please set texprad=0 to continue with this packet.");
  }
  else if ( CCTK_EQUALS( wave_dir[w_ind], "ingoing" ) || CCTK_EQUALS( wave_dir[w_ind], "outgoing" ) ) {
	CCTK_WARN(0,"Taylor expansion is only valid for time-symmetric initial data.  Please set texprad=0 to continue with this packet.");
  }

  if ( CCTK_EQUALS( ffshape[w_ind], "eppley" ) ) {

        // Evaluate the Taylor Expansions of the coefficients to 8th order
        cofs[0] = -((8*amp[w_ind]*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*ex2s2*sigp8*sigp2) -
                  (8*amp[w_ind]*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/
                  (35.*ex2s2*sigp8*sigp6) +
                  (4*amp[w_ind]*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 -
                  9450*r0p2*sigp8 + 945*sigp6*sigp4))/(315.*ex2s2*sigp8*sigp8*sigp2) -
                  (4*amp[w_ind]*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +
                  207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(10395.*ex2s2*sigp8*sigp8*sigp6) +
                  (amp[w_ind]*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 -
                  2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/
                  (135135.*ex2s2*sigp8*sigp8*sigp8*sigp2));
        cofs[1] = -((8*amp[w_ind]*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*ex2s2*sigp8*sigp2) -
                  (8*amp[w_ind]*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/
                  (21.*ex2s2*sigp8*sigp6) +
                  (4*amp[w_ind]*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 -
                  9450*r0p2*sigp8 + 945*sigp6*sigp4))/(135.*ex2s2*sigp8*sigp8*sigp2) -
                  (4*amp[w_ind]*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +
                  207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(3465.*ex2s2*sigp8*sigp8*sigp6) +
                  (amp[w_ind]*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 -
                  2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/
                  (36855.*ex2s2*sigp8*sigp8*sigp8*sigp2));
        cofs[2] = -((8*amp[w_ind]*(-8*r0p6 + 60*r0p4*sigp2 - 90*r0p2*sigp4 + 15*sigp6))/(5.*ex2s2*sigp8*sigp2) -
                  (8*amp[w_ind]*rad2*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/
                  (15.*ex2s2*sigp8*sigp6) +
                  (52*amp[w_ind]*rad4*(-32*r0p6*r0p4 + 720*r0p8*sigp2 - 5040*r0p6*sigp4 + 12600*r0p4*sigp6 -
                  9450*r0p2*sigp8 + 945*sigp6*sigp4))/(945.*ex2s2*sigp8*sigp8*sigp2) -
                  (4*amp[w_ind]*rad3*rad3*(64*r0p6*r0p6 - 2112*r0p8*r0p2*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +
                  207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp6*sigp6))/(1485.*ex2s2*sigp8*sigp8*sigp6) +
                  (31.*amp[w_ind]*rad4*rad4*(-128*r0p8*r0p6 + 5824*r0p6*r0p6*sigp2 - 96096*r0p8*r0p2*sigp4 + 720720*r0p8*sigp6 -
                  2522520*r0p6*sigp8 + 3783780*r0p4*sigp6*sigp4 - 1891890*r0p2*sigp6*sigp6 + 135135*sigp8*sigp6))/
                  (405405.*ex2s2*sigp8*sigp8*sigp8*sigp2));
        cofs[3] = (32*amp[w_ind]*rad*(8*r0p6 - 60*r0p4*sigp2 + 90*r0p2*sigp4 - 15*sigp6))/(15.*ex2s2*sigp6*sigp4) +
                  (16*amp[w_ind]*rad3*(16*r0p8 - 224*r0p6*sigp2 + 840*r0p4*sigp4 - 840*r0p2*sigp6 + 105*sigp8))/
                  (35.*ex2s2*sigp8*sigp6) +
                  (32*amp[w_ind]*rad4*rad*(32*r0p8*r0p2 - 720*r0p8*sigp2 + 5040*r0p6*sigp4 - 12600*r0p4*sigp6 + 9450*r0p2*sigp8 -
                  945*sigp4*sigp6))/(945.*ex2s2*sigp8*sigp8*sigp2) +
                  (8*amp[w_ind]*rad3*rad4*(64*r0p8*r0p4 - 2112*r0p6*r0p4*sigp2 + 23760*r0p8*sigp4 - 110880*r0p6*sigp6 +
                  207900*r0p4*sigp8 - 124740*r0p2*sigp6*sigp4 + 10395*sigp8*sigp4))/(6237.*ex2s2*sigp8*sigp8*sigp6);
        cofs[4] = rad*cofs[0]/3.;

  }
  else if ( CCTK_EQUALS( ffshape[w_ind], "eppleysmod" ) ) {

 	// Evaluate the Taylor expansions of the coefficients to 6th order

	cofs[0] = (2*amp[w_ind]*((32*r0p6 - 80*r0p4*sigp2*(3 + kp2*sigp2) - 5*sigp6*(12 + 12*kp2*sigp2 + kp4*sigp4) + 
          	  10*r0p2*sigp4*(36 + 24*kp2*sigp2 + kp4*sigp4))*ckr + 
	          k*r0*sigp2*(80*r0p4 - 40*r0p2*sigp2*(10 + kp2*sigp2) + sigp4*(300 + 60*kp2*sigp2 + kp4*sigp4))*skr))/
	   	  (5.*ex2s2*sigp8*sigp2)
	  	  - rad*(amp[w_ind]*(r0*(-64*r0p6 + 48*r0p4*sigp2*(14 + 5*kp2*sigp2) - 60*r0p2*sigp4*(28 + 20*kp2*sigp2 + kp4*sigp4) + 
           	  sigp6*(840 + 900*kp2*sigp2 + 90*kp4*sigp4 + kp6*sigp6))*ckr + 
        	  2*k*sigp2*(-96*r0p6 + 80*r0p4*sigp2*(9 + kp2*sigp2) + 3*sigp6*(60 + 20*kp2*sigp2 + kp4*sigp4) - 
           	  6*r0p2*sigp4*(180 + 40*kp2*sigp2 + kp4*sigp4))*skr))/(8.*ex2s2*sigp8*sigp4)
		  + rad2*(amp[w_ind]*((128*r0p8 - 224*r0p6*sigp2*(8 + 3*kp2*sigp2) + 280*r0p4*sigp4*(24 + 18*kp2*sigp2 + kp4*sigp4) + 
          	  7*sigp8*(120 + 180*kp2*sigp2 + 30*kp4*sigp4 + kp6*sigp6) - 
          	  14*r0p2*sigp6*(480 + 540*kp2*sigp2 + 60*kp4*sigp4 + kp6*sigp6))*ckr - 
       		  k*r0*sigp2*(-448*r0p6 + 112*r0p4*sigp2*(42 + 5*kp2*sigp2) - 28*r0p2*sigp4*(420 + 100*kp2*sigp2 + 3*kp4*sigp4) + 
          	  sigp6*(5880 + 2100*kp2*sigp2 + 126*kp4*sigp4 + kp6*sigp6))*skr))/(35.*ex2s2*sigp8*sigp6)
		  + rad3* (amp[w_ind]*(r0*(256*r0p8 - 256*r0p6*sigp2*(18 + 7*kp2*sigp2) + 224*r0p4*sigp4*(108 + 84*kp2*sigp2 + 5*kp4*sigp4) - 
          	  112*r0p2*sigp6*(360 + 420*kp2*sigp2 + 50*kp4*sigp4 + kp6*sigp6) + 
          	  sigp8*(15120 + 23520*kp2*sigp2 + 4200*kp4*sigp4 + 168*kp6*sigp6 + kp8*sigp8))*ckr + 
       		  8*k*sigp2*(128*r0p8 - 224*r0p6*sigp2*(8 + kp2*sigp2) + 56*r0p4*sigp4*(120 + 30*kp2*sigp2 + kp4*sigp4) + 
          	  sigp8*(840 + 420*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
          	  2*r0p2*sigp6*(3360 + 1260*kp2*sigp2 + 84*kp4*sigp4 + kp6*sigp6))*skr))/(192.*ex2s2*sigp8*sigp8)
		  + rad4*(amp[w_ind]*((512*r0p8*r0p2 - 2304*r0p8*sigp2*(5 + 2*kp2*sigp2) + 4032*r0p6*sigp4*(20 + 16*kp2*sigp2 + kp4*sigp4) - 
          	  672*r0p4*sigp6*(300 + 360*kp2*sigp2 + 45*kp4*sigp4 + kp6*sigp6) - 
          	  9*sigp8*sigp2*(1680 + 3360*kp2*sigp2 + 840*kp4*sigp4 + 56*kp6*sigp6 + kp8*sigp8) + 
          	  18*r0p2*sigp8*(8400 + 13440*kp2*sigp2 + 2520*kp4*sigp4 + 112*kp6*sigp6 + kp8*sigp8))*ckr + 
       		  k*r0*sigp2*(2304*r0p8 - 768*r0p6*sigp2*(54 + 7*kp2*sigp2) + 2016*r0p4*sigp4*(108 + 28*kp2*sigp2 + kp4*sigp4) - 
          	  144*r0p2*sigp6*(2520 + 980*kp2*sigp2 + 70*kp4*sigp4 + kp6*sigp6) + 
          	  sigp8*(136080 + 70560*kp2*sigp2 + 7560*kp4*sigp4 + 216*kp6*sigp6 + kp8*sigp8))*skr))/(1260.*ex2s2*sigp8*sigp8*sigp2)
		  - rad5*(amp[w_ind]*(r0*(-1024*r0p8*r0p2 + 1280*r0p8*sigp2*(22 + 9*kp2*sigp2) - 1920*r0p6*sigp4*(132 + 108*kp2*sigp2 + 7*kp4*sigp4) + 
           	  3360*r0p4*sigp6*(264 + 324*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
           	  60*r0p2*sigp8*(18480 + 30240*kp2*sigp2 + 5880*kp4*sigp4 + 280*kp6*sigp6 + 3*kp8*sigp8) + 
           	  sigp8*sigp2*(332640 + 680400*kp2*sigp2 + 176400*kp4*sigp4 + 12600*kp6*sigp6 + 270*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr + 
        	  2*k*sigp2*(-2560*r0p8*r0p2 + 3840*r0p8*sigp2*(15 + 2*kp2*sigp2) - 1344*r0p6*sigp4*(300 + 80*kp2*sigp2 + 3*kp4*sigp4) + 
           	  480*r0p4*sigp6*(2100 + 840*kp2*sigp2 + 63*kp4*sigp4 + kp6*sigp6) + 
           	  5*sigp8*sigp2*(15120 + 10080*kp2*sigp2 + 1512*kp4*sigp4 + 72*kp6*sigp6 + kp8*sigp8) - 
           	  10*r0p2*sigp8*(75600 + 40320*kp2*sigp2 + 4536*kp4*sigp4 + 144*kp6*sigp6 + kp8*sigp8))*skr))/(9600.*ex2s2*sigp8*sigp8*sigp4)
		  + rad5*rad*(amp[w_ind]*((2048*r0p6*r0p6 - 5632*r0p8*r0p2*sigp2*(12 + 5*kp2*sigp2) + 42240*r0p8*sigp4*(18 + 15*kp2*sigp2 + kp4*sigp4) - 
          	  14784*r0p6*sigp6*(240 + 300*kp2*sigp2 + 40*kp4*sigp4 + kp6*sigp6) + 
          	  1320*r0p4*sigp8*(5040 + 8400*kp2*sigp2 + 1680*kp4*sigp4 + 84*kp6*sigp6 + kp8*sigp8) + 
          	  11*sigp8*sigp4*(30240 + 75600*kp2*sigp2 + 25200*kp4*sigp4 + 2520*kp6*sigp6 + 90*kp8*sigp8 + kp8*kp2*sigp8*sigp2) - 
          	  22*r0p2*sigp8*sigp2*(181440 + 378000*kp2*sigp2 + 100800*kp4*sigp4 + 7560*kp6*sigp6 + 180*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr - 
       		  k*r0*sigp2*(-11264*r0p8*r0p2 + 14080*r0p8*sigp2*(22 + 3*kp2*sigp2) - 4224*r0p6*sigp4*(660 + 180*kp2*sigp2 + 7*kp4*sigp4) + 
          	  1056*r0p4*sigp6*(9240 + 3780*kp2*sigp2 + 294*kp4*sigp4 + 5*kp6*sigp6) - 
          	  220*r0p2*sigp8*(55440 + 30240*kp2*sigp2 + 3528*kp4*sigp4 + 120*kp6*sigp6 + kp8*sigp8) + 
          	  sigp8*sigp2*(3659040 + 2494800*kp2*sigp2 + 388080*kp4*sigp4 + 19800*kp6*sigp6 + 330*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*skr))/
   		  (83160.*ex2s2*sigp8*sigp8*sigp6);

	cofs[1] = 2*amp[w_ind]*((32*r0p6 - 80*r0p4*sigp2*(3 + kp2*sigp2) - 5*sigp6*(12 + 12*kp2*sigp2 + kp4*sigp4) + 
          	  10*r0p2*sigp4*(36 + 24*kp2*sigp2 + kp4*sigp4))*ckr + 
       		  k*r0*sigp2*(80*r0p4 - 40*r0p2*sigp2*(10 + kp2*sigp2) + sigp4*(300 + 60*kp2*sigp2 + kp4*sigp4))*skr)/
   		  (5.*ex2s2*sigp8*sigp2)
		  - rad*(amp[w_ind]*(r0*(-64*r0p6 + 48*r0p4*sigp2*(14 + 5*kp2*sigp2) - 60*r0p2*sigp4*(28 + 20*kp2*sigp2 + kp4*sigp4) + 
           	  sigp6*(840 + 900*kp2*sigp2 + 90*kp4*sigp4 + kp6*sigp6))*ckr + 
        	  2*k*sigp2*(-96*r0p6 + 80*r0p4*sigp2*(9 + kp2*sigp2) + 3*sigp6*(60 + 20*kp2*sigp2 + kp4*sigp4) - 
           	  6*r0p2*sigp4*(180 + 40*kp2*sigp2 + kp4*sigp4))*skr))/(6.*ex2s2*sigp8*sigp4)
		  + rad2*(amp[w_ind]*((128*r0p8 - 224*r0p6*sigp2*(8 + 3*kp2*sigp2) + 280*r0p4*sigp4*(24 + 18*kp2*sigp2 + kp4*sigp4) + 
          	  7*sigp8*(120 + 180*kp2*sigp2 + 30*kp4*sigp4 + kp6*sigp6) - 
          	  14*r0p2*sigp6*(480 + 540*kp2*sigp2 + 60*kp4*sigp4 + kp6*sigp6))*ckr - 
       		  k*r0*sigp2*(-448*r0p6 + 112*r0p4*sigp2*(42 + 5*kp2*sigp2) - 28*r0p2*sigp4*(420 + 100*kp2*sigp2 + 3*kp4*sigp4) + 
          	  sigp6*(5880 + 2100*kp2*sigp2 + 126*kp4*sigp4 + kp6*sigp6))*skr))/(21.*ex2s2*sigp8*sigp6)
		  + rad3*(amp[w_ind]*(r0*(256*r0p8 - 256*r0p6*sigp2*(18 + 7*kp2*sigp2) + 224*r0p4*sigp4*(108 + 84*kp2*sigp2 + 5*kp4*sigp4) - 
          	  112*r0p2*sigp6*(360 + 420*kp2*sigp2 + 50*kp4*sigp4 + kp6*sigp6) + 
          	  sigp8*(15120 + 23520*kp2*sigp2 + 4200*kp4*sigp4 + 168*kp6*sigp6 + kp8*sigp8))*ckr + 
       		  8*k*sigp2*(128*r0p8 - 224*r0p6*sigp2*(8 + kp2*sigp2) + 56*r0p4*sigp4*(120 + 30*kp2*sigp2 + kp4*sigp4) + 
          	  sigp8*(840 + 420*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
          	  2*r0p2*sigp6*(3360 + 1260*kp2*sigp2 + 84*kp4*sigp4 + kp6*sigp6))*skr))/(96.*ex2s2*sigp8*sigp8)
		  + rad4*(amp[w_ind]*((512*r0p8*r0p2 - 2304*r0p8*sigp2*(5 + 2*kp2*sigp2) + 4032*r0p6*sigp4*(20 + 16*kp2*sigp2 + kp4*sigp4) - 
          	  672*r0p4*sigp6*(300 + 360*kp2*sigp2 + 45*kp4*sigp4 + kp6*sigp6) - 
          	  9*sigp8*sigp2*(1680 + 3360*kp2*sigp2 + 840*kp4*sigp4 + 56*kp6*sigp6 + kp8*sigp8) + 
          	  18*r0p2*sigp8*(8400 + 13440*kp2*sigp2 + 2520*kp4*sigp4 + 112*kp6*sigp6 + kp8*sigp8))*ckr + 
       		  k*r0*sigp2*(2304*r0p8 - 768*r0p6*sigp2*(54 + 7*kp2*sigp2) + 2016*r0p4*sigp4*(108 + 28*kp2*sigp2 + kp4*sigp4) - 
         	  144*r0p2*sigp6*(2520 + 980*kp2*sigp2 + 70*kp4*sigp4 + kp6*sigp6) + 
        	  sigp8*(136080 + 70560*kp2*sigp2 + 7560*kp4*sigp4 + 216*kp6*sigp6 + kp8*sigp8))*skr))/(540.*ex2s2*sigp8*sigp8*sigp2)
		  - rad5*(amp[w_ind]*(r0*(-1024*r0p8*r0p2 + 1280*r0p8*sigp2*(22 + 9*kp2*sigp2) - 1920*r0p6*sigp4*(132 + 108*kp2*sigp2 + 7*kp4*sigp4) + 
         	  3360*r0p4*sigp6*(264 + 324*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
         	  60*r0p2*sigp8*(18480 + 30240*kp2*sigp2 + 5880*kp4*sigp4 + 280*kp6*sigp6 + 3*kp8*sigp8) + 
         	  sigp8*sigp2*(332640 + 680400*kp2*sigp2 + 176400*kp4*sigp4 + 12600*kp6*sigp6 + 270*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr + 
      		  2*k*sigp2*(-2560*r0p8*r0p2 + 3840*r0p8*sigp2*(15 + 2*kp2*sigp2) - 1344*r0p6*sigp4*(300 + 80*kp2*sigp2 + 3*kp4*sigp4) + 
         	  480*r0p4*sigp6*(2100 + 840*kp2*sigp2 + 63*kp4*sigp4 + kp6*sigp6) + 
         	  5*sigp8*sigp2*(15120 + 10080*kp2*sigp2 + 1512*kp4*sigp4 + 72*kp6*sigp6 + kp8*sigp8) - 
         	  10*r0p2*sigp8*(75600 + 40320*kp2*sigp2 + 4536*kp4*sigp4 + 144*kp6*sigp6 + kp8*sigp8))*skr))/
 		  (3600.*ex2s2*sigp8*sigp8*sigp4)
		  + rad5*rad*(amp[w_ind]*((2048*r0p6*r0p6 - 5632*r0p8*r0p2*sigp2*(12 + 5*kp2*sigp2) + 42240*r0p8*sigp4*(18 + 15*kp2*sigp2 + kp4*sigp4) - 
        	  14784*r0p6*sigp6*(240 + 300*kp2*sigp2 + 40*kp4*sigp4 + kp6*sigp6) + 
        	  1320*r0p4*sigp8*(5040 + 8400*kp2*sigp2 + 1680*kp4*sigp4 + 84*kp6*sigp6 + kp8*sigp8) + 
        	  11*sigp8*sigp4*(30240 + 75600*kp2*sigp2 + 25200*kp4*sigp4 + 2520*kp6*sigp6 + 90*kp8*sigp8 + kp8*kp2*sigp8*sigp2) - 
        	  22*r0p2*sigp8*sigp2*(181440 + 378000*kp2*sigp2 + 100800*kp4*sigp4 + 7560*kp6*sigp6 + 180*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr - 
       		  k*r0*sigp2*(-11264*r0p8*r0p2 + 14080*r0p8*sigp2*(22 + 3*kp2*sigp2) - 4224*r0p6*sigp4*(660 + 180*kp2*sigp2 + 7*kp4*sigp4) + 
        	  1056*r0p4*sigp6*(9240 + 3780*kp2*sigp2 + 294*kp4*sigp4 + 5*kp6*sigp6) - 
        	  220*r0p2*sigp8*(55440 + 30240*kp2*sigp2 + 3528*kp4*sigp4 + 120*kp6*sigp6 + kp8*sigp8) + 
        	  sigp8*sigp2*(3659040 + 2494800*kp2*sigp2 + 388080*kp4*sigp4 + 19800*kp6*sigp6 + 330*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*skr))/
   		  (27720.*ex2s2*sigp8*sigp8*sigp6);

	cofs[2] = (2*amp[w_ind]*((32*r0p6 - 80*r0p4*sigp2*(3 + kp2*sigp2) - 5*sigp6*(12 + 12*kp2*sigp2 + kp4*sigp4) + 
        	  10*r0p2*sigp4*(36 + 24*kp2*sigp2 + kp4*sigp4))*ckr + 
     		  k*r0*sigp2*(80*r0p4 - 40*r0p2*sigp2*(10 + kp2*sigp2) + sigp4*(300 + 60*kp2*sigp2 + kp4*sigp4))*skr))/
 		  (5.*ex2s2*sigp8*sigp2)
		  + rad*(-19*amp[w_ind]*(r0*(-64*r0p6 + 48*r0p4*sigp2*(14 + 5*kp2*sigp2) - 60*r0p2*sigp4*(28 + 20*kp2*sigp2 + kp4*sigp4) + 
        	  sigp6*(840 + 900*kp2*sigp2 + 90*kp4*sigp4 + kp6*sigp6))*ckr + 
     		  2*k*sigp2*(-96*r0p6 + 80*r0p4*sigp2*(9 + kp2*sigp2) + 3*sigp6*(60 + 20*kp2*sigp2 + kp4*sigp4) - 
     		  6*r0p2*sigp4*(180 + 40*kp2*sigp2 + kp4*sigp4))*skr))/(96.*ex2s2*sigp8*sigp4)
		  + rad2*(amp[w_ind]*((128*r0p8 - 224*r0p6*sigp2*(8 + 3*kp2*sigp2) + 280*r0p4*sigp4*(24 + 18*kp2*sigp2 + kp4*sigp4) + 
        	  7*sigp8*(120 + 180*kp2*sigp2 + 30*kp4*sigp4 + kp6*sigp6) - 
        	  14*r0p2*sigp6*(480 + 540*kp2*sigp2 + 60*kp4*sigp4 + kp6*sigp6))*ckr - 
     		  k*r0*sigp2*(-448*r0p6 + 112*r0p4*sigp2*(42 + 5*kp2*sigp2) - 28*r0p2*sigp4*(420 + 100*kp2*sigp2 + 3*kp4*sigp4) + 
        	  sigp6*(5880 + 2100*kp2*sigp2 + 126*kp4*sigp4 + kp6*sigp6))*skr))/(15.*ex2s2*sigp8*sigp6)
		  + rad3*(13*amp[w_ind]*(r0*(256*r0p8 - 256*r0p6*sigp2*(18 + 7*kp2*sigp2) + 224*r0p4*sigp4*(108 + 84*kp2*sigp2 + 5*kp4*sigp4) - 
        	  112*r0p2*sigp6*(360 + 420*kp2*sigp2 + 50*kp4*sigp4 + kp6*sigp6) + 
        	  sigp8*(15120 + 23520*kp2*sigp2 + 4200*kp4*sigp4 + 168*kp6*sigp6 + kp8*sigp8))*ckr + 
     		  8*k*sigp2*(128*r0p8 - 224*r0p6*sigp2*(8 + kp2*sigp2) + 56*r0p4*sigp4*(120 + 30*kp2*sigp2 + kp4*sigp4) + 
        	  sigp8*(840 + 420*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
        	  2*r0p2*sigp6*(3360 + 1260*kp2*sigp2 + 84*kp4*sigp4 + kp6*sigp6))*skr))/(768.*ex2s2*sigp8*sigp8)
		  + rad4*(13*amp[w_ind]*((512*r0p8*r0p2 - 2304*r0p8*sigp2*(5 + 2*kp2*sigp2) + 4032*r0p6*sigp4*(20 + 16*kp2*sigp2 + kp4*sigp4) - 
        	  672*r0p4*sigp6*(300 + 360*kp2*sigp2 + 45*kp4*sigp4 + kp6*sigp6) - 
        	  9*sigp8*sigp2*(1680 + 3360*kp2*sigp2 + 840*kp4*sigp4 + 56*kp6*sigp6 + kp8*sigp8) + 
        	  18*r0p2*sigp8*(8400 + 13440*kp2*sigp2 + 2520*kp4*sigp4 + 112*kp6*sigp6 + kp8*sigp8))*ckr + 
     		  k*r0*sigp2*(2304*r0p8 - 768*r0p6*sigp2*(54 + 7*kp2*sigp2) + 2016*r0p4*sigp4*(108 + 28*kp2*sigp2 + kp4*sigp4) - 
        	  144*r0p2*sigp6*(2520 + 980*kp2*sigp2 + 70*kp4*sigp4 + kp6*sigp6) + 
        	  sigp8*(136080 + 70560*kp2*sigp2 + 7560*kp4*sigp4 + 216*kp6*sigp6 + kp8*sigp8))*skr))/(3780.*ex2s2*sigp8*sigp8*sigp2)
		  + rad5*(-67*amp[w_ind]*(r0*(-1024*r0p8*r0p2 + 1280*r0p8*sigp2*(22 + 9*kp2*sigp2) - 1920*r0p6*sigp4*(132 + 108*kp2*sigp2 + 7*kp4*sigp4) + 
        	  3360*r0p4*sigp6*(264 + 324*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
        	  60*r0p2*sigp8*(18480 + 30240*kp2*sigp2 + 5880*kp4*sigp4 + 280*kp6*sigp6 + 3*kp8*sigp8) + 
        	  sigp8*sigp2*(332640 + 680400*kp2*sigp2 + 176400*kp4*sigp4 + 12600*kp6*sigp6 + 270*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr + 
     		  2*k*sigp2*(-2560*r0p8*r0p2 + 3840*r0p8*sigp2*(15 + 2*kp2*sigp2) - 1344*r0p6*sigp4*(300 + 80*kp2*sigp2 + 3*kp4*sigp4) + 
        	  480*r0p4*sigp6*(2100 + 840*kp2*sigp2 + 63*kp4*sigp4 + kp6*sigp6) + 
        	  5*sigp8*sigp2*(15120 + 10080*kp2*sigp2 + 1512*kp4*sigp4 + 72*kp6*sigp6 + kp8*sigp8) - 
        	  10*r0p2*sigp8*(75600 + 40320*kp2*sigp2 + 4536*kp4*sigp4 + 144*kp6*sigp6 + kp8*sigp8))*skr))/
   		  (115200.*ex2s2*sigp8*sigp8*sigp4)
		  + rad5*rad*(amp[w_ind]*((2048*r0p6*r0p6 - 5632*r0p8*r0p2*sigp2*(12 + 5*kp2*sigp2) + 42240*r0p8*sigp4*(18 + 15*kp2*sigp2 + kp4*sigp4) - 
        	  14784*r0p6*sigp6*(240 + 300*kp2*sigp2 + 40*kp4*sigp4 + kp6*sigp6) + 
        	  1320*r0p4*sigp8*(5040 + 8400*kp2*sigp2 + 1680*kp4*sigp4 + 84*kp6*sigp6 + kp8*sigp8) + 
        	  11*sigp8*sigp4*(30240 + 75600*kp2*sigp2 + 25200*kp4*sigp4 + 2520*kp6*sigp6 + 90*kp8*sigp8 + kp8*kp2*sigp8*sigp2) - 
        	  22*r0p2*sigp8*sigp2*(181440 + 378000*kp2*sigp2 + 100800*kp4*sigp4 + 7560*kp6*sigp6 + 180*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr - 
     		  k*r0*sigp2*(-11264*r0p8*r0p2 + 14080*r0p8*sigp2*(22 + 3*kp2*sigp2) - 4224*r0p6*sigp4*(660 + 180*kp2*sigp2 + 7*kp4*sigp4) + 
        	  1056*r0p4*sigp6*(9240 + 3780*kp2*sigp2 + 294*kp4*sigp4 + 5*kp6*sigp6) - 
        	  220*r0p2*sigp8*(55440 + 30240*kp2*sigp2 + 3528*kp4*sigp4 + 120*kp6*sigp6 + kp8*sigp8) + 
        	  sigp8*sigp2*(3659040 + 2494800*kp2*sigp2 + 388080*kp4*sigp4 + 19800*kp6*sigp6 + 330*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*skr))/
 		  (11880.*ex2s2*sigp8*sigp8*sigp6);

	cofs[3] = (-3*amp[w_ind]*(r0*(16*r0p4 - 8*r0p2*sigp2*(10 + 3*kp2*sigp2) + sigp4*(60 + 36*kp2*sigp2 + kp4*sigp4))*ckr + 
     		  4*k*sigp2*(8*r0p4 + sigp4*(6 + kp2*sigp2) - 2*r0p2*sigp2*(12 + kp2*sigp2))*skr))/(4.*ex2s2*sigp8)
		  + rad*(-8*amp[w_ind]*((32*r0p6 - 80*r0p4*sigp2*(3 + kp2*sigp2) - 5*sigp6*(12 + 12*kp2*sigp2 + kp4*sigp4) + 
        	  10*r0p2*sigp4*(36 + 24*kp2*sigp2 + kp4*sigp4))*ckr + 
     		  k*r0*sigp2*(80*r0p4 - 40*r0p2*sigp2*(10 + kp2*sigp2) + sigp4*(300 + 60*kp2*sigp2 + kp4*sigp4))*skr))/
 		  (15.*ex2s2*sigp8*sigp2)
		  + rad2*(5*amp[w_ind]*(r0*(-64*r0p6 + 48*r0p4*sigp2*(14 + 5*kp2*sigp2) - 60*r0p2*sigp4*(28 + 20*kp2*sigp2 + kp4*sigp4) + 
        	  sigp6*(840 + 900*kp2*sigp2 + 90*kp4*sigp4 + kp6*sigp6))*ckr + 
     		  2*k*sigp2*(-96*r0p6 + 80*r0p4*sigp2*(9 + kp2*sigp2) + 3*sigp6*(60 + 20*kp2*sigp2 + kp4*sigp4) - 
        	  6*r0p2*sigp4*(180 + 40*kp2*sigp2 + kp4*sigp4))*skr))/(24.*ex2s2*sigp8*sigp4)
		  + rad3*(2*amp[w_ind]*(-((128*r0p8 - 224*r0p6*sigp2*(8 + 3*kp2*sigp2) + 280*r0p4*sigp4*(24 + 18*kp2*sigp2 + kp4*sigp4) + 
        	  7*sigp8*(120 + 180*kp2*sigp2 + 30*kp4*sigp4 + kp6*sigp6) - 
       		  14*r0p2*sigp6*(480 + 540*kp2*sigp2 + 60*kp4*sigp4 + kp6*sigp6))*ckr) + 
     		  k*r0*sigp2*(-448*r0p6 + 112*r0p4*sigp2*(42 + 5*kp2*sigp2) - 28*r0p2*sigp4*(420 + 100*kp2*sigp2 + 3*kp4*sigp4) + 
        	  sigp6*(5880 + 2100*kp2*sigp2 + 126*kp4*sigp4 + kp6*sigp6))*skr))/(35.*ex2s2*sigp8*sigp6)
		  + rad4*(-7*amp[w_ind]*(r0*(256*r0p8 - 256*r0p6*sigp2*(18 + 7*kp2*sigp2) + 224*r0p4*sigp4*(108 + 84*kp2*sigp2 + 5*kp4*sigp4) - 
        	  112*r0p2*sigp6*(360 + 420*kp2*sigp2 + 50*kp4*sigp4 + kp6*sigp6) + 
        	  sigp8*(15120 + 23520*kp2*sigp2 + 4200*kp4*sigp4 + 168*kp6*sigp6 + kp8*sigp8))*ckr + 
     		  8*k*sigp2*(128*r0p8 - 224*r0p6*sigp2*(8 + kp2*sigp2) + 56*r0p4*sigp4*(120 + 30*kp2*sigp2 + kp4*sigp4) + 
        	  sigp8*(840 + 420*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
        	  2*r0p2*sigp6*(3360 + 1260*kp2*sigp2 + 84*kp4*sigp4 + kp6*sigp6))*skr))/(576.*ex2s2*sigp8*sigp8)
		  + rad5*(-2*amp[w_ind]*((512*r0p8*r0p2 - 2304*r0p8*sigp2*(5 + 2*kp2*sigp2) + 4032*r0p6*sigp4*(20 + 16*kp2*sigp2 + kp4*sigp4) - 
        	  672*r0p4*sigp6*(300 + 360*kp2*sigp2 + 45*kp4*sigp4 + kp6*sigp6) - 
        	  9*sigp8*sigp2*(1680 + 3360*kp2*sigp2 + 840*kp4*sigp4 + 56*kp6*sigp6 + kp8*sigp8) + 
        	  18*r0p2*sigp8*(8400 + 13440*kp2*sigp2 + 2520*kp4*sigp4 + 112*kp6*sigp6 + kp8*sigp8))*ckr + 
       		  k*r0*sigp2*(2304*r0p8 - 768*r0p6*sigp2*(54 + 7*kp2*sigp2) + 2016*r0p4*sigp4*(108 + 28*kp2*sigp2 + kp4*sigp4) - 
       		  144*r0p2*sigp6*(2520 + 980*kp2*sigp2 + 70*kp4*sigp4 + kp6*sigp6) + 
        	  sigp8*(136080 + 70560*kp2*sigp2 + 7560*kp4*sigp4 + 216*kp6*sigp6 + kp8*sigp8))*skr))/(945.*ex2s2*sigp8*sigp8*sigp2)
		  + rad5*rad*(amp[w_ind]*(r0*(-1024*r0p8*r0p2 + 1280*r0p8*sigp2*(22 + 9*kp2*sigp2) - 1920*r0p6*sigp4*(132 + 108*kp2*sigp2 + 7*kp4*sigp4) + 
        	  3360*r0p4*sigp6*(264 + 324*kp2*sigp2 + 42*kp4*sigp4 + kp6*sigp6) - 
        	  60*r0p2*sigp8*(18480 + 30240*kp2*sigp2 + 5880*kp4*sigp4 + 280*kp6*sigp6 + 3*kp8*sigp8) + 
        	  sigp8*sigp2*(332640 + 680400*kp2*sigp2 + 176400*kp4*sigp4 + 12600*kp6*sigp6 + 270*kp8*sigp8 + kp8*kp2*sigp8*sigp2))*ckr + 
     		  2*k*sigp2*(-2560*r0p8*r0p2 + 3840*r0p8*sigp2*(15 + 2*kp2*sigp2) - 1344*r0p6*sigp4*(300 + 80*kp2*sigp2 + 3*kp4*sigp4) + 
      		  480*r0p4*sigp6*(2100 + 840*kp2*sigp2 + 63*kp4*sigp4 + kp6*sigp6) + 
        	  5*sigp8*sigp2*(15120 + 10080*kp2*sigp2 + 1512*kp4*sigp4 + 72*kp6*sigp6 + kp8*sigp8) - 
        	  10*r0p2*sigp8*(75600 + 40320*kp2*sigp2 + 4536*kp4*sigp4 + 144*kp6*sigp6 + kp8*sigp8))*skr))/
   		  (3200.*ex2s2*sigp8*sigp8*sigp4);
	
	cofs[4] = rad*cofs[0]/3.;
	//CCTK_INFO("*DID* the Taylor expansion.");

  }
  else {
	CCTK_WARN(0,"Invalid wave packet for Taylor Expansion.");
  }

}
