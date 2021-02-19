// TwoPunctures:  File  "CoordTransf.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <time.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

//---------------------------------------------------------------
void
AB_To_XR (int nvar, double A, double B, double *X, double *R,
	  derivs U)
// On Entrance: U.d0[]=U[]; U.d1[] =U[]_A;  U.d2[] =U[]_B;  U.d3[] =U[]_3;  
//                          U.d11[]=U[]_AA; U.d12[]=U[]_AB; U.d13[]=U[]_A3; 
//                          U.d22[]=U[]_BB; U.d23[]=U[]_B3; U.d33[]=U[]_33; 
// At Exit:     U.d0[]=U[]; U.d1[] =U[]_X;  U.d2[] =U[]_R;  U.d3[] =U[]_3;  
//                          U.d11[]=U[]_XX; U.d12[]=U[]_XR; U.d13[]=U[]_X3; 
//                          U.d22[]=U[]_RR; U.d23[]=U[]_R3; U.d33[]=U[]_33; 
{
  DECLARE_CCTK_PARAMETERS;
  double At = 0.5 * (A + 1), A_X, A_XX, B_R, B_RR;
  int ivar;

  *X = 2 * atanh (At);
  *R = Pih + 2 * atan (B);

  A_X = 1 - At * At;
  A_XX = -At * A_X;
  B_R = 0.5 * (1 + B * B);
  B_RR = B * B_R;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    U.d11[ivar] = A_X * A_X * U.d11[ivar] + A_XX * U.d1[ivar];
    U.d12[ivar] = A_X * B_R * U.d12[ivar];
    U.d13[ivar] = A_X * U.d13[ivar];
    U.d22[ivar] = B_R * B_R * U.d22[ivar] + B_RR * U.d2[ivar];
    U.d23[ivar] = B_R * U.d23[ivar];
    U.d1[ivar] = A_X * U.d1[ivar];
    U.d2[ivar] = B_R * U.d2[ivar];
  }
}

//---------------------------------------------------------------
void
C_To_c (int nvar, double X, double R, double *x, double *r,
	derivs U)
// On Entrance: U.d0[]=U[]; U.d1[] =U[]_X;  U.d2[] =U[]_R;  U.d3[] =U[]_3;  
//                          U.d11[]=U[]_XX; U.d12[]=U[]_XR; U.d13[]=U[]_X3; 
//                          U.d22[]=U[]_RR; U.d23[]=U[]_R3; U.d33[]=U[]_33; 
// At Exit:     U.d0[]=U[]; U.d1[] =U[]_x;  U.d2[] =U[]_r;  U.d3[] =U[]_3;  
//                          U.d11[]=U[]_xx; U.d12[]=U[]_xr; U.d13[]=U[]_x3; 
//                          U.d22[]=U[]_rr; U.d23[]=U[]_r3; U.d33[]=U[]_33; 
{
  DECLARE_CCTK_PARAMETERS;
  double C_c2, U_cb, U_CB;
  dcomplex C, C_c, C_cc, c, c_C, c_CC, U_c, U_cc, U_C, U_CC, One =
    Complex (1., 0.);
  int ivar;

  C.r = X;
  C.i = R;

  c = RCmul (par_b, Ccosh (C));	// c=b*cosh(C)
  c_C = RCmul (par_b, Csinh (C));
  c_CC = c;

  C_c = Cdiv (One, c_C);
  C_cc = RCmul (-1., Cmul (Cmul (C_c, C_c), Cmul (C_c, c_CC)));
  C_c2 = C_c.r * C_c.r + C_c.i * C_c.i;

  for (ivar = 0; ivar < nvar; ivar++)
  {
    // U_C = 0.5*(U_X3-i*U_R3)
    // U_c = U_C*C_c = 0.5*(U_x3-i*U_r3)
    U_C.r = 0.5 * U.d13[ivar];
    U_C.i = -0.5 * U.d23[ivar];
    U_c = Cmul (U_C, C_c);
    U.d13[ivar] = 2. * U_c.r;
    U.d23[ivar] = -2. * U_c.i;

    // U_C = 0.5*(U_X-i*U_R)
    // U_c = U_C*C_c = 0.5*(U_x-i*U_r)
    U_C.r = 0.5 * U.d1[ivar];
    U_C.i = -0.5 * U.d2[ivar];
    U_c = Cmul (U_C, C_c);
    U.d1[ivar] = 2. * U_c.r;
    U.d2[ivar] = -2. * U_c.i;

    // U_CC = 0.25*(U_XX-U_RR-2*i*U_XR)
    // U_CB = d^2(U)/(dC*d\bar{C}) = 0.25*(U_XX+U_RR)
    U_CC.r = 0.25 * (U.d11[ivar] - U.d22[ivar]);
    U_CC.i = -0.5 * U.d12[ivar];
    U_CB = 0.25 * (U.d11[ivar] + U.d22[ivar]);

    // U_cc = C_cc*U_C+(C_c)^2*U_CC
    U_cb = U_CB * C_c2;
    U_cc = Cadd (Cmul (C_cc, U_C), Cmul (Cmul (C_c, C_c), U_CC));

    // U_xx = 2*(U_cb+Re[U_cc])
    // U_rr = 2*(U_cb-Re[U_cc])
    // U_rx = -2*Im[U_cc]
    U.d11[ivar] = 2 * (U_cb + U_cc.r);
    U.d22[ivar] = 2 * (U_cb - U_cc.r);
    U.d12[ivar] = -2 * U_cc.i;
  }

  *x = c.r;
  *r = c.i;
}

//---------------------------------------------------------------
void
rx3_To_xyz (int nvar, double x, double r, double phi,
	    double *y, double *z, derivs U)
// On Entrance: U.d0[]=U[]; U.d1[] =U[]_x;  U.d2[] =U[]_r;  U.d3[] =U[]_3;  
//                          U.d11[]=U[]_xx; U.d12[]=U[]_xr; U.d13[]=U[]_x3; 
//                          U.d22[]=U[]_rr; U.d23[]=U[]_r3; U.d33[]=U[]_33; 
// At Exit:     U.d0[]=U[]; U.d1[] =U[]_x;  U.d2[] =U[]_y;  U.dz[] =U[]_z;  
//                          U.d11[]=U[]_xx; U.d12[]=U[]_xy; U.d1z[]=U[]_xz; 
//                          U.d22[]=U[]_yy; U.d2z[]=U[]_yz; U.dzz[]=U[]_zz; 
{
  int jvar;
  double
    sin_phi = sin (phi),
    cos_phi = cos (phi),
    sin2_phi = sin_phi * sin_phi,
    cos2_phi = cos_phi * cos_phi,
    sin_2phi = 2 * sin_phi * cos_phi,
    cos_2phi = cos2_phi - sin2_phi, r_inv = 1 / r, r_inv2 = r_inv * r_inv;

  *y = r * cos_phi;
  *z = r * sin_phi;

  for (jvar = 0; jvar < nvar; jvar++)
  {
    double U_x = U.d1[jvar], U_r = U.d2[jvar], U_3 = U.d3[jvar],
      U_xx = U.d11[jvar], U_xr = U.d12[jvar], U_x3 = U.d13[jvar],
      U_rr = U.d22[jvar], U_r3 = U.d23[jvar], U_33 = U.d33[jvar];
    U.d1[jvar] = U_x;		// U_x
    U.d2[jvar] = U_r * cos_phi - U_3 * r_inv * sin_phi;	// U_y
    U.d3[jvar] = U_r * sin_phi + U_3 * r_inv * cos_phi;	// U_z
    U.d11[jvar] = U_xx;		// U_xx
    U.d12[jvar] = U_xr * cos_phi - U_x3 * r_inv * sin_phi;	// U_xy
    U.d13[jvar] = U_xr * sin_phi + U_x3 * r_inv * cos_phi;	// U_xz
    U.d22[jvar] = U_rr * cos2_phi + r_inv2 * sin2_phi * (U_33 + r * U_r)	// U_yy
      + sin_2phi * r_inv2 * (U_3 - r * U_r3);
    U.d23[jvar] = 0.5 * sin_2phi * (U_rr - r_inv * U_r - r_inv2 * U_33)	// U_yz
      - cos_2phi * r_inv2 * (U_3 - r * U_r3);
    U.d33[jvar] = U_rr * sin2_phi + r_inv2 * cos2_phi * (U_33 + r * U_r)	// U_zz
      - sin_2phi * r_inv2 * (U_3 - r * U_r3);
  }
}

//---------------------------------------------------------------
void AB3_To_xyz(double A, double B, double phi, double *x, double *y, double *z)
// strung together version of the above coordinate transformations which only 
// transforms the coordinates
{
  DECLARE_CCTK_PARAMETERS;
  double At = 0.5 * (A + 1);
  double X, R;
  double
    sin_phi = sin (phi),
    cos_phi = cos (phi);
  dcomplex C, c;
  double r;

  X = 2 * atanh (At);
  R = Pih + 2 * atan (B);

  C.r = X;
  C.i = R;

  c = RCmul (par_b, Ccosh (C));	// c=b*cosh(C)

  *x = c.r;
  r = c.i;

  *y = r * cos_phi;
  *z = r * sin_phi;
}

//---------------------------------------------------------------
