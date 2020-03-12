// TwoPunctures:  File  "Newton.c"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "cctk_Parameters.h"
#include "TP_utilities.h"
#include "TwoPunctures.h"

// functions pointers describing the phyical system to Newton.c
void (*NonLinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	              double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
void (*LinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
void (*makeInitialGuess)(double x, double y, double z, double *u, double time);
void *NewtonParam;

// local function prototypes
static void relax (double *dv, int nvar, int n1, int n2, int n3, double *rhs,
	    int *ncols, int **cols, double **JFD);

// -----------------------------------------------------------------------------------
double
norm_inf (double *F, int ntotal)
{
  CCTK_REAL dmax = -1;
#pragma omp parallel
  {
    CCTK_REAL dmax1 = -1;
#pragma omp for
    for (int j = 0; j < ntotal; j++)
      if (fabs (F[j]) > dmax1)
        dmax1 = fabs (F[j]);
#pragma omp critical
    if (dmax1 > dmax)
      dmax = dmax1;
  }
  return dmax;
}
// -----------------------------------------------------------------------------------
void
resid (double *res, int ntotal, double *dv, double *rhs,
       int *ncols, int **cols, double **JFD)
{
#pragma omp parallel for
  for (int i = 0; i < ntotal; i++)
  {
    double JFDdv_i = 0;
    for (int m = 0; m < ncols[i]; m++)
      JFDdv_i += JFD[i][m] * dv[cols[i][m]];
    res[i] = rhs[i] - JFDdv_i;
  }
}

// -----------------------------------------------------------------------------------
void
LineRelax_al (double *dv, int j, int k, int nvar, int n1, int n2, int n3,
	      double *rhs, int *ncols, int **cols, double **JFD)
{
  int i, m, Ic, Ip, Im, col, ivar;
  double *a, *b, *c, *r, *u;
  a = dvector (0, n1 - 1);
  b = dvector (0, n1 - 1);
  c = dvector (0, n1 - 1);
  r = dvector (0, n1 - 1);
  u = dvector (0, n1 - 1);

  for (ivar = 0; ivar < nvar; ivar++)
  {
    for (i = 0; i < n1; i++)
    {
      Ip = Index (ivar, i + 1, j, k, nvar, n1, n2, n3);
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index (ivar, i - 1, j, k, nvar, n1, n2, n3);
      a[i] = 0;
      b[i] = 0;
      c[i] = 0;
      r[i] = rhs[Ic];
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
	  r[i] -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im)
	    a[i] = JFD[Ic][m];
	  if (col == Ic)
	    b[i] = JFD[Ic][m];
	  if (col == Ip)
	    c[i] = JFD[Ic][m];
	}
      }
    }
    tridag (a, b, c, r, u, n1);
    for (i = 0; i < n1; i++)
    {
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = u[i];
    }
  }
  free_dvector (a, 0, n1 - 1);
  free_dvector (b, 0, n1 - 1);
  free_dvector (c, 0, n1 - 1);
  free_dvector (r, 0, n1 - 1);
  free_dvector (u, 0, n1 - 1);
}

// -----------------------------------------------------------------------------------
void
LineRelax_be (double *dv, int i, int k, int nvar, int n1, int n2, int n3,
	      double *rhs, int *ncols, int **cols, double **JFD)
{
  int j, m, Ic, Ip, Im, col, ivar;
  double *a, *b, *c, *r, *u;
  a = dvector (0, n2 - 1);
  b = dvector (0, n2 - 1);
  c = dvector (0, n2 - 1);
  r = dvector (0, n2 - 1);
  u = dvector (0, n2 - 1);

  for (ivar = 0; ivar < nvar; ivar++)
  {
    for (j = 0; j < n2; j++)
    {
      Ip = Index (ivar, i, j + 1, k, nvar, n1, n2, n3);
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      Im = Index (ivar, i, j - 1, k, nvar, n1, n2, n3);
      a[j] = 0;
      b[j] = 0;
      c[j] = 0;
      r[j] = rhs[Ic];
      for (m = 0; m < ncols[Ic]; m++)
      {
	col = cols[Ic][m];
	if (col != Ip && col != Ic && col != Im)
	  r[j] -= JFD[Ic][m] * dv[col];
	else
	{
	  if (col == Im)
	    a[j] = JFD[Ic][m];
	  if (col == Ic)
	    b[j] = JFD[Ic][m];
	  if (col == Ip)
	    c[j] = JFD[Ic][m];
	}
      }
    }
    tridag (a, b, c, r, u, n2);
    for (j = 0; j < n2; j++)
    {
      Ic = Index (ivar, i, j, k, nvar, n1, n2, n3);
      dv[Ic] = u[j];
    }
  }
  free_dvector (a, 0, n2 - 1);
  free_dvector (b, 0, n2 - 1);
  free_dvector (c, 0, n2 - 1);
  free_dvector (r, 0, n2 - 1);
  free_dvector (u, 0, n2 - 1);
}

// -----------------------------------------------------------------------------------
static void
relax (double *dv, int nvar, int n1, int n2, int n3,
       double *rhs, int *ncols, int **cols, double **JFD)
{
  for (int k = 0; k < n3; k = k + 2)
  {
    for (int n = 0; n < N_PlaneRelax; n++)
    {
#pragma omp parallel for schedule(dynamic)
      for (int i = 2; i < n1; i = i + 2)
	LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int i = 1; i < n1; i = i + 2)
	LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int j = 1; j < n2; j = j + 2)
	LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < n2; j = j + 2)
	LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
  for (int k = 1; k < n3; k = k + 2)
  {
    for (int n = 0; n < N_PlaneRelax; n++)
    {
#pragma omp parallel for schedule(dynamic)
      for (int i = 0; i < n1; i = i + 2)
	LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int i = 1; i < n1; i = i + 2)
	LineRelax_be (dv, i, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int j = 1; j < n2; j = j + 2)
	LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < n2; j = j + 2)
	LineRelax_al (dv, j, k, nvar, n1, n2, n3, rhs, ncols, cols, JFD);
    }
  }
}

// -----------------------------------------------------------------------------------
void
TestRelax (int nvar, int n1, int n2, int n3, derivs v,
	   double *dv, double time)
{
  DECLARE_CCTK_PARAMETERS;
  int ntotal = n1 * n2 * n3 * nvar, **cols, *ncols, maxcol =
    StencilSize * nvar, j;
  double *F, *res, **JFD;
  derivs u;

  F = dvector (0, ntotal - 1);
  res = dvector (0, ntotal - 1);
  allocate_derivs (&u, ntotal);

  JFD = dmatrix (0, ntotal - 1, 0, maxcol - 1);
  cols = imatrix (0, ntotal - 1, 0, maxcol - 1);
  ncols = ivector (0, ntotal - 1);

  F_of_v (nvar, n1, n2, n3, v, F, u, time);

  SetMatrix_JFD (nvar, n1, n2, n3, u, ncols, cols, JFD,time);

  for (j = 0; j < ntotal; j++)
    dv[j] = 0;
  resid (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("Before: |F|=%20.15e\n", norm1 (res, ntotal));
  fflush(stdout);
  for (j = 0; j < NRELAX; j++)
  {
    relax (dv, nvar, n1, n2, n3, F, ncols, cols, JFD);	// solves JFD*sh = s
    if (j % Step_Relax == 0)
    {
      resid (res, ntotal, dv, F, ncols, cols, JFD);
      printf ("j=%d\t |F|=%20.15e\n", j, norm1 (res, ntotal));
      fflush(stdout);
    }
  }

  resid (res, ntotal, dv, F, ncols, cols, JFD);
  printf ("After: |F|=%20.15e\n", norm1 (res, ntotal));
  fflush(stdout);

  free_dvector (F, 0, ntotal - 1);
  free_dvector (res, 0, ntotal - 1);
  free_derivs (&u, ntotal);

  free_dmatrix (JFD, 0, ntotal - 1, 0, maxcol - 1);
  free_imatrix (cols, 0, ntotal - 1, 0, maxcol - 1);
  free_ivector (ncols, 0, ntotal - 1);
}

// -----------------------------------------------------------------------------------
int
bicgstab (int nvar, int n1, int n2, int n3, derivs v,
	  derivs dv, int output, int itmax, double tol, double *normres, double time)
{
  DECLARE_CCTK_PARAMETERS;
  int ntotal = n1 * n2 * n3 * nvar, ii;
  double alpha = 0, beta = 0;
  double rho = 0, rho1 = 1, rhotol = 1e-50;
  double omega = 0, omegatol = 1e-50;
  double *p, *rt, *s, *t, *r, *vv;
  double **JFD;
  int **cols, *ncols, maxcol = StencilSize * nvar;
  double *F;
  derivs u, ph, sh;

  F = dvector (0, ntotal - 1);
  allocate_derivs (&u, ntotal);

  JFD = dmatrix (0, ntotal - 1, 0, maxcol - 1);
  cols = imatrix (0, ntotal - 1, 0, maxcol - 1);
  ncols = ivector (0, ntotal - 1);

  F_of_v (nvar, n1, n2, n3, v, F, u, time);
  SetMatrix_JFD (nvar, n1, n2, n3, u, ncols, cols, JFD, time);

  /* temporary storage */
  r = dvector (0, ntotal - 1);
  p = dvector (0, ntotal - 1);
  allocate_derivs (&ph, ntotal);
//      ph  = dvector(0, ntotal-1);
  rt = dvector (0, ntotal - 1);
  s = dvector (0, ntotal - 1);
  allocate_derivs (&sh, ntotal);
//      sh  = dvector(0, ntotal-1);
  t = dvector (0, ntotal - 1);
  vv = dvector (0, ntotal - 1);

  /* check */
  if (output == 1) {
    printf ("bicgstab:  itmax %d, tol %e\n", itmax, tol);
    fflush(stdout);
  }

  /* compute initial residual rt = r = F - J*dv */
  J_times_dv (nvar, n1, n2, n3, dv, r, u,time);
#pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < ntotal; j++)
    rt[j] = r[j] = F[j] - r[j];

  *normres = norm2 (r, ntotal);
  if (output == 1) {
    printf ("bicgstab: %5d  %10.3e\n", 0, *normres);
    fflush(stdout);
  }

  if (*normres <= tol)
    return 0;

  /* cgs iteration */
  for (ii = 0; ii < itmax; ii++)
  {
    rho = scalarproduct (rt, r, ntotal);
    if (fabs (rho) < rhotol)
      break;

    /* compute direction vector p */
    if (ii == 0)
#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < ntotal; j++)
	p[j] = r[j];
    else
    {
      beta = (rho / rho1) * (alpha / omega);
#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < ntotal; j++)
	p[j] = r[j] + beta * (p[j] - omega * vv[j]);
    }

    /* compute direction adjusting vector ph and scalar alpha */
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < ntotal; j++)
      ph.d0[j] = 0;
    for (int j = 0; j < NRELAX; j++)	// solves JFD*ph = p by relaxation
      relax (ph.d0, nvar, n1, n2, n3, p, ncols, cols, JFD);

    J_times_dv (nvar, n1, n2, n3, ph, vv, u, time);	// vv=J*ph
    alpha = rho / scalarproduct (rt, vv, ntotal);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < ntotal; j++)
      s[j] = r[j] - alpha * vv[j];

    /* early check of tolerance */
    *normres = norm2 (s, ntotal);
    if (*normres <= tol)
    {
#pragma omp parallel for schedule(dynamic)
      for (int j = 0; j < ntotal; j++)
	dv.d0[j] += alpha * ph.d0[j];
      if (output == 1) {
	printf ("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n",
		ii + 1, *normres, alpha, beta, omega);
        fflush(stdout);
      }
      break;
    }

    /* compute stabilizer vector sh and scalar omega */
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < ntotal; j++)
      sh.d0[j] = 0;
    for (int j = 0; j < NRELAX; j++)	// solves JFD*sh = s by relaxation
      relax (sh.d0, nvar, n1, n2, n3, s, ncols, cols, JFD);

    J_times_dv (nvar, n1, n2, n3, sh, t, u, time);	// t=J*sh
    omega = scalarproduct (t, s, ntotal) / scalarproduct (t, t, ntotal);

    /* compute new solution approximation */
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < ntotal; j++)
    {
      dv.d0[j] += alpha * ph.d0[j] + omega * sh.d0[j];
      r[j] = s[j] - omega * t[j];
    }
    /* are we done? */
    *normres = norm2 (r, ntotal);
    if (output == 1) {
      printf ("bicgstab: %5d  %10.3e  %10.3e  %10.3e  %10.3e\n",
	      ii + 1, *normres, alpha, beta, omega);
      fflush(stdout);
    }
    if (*normres <= tol)
      break;
    rho1 = rho;
    if (fabs (omega) < omegatol)
      break;

  }

  /* free temporary storage */
  free_dvector (r, 0, ntotal - 1);
  free_dvector (p, 0, ntotal - 1);
//      free_dvector(ph,  0, ntotal-1);
  free_derivs (&ph, ntotal);
  free_dvector (rt, 0, ntotal - 1);
  free_dvector (s, 0, ntotal - 1);
//      free_dvector(sh,  0, ntotal-1);
  free_derivs (&sh, ntotal);
  free_dvector (t, 0, ntotal - 1);
  free_dvector (vv, 0, ntotal - 1);

  free_dvector (F, 0, ntotal - 1);
  free_derivs (&u, ntotal);

  free_dmatrix (JFD, 0, ntotal - 1, 0, maxcol - 1);
  free_imatrix (cols, 0, ntotal - 1, 0, maxcol - 1);
  free_ivector (ncols, 0, ntotal - 1);

  /* iteration failed */
  if (ii > itmax)
    return -1;

  /* breakdown */
  if (fabs (rho) < rhotol)
    return -10;
  if (fabs (omega) < omegatol)
    return -11;

  /* success! */
  return ii + 1;
}


// -------------------------------------------------------------------
int
Newton (int nvar, int n1, int n2, int n3,
	derivs v, double tol, int itmax, double time)
{
  DECLARE_CCTK_PARAMETERS;
  int ntotal = n1 * n2 * n3 * nvar, ii, it;
  double *F, dmax, normres;
  derivs u, dv;

  F = dvector (0, ntotal - 1);
  allocate_derivs (&dv, ntotal);
  allocate_derivs (&u, ntotal);

/* 	TestRelax(nvar, n1, n2, n3, v, dv.d0,time); */
  getInitialGuess(nvar, n1, n2, n3, v, time);  //get initial guess for v (U in paper) = (Psi - Psi_guess)/(A-1)

  it = 0;
  dmax = 1;
  while (dmax > tol && it < itmax)
  {
    if (it == 0)
    {
      F_of_v (nvar, n1, n2, n3, v, F, u, time); //Compute the cartesian derivatives of u using spectral erivatives of v  and store hamiltonian constraint violations in F (with some additional terms)
      
      dmax = norm_inf (F, ntotal); 
    }
#pragma omp parallel for
    for (int j = 0; j < ntotal; j++)
      dv.d0[j] = 0;

    printf ("Newton: it=%d \t |F|=%e\n", it, dmax);
    fflush(stdout);
    ii =
      bicgstab (nvar, n1, n2, n3, v, dv, verbose, 100, dmax * 1.e-3, &normres, time);
#pragma omp parallel for schedule(dynamic)
    for (int j = 0; j < ntotal; j++)
      v.d0[j] -= dv.d0[j];
    F_of_v (nvar, n1, n2, n3, v, F, u, time);
    dmax = norm_inf(F, ntotal);
    it += 1;
  }
  printf ("Newton: it=%d \t |F|=%e \n", it, dmax);
  fflush(stdout);

  free_dvector (F, 0, ntotal - 1);
  free_derivs (&dv, ntotal);
  free_derivs (&u, ntotal);

  return dmax <= tol;
}

// -------------------------------------------------------------------
void getInitialGuess(int nvar, int n1, int n2, int n3, derivs v, double time)
{
  int i, j, k, ivar, indx;
  double al, be, A, B, x, phi, y, z, Am1;
  double *U;

  U = dvector (0, nvar - 1);

  // loop over all points and variables
  for (i = 0; i < n1; i++)
  {
    for (j = 0; j < n2; j++)
    {
      for (k = 0; k < n3; k++)
      {

	al = Pih * (2 * i + 1) / n1;
	A = -cos (al);   		//BK: Why is A cosine function instead of sin^2 function?
	be = Pih * (2 * j + 1) / n2;
	B = -cos (be);
	phi = 2. * Pi * k / n3;

	// translate to isotropic coordinates
	AB3_To_xyz(A, B, phi, &x, &y, &z);
    	  
	// make our guess
        makeInitialGuess(x, y, z, U, time); //BK: Psi-Psi_guess 

	Am1 = A - 1;
	for (ivar = 0; ivar < nvar; ivar++)
	{
	  indx = Index (ivar, i, j, k, nvar, n1, n2, n3);

	  // U needs no further translation (but U_{,i} etc. would)
          v.d0[indx] = U[ivar] / Am1;	 //BK: Notation user here - v.d0 -> U (of paper), U[ivar] -> u (of paper)
	}
      }
    }
  }

  free_dvector (U, 0, nvar - 1);
}

// -------------------------------------------------------------------
