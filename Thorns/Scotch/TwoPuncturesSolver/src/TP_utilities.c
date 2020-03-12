// TwoPunctures:  File  "utilities.c"

#include <math.h>
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <assert.h>
#include "TP_utilities.h"

//-----------------------------------------------------------------------------
void
nrerror (char error_text[])
/* Numerical Recipes standard error handler */
{
  fprintf (stderr, "Numerical Recipes run-time error...\n");
  fprintf (stderr, "%s\n", error_text);
  fprintf (stderr, "...now exiting to system...\n");
  exit (1);
}

//-----------------------------------------------------------------------------
int *
ivector (long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
  int *v;

  v = (int *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (int)));
  if (!v)
    nrerror ("allocation failure in ivector()");
  return v - nl + NR_END;
}

//-----------------------------------------------------------------------------
double *
dvector (long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
  double *v;

  v = (double *) malloc ((size_t) ((nh - nl + 1 + NR_END) * sizeof (double)));
  if (!v)
    nrerror ("allocation failure in dvector()");
  return v - nl + NR_END;
}

//-----------------------------------------------------------------------------
int **
imatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  int **m;

  /* allocate pointers to rows */
  m = (int **) malloc ((size_t) ((nrow + NR_END) * sizeof (int *)));
  if (!m)
    nrerror ("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;


  /* allocate rows and set pointers to them */
  m[nrl] = (int *) malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (int)));
  if (!m[nrl])
    nrerror ("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

//-----------------------------------------------------------------------------
double **
dmatrix (long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
  long i, nrow = nrh - nrl + 1, ncol = nch - ncl + 1;
  double **m;

  /* allocate pointers to rows */
  m = (double **) malloc ((size_t) ((nrow + NR_END) * sizeof (double *)));
  if (!m)
    nrerror ("allocation failure 1 in matrix()");
  m += NR_END;
  m -= nrl;

  /* allocate rows and set pointers to them */
  m[nrl] =
    (double *) malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (double)));
  if (!m[nrl])
    nrerror ("allocation failure 2 in matrix()");
  m[nrl] += NR_END;
  m[nrl] -= ncl;

  for (i = nrl + 1; i <= nrh; i++)
    m[i] = m[i - 1] + ncol;

  /* return pointer to array of pointers to rows */
  return m;
}

//-----------------------------------------------------------------------------
double ***
d3tensor (long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
  long i, j, nrow = nrh - nrl + 1, ncol = nch - ncl + 1, ndep = ndh - ndl + 1;
  double ***t;

  /* allocate pointers to pointers to rows */
  t = (double ***) malloc ((size_t) ((nrow + NR_END) * sizeof (double **)));
  if (!t)
    nrerror ("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;

  /* allocate pointers to rows and set pointers to them */
  t[nrl] =
    (double **)
    malloc ((size_t) ((nrow * ncol + NR_END) * sizeof (double *)));
  if (!t[nrl])
    nrerror ("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;

  /* allocate rows and set pointers to them */
  t[nrl][ncl] =
    (double *)
    malloc ((size_t) ((nrow * ncol * ndep + NR_END) * sizeof (double)));
  if (!t[nrl][ncl])
    nrerror ("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;

  for (j = ncl + 1; j <= nch; j++)
    t[nrl][j] = t[nrl][j - 1] + ndep;
  for (i = nrl + 1; i <= nrh; i++)
  {
    t[i] = t[i - 1] + ncol;
    t[i][ncl] = t[i - 1][ncl] + ncol * ndep;
    for (j = ncl + 1; j <= nch; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }

  /* return pointer to array of pointers to rows */
  return t;
}

//-----------------------------------------------------------------------------
void
free_ivector (int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
}

//-----------------------------------------------------------------------------
void
free_dvector (double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
  free ((FREE_ARG) (v + nl - NR_END));
}

//-----------------------------------------------------------------------------
void
free_imatrix (int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
  free ((FREE_ARG) (m[nrl] + ncl - NR_END));
  free ((FREE_ARG) (m + nrl - NR_END));
}

//-----------------------------------------------------------------------------
void
free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
  free ((FREE_ARG) (m[nrl] + ncl - NR_END));
  free ((FREE_ARG) (m + nrl - NR_END));
}

//-----------------------------------------------------------------------------
void
free_d3tensor (double ***t, long nrl, long nrh, long ncl, long nch,
	       long ndl, long ndh)
/* free a double f3tensor allocated by f3tensor() */
{
  free ((FREE_ARG) (t[nrl][ncl] + ndl - NR_END));
  free ((FREE_ARG) (t[nrl] + ncl - NR_END));
  free ((FREE_ARG) (t + nrl - NR_END));
}

//-----------------------------------------------------------------------------
int
minimum2 (int i, int j)
{
  int result = i;
  if (j < result)
    result = j;
  return result;
}

//-----------------------------------------------------------------------------
int
minimum3 (int i, int j, int k)
{
  int result = i;
  if (j < result)
    result = j;
  if (k < result)
    result = k;
  return result;
}

//-----------------------------------------------------------------------------
int
maximum2 (int i, int j)
{
  int result = i;
  if (j > result)
    result = j;
  return result;
}

//-----------------------------------------------------------------------------
int
maximum3 (int i, int j, int k)
{
  int result = i;
  if (j > result)
    result = j;
  if (k > result)
    result = k;
  return result;
}

//-----------------------------------------------------------------------------
int
pow_int (int mantisse, int exponent)
{
  int i, result = 1;

  for (i = 1; i <= exponent; i++)
    result *= mantisse;

  return result;
}

//-----------------------------------------------------------------------------
double
atanh (double x)
{
  return 0.5 * log ((1 + x) / (1 - x));
}

//-----------------------------------------------------------------------------
double
asinh (double x)
{
  return log (x + sqrt (1 + x * x));
}

//-----------------------------------------------------------------------------
double
acosh (double x)
{
  return log (x + sqrt (x * x - 1));
}

//-----------------------------------------------------------------------------
dcomplex
Cadd (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r + b.r;
  c.i = a.i + b.i;
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Csub (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r - b.r;
  c.i = a.i - b.i;
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Cmul (dcomplex a, dcomplex b)
{
  dcomplex c;
  c.r = a.r * b.r - a.i * b.i;
  c.i = a.i * b.r + a.r * b.i;
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
RCmul (double x, dcomplex a)
{
  dcomplex c;
  c.r = x * a.r;
  c.i = x * a.i;
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Cdiv (dcomplex a, dcomplex b)
{
  dcomplex c;
  double r, den;
  if (fabs (b.r) >= fabs (b.i))
  {
    r = b.i / b.r;
    den = b.r + r * b.i;
    c.r = (a.r + r * a.i) / den;
    c.i = (a.i - r * a.r) / den;
  }
  else
  {
    r = b.r / b.i;
    den = b.i + r * b.r;
    c.r = (a.r * r + a.i) / den;
    c.i = (a.i * r - a.r) / den;
  }
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Complex (double re, double im)
{
  dcomplex c;
  c.r = re;
  c.i = im;
  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Conjg (dcomplex z)
{
  dcomplex c;
  c.r = z.r;
  c.i = -z.i;
  return c;
}

//-----------------------------------------------------------------------------
double
Cabs (dcomplex z)
{
  double x, y, ans, temp;
  x = fabs (z.r);
  y = fabs (z.i);
  if (x == 0.0)
    ans = y;
  else if (y == 0.0)
    ans = x;
  else if (x > y)
  {
    temp = y / x;
    ans = x * sqrt (1.0 + temp * temp);
  }
  else
  {
    temp = x / y;
    ans = y * sqrt (1.0 + temp * temp);
  }
  return ans;
}

//-----------------------------------------------------------------------------
dcomplex
Csqrt (dcomplex z)
{
  dcomplex c;
  double x, y, w, r;
  if ((z.r == 0.0) && (z.i == 0.0))
  {
    c.r = 0.0;
    c.i = 0.0;
    return c;
  }
  else
  {
    x = fabs (z.r);
    y = fabs (z.i);
    if (x >= y)
    {
      r = y / x;
      w = sqrt (x) * sqrt (0.5 * (1.0 + sqrt (1.0 + r * r)));
    }
    else
    {
      r = x / y;
      w = sqrt (y) * sqrt (0.5 * (r + sqrt (1.0 + r * r)));
    }
    if (z.r >= 0.0)
    {
      c.r = w;
      c.i = z.i / (2.0 * w);
    }
    else
    {
      c.i = (z.i >= 0) ? w : -w;
      c.r = z.i / (2.0 * c.i);
    }
    return c;
  }
}

//-----------------------------------------------------------------------------
dcomplex
Cexp (dcomplex z)
{
  dcomplex c;
  double exp_r = exp (z.r);

  c.r = exp_r * cos (z.i);
  c.i = exp_r * sin (z.i);

  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Clog (dcomplex z)
{
  dcomplex c;

  c.r = 0.5 * log (z.r * z.r + z.i * z.i);
  c.i = atan2 (z.i, z.r);

  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Csin (dcomplex z)
{
  dcomplex c;

  c.r = sin (z.r) * cosh (z.i);
  c.i = cos (z.r) * sinh (z.i);

  return c;
}				//-----------------------------------------------------------------------------

dcomplex
Ccos (dcomplex z)
{
  dcomplex c;

  c.r = cos (z.r) * cosh (z.i);
  c.i = -sin (z.r) * sinh (z.i);

  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Ctan (dcomplex z)
{
  return Cdiv (Csin (z), Ccos (z));
}

//-----------------------------------------------------------------------------
dcomplex
Ccot (dcomplex z)
{
  return Cdiv (Ccos (z), Csin (z));
}

//-----------------------------------------------------------------------------
dcomplex
Csinh (dcomplex z)
{
  dcomplex c;

  c.r = sinh (z.r) * cos (z.i);
  c.i = cosh (z.r) * sin (z.i);

  return c;
}				//-----------------------------------------------------------------------------

dcomplex
Ccosh (dcomplex z)
{
  dcomplex c;

  c.r = cosh (z.r) * cos (z.i);
  c.i = sinh (z.r) * sin (z.i);

  return c;
}

//-----------------------------------------------------------------------------
dcomplex
Ctanh (dcomplex z)
{
  return Cdiv (Csinh (z), Ccosh (z));
}

//-----------------------------------------------------------------------------
dcomplex
Ccoth (dcomplex z)
{
  return Cdiv (Ccosh (z), Csinh (z));
}

//----------------------------------------------------------------------------
#define COSCACHE_N_MAX 256
static double *
chebft_Zeros_cache (int n)
  /* caches values of cos (Pi/n * j * (k + 0.5)) */
{
  int j, k;
  double Pion;
  static double *cache[COSCACHE_N_MAX];
  assert(n <= COSCACHE_N_MAX);

  if(cache[n-1] == NULL)
  {
#pragma omp critical
    if(cache[n-1] == NULL) /* check again in case some other thread allocated storage in omp critical sync */
    {
      /* set up new cache */
      cache[n-1] = dvector(0, n*n-1);
      /* TODO: add method to free allocated memory */
      Pion = Pi/n;
      for (j = 0; j < n; j++)
        for (k = 0; k < n; k++)
          cache[n-1][j+k*n] = cos (Pion * j * (k + 0.5));
    }
  }

  assert(cache[n-1]);
  return cache[n-1];
}


//-----------------------------------------------------------------------------
void
chebft_Zeros (double u[], int n, int inv)
{
  int k, j, isignum;
  double fac, sum, Pion, *c;
  double *Cos = chebft_Zeros_cache(n);

  c = dvector (0, n);
  Pion = Pi / n;
  if (inv == 0)
  {
    fac = 2.0 / n;
    isignum = 1;
    for (j = 0; j < n; j++)
    {
      sum = 0.0;
      for (k = 0; k < n; k++)
	sum += u[k] * Cos[j+k*n]; /*cos (Pion * j * (k + 0.5))*/
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * Cos[k+j*n] * isignum; /*cos (Pion * (j + 0.5) * k)*/
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
    u[j] = c[j];
  free_dvector (c, 0, n);
}

// -------------------------------------------------------------------------------

void
chebft_Extremes (double u[], int n, int inv)
{
  int k, j, isignum, N = n - 1;
  double fac, sum, PioN, *c;

  c = dvector (0, N);
  PioN = Pi / N;
  if (inv == 0)
  {
    fac = 2.0 / N;
    isignum = 1;
    for (j = 0; j < n; j++)
    {
      sum = 0.5 * (u[0] + u[N] * isignum);
      for (k = 1; k < N; k++)
	sum += u[k] * cos (PioN * j * k);
      c[j] = fac * sum * isignum;
      isignum = -isignum;
    }
    c[N] = 0.5 * c[N];
  }
  else
  {
    for (j = 0; j < n; j++)
    {
      sum = -0.5 * u[0];
      isignum = 1;
      for (k = 0; k < n; k++)
      {
	sum += u[k] * cos (PioN * j * k) * isignum;
	isignum = -isignum;
      }
      c[j] = sum;
    }
  }
  for (j = 0; j < n; j++)
    u[j] = c[j];
  free_dvector (c, 0, N);
}

// -------------------------------------------------------------------------------

void
chder (double *c, double *cder, int n)
{
  int j;

  cder[n] = 0.0;
  cder[n - 1] = 0.0;
  for (j = n - 2; j >= 0; j--)
    cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];
}

// -------------------------------------------------------------------------------
double
chebev (double a, double b, double c[], int m, double x)
{
  double d = 0.0, dd = 0.0, sv, y, y2;
  int j;

  y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a));
  for (j = m - 1; j >= 1; j--)
  {
    sv = d;
    d = y2 * d - dd + c[j];
    dd = sv;
  }
  return y * d - dd + 0.5 * c[0];
}

// -------------------------------------------------------------------------------
void
fourft (double *u, int N, int inv)
{
  int l, k, iy, M;
  double x, x1, fac, Pi_fac, *a, *b;

  M = N / 2;
  a = dvector (0, M);
  b = dvector (1, M);		// Actually: b=vector(1,M-1) but this is problematic if M=1
  fac = 1. / M;
  Pi_fac = Pi * fac;
  if (inv == 0)
  {
    for (l = 0; l <= M; l++)
    {
      a[l] = 0;
      if (l > 0 && l < M)
	b[l] = 0;
      x1 = Pi_fac * l;
      for (k = 0; k < N; k++)
      {
	x = x1 * k;
	a[l] += fac * u[k] * cos (x);
	if (l > 0 && l < M)
	  b[l] += fac * u[k] * sin (x);
      }
    }
    u[0] = a[0];
    u[M] = a[M];
    for (l = 1; l < M; l++)
    {
      u[l] = a[l];
      u[l + M] = b[l];
    }
  }
  else
  {
    a[0] = u[0];
    a[M] = u[M];
    for (l = 1; l < M; l++)
    {
      a[l] = u[l];
      b[l] = u[M + l];
    }
    iy = 1;
    for (k = 0; k < N; k++)
    {
      u[k] = 0.5 * (a[0] + a[M] * iy);
      x1 = Pi_fac * k;
      for (l = 1; l < M; l++)
      {
	x = x1 * l;
	u[k] += a[l] * cos (x) + b[l] * sin (x);
      }
      iy = -iy;
    }
  }
  free_dvector (a, 0, M);
  free_dvector (b, 1, M);
}

/* -----------------------------------------*/
void
fourder (double u[], double du[], int N)
{
  int l, M, lpM;

  M = N / 2;
  du[0] = 0.;
  du[M] = 0.;
  for (l = 1; l < M; l++)
  {
    lpM = l + M;
    du[l] = u[lpM] * l;
    du[lpM] = -u[l] * l;
  }
}

/* -----------------------------------------*/
void
fourder2 (double u[], double d2u[], int N)
{
  int l, l2, M, lpM;

  d2u[0] = 0.;
  M = N / 2;
  for (l = 1; l <= M; l++)
  {
    l2 = l * l;
    lpM = l + M;
    d2u[l] = -u[l] * l2;
    if (l < M)
      d2u[lpM] = -u[lpM] * l2;
  }
}

/* ----------------------------------------- */
double
fourev (double *u, int N, double x)
{
  int l, M = N / 2;
  double xl, result;

  result = 0.5 * (u[0] + u[M] * cos (x * M));
  for (l = 1; l < M; l++)
  {
    xl = x * l;
    result += u[l] * cos (xl) + u[M + l] * sin (xl);
  }
  return result;
}

// -------------------------------------------------------------------------------
void
ludcmp (double **a, int n, int *indx, double *d)
{				// Version of 'ludcmp' of the numerical recipes for
  // matrices a[0:n-1][0:n-1]
  int i, imax, j, k;
  double big, dum, sum, temp;
  double *vv;

  vv = dvector (0, n - 1);
  *d = 1.0;
  for (i = 0; i < n; i++)
  {
    big = 0.0;
    for (j = 0; j < n; j++)
      if ((temp = fabs (a[i][j])) > big)
	big = temp;
    if (big == 0.0)
      nrerror ("Singular matrix in routine ludcmp");
    vv[i] = 1.0 / big;
  }
  for (j = 0; j < n; j++)
  {
    for (i = 0; i < j; i++)
    {
      sum = a[i][j];
      for (k = 0; k < i; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
    }
    big = 0.0;
    for (i = j; i < n; i++)
    {
      sum = a[i][j];
      for (k = 0; k < j; k++)
	sum -= a[i][k] * a[k][j];
      a[i][j] = sum;
      if ((dum = vv[i] * fabs (sum)) >= big)
      {
	big = dum;
	imax = i;
      }
    }
    if (j != imax)
    {
      for (k = 0; k < n; k++)
      {
	dum = a[imax][k];
	a[imax][k] = a[j][k];
	a[j][k] = dum;
      }
      *d = -(*d);
      vv[imax] = vv[j];
    }
    indx[j] = imax;
    if (a[j][j] == 0.0)
      a[j][j] = TINY;
    if (j != n)
    {
      dum = 1.0 / (a[j][j]);
      for (i = j + 1; i < n; i++)
	a[i][j] *= dum;
    }
  }
  free_dvector (vv, 0, n - 1);
}

// -------------------------------------------------------------------------------
void
lubksb (double **a, int n, int *indx, double b[])
{				// Version of 'lubksb' of the numerical recipes for
  // matrices a[0:n-1][0:n-1] and vectors b[0:n-1]

  int i, ii = 0, ip, j;
  double sum;

  for (i = 0; i < n; i++)
  {
    ip = indx[i];
    sum = b[ip];
    b[ip] = b[i];
    if (ii)
      for (j = ii; j <= i - 1; j++)
	sum -= a[i][j] * b[j];
    else if (sum)
      ii = i;
    b[i] = sum;
  }
  for (i = n - 1; i >= 0; i--)
  {
    sum = b[i];
    for (j = i + 1; j <= n; j++)
      sum -= a[i][j] * b[j];
    b[i] = sum / a[i][i];
  }
}

// -----------------------------------------------------------------------------------
void
tridag (double a[], double b[], double c[], double r[], double u[], int n)
{				// Version of 'tridag' of the numerical recipes for
  // vectors a, b, c, r, u with indices in the range [0:n-1]
  int j;
  double bet, *gam;

  gam = dvector (0, n - 1);
  if (b[0] == 0.0)
    nrerror ("Error 1 in tridag");
  u[0] = r[0] / (bet = b[0]);
  for (j = 1; j < n; j++)
  {
    gam[j] = c[j - 1] / bet;
    bet = b[j] - a[j] * gam[j];
    if (bet == 0.0)
      nrerror ("Error 2 in tridag");
    u[j] = (r[j] - a[j] * u[j - 1]) / bet;
  }
  for (j = (n - 2); j >= 0; j--)
    u[j] -= gam[j + 1] * u[j + 1];
  free_dvector (gam, 0, n - 1);
}

// -----------------------------------------------------------------------------------
double
norm1 (double *v, int n)
{
  int i;
  double result = -1;

  for (i = 0; i < n; i++)
    if (fabs (v[i]) > result)
      result = fabs (v[i]);

  return result;
}

// -----------------------------------------------------------------------------------
double
norm2 (double *v, int n)
{
  int i;
  double result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * v[i];

  return sqrt (result);
}

// -----------------------------------------------------------------------------------
double
scalarproduct (double *v, double *w, int n)
{
  int i;
  double result = 0;

  for (i = 0; i < n; i++)
    result += v[i] * w[i];

  return result;
}

// -----------------------------------------------------------------------------------
double
plgndr (int l, int m, double x)
{
  void nrerror (char error_text[]);
  double fact, pll, pmm, pmmp1, somx2;
  int i, ll;

  if (m < 0 || m > l || fabs (x) > 1.0)
    nrerror ("Bad arguments in routine plgndr");
  pmm = 1.0;
  if (m > 0)
  {
    somx2 = sqrt ((1.0 - x) * (1.0 + x));
    fact = 1.0;
    for (i = 1; i <= m; i++)
    {
      pmm *= -fact * somx2;
      fact += 2.0;
    }
  }
  if (l == m)
    return pmm;
  else
  {
    pmmp1 = x * (2 * m + 1) * pmm;
    if (l == (m + 1))
      return pmmp1;
    else
    {
      for (ll = m + 2; ll <= l; ll++)
      {
	pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm) / (ll - m);
	pmm = pmmp1;
	pmmp1 = pll;
      }
      return pll;
    }
  }
}

// -----------------------------------------------------------------------------------

void TP_FullTransform(const double *data, double *coeffs, int nvar, int n1, int n2, int n3)
{
    int N;
    double *p;

    // temp buffer for 1d transforms
    N = maximum3(n1,n2,n3);
    p = malloc(sizeof(*p)*N*nvar);
    assert(p);

    // transform from f(A_i,B_j,phi_k) to f_k(A_i,B_j)
    // copy from data to coeffs array
    for(int j = 0 ; j < n2 ; j++)
    {
        for(int i = 0 ; i < n1 ; i++)
        {
            for(int ivar = 0 ; ivar < nvar ; ivar ++)
            {
                for(int k = 0 ; k < n3 ; k++)
                    p[k] = data[ivar + nvar * (i + n1 * (j + n2 * k))];
                fourft (p, n3, 0);
                for(int k = 0 ; k < n3 ; k++)
                    coeffs[ivar + nvar * (i + n1 * (j + n2 * k))] = p[k];
            }
        }
    }

    // transform from f_k(A_i,B_j) to f_jk(A_i)
    for(int i = 0 ; i < n1 ; i++)
    {
        for(int k = 0 ; k < n3 ; k++)
        {
            for(int ivar = 0 ; ivar < nvar ; ivar ++)
            {
                for(int j = 0 ; j < n2 ; j++)
                    p[j] = coeffs[ivar + nvar * (i + n1 * (j + n2 * k))];
                chebft_Zeros (p, n2, 0);
                for(int j = 0 ; j < n2 ; j++)
                    coeffs[ivar + nvar * (i + n1 * (j + n2 * k))] = p[j];
            }
        }
    }

    // transform from f_jk(A_i) to f_ijk
    for(int k = 0 ; k < n3 ; k++)
    {
        for(int j = 0 ; j < n2 ; j++)
        {
            for(int ivar = 0 ; ivar < nvar ; ivar ++)
            {
                for(int i = 0 ; i < n1 ; i++)
                    p[i] = coeffs[ivar + nvar * (i + n1 * (j + n2 * k))];
                chebft_Zeros (p, n1, 0);
                for(int i = 0 ; i < n1 ; i++)
                    coeffs[ivar + nvar * (i + n1 * (j + n2 * k))] = p[i];
            }
        }
    }

    // get rid of temp. memory
    free(p);
}

double traceSquare(double Xij[6])
{
	double trace=0.;
	double Fij[3][3];
	unflattenAij(Xij, Fij);
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j< 3; j++)
		{
			trace += Fij[i][j]*Fij[i][j];
		}
	}
	return trace;
}

void unflattenAij(double flat[3], double unflat[3][3])
{
	unflat[0][0] = flat[0];
	unflat[0][1] = flat[1];
	unflat[1][0] = flat[1];
	unflat[0][2] = flat[2];
	unflat[2][0] = flat[2];
	unflat[1][1] = flat[3];
	unflat[1][2] = flat[4];
	unflat[2][1] = flat[4];
	unflat[2][2] = flat[5];
}

