// TwoPunctures:  File  "utilities.h"

#define Pi  3.14159265358979323846264338328
#define Pih 1.57079632679489661923132169164	// Pi/2
#define Piq 0.78539816339744830961566084582	// Pi/4

#define TINY 1.0e-20
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 1
#define FREE_ARG char*

#ifdef __cplusplus
extern "C" {
#endif

typedef struct DCOMPLEX
{
  double r, i;
} dcomplex;

#define nrerror TP_nrerror
#define ivector TP_ivector
#define dvector TP_dvector
#define imatrix TP_imatrix
#define dmatrix TP_dmatrix
#define d3tensor TP_d3tensor
#define free_ivector TP_free_ivector
#define free_dvector TP_free_dvector
#define free_imatrix TP_free_imatrix
#define free_dmatrix TP_free_dmatrix
#define free_d3tensor TP_free_d3tensor

void nrerror (char error_text[]);
int *ivector (long nl, long nh);
double *dvector (long nl, long nh);
int **imatrix (long nrl, long nrh, long ncl, long nch);
double **dmatrix (long nrl, long nrh, long ncl, long nch);
double ***d3tensor (long nrl, long nrh, long ncl, long nch, long ndl,
		    long ndh);
void free_ivector (int *v, long nl, long nh);
void free_dvector (double *v, long nl, long nh);
void free_imatrix (int **m, long nrl, long nrh, long ncl, long nch);
void free_dmatrix (double **m, long nrl, long nrh, long ncl, long nch);
void free_d3tensor (double ***t, long nrl, long nrh, long ncl, long nch,
		    long ndl, long ndh);

int minimum2 (int i, int j);
int minimum3 (int i, int j, int k);
int maximum2 (int i, int j);
int maximum3 (int i, int j, int k);
int pow_int (int mantisse, int exponent);
double atanh (double x);
double asinh (double x);
double acosh (double x);

dcomplex Cadd (dcomplex a, dcomplex b);
dcomplex Csub (dcomplex a, dcomplex b);
dcomplex Cmul (dcomplex a, dcomplex b);
dcomplex RCmul (double x, dcomplex a);
dcomplex Cdiv (dcomplex a, dcomplex b);
dcomplex Complex (double re, double im);
dcomplex Conjg (dcomplex z);
double Cabs (dcomplex z);

dcomplex Csqrt (dcomplex z);
dcomplex Cexp (dcomplex z);
dcomplex Clog (dcomplex z);
dcomplex Csin (dcomplex z);
dcomplex Ccos (dcomplex z);
dcomplex Ctan (dcomplex z);
dcomplex Ccot (dcomplex z);
dcomplex Csinh (dcomplex z);
dcomplex Ccosh (dcomplex z);
dcomplex Ctanh (dcomplex z);
dcomplex Ccoth (dcomplex z);

void chebft_Zeros (double u[], int n, int inv);
void chebft_Extremes (double u[], int n, int inv);
void chder (double *c, double *cder, int n);
double chebev (double a, double b, double c[], int m, double x);
void fourft (double *u, int N, int inv);
void fourder (double u[], double du[], int N);
void fourder2 (double u[], double d2u[], int N);
double fourev (double *u, int N, double x);


void ludcmp (double **a, int n, int *indx, double *d);
void lubksb (double **a, int n, int *indx, double b[]);
void tridag (double a[], double b[], double c[], double r[], double u[],
	     int n);
double norm1 (double *v, int n);
double norm2 (double *v, int n);
double scalarproduct (double *v, double *w, int n);

void TP_FullTransform(const double *data, double *coeffs, int nvar, int n1, int n2, int n3);

double plgndr (int l, int m, double x);

double traceSquare(double Xij[6]);
void unflattenAij(double flat[3], double unflat[3][3]);

#ifdef __cplusplus
}
#endif
