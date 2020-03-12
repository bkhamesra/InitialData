// TwoPunctures:  File  "TwoPunctures.h"

#include "cctk.h"
#include "cctk_Arguments.h"
#define StencilSize 19
#define N_PlaneRelax 1
#define NRELAX 200
#define Step_Relax 1

typedef struct DERIVS
{
  double *d0, *d1, *d2, *d3, *d11, *d12, *d13, *d22, *d23, *d33;
} derivs;
  
enum GRID_SETUP_METHOD { GSM_Taylor_expansion, GSM_evaluation };
typedef enum GRID_SETUP_METHOD gsm_type;

/*
Files of "TwoPunctures":
	TwoPunctures.c
	FuncAndJacobian.c
	CoordTransf.c
	Equations.c
	Newton.c
	utilities.c (see utilities.h)
**************************
*/

// Routines in  "TwoPunctures.c"
double TestSolution (double A, double B, double X, double R, double phi);
void TestVector_w (double *par, int nvar, int n1, int n2, int n3, double *w);
CCTK_REAL TwoPuncturesSolver_GetU(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
CCTK_REAL TwoPuncturesSolver_GetConformalFactor(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
void TwoPuncturesSolver_ForceResolve(CCTK_INT resolve);

// Routines in  "FuncAndJacobian.c"
int Index (int ivar, int i, int j, int k, int nvar, int n1, int n2, int n3);
void allocate_derivs (derivs * v, int n);
void free_derivs (derivs * v, int n);
void Derivatives_AB3 (int nvar, int n1, int n2, int n3, derivs v);
void F_of_v (int nvar, int n1, int n2, int n3, derivs v,
	     double *F, derivs u, double time);
void J_times_dv (int nvar, int n1, int n2, int n3, derivs dv,
		 double *Jdv, derivs u, double time);
void JFD_times_dv (int i, int j, int k, int nvar, int n1,
		   int n2, int n3, derivs dv, derivs u, double *values, double time);
void SetMatrix_JFD (int nvar, int n1, int n2, int n3,
		    derivs u, int *ncols, int **cols, double **Matrix, double time);
double PunctEvalAtArbitPosition (double *v, int ivar, double A, double B, double phi,
				 int nvar, int n1, int n2, int n3);
double PunctEvalCoeffsAtArbitPosition (double *v, int ivar, double A, double B, double phi,
                                 int nvar, int n1, int n2, int n3);
void calculate_derivs (int i, int j, int k, int ivar, int nvar, int n1,
		       int n2, int n3, derivs v, derivs vv);
double interpol (double a, double b, double c, derivs v);
double PunctTaylorExpandAtArbitPosition (int ivar, int nvar, int n1,
                                         int n2, int n3, derivs v, double x,
                                         double y, double z);
double PunctIntPolAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, double x,
				   double y, double z);
double PunctIntPolCoeffsAtArbitPosition (int ivar, int nvar, int n1,
				   int n2, int n3, derivs v, double x,
				   double y, double z);
//wraps the gsm switch
gsm_type getGSM();
double ValueAtArbitPosition(int ivar, int nvar, int n1, int n2, int n3, derivs v, derivs v_coeffs, double xx, double yy, double zz);

//allows evaluation of arbitrary (non-deriv) double arrays
double PunctIntPolForArbitVector (int ivar, int nvar, int n1, int n2, int n3, double *v, double x, double y, double z);

//

// Routines in  "CoordTransf.c"
void AB_To_XR (int nvar, double A, double B, double *X,
	       double *R, derivs U);
void C_To_c (int nvar, double X, double R, double *x,
	     double *r, derivs U);
void rx3_To_xyz (int nvar, double x, double r, double phi, double *y,
		 double *z, derivs U);
void AB3_To_xyz(double A, double B, double phi, double *x, double *y, 
                 double *z);

// Routines in  "Equations.c"
void NonLinEquations_Psi (int i, int j, int k, int n1, int n2, int n3,
                      double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
void LinEquations_Psi (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
void makeInitialGuess_Psi(double x, double y, double z, double *u, double time);

// Routines in  "Newton.c"
double
norm_inf (double *F, int ntotal);
void resid (double *res, int ntotal, double *dv, double *rhs,
	    int *ncols, int **cols, double **JFD);
void LineRelax_al (double *dv, int j, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void LineRelax_be (double *dv, int i, int k, int nvar, int n1, int n2, int n3,
		   double *rhs, int *ncols, int **cols, double **JFD);
void TestRelax (int nvar, int n1, int n2, int n3, derivs v,
		double *dv, double time);
int bicgstab (int nvar, int n1, int n2, int n3, derivs v,
	      derivs dv, int output, int itmax, double tol, double *normres, double time);
int Newton (int nvar, int n1, int n2, int n3, derivs v,
	     double tol, int itmax, double time);
void getInitialGuess(int nvar, int n1, int n2, int n3, derivs v, double time);

// function pointers describing the phyical system to Newton.c
extern void (*NonLinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	              double A, double B, double X, double R,
		      double x, double r, double phi,
		      double y, double z, derivs U, double *values, double time);
extern void (*LinEquations) (int i, int j, int k, int n1, int n2, int n3, 
	           double A, double B, double X, double R,
		   double x, double r, double phi,
		   double y, double z, derivs dU, derivs U, double *values, double time);
extern void (*makeInitialGuess)(double x, double y, double z, double *u, double time);
extern void *NewtonParam; // used to pass parameters to (Non)LinEquations and makeInitialGuess


/* 
 27: -1.325691774825335e-03
 37: -1.325691778944117e-03
 47: -1.325691778942711e-03
 
 17: -1.510625972641537e-03
 21: -1.511443006977708e-03
 27: -1.511440785153687e-03
 37: -1.511440809549005e-03
 39: -1.511440809597588e-03
 */
