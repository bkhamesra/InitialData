#include "cctk.h"
#include "cctk_Arguments.h"

#ifdef __cplusplus
extern "C"
{
#endif
	void bowenid_getSourceInfo(double xx, double yy, double zz, double *rhohat, double jhat[3], double *psi, double *psiguess, double Aij[3][3]);
	double bowenid_getSourceMass(int index);
	double bowenid_getRestMass(int index);
	void bowenid_constructSources(CCTK_ARGUMENTS);
	CCTK_REAL bowenid_getExternalPsiGuess(CCTK_INT idx, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bowenid_getLorentzFactor(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bowenid_getConfDensity(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bowenid_getConfPress(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
	CCTK_REAL bowenid_getConfPolyConstant(CCTK_REAL x, CCTK_REAL y, CCTK_REAL z);
#ifdef __cplusplus
}
#endif
