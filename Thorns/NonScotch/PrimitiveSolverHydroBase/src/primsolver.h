// MatterLibs: matterlibs.h

#include "cctk.h"
#include "cctk_Arguments.h"

#ifdef __cplusplus
extern "C" {
#endif

// this one is actually kind of private (in E_to_rho.c)
int prims_E_to_rho(int eos_handle, double dE, double dJx, double dJy,
    double dJz, double psi, double *drho, double *dlorentz);
CCTK_REAL PrimitiveSolver_getRho0W(CCTK_INT polyhandle, CCTK_REAL x, CCTK_REAL y, CCTK_REAL z, CCTK_REAL psiL);

extern CCTK_REAL primsolver_poison;

#ifdef __cplusplus
}
#endif
