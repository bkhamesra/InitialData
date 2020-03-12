//
// startup routines allocating eg. memory
//

#include "primsolver.h"
#include "string.h"
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

CCTK_REAL primsolver_poison;

void PrimitiveSolver_Startup(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;

    // generate poison value
    memset(&primsolver_poison, poison_value, sizeof(primsolver_poison));
}
