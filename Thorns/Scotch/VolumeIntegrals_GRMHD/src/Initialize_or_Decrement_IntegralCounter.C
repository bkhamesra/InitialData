#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"
#include "util_Table.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>

void VI_GRMHD_InitializeIntegralCounterToZero(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  *IntegralCounter = 0;

  if(verbose==2) printf("VolumeIntegrals_GRMHD: Just set IntegralCounter to %d\n",*IntegralCounter);
}

void VI_GRMHD_InitializeIntegralCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  if(cctk_iteration%VolIntegral_out_every==0) {
    *IntegralCounter = NumIntegrals;
    if(verbose==2) printf("VolumeIntegrals_GRMHD: Just set IntegralCounter to %d == NumIntegrals\n",*IntegralCounter);
  }
}

void VI_GRMHD_DecrementIntegralCounter(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

    (*IntegralCounter)--;
    if(verbose==2) printf("VolumeIntegrals_GRMHD: Just decremented IntegralCounter to %d\n",*IntegralCounter);
}
