 /*@@
   @file      PrimitiveSolver_EOS.c
   @date      Wed March  8 2017
   @author    
   @desc 
   Sets the EOS handle number from CactusEOS for use by
any hydro evolution code. An improvement is needed in order to be compatible with the EOS_Omni Thorn, this could be a good student project.
   @enddesc 
 @@*/

#include "cctk.h"
#include "cctk_Parameters.h"
#include "cctk_Arguments.h"

#include <stdio.h>
#include <stdlib.h>

#include "EOS_Base.h"

 /*@@
   @routine    PrimitiveSolver_EOSHandle
   @date       Wed Feb  6 18:28:01 2002
   @author     Ian Hawke
   @desc 
   Sets the EOS handle number
   @enddesc 
   @calls     
   @calledby   
   @history 
 
   @endhistory 

@@*/

void PrimitiveSolver_EOSHandle(CCTK_ARGUMENTS)
{

  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;

  int ierr = 0;
  char *infoline;
  
  *PrimitiveSolver_eos_handle = EOS_Handle(PrimitiveSolver_eos_table);
  
  if (*PrimitiveSolver_eos_handle < 0) {
    ierr = *PrimitiveSolver_eos_handle;
    CCTK_WARN(0, "Table handle is incorrect");
  }

  *PrimitiveSolver_polytrope_handle = EOS_Handle("2D_Polytrope");
  
  if (*PrimitiveSolver_polytrope_handle < 0) {
    ierr = *PrimitiveSolver_polytrope_handle;
    CCTK_WARN(0, "PrimitiveSolver requires the 2D_Polytrope thorn to be active for the atmosphere");
  }

  infoline = (char *)malloc(100*sizeof(char));
  if (!infoline) {
    CCTK_WARN(0, "Failed to malloc an array of 100 chars in PrimitiveSolver_EOS.");
  }
  else {
    sprintf(infoline,"PrimitiveSolver will use the %s equation of state.",
            PrimitiveSolver_eos_table);
    CCTK_INFO(infoline);
  }
  free(infoline);
  infoline = NULL;

}
