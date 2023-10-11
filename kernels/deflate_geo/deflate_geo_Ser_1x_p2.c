#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_1x_Ser_p2(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = (-0.5590169943749475*fld[9])-0.5590169943749475*fld[8]+0.5*fld[0]; 
  fld_deflated[1] = (-0.5590169943749476*fld[15])-0.5590169943749476*fld[12]+0.5*fld[1]; 
  fld_deflated[2] = 0.5*fld[7]; 
 
}
