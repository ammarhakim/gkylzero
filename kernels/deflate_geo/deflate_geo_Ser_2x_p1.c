#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_2x_Ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]; 
  fld_deflated[2] = 0.7071067811865475*fld[2]; 
  fld_deflated[3] = 0.7071067811865475*fld[4]; 
 
}
