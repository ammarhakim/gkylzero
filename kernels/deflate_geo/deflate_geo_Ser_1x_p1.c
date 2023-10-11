#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_1x_Ser_p1(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.5*fld[0]; 
  fld_deflated[1] = 0.5*fld[1]; 
 
}
