#include "gkyl_deflate_geo_kernels.h"

GKYL_CU_DH void deflate_geo_2x_Ser_p2(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-0.7905694150420947*fld[9]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]-0.7905694150420948*fld[15]; 
  fld_deflated[2] = 0.7071067811865475*fld[2]-0.7905694150420948*fld[16]; 
  fld_deflated[3] = 0.7071067811865475*fld[4]-0.7905694150420947*fld[19]; 
  fld_deflated[4] = 0.7071067811865475*fld[7]; 
  fld_deflated[5] = 0.7071067811865475*fld[8]; 
  fld_deflated[6] = 0.7071067811865475*fld[11]; 
  fld_deflated[7] = 0.7071067811865475*fld[12]; 
 
}
