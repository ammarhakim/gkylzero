#include <gkyl_deflate_zsurf_kernels.h>

GKYL_CU_DH void deflate_zsurf_lo_2x_ser_p1_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-1.224744871391589*fld[2]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]-1.224744871391589*fld[3]; 
 
}
GKYL_CU_DH void deflate_zsurf_up_2x_ser_p1_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.224744871391589*fld[2]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.224744871391589*fld[3]+0.7071067811865475*fld[1]; 
 
}
