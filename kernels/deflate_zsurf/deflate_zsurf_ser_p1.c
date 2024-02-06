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
GKYL_CU_DH void deflate_zsurf_lo_3x_ser_p1_remz(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 0.7071067811865475*fld[0]-1.224744871391589*fld[3]; 
  fld_deflated[1] = 0.7071067811865475*fld[1]-1.224744871391589*fld[5]; 
  fld_deflated[2] = 0.7071067811865475*fld[2]-1.224744871391589*fld[6]; 
  fld_deflated[3] = 0.7071067811865475*fld[4]-1.224744871391589*fld[7]; 
 
}
GKYL_CU_DH void deflate_zsurf_up_3x_ser_p1_remz(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.224744871391589*fld[3]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.224744871391589*fld[5]+0.7071067811865475*fld[1]; 
  fld_deflated[2] = 1.224744871391589*fld[6]+0.7071067811865475*fld[2]; 
  fld_deflated[3] = 1.224744871391589*fld[7]+0.7071067811865475*fld[4]; 
 
}
