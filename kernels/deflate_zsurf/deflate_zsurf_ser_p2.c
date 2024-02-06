#include <gkyl_deflate_zsurf_kernels.h>

GKYL_CU_DH void deflate_zsurf_lo_2x_ser_p2_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.581138830084189*fld[5]-1.224744871391589*fld[2]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.581138830084189*fld[7]-1.224744871391589*fld[3]+0.7071067811865475*fld[1]; 
  fld_deflated[2] = 0.7071067811865475*fld[4]-1.224744871391589*fld[6]; 
 
}
GKYL_CU_DH void deflate_zsurf_up_2x_ser_p2_remy(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.581138830084189*fld[5]+1.224744871391589*fld[2]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.581138830084189*fld[7]+1.224744871391589*fld[3]+0.7071067811865475*fld[1]; 
  fld_deflated[2] = 1.224744871391589*fld[6]+0.7071067811865475*fld[4]; 
 
}
GKYL_CU_DH void deflate_zsurf_lo_3x_ser_p2_remz(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.581138830084189*fld[9]-1.224744871391589*fld[3]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.581138830084189*fld[15]-1.224744871391589*fld[5]+0.7071067811865475*fld[1]; 
  fld_deflated[2] = 1.581138830084189*fld[16]-1.224744871391589*fld[6]+0.7071067811865475*fld[2]; 
  fld_deflated[3] = 1.581138830084189*fld[19]-1.224744871391589*fld[10]+0.7071067811865475*fld[4]; 
  fld_deflated[4] = 0.7071067811865475*fld[7]-1.224744871391589*fld[13]; 
  fld_deflated[5] = 0.7071067811865475*fld[8]-1.224744871391589*fld[14]; 
  fld_deflated[6] = 0.7071067811865475*fld[11]-1.224744871391589*fld[17]; 
  fld_deflated[7] = 0.7071067811865475*fld[12]-1.224744871391589*fld[18]; 
 
}
GKYL_CU_DH void deflate_zsurf_up_3x_ser_p2_remz(const double *fld, double *fld_deflated) 
{ 
  fld_deflated[0] = 1.581138830084189*fld[9]+1.224744871391589*fld[3]+0.7071067811865475*fld[0]; 
  fld_deflated[1] = 1.581138830084189*fld[15]+1.224744871391589*fld[5]+0.7071067811865475*fld[1]; 
  fld_deflated[2] = 1.581138830084189*fld[16]+1.224744871391589*fld[6]+0.7071067811865475*fld[2]; 
  fld_deflated[3] = 1.581138830084189*fld[19]+1.224744871391589*fld[10]+0.7071067811865475*fld[4]; 
  fld_deflated[4] = 1.224744871391589*fld[13]+0.7071067811865475*fld[7]; 
  fld_deflated[5] = 1.224744871391589*fld[14]+0.7071067811865475*fld[8]; 
  fld_deflated[6] = 1.224744871391589*fld[17]+0.7071067811865475*fld[11]; 
  fld_deflated[7] = 1.224744871391589*fld[18]+0.7071067811865475*fld[12]; 
 
}
