#include <gkyl_inflate_surf_kernels.h>

GKYL_CU_DH void inflate_surfx_3x_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[2] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[6] = 1.4142135623730951*fld_deflated[3]; 
}
 
GKYL_CU_DH void inflate_surfy_3x_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[3] = 1.4142135623730951*fld_deflated[2]; 
  fld[5] = 1.4142135623730951*fld_deflated[3]; 
}
 
GKYL_CU_DH void inflate_surfz_3x_ser_p1(const double *fld_deflated, double *fld) 
{ 
  fld[0] = 1.4142135623730951*fld_deflated[0]; 
  fld[1] = 1.4142135623730951*fld_deflated[1]; 
  fld[2] = 1.4142135623730951*fld_deflated[2]; 
  fld[4] = 1.4142135623730951*fld_deflated[3]; 
}
 
