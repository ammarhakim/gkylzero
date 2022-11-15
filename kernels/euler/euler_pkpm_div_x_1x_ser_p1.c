#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_div_x_1x_ser_p1(const double *dxv, const double *Al, const double *Ac, const double *Ar, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *A_l = &Al[0]; 
  const double *A_c = &Ac[0]; 
  const double *A_r = &Ar[0]; 

  out[0] += ((-0.2886751345948129*(A_r[1]+A_l[1]))+0.5773502691896258*A_c[1]+0.25*A_r[0]-0.25*A_l[0])*dx1; 
  out[1] += ((-0.5*A_r[1])+0.5*A_l[1]+0.4330127018922193*(A_r[0]+A_l[0])-0.8660254037844386*A_c[0])*dx1; 
} 
