#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_div_x_1x_ser_p2(const double *dxv, const double *Al, const double *Ac, const double *Ar, double* GKYL_RESTRICT out) 
{ 
  // dxv[NDIM]: Cell spacing.
  // Al/Ac/Ar:  Inpute vector in left/center/right cells.
  // out:       Increment to volume expansion of div(A) in one direction.

  const double dx1 = 2.0/dxv[0]; 

  const double *A_l = &Al[0]; 
  const double *A_c = &Ac[0]; 
  const double *A_r = &Ar[0]; 

  out[0] += (0.2445699350390395*A_r[2]-0.2445699350390395*A_l[2]-0.3518228202874282*(A_r[1]+A_l[1])+0.7036456405748563*A_c[1]+0.25*A_r[0]-0.25*A_l[0])*dx1; 
  out[1] += (0.4236075534914363*(A_r[2]+A_l[2])+0.8472151069828725*A_c[2]-0.609375*A_r[1]+0.609375*A_l[1]+0.4330127018922193*(A_r[0]+A_l[0])-0.8660254037844386*A_c[0])*dx1; 
  out[2] += (0.546875*A_r[2]-0.546875*A_l[2]-0.7866997421983816*(A_r[1]+A_l[1])-2.299583861810654*A_c[1]+0.5590169943749475*A_r[0]-0.5590169943749475*A_l[0])*dx1; 
} 
