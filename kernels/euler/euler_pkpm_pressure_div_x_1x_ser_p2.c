#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_div_x_1x_ser_p2(const double *dxv, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* GKYL_RESTRICT div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // div_p: Volume expansion of div(p).

  const double dx1 = 2.0/dxv[0]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[3]; 
  const double *Pxz_l = &p_ijl[6]; 
  const double *Pyy_l = &p_ijl[9]; 
  const double *Pyz_l = &p_ijl[12]; 
  const double *Pzz_l = &p_ijl[15]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[3]; 
  const double *Pxz_c = &p_ijc[6]; 
  const double *Pyy_c = &p_ijc[9]; 
  const double *Pyz_c = &p_ijc[12]; 
  const double *Pzz_c = &p_ijc[15]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[3]; 
  const double *Pxz_r = &p_ijr[6]; 
  const double *Pyy_r = &p_ijr[9]; 
  const double *Pyz_r = &p_ijr[12]; 
  const double *Pzz_r = &p_ijr[15]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[3]; 
  double *div_p_z = &div_p[6]; 
  div_p_x[0] = 0.2445699350390395*Pxx_r[2]*dx1-0.2445699350390395*Pxx_l[2]*dx1-0.3518228202874282*Pxx_r[1]*dx1-0.3518228202874282*Pxx_l[1]*dx1+0.7036456405748563*Pxx_c[1]*dx1+0.25*Pxx_r[0]*dx1-0.25*Pxx_l[0]*dx1; 
  div_p_x[1] = 0.4236075534914363*Pxx_r[2]*dx1+0.4236075534914363*Pxx_l[2]*dx1+0.8472151069828725*Pxx_c[2]*dx1-0.609375*Pxx_r[1]*dx1+0.609375*Pxx_l[1]*dx1+0.4330127018922193*Pxx_r[0]*dx1+0.4330127018922193*Pxx_l[0]*dx1-0.8660254037844386*Pxx_c[0]*dx1; 
  div_p_x[2] = 0.546875*Pxx_r[2]*dx1-0.546875*Pxx_l[2]*dx1-0.7866997421983816*Pxx_r[1]*dx1-0.7866997421983816*Pxx_l[1]*dx1-2.299583861810654*Pxx_c[1]*dx1+0.5590169943749475*Pxx_r[0]*dx1-0.5590169943749475*Pxx_l[0]*dx1; 

  div_p_y[0] = 0.2445699350390395*Pxy_r[2]*dx1-0.2445699350390395*Pxy_l[2]*dx1-0.3518228202874282*Pxy_r[1]*dx1-0.3518228202874282*Pxy_l[1]*dx1+0.7036456405748563*Pxy_c[1]*dx1+0.25*Pxy_r[0]*dx1-0.25*Pxy_l[0]*dx1; 
  div_p_y[1] = 0.4236075534914363*Pxy_r[2]*dx1+0.4236075534914363*Pxy_l[2]*dx1+0.8472151069828725*Pxy_c[2]*dx1-0.609375*Pxy_r[1]*dx1+0.609375*Pxy_l[1]*dx1+0.4330127018922193*Pxy_r[0]*dx1+0.4330127018922193*Pxy_l[0]*dx1-0.8660254037844386*Pxy_c[0]*dx1; 
  div_p_y[2] = 0.546875*Pxy_r[2]*dx1-0.546875*Pxy_l[2]*dx1-0.7866997421983816*Pxy_r[1]*dx1-0.7866997421983816*Pxy_l[1]*dx1-2.299583861810654*Pxy_c[1]*dx1+0.5590169943749475*Pxy_r[0]*dx1-0.5590169943749475*Pxy_l[0]*dx1; 

  div_p_z[0] = 0.2445699350390395*Pxz_r[2]*dx1-0.2445699350390395*Pxz_l[2]*dx1-0.3518228202874282*Pxz_r[1]*dx1-0.3518228202874282*Pxz_l[1]*dx1+0.7036456405748563*Pxz_c[1]*dx1+0.25*Pxz_r[0]*dx1-0.25*Pxz_l[0]*dx1; 
  div_p_z[1] = 0.4236075534914363*Pxz_r[2]*dx1+0.4236075534914363*Pxz_l[2]*dx1+0.8472151069828725*Pxz_c[2]*dx1-0.609375*Pxz_r[1]*dx1+0.609375*Pxz_l[1]*dx1+0.4330127018922193*Pxz_r[0]*dx1+0.4330127018922193*Pxz_l[0]*dx1-0.8660254037844386*Pxz_c[0]*dx1; 
  div_p_z[2] = 0.546875*Pxz_r[2]*dx1-0.546875*Pxz_l[2]*dx1-0.7866997421983816*Pxz_r[1]*dx1-0.7866997421983816*Pxz_l[1]*dx1-2.299583861810654*Pxz_c[1]*dx1+0.5590169943749475*Pxz_r[0]*dx1-0.5590169943749475*Pxz_l[0]*dx1; 

} 
