#include <gkyl_euler_kernels.h> 
GKYL_CU_DH void euler_pkpm_pressure_div_x_1x_ser_p3(const double *dxv, const double *p_ijl, const double *p_ijc, const double *p_ijr, double* GKYL_RESTRICT div_p) 
{ 
  // dxv[NDIM]: Cell spacing.
  // p_ijl/p_ijc/p_ijr: Pressure tensor in left/center/right cells.
  // div_p: Volume expansion of div(p).

  const double dx1 = 2.0/dxv[0]; 

  const double *Pxx_l = &p_ijl[0]; 
  const double *Pxy_l = &p_ijl[4]; 
  const double *Pxz_l = &p_ijl[8]; 
  const double *Pyy_l = &p_ijl[12]; 
  const double *Pyz_l = &p_ijl[16]; 
  const double *Pzz_l = &p_ijl[20]; 

  const double *Pxx_c = &p_ijc[0]; 
  const double *Pxy_c = &p_ijc[4]; 
  const double *Pxz_c = &p_ijc[8]; 
  const double *Pyy_c = &p_ijc[12]; 
  const double *Pyz_c = &p_ijc[16]; 
  const double *Pzz_c = &p_ijc[20]; 

  const double *Pxx_r = &p_ijr[0]; 
  const double *Pxy_r = &p_ijr[4]; 
  const double *Pxz_r = &p_ijr[8]; 
  const double *Pyy_r = &p_ijr[12]; 
  const double *Pyz_r = &p_ijr[16]; 
  const double *Pzz_r = &p_ijr[20]; 

  double *div_p_x = &div_p[0]; 
  double *div_p_y = &div_p[4]; 
  double *div_p_z = &div_p[8]; 
  div_p_x[0] = (-0.1889822365046136*Pxx_r[3]*dx1)-0.1889822365046136*Pxx_l[3]*dx1+0.3779644730092272*Pxx_c[3]*dx1+0.3493856214843422*Pxx_r[2]*dx1-0.3493856214843422*Pxx_l[2]*dx1-0.3788861141556919*Pxx_r[1]*dx1-0.3788861141556919*Pxx_l[1]*dx1+0.7577722283113838*Pxx_c[1]*dx1+0.25*Pxx_r[0]*dx1-0.25*Pxx_l[0]*dx1; 
  div_p_x[1] = (-0.3273268353539885*Pxx_r[3]*dx1)+0.3273268353539885*Pxx_l[3]*dx1+0.605153647844909*Pxx_r[2]*dx1+0.605153647844909*Pxx_l[2]*dx1+1.210307295689818*Pxx_c[2]*dx1-0.65625*Pxx_r[1]*dx1+0.65625*Pxx_l[1]*dx1+0.4330127018922193*Pxx_r[0]*dx1+0.4330127018922193*Pxx_l[0]*dx1-0.8660254037844386*Pxx_c[0]*dx1; 
  div_p_x[2] = (-0.4225771273642583*Pxx_r[3]*dx1)-0.4225771273642583*Pxx_l[3]*dx1+0.8451542547285166*Pxx_c[3]*dx1+0.78125*Pxx_r[2]*dx1-0.78125*Pxx_l[2]*dx1-0.8472151069828725*Pxx_r[1]*dx1-0.8472151069828725*Pxx_l[1]*dx1-2.178553132241671*Pxx_c[1]*dx1+0.5590169943749475*Pxx_r[0]*dx1-0.5590169943749475*Pxx_l[0]*dx1; 
  div_p_x[3] = (-0.5*Pxx_r[3]*dx1)+0.5*Pxx_l[3]*dx1+0.9243874661093152*Pxx_r[2]*dx1+0.9243874661093152*Pxx_l[2]*dx1-4.067304850880986*Pxx_c[2]*dx1-1.00243843327159*Pxx_r[1]*dx1+1.00243843327159*Pxx_l[1]*dx1+0.6614378277661477*Pxx_r[0]*dx1+0.6614378277661477*Pxx_l[0]*dx1-1.322875655532295*Pxx_c[0]*dx1; 

  div_p_y[0] = (-0.1889822365046136*Pxy_r[3]*dx1)-0.1889822365046136*Pxy_l[3]*dx1+0.3779644730092272*Pxy_c[3]*dx1+0.3493856214843422*Pxy_r[2]*dx1-0.3493856214843422*Pxy_l[2]*dx1-0.3788861141556919*Pxy_r[1]*dx1-0.3788861141556919*Pxy_l[1]*dx1+0.7577722283113838*Pxy_c[1]*dx1+0.25*Pxy_r[0]*dx1-0.25*Pxy_l[0]*dx1; 
  div_p_y[1] = (-0.3273268353539885*Pxy_r[3]*dx1)+0.3273268353539885*Pxy_l[3]*dx1+0.605153647844909*Pxy_r[2]*dx1+0.605153647844909*Pxy_l[2]*dx1+1.210307295689818*Pxy_c[2]*dx1-0.65625*Pxy_r[1]*dx1+0.65625*Pxy_l[1]*dx1+0.4330127018922193*Pxy_r[0]*dx1+0.4330127018922193*Pxy_l[0]*dx1-0.8660254037844386*Pxy_c[0]*dx1; 
  div_p_y[2] = (-0.4225771273642583*Pxy_r[3]*dx1)-0.4225771273642583*Pxy_l[3]*dx1+0.8451542547285166*Pxy_c[3]*dx1+0.78125*Pxy_r[2]*dx1-0.78125*Pxy_l[2]*dx1-0.8472151069828725*Pxy_r[1]*dx1-0.8472151069828725*Pxy_l[1]*dx1-2.178553132241671*Pxy_c[1]*dx1+0.5590169943749475*Pxy_r[0]*dx1-0.5590169943749475*Pxy_l[0]*dx1; 
  div_p_y[3] = (-0.5*Pxy_r[3]*dx1)+0.5*Pxy_l[3]*dx1+0.9243874661093152*Pxy_r[2]*dx1+0.9243874661093152*Pxy_l[2]*dx1-4.067304850880986*Pxy_c[2]*dx1-1.00243843327159*Pxy_r[1]*dx1+1.00243843327159*Pxy_l[1]*dx1+0.6614378277661477*Pxy_r[0]*dx1+0.6614378277661477*Pxy_l[0]*dx1-1.322875655532295*Pxy_c[0]*dx1; 

  div_p_z[0] = (-0.1889822365046136*Pxz_r[3]*dx1)-0.1889822365046136*Pxz_l[3]*dx1+0.3779644730092272*Pxz_c[3]*dx1+0.3493856214843422*Pxz_r[2]*dx1-0.3493856214843422*Pxz_l[2]*dx1-0.3788861141556919*Pxz_r[1]*dx1-0.3788861141556919*Pxz_l[1]*dx1+0.7577722283113838*Pxz_c[1]*dx1+0.25*Pxz_r[0]*dx1-0.25*Pxz_l[0]*dx1; 
  div_p_z[1] = (-0.3273268353539885*Pxz_r[3]*dx1)+0.3273268353539885*Pxz_l[3]*dx1+0.605153647844909*Pxz_r[2]*dx1+0.605153647844909*Pxz_l[2]*dx1+1.210307295689818*Pxz_c[2]*dx1-0.65625*Pxz_r[1]*dx1+0.65625*Pxz_l[1]*dx1+0.4330127018922193*Pxz_r[0]*dx1+0.4330127018922193*Pxz_l[0]*dx1-0.8660254037844386*Pxz_c[0]*dx1; 
  div_p_z[2] = (-0.4225771273642583*Pxz_r[3]*dx1)-0.4225771273642583*Pxz_l[3]*dx1+0.8451542547285166*Pxz_c[3]*dx1+0.78125*Pxz_r[2]*dx1-0.78125*Pxz_l[2]*dx1-0.8472151069828725*Pxz_r[1]*dx1-0.8472151069828725*Pxz_l[1]*dx1-2.178553132241671*Pxz_c[1]*dx1+0.5590169943749475*Pxz_r[0]*dx1-0.5590169943749475*Pxz_l[0]*dx1; 
  div_p_z[3] = (-0.5*Pxz_r[3]*dx1)+0.5*Pxz_l[3]*dx1+0.9243874661093152*Pxz_r[2]*dx1+0.9243874661093152*Pxz_l[2]*dx1-4.067304850880986*Pxz_c[2]*dx1-1.00243843327159*Pxz_r[1]*dx1+1.00243843327159*Pxz_l[1]*dx1+0.6614378277661477*Pxz_r[0]*dx1+0.6614378277661477*Pxz_l[0]*dx1-1.322875655532295*Pxz_c[0]*dx1; 

} 
