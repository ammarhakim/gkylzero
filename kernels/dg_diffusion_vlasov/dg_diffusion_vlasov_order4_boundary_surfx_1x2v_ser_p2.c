#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[7]-9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]+8.11898816047911*coeff[0]*fSkin[4]+8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[13]-6.708203932499369*coeff[0]*fEdge[13]+8.11898816047911*coeff[0]*fSkin[5]+8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 14.16059535957087*coeff[0]*fSkin[11]-9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 14.16059535957087*coeff[0]*fSkin[13]-9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]+8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]+8.11898816047911*coeff[0]*fSkin[10]+8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 8.118988160479114*coeff[0]*fSkin[12]+8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 8.118988160479114*coeff[0]*fSkin[15]+8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 14.16059535957086*coeff[0]*fSkin[17]-9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]+8.11898816047911*coeff[0]*fSkin[6]-8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]+15.61296411439865*coeff[0]*fSkin[4]+4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]+8.118988160479114*coeff[0]*fSkin[8]-8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]+15.61296411439865*coeff[0]*fSkin[5]+4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*fSkin[18]+8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]+8.118988160479114*coeff[0]*fSkin[9]-8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 8.118988160479114*coeff[0]*fSkin[19]+8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]+15.61296411439865*coeff[0]*fSkin[10]+4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]+8.118988160479114*coeff[0]*fSkin[14]-8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]+8.118988160479114*coeff[0]*fSkin[16]-8.118988160479114*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[7]; 
  boundSurf_incr[4] = 2.8125*coeff[0]*fSkin[4]-5.083290641897235*coeff[0]*fSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*fSkin[5]-5.083290641897235*coeff[0]*fSkin[13]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]-10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 2.8125*coeff[0]*fSkin[10]-5.083290641897234*coeff[0]*fSkin[17]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]-10.89276566120836*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]-10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]-10.89276566120836*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[7]-6.708203932499369*coeff[0]*fEdge[7]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[7])+9.077304717673634*coeff[0]*fEdge[7]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]-8.11898816047911*coeff[0]*fSkin[4]-8.11898816047911*coeff[0]*fEdge[4]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[13]-6.708203932499369*coeff[0]*fEdge[13]-8.11898816047911*coeff[0]*fSkin[5]-8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-14.16059535957087*coeff[0]*fSkin[11])+9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[4]+12.65625*coeff[0]*fEdge[4]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-14.16059535957087*coeff[0]*fSkin[13])+9.077304717673634*coeff[0]*fEdge[13]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]-8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 6.708203932499369*coeff[0]*fSkin[17]-6.708203932499369*coeff[0]*fEdge[17]-8.11898816047911*coeff[0]*fSkin[10]-8.11898816047911*coeff[0]*fEdge[10]+4.6875*coeff[0]*fSkin[6]-4.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 20.34375*coeff[0]*fSkin[7]-0.65625*coeff[0]*fEdge[7]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-8.118988160479114*coeff[0]*fSkin[12])-8.118988160479114*coeff[0]*fEdge[12]+4.6875*coeff[0]*fSkin[8]-4.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-8.118988160479114*coeff[0]*fSkin[15])-8.118988160479114*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-14.16059535957086*coeff[0]*fSkin[17])+9.077304717673634*coeff[0]*fEdge[17]+15.46875*coeff[0]*fSkin[10]+12.65625*coeff[0]*fEdge[10]-8.11898816047911*coeff[0]*fSkin[6]+8.11898816047911*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]-15.61296411439865*coeff[0]*fSkin[4]-4.72019845319029*coeff[0]*fEdge[4]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 15.46875*coeff[0]*fSkin[12]+12.65625*coeff[0]*fEdge[12]-8.118988160479114*coeff[0]*fSkin[8]+8.118988160479114*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 20.34375*coeff[0]*fSkin[13]-0.65625*coeff[0]*fEdge[13]-15.61296411439865*coeff[0]*fSkin[5]-4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-8.118988160479114*coeff[0]*fSkin[18])-8.118988160479114*coeff[0]*fEdge[18]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]-8.118988160479114*coeff[0]*fSkin[9]+8.118988160479114*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-8.118988160479114*coeff[0]*fSkin[19])-8.118988160479114*coeff[0]*fEdge[19]+4.6875*coeff[0]*fSkin[16]-4.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 20.34375*coeff[0]*fSkin[17]-0.65625*coeff[0]*fEdge[17]-15.61296411439865*coeff[0]*fSkin[10]-4.720198453190292*coeff[0]*fEdge[10]+4.192627457812107*coeff[0]*fSkin[6]-4.192627457812107*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 15.46875*coeff[0]*fSkin[18]+12.65625*coeff[0]*fEdge[18]-8.118988160479114*coeff[0]*fSkin[14]+8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 15.46875*coeff[0]*fSkin[19]+12.65625*coeff[0]*fEdge[19]-8.118988160479114*coeff[0]*fSkin[16]+8.118988160479114*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[7]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 5.083290641897235*coeff[0]*fSkin[11]+2.8125*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*fSkin[13]+2.8125*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 19.6875*coeff[0]*fSkin[7]+10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = 5.083290641897234*coeff[0]*fSkin[17]+2.8125*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]+10.89276566120836*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 2.8125*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.6875*coeff[0]*fSkin[13]+10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.6875*coeff[0]*fSkin[17]+10.89276566120836*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 2.8125*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 2.8125*coeff[0]*fSkin[19]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += -1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += -1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += -1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += -1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += -1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += -1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += -1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += -1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += -1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += -1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += -1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += -1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[20] = {0.0}; 
  vol_incr[7] = -31.81980515339464*fSkin[0]*coeff[2]; 
  vol_incr[11] = -31.81980515339463*coeff[2]*fSkin[2]; 
  vol_incr[13] = -31.81980515339463*coeff[2]*fSkin[3]; 
  vol_incr[17] = -31.81980515339464*coeff[2]*fSkin[6]; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[7])-1.797214641813825*coeff[1]*fSkin[7]+4.743416490252569*coeff[0]*fSkin[7]+5.303300858899105*coeff[2]*fEdge[7]-1.797214641813825*coeff[1]*fEdge[7]-4.743416490252569*coeff[0]*fEdge[7]-6.418623720763661*fSkin[1]*coeff[2]-6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]+5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]+5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-11.19493359006374*coeff[2]*fSkin[7])-3.112867071728247*coeff[1]*fSkin[7]+10.01305300439132*coeff[0]*fSkin[7]+7.176239480810091*coeff[2]*fEdge[7]-3.112867071728247*coeff[1]*fEdge[7]-6.418623720763665*coeff[0]*fEdge[7]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]-6.418623720763661*fSkin[0]*coeff[2]+6.418623720763661*fEdge[0]*coeff[2]-1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]+1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]+5.74099158464807*coeff[0]*fSkin[0]-5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-5.303300858899106*coeff[2]*fSkin[11])-1.797214641813825*coeff[1]*fSkin[11]+4.743416490252569*coeff[0]*fSkin[11]+5.303300858899106*coeff[2]*fEdge[11]-1.797214641813825*coeff[1]*fEdge[11]-4.743416490252569*coeff[0]*fEdge[11]-6.418623720763661*coeff[2]*fSkin[4]-0.9943689110435817*coeff[1]*fSkin[4]+5.74099158464807*coeff[0]*fSkin[4]-6.418623720763661*coeff[2]*fEdge[4]+0.9943689110435817*coeff[1]*fEdge[4]+5.74099158464807*coeff[0]*fEdge[4]-3.705794133009818*coeff[2]*fSkin[2]+3.31456303681194*coeff[0]*fSkin[2]+3.705794133009818*coeff[2]*fEdge[2]-3.31456303681194*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-5.303300858899106*coeff[2]*fSkin[13])-1.797214641813825*coeff[1]*fSkin[13]+4.743416490252569*coeff[0]*fSkin[13]+5.303300858899106*coeff[2]*fEdge[13]-1.797214641813825*coeff[1]*fEdge[13]-4.743416490252569*coeff[0]*fEdge[13]-6.418623720763661*coeff[2]*fSkin[5]-0.9943689110435817*coeff[1]*fSkin[5]+5.74099158464807*coeff[0]*fSkin[5]-6.418623720763661*coeff[2]*fEdge[5]+0.9943689110435817*coeff[1]*fEdge[5]+5.74099158464807*coeff[0]*fEdge[5]-3.705794133009818*coeff[2]*fSkin[3]+3.31456303681194*coeff[0]*fSkin[3]+3.705794133009818*coeff[2]*fEdge[3]-3.31456303681194*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-11.19493359006374*coeff[2]*fSkin[11])-3.112867071728246*coeff[1]*fSkin[11]+10.01305300439132*coeff[0]*fSkin[11]+7.176239480810093*coeff[2]*fEdge[11]-3.112867071728246*coeff[1]*fEdge[11]-6.418623720763666*coeff[0]*fEdge[11]-12.2291206389324*coeff[2]*fSkin[4]-1.72229747539442*coeff[1]*fSkin[4]+10.9380580214794*coeff[0]*fSkin[4]-10.00564415912651*coeff[2]*fEdge[4]+1.72229747539442*coeff[1]*fEdge[4]+8.949320199392238*coeff[0]*fEdge[4]-6.418623720763661*coeff[2]*fSkin[2]+5.74099158464807*coeff[0]*fSkin[2]+6.418623720763661*coeff[2]*fEdge[2]-5.74099158464807*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-11.19493359006374*coeff[2]*fSkin[13])-3.112867071728246*coeff[1]*fSkin[13]+10.01305300439132*coeff[0]*fSkin[13]+7.176239480810093*coeff[2]*fEdge[13]-3.112867071728246*coeff[1]*fEdge[13]-6.418623720763666*coeff[0]*fEdge[13]-12.2291206389324*coeff[2]*fSkin[5]-1.72229747539442*coeff[1]*fSkin[5]+10.9380580214794*coeff[0]*fSkin[5]-10.00564415912651*coeff[2]*fEdge[5]+1.72229747539442*coeff[1]*fEdge[5]+8.949320199392238*coeff[0]*fEdge[5]-6.418623720763661*coeff[2]*fSkin[3]+5.74099158464807*coeff[0]*fSkin[3]+6.418623720763661*coeff[2]*fEdge[3]-5.74099158464807*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-5.303300858899105*coeff[2]*fSkin[17])-1.797214641813825*coeff[1]*fSkin[17]+4.743416490252569*coeff[0]*fSkin[17]+5.303300858899105*coeff[2]*fEdge[17]-1.797214641813825*coeff[1]*fEdge[17]-4.743416490252569*coeff[0]*fEdge[17]-6.418623720763661*coeff[2]*fSkin[10]-0.9943689110435817*coeff[1]*fSkin[10]+5.74099158464807*coeff[0]*fSkin[10]-6.418623720763661*coeff[2]*fEdge[10]+0.9943689110435817*coeff[1]*fEdge[10]+5.74099158464807*coeff[0]*fEdge[10]-3.705794133009818*coeff[2]*fSkin[6]+3.31456303681194*coeff[0]*fSkin[6]+3.705794133009818*coeff[2]*fEdge[6]-3.31456303681194*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-34.53800131965151*coeff[2]*fSkin[7])-11.53939308514262*coeff[1]*fSkin[7]+14.38520357976382*coeff[0]*fSkin[7]+3.409330602369041*coeff[2]*fEdge[7]-0.516689242618324*coeff[1]*fEdge[7]-0.4640388251536773*coeff[0]*fEdge[7]-42.4833377263957*fSkin[1]*coeff[2]-11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]+11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]+3.337684334797107*coeff[0]*fEdge[1]-9.756308055560766*fSkin[0]*coeff[1]+5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-6.418623720763661*coeff[2]*fSkin[12])-0.9943689110435818*coeff[1]*fSkin[12]+5.740991584648071*coeff[0]*fSkin[12]-6.418623720763661*coeff[2]*fEdge[12]+0.9943689110435818*coeff[1]*fEdge[12]+5.740991584648071*coeff[0]*fEdge[12]-3.705794133009818*coeff[2]*fSkin[8]+3.31456303681194*coeff[0]*fSkin[8]+3.705794133009818*coeff[2]*fEdge[8]-3.31456303681194*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-6.418623720763661*coeff[2]*fSkin[15])-0.9943689110435818*coeff[1]*fSkin[15]+5.740991584648071*coeff[0]*fSkin[15]-6.418623720763661*coeff[2]*fEdge[15]+0.9943689110435818*coeff[1]*fEdge[15]+5.740991584648071*coeff[0]*fEdge[15]-3.705794133009818*coeff[2]*fSkin[9]+3.31456303681194*coeff[0]*fSkin[9]+3.705794133009818*coeff[2]*fEdge[9]-3.31456303681194*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-11.19493359006374*coeff[2]*fSkin[17])-3.112867071728247*coeff[1]*fSkin[17]+10.01305300439132*coeff[0]*fSkin[17]+7.176239480810091*coeff[2]*fEdge[17]-3.112867071728247*coeff[1]*fEdge[17]-6.418623720763665*coeff[0]*fEdge[17]-12.2291206389324*coeff[2]*fSkin[10]-1.72229747539442*coeff[1]*fSkin[10]+10.9380580214794*coeff[0]*fSkin[10]-10.00564415912651*coeff[2]*fEdge[10]+1.72229747539442*coeff[1]*fEdge[10]+8.949320199392238*coeff[0]*fEdge[10]-6.418623720763661*coeff[2]*fSkin[6]+5.74099158464807*coeff[0]*fSkin[6]+6.418623720763661*coeff[2]*fEdge[6]-5.74099158464807*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-34.53800131965151*coeff[2]*fSkin[11])-11.53939308514262*coeff[1]*fSkin[11]+14.38520357976382*coeff[0]*fSkin[11]+3.409330602369041*coeff[2]*fEdge[11]-0.516689242618324*coeff[1]*fEdge[11]-0.4640388251536773*coeff[0]*fEdge[11]-42.48333772639572*coeff[2]*fSkin[4]-14.89729241469946*coeff[1]*fSkin[4]+11.0400327997135*coeff[0]*fSkin[4]-11.48198316929615*coeff[2]*fEdge[4]-4.669300607592371*coeff[1]*fEdge[4]+3.337684334797105*coeff[0]*fEdge[4]-26.18504799081433*coeff[2]*fSkin[2]-9.756308055560769*coeff[1]*fSkin[2]+2.964635306407854*coeff[0]*fSkin[2]+10.27514541411701*coeff[2]*fEdge[2]+5.648388874272023*coeff[1]*fEdge[2]-2.964635306407854*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-12.2291206389324*coeff[2]*fSkin[12])-1.72229747539442*coeff[1]*fSkin[12]+10.9380580214794*coeff[0]*fSkin[12]-10.00564415912651*coeff[2]*fEdge[12]+1.72229747539442*coeff[1]*fEdge[12]+8.949320199392238*coeff[0]*fEdge[12]-6.418623720763661*coeff[2]*fSkin[8]+5.740991584648071*coeff[0]*fSkin[8]+6.418623720763661*coeff[2]*fEdge[8]-5.740991584648071*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-34.53800131965151*coeff[2]*fSkin[13])-11.53939308514262*coeff[1]*fSkin[13]+14.38520357976382*coeff[0]*fSkin[13]+3.409330602369041*coeff[2]*fEdge[13]-0.516689242618324*coeff[1]*fEdge[13]-0.4640388251536773*coeff[0]*fEdge[13]-42.48333772639572*coeff[2]*fSkin[5]-14.89729241469946*coeff[1]*fSkin[5]+11.0400327997135*coeff[0]*fSkin[5]-11.48198316929615*coeff[2]*fEdge[5]-4.669300607592371*coeff[1]*fEdge[5]+3.337684334797105*coeff[0]*fEdge[5]-26.18504799081433*coeff[2]*fSkin[3]-9.756308055560769*coeff[1]*fSkin[3]+2.964635306407854*coeff[0]*fSkin[3]+10.27514541411701*coeff[2]*fEdge[3]+5.648388874272023*coeff[1]*fEdge[3]-2.964635306407854*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-6.418623720763661*coeff[2]*fSkin[18])-0.9943689110435818*coeff[1]*fSkin[18]+5.740991584648071*coeff[0]*fSkin[18]-6.418623720763661*coeff[2]*fEdge[18]+0.9943689110435818*coeff[1]*fEdge[18]+5.740991584648071*coeff[0]*fEdge[18]-3.705794133009818*coeff[2]*fSkin[14]+3.31456303681194*coeff[0]*fSkin[14]+3.705794133009818*coeff[2]*fEdge[14]-3.31456303681194*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-12.2291206389324*coeff[2]*fSkin[15])-1.72229747539442*coeff[1]*fSkin[15]+10.9380580214794*coeff[0]*fSkin[15]-10.00564415912651*coeff[2]*fEdge[15]+1.72229747539442*coeff[1]*fEdge[15]+8.949320199392238*coeff[0]*fEdge[15]-6.418623720763661*coeff[2]*fSkin[9]+5.740991584648071*coeff[0]*fSkin[9]+6.418623720763661*coeff[2]*fEdge[9]-5.740991584648071*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-6.418623720763661*coeff[2]*fSkin[19])-0.9943689110435818*coeff[1]*fSkin[19]+5.740991584648071*coeff[0]*fSkin[19]-6.418623720763661*coeff[2]*fEdge[19]+0.9943689110435818*coeff[1]*fEdge[19]+5.740991584648071*coeff[0]*fEdge[19]-3.705794133009818*coeff[2]*fSkin[16]+3.31456303681194*coeff[0]*fSkin[16]+3.705794133009818*coeff[2]*fEdge[16]-3.31456303681194*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-34.53800131965151*coeff[2]*fSkin[17])-11.53939308514262*coeff[1]*fSkin[17]+14.38520357976382*coeff[0]*fSkin[17]+3.409330602369041*coeff[2]*fEdge[17]-0.516689242618324*coeff[1]*fEdge[17]-0.4640388251536773*coeff[0]*fEdge[17]-42.4833377263957*coeff[2]*fSkin[10]-14.89729241469947*coeff[1]*fSkin[10]+11.0400327997135*coeff[0]*fSkin[10]-11.48198316929614*coeff[2]*fEdge[10]-4.669300607592371*coeff[1]*fEdge[10]+3.337684334797107*coeff[0]*fEdge[10]-26.18504799081432*coeff[2]*fSkin[6]-9.756308055560766*coeff[1]*fSkin[6]+2.964635306407856*coeff[0]*fSkin[6]+10.27514541411701*coeff[2]*fEdge[6]+5.648388874272021*coeff[1]*fEdge[6]-2.964635306407856*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-12.2291206389324*coeff[2]*fSkin[18])-1.72229747539442*coeff[1]*fSkin[18]+10.9380580214794*coeff[0]*fSkin[18]-10.00564415912651*coeff[2]*fEdge[18]+1.72229747539442*coeff[1]*fEdge[18]+8.949320199392238*coeff[0]*fEdge[18]-6.418623720763661*coeff[2]*fSkin[14]+5.740991584648071*coeff[0]*fSkin[14]+6.418623720763661*coeff[2]*fEdge[14]-5.740991584648071*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-12.2291206389324*coeff[2]*fSkin[19])-1.72229747539442*coeff[1]*fSkin[19]+10.9380580214794*coeff[0]*fSkin[19]-10.00564415912651*coeff[2]*fEdge[19]+1.72229747539442*coeff[1]*fEdge[19]+8.949320199392238*coeff[0]*fEdge[19]-6.418623720763661*coeff[2]*fSkin[16]+5.740991584648071*coeff[0]*fSkin[16]+6.418623720763661*coeff[2]*fEdge[16]-5.740991584648071*coeff[0]*fEdge[16]; 

  boundSurf_incr[0] = 3.59442928362765*coeff[1]*fSkin[7]-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 4.018694109253648*coeff[2]*fSkin[7]-6.225734143456494*coeff[1]*fSkin[7]-3.59442928362765*coeff[0]*fSkin[7]-2.223476479805891*fSkin[1]*coeff[2]+3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 3.594429283627651*coeff[1]*fSkin[11]-1.988737822087164*coeff[1]*fSkin[4]; 
  boundSurf_incr[3] = 3.594429283627651*coeff[1]*fSkin[13]-1.988737822087164*coeff[1]*fSkin[5]; 
  boundSurf_incr[4] = 4.018694109253649*coeff[2]*fSkin[11]-6.225734143456493*coeff[1]*fSkin[11]-3.594429283627651*coeff[0]*fSkin[11]-2.223476479805891*coeff[2]*fSkin[4]+3.444594950788841*coeff[1]*fSkin[4]+1.988737822087164*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 4.018694109253649*coeff[2]*fSkin[13]-6.225734143456493*coeff[1]*fSkin[13]-3.594429283627651*coeff[0]*fSkin[13]-2.223476479805891*coeff[2]*fSkin[5]+3.444594950788841*coeff[1]*fSkin[5]+1.988737822087164*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 3.59442928362765*coeff[1]*fSkin[17]-1.988737822087164*coeff[1]*fSkin[10]; 
  boundSurf_incr[7] = (-31.12867071728247*coeff[2]*fSkin[7])+12.05608232776095*coeff[1]*fSkin[7]+13.92116475461015*coeff[0]*fSkin[7]+31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]-7.702348464916393*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]; 
  boundSurf_incr[8] = -1.988737822087164*coeff[1]*fSkin[12]; 
  boundSurf_incr[9] = -1.988737822087164*coeff[1]*fSkin[15]; 
  boundSurf_incr[10] = 4.018694109253648*coeff[2]*fSkin[17]-6.225734143456494*coeff[1]*fSkin[17]-3.59442928362765*coeff[0]*fSkin[17]-2.223476479805891*coeff[2]*fSkin[10]+3.444594950788841*coeff[1]*fSkin[10]+1.988737822087164*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = (-31.12867071728247*coeff[2]*fSkin[11])+12.05608232776095*coeff[1]*fSkin[11]+13.92116475461015*coeff[0]*fSkin[11]+31.00135455709958*coeff[2]*fSkin[4]-10.22799180710709*coeff[1]*fSkin[4]-7.702348464916396*coeff[0]*fSkin[4]-15.90990257669732*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]; 
  boundSurf_incr[12] = (-2.223476479805891*coeff[2]*fSkin[12])+3.444594950788841*coeff[1]*fSkin[12]+1.988737822087164*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = (-31.12867071728247*coeff[2]*fSkin[13])+12.05608232776095*coeff[1]*fSkin[13]+13.92116475461015*coeff[0]*fSkin[13]+31.00135455709958*coeff[2]*fSkin[5]-10.22799180710709*coeff[1]*fSkin[5]-7.702348464916396*coeff[0]*fSkin[5]-15.90990257669732*coeff[2]*fSkin[3]+4.107919181288745*coeff[1]*fSkin[3]; 
  boundSurf_incr[14] = -1.988737822087164*coeff[1]*fSkin[18]; 
  boundSurf_incr[15] = (-2.223476479805891*coeff[2]*fSkin[15])+3.444594950788841*coeff[1]*fSkin[15]+1.988737822087164*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = -1.988737822087164*coeff[1]*fSkin[19]; 
  boundSurf_incr[17] = (-31.12867071728247*coeff[2]*fSkin[17])+12.05608232776095*coeff[1]*fSkin[17]+13.92116475461015*coeff[0]*fSkin[17]+31.00135455709956*coeff[2]*fSkin[10]-10.2279918071071*coeff[1]*fSkin[10]-7.702348464916393*coeff[0]*fSkin[10]-15.90990257669731*coeff[2]*fSkin[6]+4.107919181288745*coeff[1]*fSkin[6]; 
  boundSurf_incr[18] = (-2.223476479805891*coeff[2]*fSkin[18])+3.444594950788841*coeff[1]*fSkin[18]+1.988737822087164*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = (-2.223476479805891*coeff[2]*fSkin[19])+3.444594950788841*coeff[1]*fSkin[19]+1.988737822087164*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[7])+1.797214641813825*coeff[1]*fSkin[7]+4.743416490252569*coeff[0]*fSkin[7]+5.303300858899105*coeff[2]*fEdge[7]+1.797214641813825*coeff[1]*fEdge[7]-4.743416490252569*coeff[0]*fEdge[7]+6.418623720763661*fSkin[1]*coeff[2]+6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]-5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]-5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 11.19493359006374*coeff[2]*fSkin[7]-3.112867071728247*coeff[1]*fSkin[7]-10.01305300439132*coeff[0]*fSkin[7]-7.176239480810091*coeff[2]*fEdge[7]-3.112867071728247*coeff[1]*fEdge[7]+6.418623720763665*coeff[0]*fEdge[7]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]+6.418623720763661*fSkin[0]*coeff[2]-6.418623720763661*fEdge[0]*coeff[2]+1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]-1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]-5.74099158464807*coeff[0]*fSkin[0]+5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-5.303300858899106*coeff[2]*fSkin[11])+1.797214641813825*coeff[1]*fSkin[11]+4.743416490252569*coeff[0]*fSkin[11]+5.303300858899106*coeff[2]*fEdge[11]+1.797214641813825*coeff[1]*fEdge[11]-4.743416490252569*coeff[0]*fEdge[11]+6.418623720763661*coeff[2]*fSkin[4]-0.9943689110435817*coeff[1]*fSkin[4]-5.74099158464807*coeff[0]*fSkin[4]+6.418623720763661*coeff[2]*fEdge[4]+0.9943689110435817*coeff[1]*fEdge[4]-5.74099158464807*coeff[0]*fEdge[4]-3.705794133009818*coeff[2]*fSkin[2]+3.31456303681194*coeff[0]*fSkin[2]+3.705794133009818*coeff[2]*fEdge[2]-3.31456303681194*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-5.303300858899106*coeff[2]*fSkin[13])+1.797214641813825*coeff[1]*fSkin[13]+4.743416490252569*coeff[0]*fSkin[13]+5.303300858899106*coeff[2]*fEdge[13]+1.797214641813825*coeff[1]*fEdge[13]-4.743416490252569*coeff[0]*fEdge[13]+6.418623720763661*coeff[2]*fSkin[5]-0.9943689110435817*coeff[1]*fSkin[5]-5.74099158464807*coeff[0]*fSkin[5]+6.418623720763661*coeff[2]*fEdge[5]+0.9943689110435817*coeff[1]*fEdge[5]-5.74099158464807*coeff[0]*fEdge[5]-3.705794133009818*coeff[2]*fSkin[3]+3.31456303681194*coeff[0]*fSkin[3]+3.705794133009818*coeff[2]*fEdge[3]-3.31456303681194*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 11.19493359006374*coeff[2]*fSkin[11]-3.112867071728246*coeff[1]*fSkin[11]-10.01305300439132*coeff[0]*fSkin[11]-7.176239480810093*coeff[2]*fEdge[11]-3.112867071728246*coeff[1]*fEdge[11]+6.418623720763666*coeff[0]*fEdge[11]-12.2291206389324*coeff[2]*fSkin[4]+1.72229747539442*coeff[1]*fSkin[4]+10.9380580214794*coeff[0]*fSkin[4]-10.00564415912651*coeff[2]*fEdge[4]-1.72229747539442*coeff[1]*fEdge[4]+8.949320199392238*coeff[0]*fEdge[4]+6.418623720763661*coeff[2]*fSkin[2]-5.74099158464807*coeff[0]*fSkin[2]-6.418623720763661*coeff[2]*fEdge[2]+5.74099158464807*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 11.19493359006374*coeff[2]*fSkin[13]-3.112867071728246*coeff[1]*fSkin[13]-10.01305300439132*coeff[0]*fSkin[13]-7.176239480810093*coeff[2]*fEdge[13]-3.112867071728246*coeff[1]*fEdge[13]+6.418623720763666*coeff[0]*fEdge[13]-12.2291206389324*coeff[2]*fSkin[5]+1.72229747539442*coeff[1]*fSkin[5]+10.9380580214794*coeff[0]*fSkin[5]-10.00564415912651*coeff[2]*fEdge[5]-1.72229747539442*coeff[1]*fEdge[5]+8.949320199392238*coeff[0]*fEdge[5]+6.418623720763661*coeff[2]*fSkin[3]-5.74099158464807*coeff[0]*fSkin[3]-6.418623720763661*coeff[2]*fEdge[3]+5.74099158464807*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-5.303300858899105*coeff[2]*fSkin[17])+1.797214641813825*coeff[1]*fSkin[17]+4.743416490252569*coeff[0]*fSkin[17]+5.303300858899105*coeff[2]*fEdge[17]+1.797214641813825*coeff[1]*fEdge[17]-4.743416490252569*coeff[0]*fEdge[17]+6.418623720763661*coeff[2]*fSkin[10]-0.9943689110435817*coeff[1]*fSkin[10]-5.74099158464807*coeff[0]*fSkin[10]+6.418623720763661*coeff[2]*fEdge[10]+0.9943689110435817*coeff[1]*fEdge[10]-5.74099158464807*coeff[0]*fEdge[10]-3.705794133009818*coeff[2]*fSkin[6]+3.31456303681194*coeff[0]*fSkin[6]+3.705794133009818*coeff[2]*fEdge[6]-3.31456303681194*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-34.53800131965151*coeff[2]*fSkin[7])+11.53939308514262*coeff[1]*fSkin[7]+14.38520357976382*coeff[0]*fSkin[7]+3.409330602369041*coeff[2]*fEdge[7]+0.516689242618324*coeff[1]*fEdge[7]-0.4640388251536773*coeff[0]*fEdge[7]+42.4833377263957*fSkin[1]*coeff[2]+11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]-11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]-3.337684334797107*coeff[0]*fEdge[1]+9.756308055560766*fSkin[0]*coeff[1]-5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 6.418623720763661*coeff[2]*fSkin[12]-0.9943689110435818*coeff[1]*fSkin[12]-5.740991584648071*coeff[0]*fSkin[12]+6.418623720763661*coeff[2]*fEdge[12]+0.9943689110435818*coeff[1]*fEdge[12]-5.740991584648071*coeff[0]*fEdge[12]-3.705794133009818*coeff[2]*fSkin[8]+3.31456303681194*coeff[0]*fSkin[8]+3.705794133009818*coeff[2]*fEdge[8]-3.31456303681194*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 6.418623720763661*coeff[2]*fSkin[15]-0.9943689110435818*coeff[1]*fSkin[15]-5.740991584648071*coeff[0]*fSkin[15]+6.418623720763661*coeff[2]*fEdge[15]+0.9943689110435818*coeff[1]*fEdge[15]-5.740991584648071*coeff[0]*fEdge[15]-3.705794133009818*coeff[2]*fSkin[9]+3.31456303681194*coeff[0]*fSkin[9]+3.705794133009818*coeff[2]*fEdge[9]-3.31456303681194*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 11.19493359006374*coeff[2]*fSkin[17]-3.112867071728247*coeff[1]*fSkin[17]-10.01305300439132*coeff[0]*fSkin[17]-7.176239480810091*coeff[2]*fEdge[17]-3.112867071728247*coeff[1]*fEdge[17]+6.418623720763665*coeff[0]*fEdge[17]-12.2291206389324*coeff[2]*fSkin[10]+1.72229747539442*coeff[1]*fSkin[10]+10.9380580214794*coeff[0]*fSkin[10]-10.00564415912651*coeff[2]*fEdge[10]-1.72229747539442*coeff[1]*fEdge[10]+8.949320199392238*coeff[0]*fEdge[10]+6.418623720763661*coeff[2]*fSkin[6]-5.74099158464807*coeff[0]*fSkin[6]-6.418623720763661*coeff[2]*fEdge[6]+5.74099158464807*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-34.53800131965151*coeff[2]*fSkin[11])+11.53939308514262*coeff[1]*fSkin[11]+14.38520357976382*coeff[0]*fSkin[11]+3.409330602369041*coeff[2]*fEdge[11]+0.516689242618324*coeff[1]*fEdge[11]-0.4640388251536773*coeff[0]*fEdge[11]+42.48333772639572*coeff[2]*fSkin[4]-14.89729241469946*coeff[1]*fSkin[4]-11.0400327997135*coeff[0]*fSkin[4]+11.48198316929615*coeff[2]*fEdge[4]-4.669300607592371*coeff[1]*fEdge[4]-3.337684334797105*coeff[0]*fEdge[4]-26.18504799081433*coeff[2]*fSkin[2]+9.756308055560769*coeff[1]*fSkin[2]+2.964635306407854*coeff[0]*fSkin[2]+10.27514541411701*coeff[2]*fEdge[2]-5.648388874272023*coeff[1]*fEdge[2]-2.964635306407854*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-12.2291206389324*coeff[2]*fSkin[12])+1.72229747539442*coeff[1]*fSkin[12]+10.9380580214794*coeff[0]*fSkin[12]-10.00564415912651*coeff[2]*fEdge[12]-1.72229747539442*coeff[1]*fEdge[12]+8.949320199392238*coeff[0]*fEdge[12]+6.418623720763661*coeff[2]*fSkin[8]-5.740991584648071*coeff[0]*fSkin[8]-6.418623720763661*coeff[2]*fEdge[8]+5.740991584648071*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-34.53800131965151*coeff[2]*fSkin[13])+11.53939308514262*coeff[1]*fSkin[13]+14.38520357976382*coeff[0]*fSkin[13]+3.409330602369041*coeff[2]*fEdge[13]+0.516689242618324*coeff[1]*fEdge[13]-0.4640388251536773*coeff[0]*fEdge[13]+42.48333772639572*coeff[2]*fSkin[5]-14.89729241469946*coeff[1]*fSkin[5]-11.0400327997135*coeff[0]*fSkin[5]+11.48198316929615*coeff[2]*fEdge[5]-4.669300607592371*coeff[1]*fEdge[5]-3.337684334797105*coeff[0]*fEdge[5]-26.18504799081433*coeff[2]*fSkin[3]+9.756308055560769*coeff[1]*fSkin[3]+2.964635306407854*coeff[0]*fSkin[3]+10.27514541411701*coeff[2]*fEdge[3]-5.648388874272023*coeff[1]*fEdge[3]-2.964635306407854*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 6.418623720763661*coeff[2]*fSkin[18]-0.9943689110435818*coeff[1]*fSkin[18]-5.740991584648071*coeff[0]*fSkin[18]+6.418623720763661*coeff[2]*fEdge[18]+0.9943689110435818*coeff[1]*fEdge[18]-5.740991584648071*coeff[0]*fEdge[18]-3.705794133009818*coeff[2]*fSkin[14]+3.31456303681194*coeff[0]*fSkin[14]+3.705794133009818*coeff[2]*fEdge[14]-3.31456303681194*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-12.2291206389324*coeff[2]*fSkin[15])+1.72229747539442*coeff[1]*fSkin[15]+10.9380580214794*coeff[0]*fSkin[15]-10.00564415912651*coeff[2]*fEdge[15]-1.72229747539442*coeff[1]*fEdge[15]+8.949320199392238*coeff[0]*fEdge[15]+6.418623720763661*coeff[2]*fSkin[9]-5.740991584648071*coeff[0]*fSkin[9]-6.418623720763661*coeff[2]*fEdge[9]+5.740991584648071*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 6.418623720763661*coeff[2]*fSkin[19]-0.9943689110435818*coeff[1]*fSkin[19]-5.740991584648071*coeff[0]*fSkin[19]+6.418623720763661*coeff[2]*fEdge[19]+0.9943689110435818*coeff[1]*fEdge[19]-5.740991584648071*coeff[0]*fEdge[19]-3.705794133009818*coeff[2]*fSkin[16]+3.31456303681194*coeff[0]*fSkin[16]+3.705794133009818*coeff[2]*fEdge[16]-3.31456303681194*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-34.53800131965151*coeff[2]*fSkin[17])+11.53939308514262*coeff[1]*fSkin[17]+14.38520357976382*coeff[0]*fSkin[17]+3.409330602369041*coeff[2]*fEdge[17]+0.516689242618324*coeff[1]*fEdge[17]-0.4640388251536773*coeff[0]*fEdge[17]+42.4833377263957*coeff[2]*fSkin[10]-14.89729241469947*coeff[1]*fSkin[10]-11.0400327997135*coeff[0]*fSkin[10]+11.48198316929614*coeff[2]*fEdge[10]-4.669300607592371*coeff[1]*fEdge[10]-3.337684334797107*coeff[0]*fEdge[10]-26.18504799081432*coeff[2]*fSkin[6]+9.756308055560766*coeff[1]*fSkin[6]+2.964635306407856*coeff[0]*fSkin[6]+10.27514541411701*coeff[2]*fEdge[6]-5.648388874272021*coeff[1]*fEdge[6]-2.964635306407856*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-12.2291206389324*coeff[2]*fSkin[18])+1.72229747539442*coeff[1]*fSkin[18]+10.9380580214794*coeff[0]*fSkin[18]-10.00564415912651*coeff[2]*fEdge[18]-1.72229747539442*coeff[1]*fEdge[18]+8.949320199392238*coeff[0]*fEdge[18]+6.418623720763661*coeff[2]*fSkin[14]-5.740991584648071*coeff[0]*fSkin[14]-6.418623720763661*coeff[2]*fEdge[14]+5.740991584648071*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-12.2291206389324*coeff[2]*fSkin[19])+1.72229747539442*coeff[1]*fSkin[19]+10.9380580214794*coeff[0]*fSkin[19]-10.00564415912651*coeff[2]*fEdge[19]-1.72229747539442*coeff[1]*fEdge[19]+8.949320199392238*coeff[0]*fEdge[19]+6.418623720763661*coeff[2]*fSkin[16]-5.740991584648071*coeff[0]*fSkin[16]-6.418623720763661*coeff[2]*fEdge[16]+5.740991584648071*coeff[0]*fEdge[16]; 

  boundSurf_incr[0] = (-3.59442928362765*coeff[1]*fSkin[7])-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-4.018694109253648*coeff[2]*fSkin[7])-6.225734143456494*coeff[1]*fSkin[7]+3.59442928362765*coeff[0]*fSkin[7]-2.223476479805891*fSkin[1]*coeff[2]-3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-3.594429283627651*coeff[1]*fSkin[11])-1.988737822087164*coeff[1]*fSkin[4]; 
  boundSurf_incr[3] = (-3.594429283627651*coeff[1]*fSkin[13])-1.988737822087164*coeff[1]*fSkin[5]; 
  boundSurf_incr[4] = (-4.018694109253649*coeff[2]*fSkin[11])-6.225734143456493*coeff[1]*fSkin[11]+3.594429283627651*coeff[0]*fSkin[11]-2.223476479805891*coeff[2]*fSkin[4]-3.444594950788841*coeff[1]*fSkin[4]+1.988737822087164*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = (-4.018694109253649*coeff[2]*fSkin[13])-6.225734143456493*coeff[1]*fSkin[13]+3.594429283627651*coeff[0]*fSkin[13]-2.223476479805891*coeff[2]*fSkin[5]-3.444594950788841*coeff[1]*fSkin[5]+1.988737822087164*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = (-3.59442928362765*coeff[1]*fSkin[17])-1.988737822087164*coeff[1]*fSkin[10]; 
  boundSurf_incr[7] = (-31.12867071728247*coeff[2]*fSkin[7])-12.05608232776095*coeff[1]*fSkin[7]+13.92116475461015*coeff[0]*fSkin[7]-31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]+7.702348464916393*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]; 
  boundSurf_incr[8] = -1.988737822087164*coeff[1]*fSkin[12]; 
  boundSurf_incr[9] = -1.988737822087164*coeff[1]*fSkin[15]; 
  boundSurf_incr[10] = (-4.018694109253648*coeff[2]*fSkin[17])-6.225734143456494*coeff[1]*fSkin[17]+3.59442928362765*coeff[0]*fSkin[17]-2.223476479805891*coeff[2]*fSkin[10]-3.444594950788841*coeff[1]*fSkin[10]+1.988737822087164*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = (-31.12867071728247*coeff[2]*fSkin[11])-12.05608232776095*coeff[1]*fSkin[11]+13.92116475461015*coeff[0]*fSkin[11]-31.00135455709958*coeff[2]*fSkin[4]-10.22799180710709*coeff[1]*fSkin[4]+7.702348464916396*coeff[0]*fSkin[4]-15.90990257669732*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]; 
  boundSurf_incr[12] = (-2.223476479805891*coeff[2]*fSkin[12])-3.444594950788841*coeff[1]*fSkin[12]+1.988737822087164*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = (-31.12867071728247*coeff[2]*fSkin[13])-12.05608232776095*coeff[1]*fSkin[13]+13.92116475461015*coeff[0]*fSkin[13]-31.00135455709958*coeff[2]*fSkin[5]-10.22799180710709*coeff[1]*fSkin[5]+7.702348464916396*coeff[0]*fSkin[5]-15.90990257669732*coeff[2]*fSkin[3]-4.107919181288745*coeff[1]*fSkin[3]; 
  boundSurf_incr[14] = -1.988737822087164*coeff[1]*fSkin[18]; 
  boundSurf_incr[15] = (-2.223476479805891*coeff[2]*fSkin[15])-3.444594950788841*coeff[1]*fSkin[15]+1.988737822087164*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = -1.988737822087164*coeff[1]*fSkin[19]; 
  boundSurf_incr[17] = (-31.12867071728247*coeff[2]*fSkin[17])-12.05608232776095*coeff[1]*fSkin[17]+13.92116475461015*coeff[0]*fSkin[17]-31.00135455709956*coeff[2]*fSkin[10]-10.2279918071071*coeff[1]*fSkin[10]+7.702348464916393*coeff[0]*fSkin[10]-15.90990257669731*coeff[2]*fSkin[6]-4.107919181288745*coeff[1]*fSkin[6]; 
  boundSurf_incr[18] = (-2.223476479805891*coeff[2]*fSkin[18])-3.444594950788841*coeff[1]*fSkin[18]+1.988737822087164*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = (-2.223476479805891*coeff[2]*fSkin[19])-3.444594950788841*coeff[1]*fSkin[19]+1.988737822087164*coeff[0]*fSkin[19]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += -1.0*(vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += -1.0*(vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += -1.0*(vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += -1.0*(vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += -1.0*(vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += -1.0*(vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += -1.0*(vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += -1.0*(vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += -1.0*(vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += -1.0*(vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += -1.0*(vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += -1.0*(vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 

  }

  return 0.;
}

