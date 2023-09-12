#include <gkyl_dg_diffusion_vlasov_kernels.h>

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_1x3v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]+8.11898816047911*coeff[0]*fSkin[1]+8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 14.16059535957086*coeff[0]*fSkin[11]-9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]+8.11898816047911*coeff[0]*fSkin[0]-8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[19]-6.708203932499369*coeff[0]*fEdge[19]+8.11898816047911*coeff[0]*fSkin[5]+8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[21]-6.708203932499369*coeff[0]*fEdge[21]+8.11898816047911*coeff[0]*fSkin[6]+8.11898816047911*coeff[0]*fEdge[6]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[0]*fSkin[25]-6.708203932499369*coeff[0]*fEdge[25]+8.11898816047911*coeff[0]*fSkin[8]+8.11898816047911*coeff[0]*fEdge[8]+4.6875*coeff[0]*fSkin[4]-4.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 14.16059535957087*coeff[0]*fSkin[19]-9.077304717673634*coeff[0]*fEdge[19]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]+8.11898816047911*coeff[0]*fSkin[2]-8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 14.16059535957087*coeff[0]*fSkin[21]-9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[6]+12.65625*coeff[0]*fEdge[6]+8.11898816047911*coeff[0]*fSkin[3]-8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*fSkin[32]-6.708203932499369*coeff[0]*fEdge[32]+8.11898816047911*coeff[0]*fSkin[15]+8.11898816047911*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[7]-4.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 14.16059535957087*coeff[0]*fSkin[25]-9.077304717673634*coeff[0]*fEdge[25]+15.46875*coeff[0]*fSkin[8]+12.65625*coeff[0]*fEdge[8]+8.11898816047911*coeff[0]*fSkin[4]-8.11898816047911*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[35]-6.708203932499369*coeff[0]*fEdge[35]+8.11898816047911*coeff[0]*fSkin[16]+8.11898816047911*coeff[0]*fEdge[16]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*fSkin[37]-6.708203932499369*coeff[0]*fEdge[37]+8.11898816047911*coeff[0]*fSkin[17]+8.11898816047911*coeff[0]*fEdge[17]+4.6875*coeff[0]*fSkin[10]-4.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]+15.61296411439865*coeff[0]*fSkin[1]+4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 8.118988160479114*coeff[0]*fSkin[20]+8.118988160479114*coeff[0]*fEdge[20]+4.6875*coeff[0]*fSkin[12]-4.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 8.118988160479114*coeff[0]*fSkin[23]+8.118988160479114*coeff[0]*fEdge[23]+4.6875*coeff[0]*fSkin[13]-4.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 8.118988160479114*coeff[0]*fSkin[28]+8.118988160479114*coeff[0]*fEdge[28]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 14.16059535957086*coeff[0]*fSkin[32]-9.077304717673634*coeff[0]*fEdge[32]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]+8.11898816047911*coeff[0]*fSkin[7]-8.11898816047911*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 14.16059535957086*coeff[0]*fSkin[35]-9.077304717673634*coeff[0]*fEdge[35]+15.46875*coeff[0]*fSkin[16]+12.65625*coeff[0]*fEdge[16]+8.11898816047911*coeff[0]*fSkin[9]-8.11898816047911*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 14.16059535957086*coeff[0]*fSkin[37]-9.077304717673634*coeff[0]*fEdge[37]+15.46875*coeff[0]*fSkin[17]+12.65625*coeff[0]*fEdge[17]+8.11898816047911*coeff[0]*fSkin[10]-8.11898816047911*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 6.708203932499369*coeff[0]*fSkin[44]-6.708203932499369*coeff[0]*fEdge[44]+8.11898816047911*coeff[0]*fSkin[31]+8.11898816047911*coeff[0]*fEdge[31]+4.6875*coeff[0]*fSkin[18]-4.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*fSkin[19]-0.65625*coeff[0]*fEdge[19]+15.61296411439865*coeff[0]*fSkin[5]+4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*fSkin[20]+12.65625*coeff[0]*fEdge[20]+8.118988160479114*coeff[0]*fSkin[12]-8.118988160479114*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]+15.61296411439865*coeff[0]*fSkin[6]+4.72019845319029*coeff[0]*fEdge[6]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 8.118988160479114*coeff[0]*fSkin[33]+8.118988160479114*coeff[0]*fEdge[33]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*fSkin[23]+12.65625*coeff[0]*fEdge[23]+8.118988160479114*coeff[0]*fSkin[13]-8.118988160479114*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 8.118988160479114*coeff[0]*fSkin[34]+8.118988160479114*coeff[0]*fEdge[34]+4.6875*coeff[0]*fSkin[24]-4.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*fSkin[25]-0.65625*coeff[0]*fEdge[25]+15.61296411439865*coeff[0]*fSkin[8]+4.72019845319029*coeff[0]*fEdge[8]+4.192627457812107*coeff[0]*fSkin[4]-4.192627457812107*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 8.118988160479114*coeff[0]*fSkin[36]+8.118988160479114*coeff[0]*fEdge[36]+4.6875*coeff[0]*fSkin[26]-4.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 8.118988160479114*coeff[0]*fSkin[39]+8.118988160479114*coeff[0]*fEdge[39]+4.6875*coeff[0]*fSkin[27]-4.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*fSkin[28]+12.65625*coeff[0]*fEdge[28]+8.118988160479114*coeff[0]*fSkin[14]-8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 8.118988160479114*coeff[0]*fSkin[41]+8.118988160479114*coeff[0]*fEdge[41]+4.6875*coeff[0]*fSkin[29]-4.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 8.118988160479114*coeff[0]*fSkin[42]+8.118988160479114*coeff[0]*fEdge[42]+4.6875*coeff[0]*fSkin[30]-4.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 14.16059535957087*coeff[0]*fSkin[44]-9.077304717673634*coeff[0]*fEdge[44]+15.46875*coeff[0]*fSkin[31]+12.65625*coeff[0]*fEdge[31]+8.11898816047911*coeff[0]*fSkin[18]-8.11898816047911*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*fSkin[32]-0.65625*coeff[0]*fEdge[32]+15.61296411439865*coeff[0]*fSkin[15]+4.720198453190292*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[7]-4.192627457812107*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*fSkin[33]+12.65625*coeff[0]*fEdge[33]+8.118988160479114*coeff[0]*fSkin[22]-8.118988160479114*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*fSkin[34]+12.65625*coeff[0]*fEdge[34]+8.118988160479114*coeff[0]*fSkin[24]-8.118988160479114*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*fSkin[35]-0.65625*coeff[0]*fEdge[35]+15.61296411439865*coeff[0]*fSkin[16]+4.720198453190292*coeff[0]*fEdge[16]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*fSkin[36]+12.65625*coeff[0]*fEdge[36]+8.118988160479114*coeff[0]*fSkin[26]-8.118988160479114*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*fSkin[37]-0.65625*coeff[0]*fEdge[37]+15.61296411439865*coeff[0]*fSkin[17]+4.720198453190292*coeff[0]*fEdge[17]+4.192627457812107*coeff[0]*fSkin[10]-4.192627457812107*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 8.118988160479114*coeff[0]*fSkin[45]+8.118988160479114*coeff[0]*fEdge[45]+4.6875*coeff[0]*fSkin[38]-4.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*fSkin[39]+12.65625*coeff[0]*fEdge[39]+8.118988160479114*coeff[0]*fSkin[27]-8.118988160479114*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 8.118988160479114*coeff[0]*fSkin[46]+8.118988160479114*coeff[0]*fEdge[46]+4.6875*coeff[0]*fSkin[40]-4.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*fSkin[41]+12.65625*coeff[0]*fEdge[41]+8.118988160479114*coeff[0]*fSkin[29]-8.118988160479114*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*fSkin[42]+12.65625*coeff[0]*fEdge[42]+8.118988160479114*coeff[0]*fSkin[30]-8.118988160479114*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 8.118988160479114*coeff[0]*fSkin[47]+8.118988160479114*coeff[0]*fEdge[47]+4.6875*coeff[0]*fSkin[43]-4.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*fSkin[44]-0.65625*coeff[0]*fEdge[44]+15.61296411439865*coeff[0]*fSkin[31]+4.72019845319029*coeff[0]*fEdge[31]+4.192627457812107*coeff[0]*fSkin[18]-4.192627457812107*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*fSkin[45]+12.65625*coeff[0]*fEdge[45]+8.118988160479114*coeff[0]*fSkin[38]-8.118988160479114*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*fSkin[46]+12.65625*coeff[0]*fEdge[46]+8.118988160479114*coeff[0]*fSkin[40]-8.118988160479114*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*fSkin[47]+12.65625*coeff[0]*fEdge[47]+8.118988160479114*coeff[0]*fSkin[43]-8.118988160479114*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 2.8125*coeff[0]*fSkin[1]-5.083290641897234*coeff[0]*fSkin[11]; 
  boundSurf_incr[5] = 2.8125*coeff[0]*fSkin[5]-5.083290641897235*coeff[0]*fSkin[19]; 
  boundSurf_incr[6] = 2.8125*coeff[0]*fSkin[6]-5.083290641897235*coeff[0]*fSkin[21]; 
  boundSurf_incr[8] = 2.8125*coeff[0]*fSkin[8]-5.083290641897235*coeff[0]*fSkin[25]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]-10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = 2.8125*coeff[0]*fSkin[15]-5.083290641897234*coeff[0]*fSkin[32]; 
  boundSurf_incr[16] = 2.8125*coeff[0]*fSkin[16]-5.083290641897234*coeff[0]*fSkin[35]; 
  boundSurf_incr[17] = 2.8125*coeff[0]*fSkin[17]-5.083290641897234*coeff[0]*fSkin[37]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*fSkin[19]-10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]-10.89276566120836*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*fSkin[25]-10.89276566120836*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 2.8125*coeff[0]*fSkin[31]-5.083290641897235*coeff[0]*fSkin[44]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*fSkin[32]-10.89276566120836*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*fSkin[35]-10.89276566120836*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*fSkin[37]-10.89276566120836*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*fSkin[44]-10.89276566120836*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*fSkin[47]; 

  } else { 

  edgeSurf_incr[0] = 6.708203932499369*coeff[0]*fSkin[11]-6.708203932499369*coeff[0]*fEdge[11]-8.11898816047911*coeff[0]*fSkin[1]-8.11898816047911*coeff[0]*fEdge[1]+4.6875*coeff[0]*fSkin[0]-4.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-14.16059535957086*coeff[0]*fSkin[11])+9.077304717673634*coeff[0]*fEdge[11]+15.46875*coeff[0]*fSkin[1]+12.65625*coeff[0]*fEdge[1]-8.11898816047911*coeff[0]*fSkin[0]+8.11898816047911*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 6.708203932499369*coeff[0]*fSkin[19]-6.708203932499369*coeff[0]*fEdge[19]-8.11898816047911*coeff[0]*fSkin[5]-8.11898816047911*coeff[0]*fEdge[5]+4.6875*coeff[0]*fSkin[2]-4.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 6.708203932499369*coeff[0]*fSkin[21]-6.708203932499369*coeff[0]*fEdge[21]-8.11898816047911*coeff[0]*fSkin[6]-8.11898816047911*coeff[0]*fEdge[6]+4.6875*coeff[0]*fSkin[3]-4.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 6.708203932499369*coeff[0]*fSkin[25]-6.708203932499369*coeff[0]*fEdge[25]-8.11898816047911*coeff[0]*fSkin[8]-8.11898816047911*coeff[0]*fEdge[8]+4.6875*coeff[0]*fSkin[4]-4.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-14.16059535957087*coeff[0]*fSkin[19])+9.077304717673634*coeff[0]*fEdge[19]+15.46875*coeff[0]*fSkin[5]+12.65625*coeff[0]*fEdge[5]-8.11898816047911*coeff[0]*fSkin[2]+8.11898816047911*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-14.16059535957087*coeff[0]*fSkin[21])+9.077304717673634*coeff[0]*fEdge[21]+15.46875*coeff[0]*fSkin[6]+12.65625*coeff[0]*fEdge[6]-8.11898816047911*coeff[0]*fSkin[3]+8.11898816047911*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 6.708203932499369*coeff[0]*fSkin[32]-6.708203932499369*coeff[0]*fEdge[32]-8.11898816047911*coeff[0]*fSkin[15]-8.11898816047911*coeff[0]*fEdge[15]+4.6875*coeff[0]*fSkin[7]-4.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-14.16059535957087*coeff[0]*fSkin[25])+9.077304717673634*coeff[0]*fEdge[25]+15.46875*coeff[0]*fSkin[8]+12.65625*coeff[0]*fEdge[8]-8.11898816047911*coeff[0]*fSkin[4]+8.11898816047911*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 6.708203932499369*coeff[0]*fSkin[35]-6.708203932499369*coeff[0]*fEdge[35]-8.11898816047911*coeff[0]*fSkin[16]-8.11898816047911*coeff[0]*fEdge[16]+4.6875*coeff[0]*fSkin[9]-4.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 6.708203932499369*coeff[0]*fSkin[37]-6.708203932499369*coeff[0]*fEdge[37]-8.11898816047911*coeff[0]*fSkin[17]-8.11898816047911*coeff[0]*fEdge[17]+4.6875*coeff[0]*fSkin[10]-4.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 20.34375*coeff[0]*fSkin[11]-0.65625*coeff[0]*fEdge[11]-15.61296411439865*coeff[0]*fSkin[1]-4.720198453190292*coeff[0]*fEdge[1]+4.192627457812107*coeff[0]*fSkin[0]-4.192627457812107*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = (-8.118988160479114*coeff[0]*fSkin[20])-8.118988160479114*coeff[0]*fEdge[20]+4.6875*coeff[0]*fSkin[12]-4.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = (-8.118988160479114*coeff[0]*fSkin[23])-8.118988160479114*coeff[0]*fEdge[23]+4.6875*coeff[0]*fSkin[13]-4.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = (-8.118988160479114*coeff[0]*fSkin[28])-8.118988160479114*coeff[0]*fEdge[28]+4.6875*coeff[0]*fSkin[14]-4.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-14.16059535957086*coeff[0]*fSkin[32])+9.077304717673634*coeff[0]*fEdge[32]+15.46875*coeff[0]*fSkin[15]+12.65625*coeff[0]*fEdge[15]-8.11898816047911*coeff[0]*fSkin[7]+8.11898816047911*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = (-14.16059535957086*coeff[0]*fSkin[35])+9.077304717673634*coeff[0]*fEdge[35]+15.46875*coeff[0]*fSkin[16]+12.65625*coeff[0]*fEdge[16]-8.11898816047911*coeff[0]*fSkin[9]+8.11898816047911*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = (-14.16059535957086*coeff[0]*fSkin[37])+9.077304717673634*coeff[0]*fEdge[37]+15.46875*coeff[0]*fSkin[17]+12.65625*coeff[0]*fEdge[17]-8.11898816047911*coeff[0]*fSkin[10]+8.11898816047911*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 6.708203932499369*coeff[0]*fSkin[44]-6.708203932499369*coeff[0]*fEdge[44]-8.11898816047911*coeff[0]*fSkin[31]-8.11898816047911*coeff[0]*fEdge[31]+4.6875*coeff[0]*fSkin[18]-4.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 20.34375*coeff[0]*fSkin[19]-0.65625*coeff[0]*fEdge[19]-15.61296411439865*coeff[0]*fSkin[5]-4.72019845319029*coeff[0]*fEdge[5]+4.192627457812107*coeff[0]*fSkin[2]-4.192627457812107*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 15.46875*coeff[0]*fSkin[20]+12.65625*coeff[0]*fEdge[20]-8.118988160479114*coeff[0]*fSkin[12]+8.118988160479114*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 20.34375*coeff[0]*fSkin[21]-0.65625*coeff[0]*fEdge[21]-15.61296411439865*coeff[0]*fSkin[6]-4.72019845319029*coeff[0]*fEdge[6]+4.192627457812107*coeff[0]*fSkin[3]-4.192627457812107*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = (-8.118988160479114*coeff[0]*fSkin[33])-8.118988160479114*coeff[0]*fEdge[33]+4.6875*coeff[0]*fSkin[22]-4.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 15.46875*coeff[0]*fSkin[23]+12.65625*coeff[0]*fEdge[23]-8.118988160479114*coeff[0]*fSkin[13]+8.118988160479114*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = (-8.118988160479114*coeff[0]*fSkin[34])-8.118988160479114*coeff[0]*fEdge[34]+4.6875*coeff[0]*fSkin[24]-4.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 20.34375*coeff[0]*fSkin[25]-0.65625*coeff[0]*fEdge[25]-15.61296411439865*coeff[0]*fSkin[8]-4.72019845319029*coeff[0]*fEdge[8]+4.192627457812107*coeff[0]*fSkin[4]-4.192627457812107*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = (-8.118988160479114*coeff[0]*fSkin[36])-8.118988160479114*coeff[0]*fEdge[36]+4.6875*coeff[0]*fSkin[26]-4.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = (-8.118988160479114*coeff[0]*fSkin[39])-8.118988160479114*coeff[0]*fEdge[39]+4.6875*coeff[0]*fSkin[27]-4.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 15.46875*coeff[0]*fSkin[28]+12.65625*coeff[0]*fEdge[28]-8.118988160479114*coeff[0]*fSkin[14]+8.118988160479114*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = (-8.118988160479114*coeff[0]*fSkin[41])-8.118988160479114*coeff[0]*fEdge[41]+4.6875*coeff[0]*fSkin[29]-4.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = (-8.118988160479114*coeff[0]*fSkin[42])-8.118988160479114*coeff[0]*fEdge[42]+4.6875*coeff[0]*fSkin[30]-4.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = (-14.16059535957087*coeff[0]*fSkin[44])+9.077304717673634*coeff[0]*fEdge[44]+15.46875*coeff[0]*fSkin[31]+12.65625*coeff[0]*fEdge[31]-8.11898816047911*coeff[0]*fSkin[18]+8.11898816047911*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 20.34375*coeff[0]*fSkin[32]-0.65625*coeff[0]*fEdge[32]-15.61296411439865*coeff[0]*fSkin[15]-4.720198453190292*coeff[0]*fEdge[15]+4.192627457812107*coeff[0]*fSkin[7]-4.192627457812107*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 15.46875*coeff[0]*fSkin[33]+12.65625*coeff[0]*fEdge[33]-8.118988160479114*coeff[0]*fSkin[22]+8.118988160479114*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 15.46875*coeff[0]*fSkin[34]+12.65625*coeff[0]*fEdge[34]-8.118988160479114*coeff[0]*fSkin[24]+8.118988160479114*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 20.34375*coeff[0]*fSkin[35]-0.65625*coeff[0]*fEdge[35]-15.61296411439865*coeff[0]*fSkin[16]-4.720198453190292*coeff[0]*fEdge[16]+4.192627457812107*coeff[0]*fSkin[9]-4.192627457812107*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 15.46875*coeff[0]*fSkin[36]+12.65625*coeff[0]*fEdge[36]-8.118988160479114*coeff[0]*fSkin[26]+8.118988160479114*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 20.34375*coeff[0]*fSkin[37]-0.65625*coeff[0]*fEdge[37]-15.61296411439865*coeff[0]*fSkin[17]-4.720198453190292*coeff[0]*fEdge[17]+4.192627457812107*coeff[0]*fSkin[10]-4.192627457812107*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = (-8.118988160479114*coeff[0]*fSkin[45])-8.118988160479114*coeff[0]*fEdge[45]+4.6875*coeff[0]*fSkin[38]-4.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 15.46875*coeff[0]*fSkin[39]+12.65625*coeff[0]*fEdge[39]-8.118988160479114*coeff[0]*fSkin[27]+8.118988160479114*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = (-8.118988160479114*coeff[0]*fSkin[46])-8.118988160479114*coeff[0]*fEdge[46]+4.6875*coeff[0]*fSkin[40]-4.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 15.46875*coeff[0]*fSkin[41]+12.65625*coeff[0]*fEdge[41]-8.118988160479114*coeff[0]*fSkin[29]+8.118988160479114*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 15.46875*coeff[0]*fSkin[42]+12.65625*coeff[0]*fEdge[42]-8.118988160479114*coeff[0]*fSkin[30]+8.118988160479114*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = (-8.118988160479114*coeff[0]*fSkin[47])-8.118988160479114*coeff[0]*fEdge[47]+4.6875*coeff[0]*fSkin[43]-4.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 20.34375*coeff[0]*fSkin[44]-0.65625*coeff[0]*fEdge[44]-15.61296411439865*coeff[0]*fSkin[31]-4.72019845319029*coeff[0]*fEdge[31]+4.192627457812107*coeff[0]*fSkin[18]-4.192627457812107*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 15.46875*coeff[0]*fSkin[45]+12.65625*coeff[0]*fEdge[45]-8.118988160479114*coeff[0]*fSkin[38]+8.118988160479114*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 15.46875*coeff[0]*fSkin[46]+12.65625*coeff[0]*fEdge[46]-8.118988160479114*coeff[0]*fSkin[40]+8.118988160479114*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 15.46875*coeff[0]*fSkin[47]+12.65625*coeff[0]*fEdge[47]-8.118988160479114*coeff[0]*fSkin[43]+8.118988160479114*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 5.083290641897234*coeff[0]*fSkin[11]+2.8125*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 5.083290641897235*coeff[0]*fSkin[19]+2.8125*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 5.083290641897235*coeff[0]*fSkin[21]+2.8125*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 5.083290641897235*coeff[0]*fSkin[25]+2.8125*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 19.6875*coeff[0]*fSkin[11]+10.89276566120836*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = 5.083290641897234*coeff[0]*fSkin[32]+2.8125*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 5.083290641897234*coeff[0]*fSkin[35]+2.8125*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 5.083290641897234*coeff[0]*fSkin[37]+2.8125*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = 19.6875*coeff[0]*fSkin[19]+10.89276566120836*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 2.8125*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 19.6875*coeff[0]*fSkin[21]+10.89276566120836*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = 2.8125*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = 19.6875*coeff[0]*fSkin[25]+10.89276566120836*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = 2.8125*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 5.083290641897235*coeff[0]*fSkin[44]+2.8125*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 19.6875*coeff[0]*fSkin[32]+10.89276566120836*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 2.8125*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 2.8125*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 19.6875*coeff[0]*fSkin[35]+10.89276566120836*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 2.8125*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 19.6875*coeff[0]*fSkin[37]+10.89276566120836*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = 2.8125*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = 2.8125*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 2.8125*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = 19.6875*coeff[0]*fSkin[44]+10.89276566120836*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 2.8125*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 2.8125*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 2.8125*coeff[0]*fSkin[47]; 

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
  out[20] += -1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += -1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += -1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += -1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += -1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += -1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += -1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += -1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += -1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += -1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += -1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += -1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += -1.0*(vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += -1.0*(vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += -1.0*(vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += -1.0*(vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += -1.0*(vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += -1.0*(vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += -1.0*(vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += -1.0*(vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += -1.0*(vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += -1.0*(vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += -1.0*(vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += -1.0*(vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += -1.0*(vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += -1.0*(vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += -1.0*(vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += -1.0*(vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_vlasov_order4_boundary_surfx_1x3v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[48] = {0.0}; 
  vol_incr[11] = -31.81980515339464*fSkin[0]*coeff[2]; 
  vol_incr[19] = -31.81980515339463*coeff[2]*fSkin[2]; 
  vol_incr[21] = -31.81980515339463*coeff[2]*fSkin[3]; 
  vol_incr[25] = -31.81980515339463*coeff[2]*fSkin[4]; 
  vol_incr[32] = -31.81980515339464*coeff[2]*fSkin[7]; 
  vol_incr[35] = -31.81980515339464*coeff[2]*fSkin[9]; 
  vol_incr[37] = -31.81980515339464*coeff[2]*fSkin[10]; 
  vol_incr[44] = -31.81980515339463*coeff[2]*fSkin[18]; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[11])-1.797214641813825*coeff[1]*fSkin[11]+4.743416490252569*coeff[0]*fSkin[11]+5.303300858899105*coeff[2]*fEdge[11]-1.797214641813825*coeff[1]*fEdge[11]-4.743416490252569*coeff[0]*fEdge[11]-6.418623720763661*fSkin[1]*coeff[2]-6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]+5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]+5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-11.19493359006374*coeff[2]*fSkin[11])-3.112867071728247*coeff[1]*fSkin[11]+10.01305300439132*coeff[0]*fSkin[11]+7.176239480810091*coeff[2]*fEdge[11]-3.112867071728247*coeff[1]*fEdge[11]-6.418623720763665*coeff[0]*fEdge[11]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]-6.418623720763661*fSkin[0]*coeff[2]+6.418623720763661*fEdge[0]*coeff[2]-1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]+1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]+5.74099158464807*coeff[0]*fSkin[0]-5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-5.303300858899106*coeff[2]*fSkin[19])-1.797214641813825*coeff[1]*fSkin[19]+4.743416490252569*coeff[0]*fSkin[19]+5.303300858899106*coeff[2]*fEdge[19]-1.797214641813825*coeff[1]*fEdge[19]-4.743416490252569*coeff[0]*fEdge[19]-6.418623720763661*coeff[2]*fSkin[5]-0.9943689110435817*coeff[1]*fSkin[5]+5.74099158464807*coeff[0]*fSkin[5]-6.418623720763661*coeff[2]*fEdge[5]+0.9943689110435817*coeff[1]*fEdge[5]+5.74099158464807*coeff[0]*fEdge[5]-3.705794133009818*coeff[2]*fSkin[2]+3.31456303681194*coeff[0]*fSkin[2]+3.705794133009818*coeff[2]*fEdge[2]-3.31456303681194*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-5.303300858899106*coeff[2]*fSkin[21])-1.797214641813825*coeff[1]*fSkin[21]+4.743416490252569*coeff[0]*fSkin[21]+5.303300858899106*coeff[2]*fEdge[21]-1.797214641813825*coeff[1]*fEdge[21]-4.743416490252569*coeff[0]*fEdge[21]-6.418623720763661*coeff[2]*fSkin[6]-0.9943689110435817*coeff[1]*fSkin[6]+5.74099158464807*coeff[0]*fSkin[6]-6.418623720763661*coeff[2]*fEdge[6]+0.9943689110435817*coeff[1]*fEdge[6]+5.74099158464807*coeff[0]*fEdge[6]-3.705794133009818*coeff[2]*fSkin[3]+3.31456303681194*coeff[0]*fSkin[3]+3.705794133009818*coeff[2]*fEdge[3]-3.31456303681194*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-5.303300858899106*coeff[2]*fSkin[25])-1.797214641813825*coeff[1]*fSkin[25]+4.743416490252569*coeff[0]*fSkin[25]+5.303300858899106*coeff[2]*fEdge[25]-1.797214641813825*coeff[1]*fEdge[25]-4.743416490252569*coeff[0]*fEdge[25]-6.418623720763661*coeff[2]*fSkin[8]-0.9943689110435817*coeff[1]*fSkin[8]+5.74099158464807*coeff[0]*fSkin[8]-6.418623720763661*coeff[2]*fEdge[8]+0.9943689110435817*coeff[1]*fEdge[8]+5.74099158464807*coeff[0]*fEdge[8]-3.705794133009818*coeff[2]*fSkin[4]+3.31456303681194*coeff[0]*fSkin[4]+3.705794133009818*coeff[2]*fEdge[4]-3.31456303681194*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-11.19493359006374*coeff[2]*fSkin[19])-3.112867071728246*coeff[1]*fSkin[19]+10.01305300439132*coeff[0]*fSkin[19]+7.176239480810093*coeff[2]*fEdge[19]-3.112867071728246*coeff[1]*fEdge[19]-6.418623720763666*coeff[0]*fEdge[19]-12.2291206389324*coeff[2]*fSkin[5]-1.72229747539442*coeff[1]*fSkin[5]+10.9380580214794*coeff[0]*fSkin[5]-10.00564415912651*coeff[2]*fEdge[5]+1.72229747539442*coeff[1]*fEdge[5]+8.949320199392238*coeff[0]*fEdge[5]-6.418623720763661*coeff[2]*fSkin[2]+5.74099158464807*coeff[0]*fSkin[2]+6.418623720763661*coeff[2]*fEdge[2]-5.74099158464807*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-11.19493359006374*coeff[2]*fSkin[21])-3.112867071728246*coeff[1]*fSkin[21]+10.01305300439132*coeff[0]*fSkin[21]+7.176239480810093*coeff[2]*fEdge[21]-3.112867071728246*coeff[1]*fEdge[21]-6.418623720763666*coeff[0]*fEdge[21]-12.2291206389324*coeff[2]*fSkin[6]-1.72229747539442*coeff[1]*fSkin[6]+10.9380580214794*coeff[0]*fSkin[6]-10.00564415912651*coeff[2]*fEdge[6]+1.72229747539442*coeff[1]*fEdge[6]+8.949320199392238*coeff[0]*fEdge[6]-6.418623720763661*coeff[2]*fSkin[3]+5.74099158464807*coeff[0]*fSkin[3]+6.418623720763661*coeff[2]*fEdge[3]-5.74099158464807*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-5.303300858899105*coeff[2]*fSkin[32])-1.797214641813825*coeff[1]*fSkin[32]+4.743416490252569*coeff[0]*fSkin[32]+5.303300858899105*coeff[2]*fEdge[32]-1.797214641813825*coeff[1]*fEdge[32]-4.743416490252569*coeff[0]*fEdge[32]-6.418623720763661*coeff[2]*fSkin[15]-0.9943689110435817*coeff[1]*fSkin[15]+5.74099158464807*coeff[0]*fSkin[15]-6.418623720763661*coeff[2]*fEdge[15]+0.9943689110435817*coeff[1]*fEdge[15]+5.74099158464807*coeff[0]*fEdge[15]-3.705794133009818*coeff[2]*fSkin[7]+3.31456303681194*coeff[0]*fSkin[7]+3.705794133009818*coeff[2]*fEdge[7]-3.31456303681194*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-11.19493359006374*coeff[2]*fSkin[25])-3.112867071728246*coeff[1]*fSkin[25]+10.01305300439132*coeff[0]*fSkin[25]+7.176239480810093*coeff[2]*fEdge[25]-3.112867071728246*coeff[1]*fEdge[25]-6.418623720763666*coeff[0]*fEdge[25]-12.2291206389324*coeff[2]*fSkin[8]-1.72229747539442*coeff[1]*fSkin[8]+10.9380580214794*coeff[0]*fSkin[8]-10.00564415912651*coeff[2]*fEdge[8]+1.72229747539442*coeff[1]*fEdge[8]+8.949320199392238*coeff[0]*fEdge[8]-6.418623720763661*coeff[2]*fSkin[4]+5.74099158464807*coeff[0]*fSkin[4]+6.418623720763661*coeff[2]*fEdge[4]-5.74099158464807*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-5.303300858899105*coeff[2]*fSkin[35])-1.797214641813825*coeff[1]*fSkin[35]+4.743416490252569*coeff[0]*fSkin[35]+5.303300858899105*coeff[2]*fEdge[35]-1.797214641813825*coeff[1]*fEdge[35]-4.743416490252569*coeff[0]*fEdge[35]-6.418623720763661*coeff[2]*fSkin[16]-0.9943689110435817*coeff[1]*fSkin[16]+5.74099158464807*coeff[0]*fSkin[16]-6.418623720763661*coeff[2]*fEdge[16]+0.9943689110435817*coeff[1]*fEdge[16]+5.74099158464807*coeff[0]*fEdge[16]-3.705794133009818*coeff[2]*fSkin[9]+3.31456303681194*coeff[0]*fSkin[9]+3.705794133009818*coeff[2]*fEdge[9]-3.31456303681194*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-5.303300858899105*coeff[2]*fSkin[37])-1.797214641813825*coeff[1]*fSkin[37]+4.743416490252569*coeff[0]*fSkin[37]+5.303300858899105*coeff[2]*fEdge[37]-1.797214641813825*coeff[1]*fEdge[37]-4.743416490252569*coeff[0]*fEdge[37]-6.418623720763661*coeff[2]*fSkin[17]-0.9943689110435817*coeff[1]*fSkin[17]+5.74099158464807*coeff[0]*fSkin[17]-6.418623720763661*coeff[2]*fEdge[17]+0.9943689110435817*coeff[1]*fEdge[17]+5.74099158464807*coeff[0]*fEdge[17]-3.705794133009818*coeff[2]*fSkin[10]+3.31456303681194*coeff[0]*fSkin[10]+3.705794133009818*coeff[2]*fEdge[10]-3.31456303681194*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-34.53800131965151*coeff[2]*fSkin[11])-11.53939308514262*coeff[1]*fSkin[11]+14.38520357976382*coeff[0]*fSkin[11]+3.409330602369041*coeff[2]*fEdge[11]-0.516689242618324*coeff[1]*fEdge[11]-0.4640388251536773*coeff[0]*fEdge[11]-42.4833377263957*fSkin[1]*coeff[2]-11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]+11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]+3.337684334797107*coeff[0]*fEdge[1]-9.756308055560766*fSkin[0]*coeff[1]+5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = (-6.418623720763661*coeff[2]*fSkin[20])-0.9943689110435818*coeff[1]*fSkin[20]+5.740991584648071*coeff[0]*fSkin[20]-6.418623720763661*coeff[2]*fEdge[20]+0.9943689110435818*coeff[1]*fEdge[20]+5.740991584648071*coeff[0]*fEdge[20]-3.705794133009818*coeff[2]*fSkin[12]+3.31456303681194*coeff[0]*fSkin[12]+3.705794133009818*coeff[2]*fEdge[12]-3.31456303681194*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = (-6.418623720763661*coeff[2]*fSkin[23])-0.9943689110435818*coeff[1]*fSkin[23]+5.740991584648071*coeff[0]*fSkin[23]-6.418623720763661*coeff[2]*fEdge[23]+0.9943689110435818*coeff[1]*fEdge[23]+5.740991584648071*coeff[0]*fEdge[23]-3.705794133009818*coeff[2]*fSkin[13]+3.31456303681194*coeff[0]*fSkin[13]+3.705794133009818*coeff[2]*fEdge[13]-3.31456303681194*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = (-6.418623720763661*coeff[2]*fSkin[28])-0.9943689110435818*coeff[1]*fSkin[28]+5.740991584648071*coeff[0]*fSkin[28]-6.418623720763661*coeff[2]*fEdge[28]+0.9943689110435818*coeff[1]*fEdge[28]+5.740991584648071*coeff[0]*fEdge[28]-3.705794133009818*coeff[2]*fSkin[14]+3.31456303681194*coeff[0]*fSkin[14]+3.705794133009818*coeff[2]*fEdge[14]-3.31456303681194*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-11.19493359006374*coeff[2]*fSkin[32])-3.112867071728247*coeff[1]*fSkin[32]+10.01305300439132*coeff[0]*fSkin[32]+7.176239480810091*coeff[2]*fEdge[32]-3.112867071728247*coeff[1]*fEdge[32]-6.418623720763665*coeff[0]*fEdge[32]-12.2291206389324*coeff[2]*fSkin[15]-1.72229747539442*coeff[1]*fSkin[15]+10.9380580214794*coeff[0]*fSkin[15]-10.00564415912651*coeff[2]*fEdge[15]+1.72229747539442*coeff[1]*fEdge[15]+8.949320199392238*coeff[0]*fEdge[15]-6.418623720763661*coeff[2]*fSkin[7]+5.74099158464807*coeff[0]*fSkin[7]+6.418623720763661*coeff[2]*fEdge[7]-5.74099158464807*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = (-11.19493359006374*coeff[2]*fSkin[35])-3.112867071728247*coeff[1]*fSkin[35]+10.01305300439132*coeff[0]*fSkin[35]+7.176239480810091*coeff[2]*fEdge[35]-3.112867071728247*coeff[1]*fEdge[35]-6.418623720763665*coeff[0]*fEdge[35]-12.2291206389324*coeff[2]*fSkin[16]-1.72229747539442*coeff[1]*fSkin[16]+10.9380580214794*coeff[0]*fSkin[16]-10.00564415912651*coeff[2]*fEdge[16]+1.72229747539442*coeff[1]*fEdge[16]+8.949320199392238*coeff[0]*fEdge[16]-6.418623720763661*coeff[2]*fSkin[9]+5.74099158464807*coeff[0]*fSkin[9]+6.418623720763661*coeff[2]*fEdge[9]-5.74099158464807*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = (-11.19493359006374*coeff[2]*fSkin[37])-3.112867071728247*coeff[1]*fSkin[37]+10.01305300439132*coeff[0]*fSkin[37]+7.176239480810091*coeff[2]*fEdge[37]-3.112867071728247*coeff[1]*fEdge[37]-6.418623720763665*coeff[0]*fEdge[37]-12.2291206389324*coeff[2]*fSkin[17]-1.72229747539442*coeff[1]*fSkin[17]+10.9380580214794*coeff[0]*fSkin[17]-10.00564415912651*coeff[2]*fEdge[17]+1.72229747539442*coeff[1]*fEdge[17]+8.949320199392238*coeff[0]*fEdge[17]-6.418623720763661*coeff[2]*fSkin[10]+5.74099158464807*coeff[0]*fSkin[10]+6.418623720763661*coeff[2]*fEdge[10]-5.74099158464807*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = (-5.303300858899106*coeff[2]*fSkin[44])-1.797214641813825*coeff[1]*fSkin[44]+4.743416490252569*coeff[0]*fSkin[44]+5.303300858899106*coeff[2]*fEdge[44]-1.797214641813825*coeff[1]*fEdge[44]-4.743416490252569*coeff[0]*fEdge[44]-6.418623720763661*coeff[2]*fSkin[31]-0.9943689110435817*coeff[1]*fSkin[31]+5.74099158464807*coeff[0]*fSkin[31]-6.418623720763661*coeff[2]*fEdge[31]+0.9943689110435817*coeff[1]*fEdge[31]+5.74099158464807*coeff[0]*fEdge[31]-3.705794133009818*coeff[2]*fSkin[18]+3.31456303681194*coeff[0]*fSkin[18]+3.705794133009818*coeff[2]*fEdge[18]-3.31456303681194*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-34.53800131965151*coeff[2]*fSkin[19])-11.53939308514262*coeff[1]*fSkin[19]+14.38520357976382*coeff[0]*fSkin[19]+3.409330602369041*coeff[2]*fEdge[19]-0.516689242618324*coeff[1]*fEdge[19]-0.4640388251536773*coeff[0]*fEdge[19]-42.48333772639572*coeff[2]*fSkin[5]-14.89729241469946*coeff[1]*fSkin[5]+11.0400327997135*coeff[0]*fSkin[5]-11.48198316929615*coeff[2]*fEdge[5]-4.669300607592371*coeff[1]*fEdge[5]+3.337684334797105*coeff[0]*fEdge[5]-26.18504799081433*coeff[2]*fSkin[2]-9.756308055560769*coeff[1]*fSkin[2]+2.964635306407854*coeff[0]*fSkin[2]+10.27514541411701*coeff[2]*fEdge[2]+5.648388874272023*coeff[1]*fEdge[2]-2.964635306407854*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = (-12.2291206389324*coeff[2]*fSkin[20])-1.72229747539442*coeff[1]*fSkin[20]+10.9380580214794*coeff[0]*fSkin[20]-10.00564415912651*coeff[2]*fEdge[20]+1.72229747539442*coeff[1]*fEdge[20]+8.949320199392238*coeff[0]*fEdge[20]-6.418623720763661*coeff[2]*fSkin[12]+5.740991584648071*coeff[0]*fSkin[12]+6.418623720763661*coeff[2]*fEdge[12]-5.740991584648071*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = (-34.53800131965151*coeff[2]*fSkin[21])-11.53939308514262*coeff[1]*fSkin[21]+14.38520357976382*coeff[0]*fSkin[21]+3.409330602369041*coeff[2]*fEdge[21]-0.516689242618324*coeff[1]*fEdge[21]-0.4640388251536773*coeff[0]*fEdge[21]-42.48333772639572*coeff[2]*fSkin[6]-14.89729241469946*coeff[1]*fSkin[6]+11.0400327997135*coeff[0]*fSkin[6]-11.48198316929615*coeff[2]*fEdge[6]-4.669300607592371*coeff[1]*fEdge[6]+3.337684334797105*coeff[0]*fEdge[6]-26.18504799081433*coeff[2]*fSkin[3]-9.756308055560769*coeff[1]*fSkin[3]+2.964635306407854*coeff[0]*fSkin[3]+10.27514541411701*coeff[2]*fEdge[3]+5.648388874272023*coeff[1]*fEdge[3]-2.964635306407854*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = (-6.418623720763661*coeff[2]*fSkin[33])-0.9943689110435818*coeff[1]*fSkin[33]+5.740991584648071*coeff[0]*fSkin[33]-6.418623720763661*coeff[2]*fEdge[33]+0.9943689110435818*coeff[1]*fEdge[33]+5.740991584648071*coeff[0]*fEdge[33]-3.705794133009818*coeff[2]*fSkin[22]+3.31456303681194*coeff[0]*fSkin[22]+3.705794133009818*coeff[2]*fEdge[22]-3.31456303681194*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-12.2291206389324*coeff[2]*fSkin[23])-1.72229747539442*coeff[1]*fSkin[23]+10.9380580214794*coeff[0]*fSkin[23]-10.00564415912651*coeff[2]*fEdge[23]+1.72229747539442*coeff[1]*fEdge[23]+8.949320199392238*coeff[0]*fEdge[23]-6.418623720763661*coeff[2]*fSkin[13]+5.740991584648071*coeff[0]*fSkin[13]+6.418623720763661*coeff[2]*fEdge[13]-5.740991584648071*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = (-6.418623720763661*coeff[2]*fSkin[34])-0.9943689110435818*coeff[1]*fSkin[34]+5.740991584648071*coeff[0]*fSkin[34]-6.418623720763661*coeff[2]*fEdge[34]+0.9943689110435818*coeff[1]*fEdge[34]+5.740991584648071*coeff[0]*fEdge[34]-3.705794133009818*coeff[2]*fSkin[24]+3.31456303681194*coeff[0]*fSkin[24]+3.705794133009818*coeff[2]*fEdge[24]-3.31456303681194*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = (-34.53800131965151*coeff[2]*fSkin[25])-11.53939308514262*coeff[1]*fSkin[25]+14.38520357976382*coeff[0]*fSkin[25]+3.409330602369041*coeff[2]*fEdge[25]-0.516689242618324*coeff[1]*fEdge[25]-0.4640388251536773*coeff[0]*fEdge[25]-42.48333772639572*coeff[2]*fSkin[8]-14.89729241469946*coeff[1]*fSkin[8]+11.0400327997135*coeff[0]*fSkin[8]-11.48198316929615*coeff[2]*fEdge[8]-4.669300607592371*coeff[1]*fEdge[8]+3.337684334797105*coeff[0]*fEdge[8]-26.18504799081433*coeff[2]*fSkin[4]-9.756308055560769*coeff[1]*fSkin[4]+2.964635306407854*coeff[0]*fSkin[4]+10.27514541411701*coeff[2]*fEdge[4]+5.648388874272023*coeff[1]*fEdge[4]-2.964635306407854*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = (-6.418623720763661*coeff[2]*fSkin[36])-0.9943689110435818*coeff[1]*fSkin[36]+5.740991584648071*coeff[0]*fSkin[36]-6.418623720763661*coeff[2]*fEdge[36]+0.9943689110435818*coeff[1]*fEdge[36]+5.740991584648071*coeff[0]*fEdge[36]-3.705794133009818*coeff[2]*fSkin[26]+3.31456303681194*coeff[0]*fSkin[26]+3.705794133009818*coeff[2]*fEdge[26]-3.31456303681194*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = (-6.418623720763661*coeff[2]*fSkin[39])-0.9943689110435818*coeff[1]*fSkin[39]+5.740991584648071*coeff[0]*fSkin[39]-6.418623720763661*coeff[2]*fEdge[39]+0.9943689110435818*coeff[1]*fEdge[39]+5.740991584648071*coeff[0]*fEdge[39]-3.705794133009818*coeff[2]*fSkin[27]+3.31456303681194*coeff[0]*fSkin[27]+3.705794133009818*coeff[2]*fEdge[27]-3.31456303681194*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = (-12.2291206389324*coeff[2]*fSkin[28])-1.72229747539442*coeff[1]*fSkin[28]+10.9380580214794*coeff[0]*fSkin[28]-10.00564415912651*coeff[2]*fEdge[28]+1.72229747539442*coeff[1]*fEdge[28]+8.949320199392238*coeff[0]*fEdge[28]-6.418623720763661*coeff[2]*fSkin[14]+5.740991584648071*coeff[0]*fSkin[14]+6.418623720763661*coeff[2]*fEdge[14]-5.740991584648071*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = (-6.418623720763661*coeff[2]*fSkin[41])-0.9943689110435818*coeff[1]*fSkin[41]+5.740991584648071*coeff[0]*fSkin[41]-6.418623720763661*coeff[2]*fEdge[41]+0.9943689110435818*coeff[1]*fEdge[41]+5.740991584648071*coeff[0]*fEdge[41]-3.705794133009818*coeff[2]*fSkin[29]+3.31456303681194*coeff[0]*fSkin[29]+3.705794133009818*coeff[2]*fEdge[29]-3.31456303681194*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = (-6.418623720763661*coeff[2]*fSkin[42])-0.9943689110435818*coeff[1]*fSkin[42]+5.740991584648071*coeff[0]*fSkin[42]-6.418623720763661*coeff[2]*fEdge[42]+0.9943689110435818*coeff[1]*fEdge[42]+5.740991584648071*coeff[0]*fEdge[42]-3.705794133009818*coeff[2]*fSkin[30]+3.31456303681194*coeff[0]*fSkin[30]+3.705794133009818*coeff[2]*fEdge[30]-3.31456303681194*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = (-11.19493359006374*coeff[2]*fSkin[44])-3.112867071728246*coeff[1]*fSkin[44]+10.01305300439132*coeff[0]*fSkin[44]+7.176239480810093*coeff[2]*fEdge[44]-3.112867071728246*coeff[1]*fEdge[44]-6.418623720763666*coeff[0]*fEdge[44]-12.2291206389324*coeff[2]*fSkin[31]-1.72229747539442*coeff[1]*fSkin[31]+10.9380580214794*coeff[0]*fSkin[31]-10.00564415912651*coeff[2]*fEdge[31]+1.72229747539442*coeff[1]*fEdge[31]+8.949320199392238*coeff[0]*fEdge[31]-6.418623720763661*coeff[2]*fSkin[18]+5.74099158464807*coeff[0]*fSkin[18]+6.418623720763661*coeff[2]*fEdge[18]-5.74099158464807*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = (-34.53800131965151*coeff[2]*fSkin[32])-11.53939308514262*coeff[1]*fSkin[32]+14.38520357976382*coeff[0]*fSkin[32]+3.409330602369041*coeff[2]*fEdge[32]-0.516689242618324*coeff[1]*fEdge[32]-0.4640388251536773*coeff[0]*fEdge[32]-42.4833377263957*coeff[2]*fSkin[15]-14.89729241469947*coeff[1]*fSkin[15]+11.0400327997135*coeff[0]*fSkin[15]-11.48198316929614*coeff[2]*fEdge[15]-4.669300607592371*coeff[1]*fEdge[15]+3.337684334797107*coeff[0]*fEdge[15]-26.18504799081432*coeff[2]*fSkin[7]-9.756308055560766*coeff[1]*fSkin[7]+2.964635306407856*coeff[0]*fSkin[7]+10.27514541411701*coeff[2]*fEdge[7]+5.648388874272021*coeff[1]*fEdge[7]-2.964635306407856*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = (-12.2291206389324*coeff[2]*fSkin[33])-1.72229747539442*coeff[1]*fSkin[33]+10.9380580214794*coeff[0]*fSkin[33]-10.00564415912651*coeff[2]*fEdge[33]+1.72229747539442*coeff[1]*fEdge[33]+8.949320199392238*coeff[0]*fEdge[33]-6.418623720763661*coeff[2]*fSkin[22]+5.740991584648071*coeff[0]*fSkin[22]+6.418623720763661*coeff[2]*fEdge[22]-5.740991584648071*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = (-12.2291206389324*coeff[2]*fSkin[34])-1.72229747539442*coeff[1]*fSkin[34]+10.9380580214794*coeff[0]*fSkin[34]-10.00564415912651*coeff[2]*fEdge[34]+1.72229747539442*coeff[1]*fEdge[34]+8.949320199392238*coeff[0]*fEdge[34]-6.418623720763661*coeff[2]*fSkin[24]+5.740991584648071*coeff[0]*fSkin[24]+6.418623720763661*coeff[2]*fEdge[24]-5.740991584648071*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = (-34.53800131965151*coeff[2]*fSkin[35])-11.53939308514262*coeff[1]*fSkin[35]+14.38520357976382*coeff[0]*fSkin[35]+3.409330602369041*coeff[2]*fEdge[35]-0.516689242618324*coeff[1]*fEdge[35]-0.4640388251536773*coeff[0]*fEdge[35]-42.4833377263957*coeff[2]*fSkin[16]-14.89729241469947*coeff[1]*fSkin[16]+11.0400327997135*coeff[0]*fSkin[16]-11.48198316929614*coeff[2]*fEdge[16]-4.669300607592371*coeff[1]*fEdge[16]+3.337684334797107*coeff[0]*fEdge[16]-26.18504799081432*coeff[2]*fSkin[9]-9.756308055560766*coeff[1]*fSkin[9]+2.964635306407856*coeff[0]*fSkin[9]+10.27514541411701*coeff[2]*fEdge[9]+5.648388874272021*coeff[1]*fEdge[9]-2.964635306407856*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = (-12.2291206389324*coeff[2]*fSkin[36])-1.72229747539442*coeff[1]*fSkin[36]+10.9380580214794*coeff[0]*fSkin[36]-10.00564415912651*coeff[2]*fEdge[36]+1.72229747539442*coeff[1]*fEdge[36]+8.949320199392238*coeff[0]*fEdge[36]-6.418623720763661*coeff[2]*fSkin[26]+5.740991584648071*coeff[0]*fSkin[26]+6.418623720763661*coeff[2]*fEdge[26]-5.740991584648071*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = (-34.53800131965151*coeff[2]*fSkin[37])-11.53939308514262*coeff[1]*fSkin[37]+14.38520357976382*coeff[0]*fSkin[37]+3.409330602369041*coeff[2]*fEdge[37]-0.516689242618324*coeff[1]*fEdge[37]-0.4640388251536773*coeff[0]*fEdge[37]-42.4833377263957*coeff[2]*fSkin[17]-14.89729241469947*coeff[1]*fSkin[17]+11.0400327997135*coeff[0]*fSkin[17]-11.48198316929614*coeff[2]*fEdge[17]-4.669300607592371*coeff[1]*fEdge[17]+3.337684334797107*coeff[0]*fEdge[17]-26.18504799081432*coeff[2]*fSkin[10]-9.756308055560766*coeff[1]*fSkin[10]+2.964635306407856*coeff[0]*fSkin[10]+10.27514541411701*coeff[2]*fEdge[10]+5.648388874272021*coeff[1]*fEdge[10]-2.964635306407856*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = (-6.418623720763661*coeff[2]*fSkin[45])-0.9943689110435818*coeff[1]*fSkin[45]+5.740991584648071*coeff[0]*fSkin[45]-6.418623720763661*coeff[2]*fEdge[45]+0.9943689110435818*coeff[1]*fEdge[45]+5.740991584648071*coeff[0]*fEdge[45]-3.705794133009818*coeff[2]*fSkin[38]+3.31456303681194*coeff[0]*fSkin[38]+3.705794133009818*coeff[2]*fEdge[38]-3.31456303681194*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = (-12.2291206389324*coeff[2]*fSkin[39])-1.72229747539442*coeff[1]*fSkin[39]+10.9380580214794*coeff[0]*fSkin[39]-10.00564415912651*coeff[2]*fEdge[39]+1.72229747539442*coeff[1]*fEdge[39]+8.949320199392238*coeff[0]*fEdge[39]-6.418623720763661*coeff[2]*fSkin[27]+5.740991584648071*coeff[0]*fSkin[27]+6.418623720763661*coeff[2]*fEdge[27]-5.740991584648071*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = (-6.418623720763661*coeff[2]*fSkin[46])-0.9943689110435818*coeff[1]*fSkin[46]+5.740991584648071*coeff[0]*fSkin[46]-6.418623720763661*coeff[2]*fEdge[46]+0.9943689110435818*coeff[1]*fEdge[46]+5.740991584648071*coeff[0]*fEdge[46]-3.705794133009818*coeff[2]*fSkin[40]+3.31456303681194*coeff[0]*fSkin[40]+3.705794133009818*coeff[2]*fEdge[40]-3.31456303681194*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = (-12.2291206389324*coeff[2]*fSkin[41])-1.72229747539442*coeff[1]*fSkin[41]+10.9380580214794*coeff[0]*fSkin[41]-10.00564415912651*coeff[2]*fEdge[41]+1.72229747539442*coeff[1]*fEdge[41]+8.949320199392238*coeff[0]*fEdge[41]-6.418623720763661*coeff[2]*fSkin[29]+5.740991584648071*coeff[0]*fSkin[29]+6.418623720763661*coeff[2]*fEdge[29]-5.740991584648071*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = (-12.2291206389324*coeff[2]*fSkin[42])-1.72229747539442*coeff[1]*fSkin[42]+10.9380580214794*coeff[0]*fSkin[42]-10.00564415912651*coeff[2]*fEdge[42]+1.72229747539442*coeff[1]*fEdge[42]+8.949320199392238*coeff[0]*fEdge[42]-6.418623720763661*coeff[2]*fSkin[30]+5.740991584648071*coeff[0]*fSkin[30]+6.418623720763661*coeff[2]*fEdge[30]-5.740991584648071*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = (-6.418623720763661*coeff[2]*fSkin[47])-0.9943689110435818*coeff[1]*fSkin[47]+5.740991584648071*coeff[0]*fSkin[47]-6.418623720763661*coeff[2]*fEdge[47]+0.9943689110435818*coeff[1]*fEdge[47]+5.740991584648071*coeff[0]*fEdge[47]-3.705794133009818*coeff[2]*fSkin[43]+3.31456303681194*coeff[0]*fSkin[43]+3.705794133009818*coeff[2]*fEdge[43]-3.31456303681194*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = (-34.53800131965151*coeff[2]*fSkin[44])-11.53939308514262*coeff[1]*fSkin[44]+14.38520357976382*coeff[0]*fSkin[44]+3.409330602369041*coeff[2]*fEdge[44]-0.516689242618324*coeff[1]*fEdge[44]-0.4640388251536773*coeff[0]*fEdge[44]-42.48333772639572*coeff[2]*fSkin[31]-14.89729241469946*coeff[1]*fSkin[31]+11.0400327997135*coeff[0]*fSkin[31]-11.48198316929615*coeff[2]*fEdge[31]-4.669300607592371*coeff[1]*fEdge[31]+3.337684334797105*coeff[0]*fEdge[31]-26.18504799081433*coeff[2]*fSkin[18]-9.756308055560769*coeff[1]*fSkin[18]+2.964635306407854*coeff[0]*fSkin[18]+10.27514541411701*coeff[2]*fEdge[18]+5.648388874272023*coeff[1]*fEdge[18]-2.964635306407854*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = (-12.2291206389324*coeff[2]*fSkin[45])-1.72229747539442*coeff[1]*fSkin[45]+10.9380580214794*coeff[0]*fSkin[45]-10.00564415912651*coeff[2]*fEdge[45]+1.72229747539442*coeff[1]*fEdge[45]+8.949320199392238*coeff[0]*fEdge[45]-6.418623720763661*coeff[2]*fSkin[38]+5.740991584648071*coeff[0]*fSkin[38]+6.418623720763661*coeff[2]*fEdge[38]-5.740991584648071*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = (-12.2291206389324*coeff[2]*fSkin[46])-1.72229747539442*coeff[1]*fSkin[46]+10.9380580214794*coeff[0]*fSkin[46]-10.00564415912651*coeff[2]*fEdge[46]+1.72229747539442*coeff[1]*fEdge[46]+8.949320199392238*coeff[0]*fEdge[46]-6.418623720763661*coeff[2]*fSkin[40]+5.740991584648071*coeff[0]*fSkin[40]+6.418623720763661*coeff[2]*fEdge[40]-5.740991584648071*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = (-12.2291206389324*coeff[2]*fSkin[47])-1.72229747539442*coeff[1]*fSkin[47]+10.9380580214794*coeff[0]*fSkin[47]-10.00564415912651*coeff[2]*fEdge[47]+1.72229747539442*coeff[1]*fEdge[47]+8.949320199392238*coeff[0]*fEdge[47]-6.418623720763661*coeff[2]*fSkin[43]+5.740991584648071*coeff[0]*fSkin[43]+6.418623720763661*coeff[2]*fEdge[43]-5.740991584648071*coeff[0]*fEdge[43]; 

  boundSurf_incr[0] = 3.59442928362765*coeff[1]*fSkin[11]-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 4.018694109253648*coeff[2]*fSkin[11]-6.225734143456494*coeff[1]*fSkin[11]-3.59442928362765*coeff[0]*fSkin[11]-2.223476479805891*fSkin[1]*coeff[2]+3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 3.594429283627651*coeff[1]*fSkin[19]-1.988737822087164*coeff[1]*fSkin[5]; 
  boundSurf_incr[3] = 3.594429283627651*coeff[1]*fSkin[21]-1.988737822087164*coeff[1]*fSkin[6]; 
  boundSurf_incr[4] = 3.594429283627651*coeff[1]*fSkin[25]-1.988737822087164*coeff[1]*fSkin[8]; 
  boundSurf_incr[5] = 4.018694109253649*coeff[2]*fSkin[19]-6.225734143456493*coeff[1]*fSkin[19]-3.594429283627651*coeff[0]*fSkin[19]-2.223476479805891*coeff[2]*fSkin[5]+3.444594950788841*coeff[1]*fSkin[5]+1.988737822087164*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 4.018694109253649*coeff[2]*fSkin[21]-6.225734143456493*coeff[1]*fSkin[21]-3.594429283627651*coeff[0]*fSkin[21]-2.223476479805891*coeff[2]*fSkin[6]+3.444594950788841*coeff[1]*fSkin[6]+1.988737822087164*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = 3.59442928362765*coeff[1]*fSkin[32]-1.988737822087164*coeff[1]*fSkin[15]; 
  boundSurf_incr[8] = 4.018694109253649*coeff[2]*fSkin[25]-6.225734143456493*coeff[1]*fSkin[25]-3.594429283627651*coeff[0]*fSkin[25]-2.223476479805891*coeff[2]*fSkin[8]+3.444594950788841*coeff[1]*fSkin[8]+1.988737822087164*coeff[0]*fSkin[8]; 
  boundSurf_incr[9] = 3.59442928362765*coeff[1]*fSkin[35]-1.988737822087164*coeff[1]*fSkin[16]; 
  boundSurf_incr[10] = 3.59442928362765*coeff[1]*fSkin[37]-1.988737822087164*coeff[1]*fSkin[17]; 
  boundSurf_incr[11] = (-31.12867071728247*coeff[2]*fSkin[11])+12.05608232776095*coeff[1]*fSkin[11]+13.92116475461015*coeff[0]*fSkin[11]+31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]-7.702348464916393*coeff[0]*fSkin[1]+4.107919181288745*fSkin[0]*coeff[1]; 
  boundSurf_incr[12] = -1.988737822087164*coeff[1]*fSkin[20]; 
  boundSurf_incr[13] = -1.988737822087164*coeff[1]*fSkin[23]; 
  boundSurf_incr[14] = -1.988737822087164*coeff[1]*fSkin[28]; 
  boundSurf_incr[15] = 4.018694109253648*coeff[2]*fSkin[32]-6.225734143456494*coeff[1]*fSkin[32]-3.59442928362765*coeff[0]*fSkin[32]-2.223476479805891*coeff[2]*fSkin[15]+3.444594950788841*coeff[1]*fSkin[15]+1.988737822087164*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 4.018694109253648*coeff[2]*fSkin[35]-6.225734143456494*coeff[1]*fSkin[35]-3.59442928362765*coeff[0]*fSkin[35]-2.223476479805891*coeff[2]*fSkin[16]+3.444594950788841*coeff[1]*fSkin[16]+1.988737822087164*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 4.018694109253648*coeff[2]*fSkin[37]-6.225734143456494*coeff[1]*fSkin[37]-3.59442928362765*coeff[0]*fSkin[37]-2.223476479805891*coeff[2]*fSkin[17]+3.444594950788841*coeff[1]*fSkin[17]+1.988737822087164*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = 3.594429283627651*coeff[1]*fSkin[44]-1.988737822087164*coeff[1]*fSkin[31]; 
  boundSurf_incr[19] = (-31.12867071728247*coeff[2]*fSkin[19])+12.05608232776095*coeff[1]*fSkin[19]+13.92116475461015*coeff[0]*fSkin[19]+31.00135455709958*coeff[2]*fSkin[5]-10.22799180710709*coeff[1]*fSkin[5]-7.702348464916396*coeff[0]*fSkin[5]-15.90990257669732*coeff[2]*fSkin[2]+4.107919181288745*coeff[1]*fSkin[2]; 
  boundSurf_incr[20] = (-2.223476479805891*coeff[2]*fSkin[20])+3.444594950788841*coeff[1]*fSkin[20]+1.988737822087164*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = (-31.12867071728247*coeff[2]*fSkin[21])+12.05608232776095*coeff[1]*fSkin[21]+13.92116475461015*coeff[0]*fSkin[21]+31.00135455709958*coeff[2]*fSkin[6]-10.22799180710709*coeff[1]*fSkin[6]-7.702348464916396*coeff[0]*fSkin[6]-15.90990257669732*coeff[2]*fSkin[3]+4.107919181288745*coeff[1]*fSkin[3]; 
  boundSurf_incr[22] = -1.988737822087164*coeff[1]*fSkin[33]; 
  boundSurf_incr[23] = (-2.223476479805891*coeff[2]*fSkin[23])+3.444594950788841*coeff[1]*fSkin[23]+1.988737822087164*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = -1.988737822087164*coeff[1]*fSkin[34]; 
  boundSurf_incr[25] = (-31.12867071728247*coeff[2]*fSkin[25])+12.05608232776095*coeff[1]*fSkin[25]+13.92116475461015*coeff[0]*fSkin[25]+31.00135455709958*coeff[2]*fSkin[8]-10.22799180710709*coeff[1]*fSkin[8]-7.702348464916396*coeff[0]*fSkin[8]-15.90990257669732*coeff[2]*fSkin[4]+4.107919181288745*coeff[1]*fSkin[4]; 
  boundSurf_incr[26] = -1.988737822087164*coeff[1]*fSkin[36]; 
  boundSurf_incr[27] = -1.988737822087164*coeff[1]*fSkin[39]; 
  boundSurf_incr[28] = (-2.223476479805891*coeff[2]*fSkin[28])+3.444594950788841*coeff[1]*fSkin[28]+1.988737822087164*coeff[0]*fSkin[28]; 
  boundSurf_incr[29] = -1.988737822087164*coeff[1]*fSkin[41]; 
  boundSurf_incr[30] = -1.988737822087164*coeff[1]*fSkin[42]; 
  boundSurf_incr[31] = 4.018694109253649*coeff[2]*fSkin[44]-6.225734143456493*coeff[1]*fSkin[44]-3.594429283627651*coeff[0]*fSkin[44]-2.223476479805891*coeff[2]*fSkin[31]+3.444594950788841*coeff[1]*fSkin[31]+1.988737822087164*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = (-31.12867071728247*coeff[2]*fSkin[32])+12.05608232776095*coeff[1]*fSkin[32]+13.92116475461015*coeff[0]*fSkin[32]+31.00135455709956*coeff[2]*fSkin[15]-10.2279918071071*coeff[1]*fSkin[15]-7.702348464916393*coeff[0]*fSkin[15]-15.90990257669731*coeff[2]*fSkin[7]+4.107919181288745*coeff[1]*fSkin[7]; 
  boundSurf_incr[33] = (-2.223476479805891*coeff[2]*fSkin[33])+3.444594950788841*coeff[1]*fSkin[33]+1.988737822087164*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = (-2.223476479805891*coeff[2]*fSkin[34])+3.444594950788841*coeff[1]*fSkin[34]+1.988737822087164*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = (-31.12867071728247*coeff[2]*fSkin[35])+12.05608232776095*coeff[1]*fSkin[35]+13.92116475461015*coeff[0]*fSkin[35]+31.00135455709956*coeff[2]*fSkin[16]-10.2279918071071*coeff[1]*fSkin[16]-7.702348464916393*coeff[0]*fSkin[16]-15.90990257669731*coeff[2]*fSkin[9]+4.107919181288745*coeff[1]*fSkin[9]; 
  boundSurf_incr[36] = (-2.223476479805891*coeff[2]*fSkin[36])+3.444594950788841*coeff[1]*fSkin[36]+1.988737822087164*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = (-31.12867071728247*coeff[2]*fSkin[37])+12.05608232776095*coeff[1]*fSkin[37]+13.92116475461015*coeff[0]*fSkin[37]+31.00135455709956*coeff[2]*fSkin[17]-10.2279918071071*coeff[1]*fSkin[17]-7.702348464916393*coeff[0]*fSkin[17]-15.90990257669731*coeff[2]*fSkin[10]+4.107919181288745*coeff[1]*fSkin[10]; 
  boundSurf_incr[38] = -1.988737822087164*coeff[1]*fSkin[45]; 
  boundSurf_incr[39] = (-2.223476479805891*coeff[2]*fSkin[39])+3.444594950788841*coeff[1]*fSkin[39]+1.988737822087164*coeff[0]*fSkin[39]; 
  boundSurf_incr[40] = -1.988737822087164*coeff[1]*fSkin[46]; 
  boundSurf_incr[41] = (-2.223476479805891*coeff[2]*fSkin[41])+3.444594950788841*coeff[1]*fSkin[41]+1.988737822087164*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = (-2.223476479805891*coeff[2]*fSkin[42])+3.444594950788841*coeff[1]*fSkin[42]+1.988737822087164*coeff[0]*fSkin[42]; 
  boundSurf_incr[43] = -1.988737822087164*coeff[1]*fSkin[47]; 
  boundSurf_incr[44] = (-31.12867071728247*coeff[2]*fSkin[44])+12.05608232776095*coeff[1]*fSkin[44]+13.92116475461015*coeff[0]*fSkin[44]+31.00135455709958*coeff[2]*fSkin[31]-10.22799180710709*coeff[1]*fSkin[31]-7.702348464916396*coeff[0]*fSkin[31]-15.90990257669732*coeff[2]*fSkin[18]+4.107919181288745*coeff[1]*fSkin[18]; 
  boundSurf_incr[45] = (-2.223476479805891*coeff[2]*fSkin[45])+3.444594950788841*coeff[1]*fSkin[45]+1.988737822087164*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = (-2.223476479805891*coeff[2]*fSkin[46])+3.444594950788841*coeff[1]*fSkin[46]+1.988737822087164*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = (-2.223476479805891*coeff[2]*fSkin[47])+3.444594950788841*coeff[1]*fSkin[47]+1.988737822087164*coeff[0]*fSkin[47]; 

  } else { 

  edgeSurf_incr[0] = (-5.303300858899105*coeff[2]*fSkin[11])+1.797214641813825*coeff[1]*fSkin[11]+4.743416490252569*coeff[0]*fSkin[11]+5.303300858899105*coeff[2]*fEdge[11]+1.797214641813825*coeff[1]*fEdge[11]-4.743416490252569*coeff[0]*fEdge[11]+6.418623720763661*fSkin[1]*coeff[2]+6.418623720763661*fEdge[1]*coeff[2]-3.705794133009818*fSkin[0]*coeff[2]+3.705794133009818*fEdge[0]*coeff[2]-0.9943689110435817*coeff[1]*fSkin[1]-5.74099158464807*coeff[0]*fSkin[1]+0.9943689110435817*coeff[1]*fEdge[1]-5.74099158464807*coeff[0]*fEdge[1]+3.31456303681194*coeff[0]*fSkin[0]-3.31456303681194*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 11.19493359006374*coeff[2]*fSkin[11]-3.112867071728247*coeff[1]*fSkin[11]-10.01305300439132*coeff[0]*fSkin[11]-7.176239480810091*coeff[2]*fEdge[11]-3.112867071728247*coeff[1]*fEdge[11]+6.418623720763665*coeff[0]*fEdge[11]-12.2291206389324*fSkin[1]*coeff[2]-10.00564415912651*fEdge[1]*coeff[2]+6.418623720763661*fSkin[0]*coeff[2]-6.418623720763661*fEdge[0]*coeff[2]+1.72229747539442*coeff[1]*fSkin[1]+10.9380580214794*coeff[0]*fSkin[1]-1.72229747539442*coeff[1]*fEdge[1]+8.949320199392238*coeff[0]*fEdge[1]-5.74099158464807*coeff[0]*fSkin[0]+5.74099158464807*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-5.303300858899106*coeff[2]*fSkin[19])+1.797214641813825*coeff[1]*fSkin[19]+4.743416490252569*coeff[0]*fSkin[19]+5.303300858899106*coeff[2]*fEdge[19]+1.797214641813825*coeff[1]*fEdge[19]-4.743416490252569*coeff[0]*fEdge[19]+6.418623720763661*coeff[2]*fSkin[5]-0.9943689110435817*coeff[1]*fSkin[5]-5.74099158464807*coeff[0]*fSkin[5]+6.418623720763661*coeff[2]*fEdge[5]+0.9943689110435817*coeff[1]*fEdge[5]-5.74099158464807*coeff[0]*fEdge[5]-3.705794133009818*coeff[2]*fSkin[2]+3.31456303681194*coeff[0]*fSkin[2]+3.705794133009818*coeff[2]*fEdge[2]-3.31456303681194*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-5.303300858899106*coeff[2]*fSkin[21])+1.797214641813825*coeff[1]*fSkin[21]+4.743416490252569*coeff[0]*fSkin[21]+5.303300858899106*coeff[2]*fEdge[21]+1.797214641813825*coeff[1]*fEdge[21]-4.743416490252569*coeff[0]*fEdge[21]+6.418623720763661*coeff[2]*fSkin[6]-0.9943689110435817*coeff[1]*fSkin[6]-5.74099158464807*coeff[0]*fSkin[6]+6.418623720763661*coeff[2]*fEdge[6]+0.9943689110435817*coeff[1]*fEdge[6]-5.74099158464807*coeff[0]*fEdge[6]-3.705794133009818*coeff[2]*fSkin[3]+3.31456303681194*coeff[0]*fSkin[3]+3.705794133009818*coeff[2]*fEdge[3]-3.31456303681194*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-5.303300858899106*coeff[2]*fSkin[25])+1.797214641813825*coeff[1]*fSkin[25]+4.743416490252569*coeff[0]*fSkin[25]+5.303300858899106*coeff[2]*fEdge[25]+1.797214641813825*coeff[1]*fEdge[25]-4.743416490252569*coeff[0]*fEdge[25]+6.418623720763661*coeff[2]*fSkin[8]-0.9943689110435817*coeff[1]*fSkin[8]-5.74099158464807*coeff[0]*fSkin[8]+6.418623720763661*coeff[2]*fEdge[8]+0.9943689110435817*coeff[1]*fEdge[8]-5.74099158464807*coeff[0]*fEdge[8]-3.705794133009818*coeff[2]*fSkin[4]+3.31456303681194*coeff[0]*fSkin[4]+3.705794133009818*coeff[2]*fEdge[4]-3.31456303681194*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 11.19493359006374*coeff[2]*fSkin[19]-3.112867071728246*coeff[1]*fSkin[19]-10.01305300439132*coeff[0]*fSkin[19]-7.176239480810093*coeff[2]*fEdge[19]-3.112867071728246*coeff[1]*fEdge[19]+6.418623720763666*coeff[0]*fEdge[19]-12.2291206389324*coeff[2]*fSkin[5]+1.72229747539442*coeff[1]*fSkin[5]+10.9380580214794*coeff[0]*fSkin[5]-10.00564415912651*coeff[2]*fEdge[5]-1.72229747539442*coeff[1]*fEdge[5]+8.949320199392238*coeff[0]*fEdge[5]+6.418623720763661*coeff[2]*fSkin[2]-5.74099158464807*coeff[0]*fSkin[2]-6.418623720763661*coeff[2]*fEdge[2]+5.74099158464807*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 11.19493359006374*coeff[2]*fSkin[21]-3.112867071728246*coeff[1]*fSkin[21]-10.01305300439132*coeff[0]*fSkin[21]-7.176239480810093*coeff[2]*fEdge[21]-3.112867071728246*coeff[1]*fEdge[21]+6.418623720763666*coeff[0]*fEdge[21]-12.2291206389324*coeff[2]*fSkin[6]+1.72229747539442*coeff[1]*fSkin[6]+10.9380580214794*coeff[0]*fSkin[6]-10.00564415912651*coeff[2]*fEdge[6]-1.72229747539442*coeff[1]*fEdge[6]+8.949320199392238*coeff[0]*fEdge[6]+6.418623720763661*coeff[2]*fSkin[3]-5.74099158464807*coeff[0]*fSkin[3]-6.418623720763661*coeff[2]*fEdge[3]+5.74099158464807*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-5.303300858899105*coeff[2]*fSkin[32])+1.797214641813825*coeff[1]*fSkin[32]+4.743416490252569*coeff[0]*fSkin[32]+5.303300858899105*coeff[2]*fEdge[32]+1.797214641813825*coeff[1]*fEdge[32]-4.743416490252569*coeff[0]*fEdge[32]+6.418623720763661*coeff[2]*fSkin[15]-0.9943689110435817*coeff[1]*fSkin[15]-5.74099158464807*coeff[0]*fSkin[15]+6.418623720763661*coeff[2]*fEdge[15]+0.9943689110435817*coeff[1]*fEdge[15]-5.74099158464807*coeff[0]*fEdge[15]-3.705794133009818*coeff[2]*fSkin[7]+3.31456303681194*coeff[0]*fSkin[7]+3.705794133009818*coeff[2]*fEdge[7]-3.31456303681194*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 11.19493359006374*coeff[2]*fSkin[25]-3.112867071728246*coeff[1]*fSkin[25]-10.01305300439132*coeff[0]*fSkin[25]-7.176239480810093*coeff[2]*fEdge[25]-3.112867071728246*coeff[1]*fEdge[25]+6.418623720763666*coeff[0]*fEdge[25]-12.2291206389324*coeff[2]*fSkin[8]+1.72229747539442*coeff[1]*fSkin[8]+10.9380580214794*coeff[0]*fSkin[8]-10.00564415912651*coeff[2]*fEdge[8]-1.72229747539442*coeff[1]*fEdge[8]+8.949320199392238*coeff[0]*fEdge[8]+6.418623720763661*coeff[2]*fSkin[4]-5.74099158464807*coeff[0]*fSkin[4]-6.418623720763661*coeff[2]*fEdge[4]+5.74099158464807*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-5.303300858899105*coeff[2]*fSkin[35])+1.797214641813825*coeff[1]*fSkin[35]+4.743416490252569*coeff[0]*fSkin[35]+5.303300858899105*coeff[2]*fEdge[35]+1.797214641813825*coeff[1]*fEdge[35]-4.743416490252569*coeff[0]*fEdge[35]+6.418623720763661*coeff[2]*fSkin[16]-0.9943689110435817*coeff[1]*fSkin[16]-5.74099158464807*coeff[0]*fSkin[16]+6.418623720763661*coeff[2]*fEdge[16]+0.9943689110435817*coeff[1]*fEdge[16]-5.74099158464807*coeff[0]*fEdge[16]-3.705794133009818*coeff[2]*fSkin[9]+3.31456303681194*coeff[0]*fSkin[9]+3.705794133009818*coeff[2]*fEdge[9]-3.31456303681194*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-5.303300858899105*coeff[2]*fSkin[37])+1.797214641813825*coeff[1]*fSkin[37]+4.743416490252569*coeff[0]*fSkin[37]+5.303300858899105*coeff[2]*fEdge[37]+1.797214641813825*coeff[1]*fEdge[37]-4.743416490252569*coeff[0]*fEdge[37]+6.418623720763661*coeff[2]*fSkin[17]-0.9943689110435817*coeff[1]*fSkin[17]-5.74099158464807*coeff[0]*fSkin[17]+6.418623720763661*coeff[2]*fEdge[17]+0.9943689110435817*coeff[1]*fEdge[17]-5.74099158464807*coeff[0]*fEdge[17]-3.705794133009818*coeff[2]*fSkin[10]+3.31456303681194*coeff[0]*fSkin[10]+3.705794133009818*coeff[2]*fEdge[10]-3.31456303681194*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-34.53800131965151*coeff[2]*fSkin[11])+11.53939308514262*coeff[1]*fSkin[11]+14.38520357976382*coeff[0]*fSkin[11]+3.409330602369041*coeff[2]*fEdge[11]+0.516689242618324*coeff[1]*fEdge[11]-0.4640388251536773*coeff[0]*fEdge[11]+42.4833377263957*fSkin[1]*coeff[2]+11.48198316929614*fEdge[1]*coeff[2]-26.18504799081432*fSkin[0]*coeff[2]+10.27514541411701*fEdge[0]*coeff[2]-14.89729241469947*coeff[1]*fSkin[1]-11.0400327997135*coeff[0]*fSkin[1]-4.669300607592371*coeff[1]*fEdge[1]-3.337684334797107*coeff[0]*fEdge[1]+9.756308055560766*fSkin[0]*coeff[1]-5.648388874272021*fEdge[0]*coeff[1]+2.964635306407856*coeff[0]*fSkin[0]-2.964635306407856*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 6.418623720763661*coeff[2]*fSkin[20]-0.9943689110435818*coeff[1]*fSkin[20]-5.740991584648071*coeff[0]*fSkin[20]+6.418623720763661*coeff[2]*fEdge[20]+0.9943689110435818*coeff[1]*fEdge[20]-5.740991584648071*coeff[0]*fEdge[20]-3.705794133009818*coeff[2]*fSkin[12]+3.31456303681194*coeff[0]*fSkin[12]+3.705794133009818*coeff[2]*fEdge[12]-3.31456303681194*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 6.418623720763661*coeff[2]*fSkin[23]-0.9943689110435818*coeff[1]*fSkin[23]-5.740991584648071*coeff[0]*fSkin[23]+6.418623720763661*coeff[2]*fEdge[23]+0.9943689110435818*coeff[1]*fEdge[23]-5.740991584648071*coeff[0]*fEdge[23]-3.705794133009818*coeff[2]*fSkin[13]+3.31456303681194*coeff[0]*fSkin[13]+3.705794133009818*coeff[2]*fEdge[13]-3.31456303681194*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 6.418623720763661*coeff[2]*fSkin[28]-0.9943689110435818*coeff[1]*fSkin[28]-5.740991584648071*coeff[0]*fSkin[28]+6.418623720763661*coeff[2]*fEdge[28]+0.9943689110435818*coeff[1]*fEdge[28]-5.740991584648071*coeff[0]*fEdge[28]-3.705794133009818*coeff[2]*fSkin[14]+3.31456303681194*coeff[0]*fSkin[14]+3.705794133009818*coeff[2]*fEdge[14]-3.31456303681194*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 11.19493359006374*coeff[2]*fSkin[32]-3.112867071728247*coeff[1]*fSkin[32]-10.01305300439132*coeff[0]*fSkin[32]-7.176239480810091*coeff[2]*fEdge[32]-3.112867071728247*coeff[1]*fEdge[32]+6.418623720763665*coeff[0]*fEdge[32]-12.2291206389324*coeff[2]*fSkin[15]+1.72229747539442*coeff[1]*fSkin[15]+10.9380580214794*coeff[0]*fSkin[15]-10.00564415912651*coeff[2]*fEdge[15]-1.72229747539442*coeff[1]*fEdge[15]+8.949320199392238*coeff[0]*fEdge[15]+6.418623720763661*coeff[2]*fSkin[7]-5.74099158464807*coeff[0]*fSkin[7]-6.418623720763661*coeff[2]*fEdge[7]+5.74099158464807*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 11.19493359006374*coeff[2]*fSkin[35]-3.112867071728247*coeff[1]*fSkin[35]-10.01305300439132*coeff[0]*fSkin[35]-7.176239480810091*coeff[2]*fEdge[35]-3.112867071728247*coeff[1]*fEdge[35]+6.418623720763665*coeff[0]*fEdge[35]-12.2291206389324*coeff[2]*fSkin[16]+1.72229747539442*coeff[1]*fSkin[16]+10.9380580214794*coeff[0]*fSkin[16]-10.00564415912651*coeff[2]*fEdge[16]-1.72229747539442*coeff[1]*fEdge[16]+8.949320199392238*coeff[0]*fEdge[16]+6.418623720763661*coeff[2]*fSkin[9]-5.74099158464807*coeff[0]*fSkin[9]-6.418623720763661*coeff[2]*fEdge[9]+5.74099158464807*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 11.19493359006374*coeff[2]*fSkin[37]-3.112867071728247*coeff[1]*fSkin[37]-10.01305300439132*coeff[0]*fSkin[37]-7.176239480810091*coeff[2]*fEdge[37]-3.112867071728247*coeff[1]*fEdge[37]+6.418623720763665*coeff[0]*fEdge[37]-12.2291206389324*coeff[2]*fSkin[17]+1.72229747539442*coeff[1]*fSkin[17]+10.9380580214794*coeff[0]*fSkin[17]-10.00564415912651*coeff[2]*fEdge[17]-1.72229747539442*coeff[1]*fEdge[17]+8.949320199392238*coeff[0]*fEdge[17]+6.418623720763661*coeff[2]*fSkin[10]-5.74099158464807*coeff[0]*fSkin[10]-6.418623720763661*coeff[2]*fEdge[10]+5.74099158464807*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = (-5.303300858899106*coeff[2]*fSkin[44])+1.797214641813825*coeff[1]*fSkin[44]+4.743416490252569*coeff[0]*fSkin[44]+5.303300858899106*coeff[2]*fEdge[44]+1.797214641813825*coeff[1]*fEdge[44]-4.743416490252569*coeff[0]*fEdge[44]+6.418623720763661*coeff[2]*fSkin[31]-0.9943689110435817*coeff[1]*fSkin[31]-5.74099158464807*coeff[0]*fSkin[31]+6.418623720763661*coeff[2]*fEdge[31]+0.9943689110435817*coeff[1]*fEdge[31]-5.74099158464807*coeff[0]*fEdge[31]-3.705794133009818*coeff[2]*fSkin[18]+3.31456303681194*coeff[0]*fSkin[18]+3.705794133009818*coeff[2]*fEdge[18]-3.31456303681194*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-34.53800131965151*coeff[2]*fSkin[19])+11.53939308514262*coeff[1]*fSkin[19]+14.38520357976382*coeff[0]*fSkin[19]+3.409330602369041*coeff[2]*fEdge[19]+0.516689242618324*coeff[1]*fEdge[19]-0.4640388251536773*coeff[0]*fEdge[19]+42.48333772639572*coeff[2]*fSkin[5]-14.89729241469946*coeff[1]*fSkin[5]-11.0400327997135*coeff[0]*fSkin[5]+11.48198316929615*coeff[2]*fEdge[5]-4.669300607592371*coeff[1]*fEdge[5]-3.337684334797105*coeff[0]*fEdge[5]-26.18504799081433*coeff[2]*fSkin[2]+9.756308055560769*coeff[1]*fSkin[2]+2.964635306407854*coeff[0]*fSkin[2]+10.27514541411701*coeff[2]*fEdge[2]-5.648388874272023*coeff[1]*fEdge[2]-2.964635306407854*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = (-12.2291206389324*coeff[2]*fSkin[20])+1.72229747539442*coeff[1]*fSkin[20]+10.9380580214794*coeff[0]*fSkin[20]-10.00564415912651*coeff[2]*fEdge[20]-1.72229747539442*coeff[1]*fEdge[20]+8.949320199392238*coeff[0]*fEdge[20]+6.418623720763661*coeff[2]*fSkin[12]-5.740991584648071*coeff[0]*fSkin[12]-6.418623720763661*coeff[2]*fEdge[12]+5.740991584648071*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = (-34.53800131965151*coeff[2]*fSkin[21])+11.53939308514262*coeff[1]*fSkin[21]+14.38520357976382*coeff[0]*fSkin[21]+3.409330602369041*coeff[2]*fEdge[21]+0.516689242618324*coeff[1]*fEdge[21]-0.4640388251536773*coeff[0]*fEdge[21]+42.48333772639572*coeff[2]*fSkin[6]-14.89729241469946*coeff[1]*fSkin[6]-11.0400327997135*coeff[0]*fSkin[6]+11.48198316929615*coeff[2]*fEdge[6]-4.669300607592371*coeff[1]*fEdge[6]-3.337684334797105*coeff[0]*fEdge[6]-26.18504799081433*coeff[2]*fSkin[3]+9.756308055560769*coeff[1]*fSkin[3]+2.964635306407854*coeff[0]*fSkin[3]+10.27514541411701*coeff[2]*fEdge[3]-5.648388874272023*coeff[1]*fEdge[3]-2.964635306407854*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 6.418623720763661*coeff[2]*fSkin[33]-0.9943689110435818*coeff[1]*fSkin[33]-5.740991584648071*coeff[0]*fSkin[33]+6.418623720763661*coeff[2]*fEdge[33]+0.9943689110435818*coeff[1]*fEdge[33]-5.740991584648071*coeff[0]*fEdge[33]-3.705794133009818*coeff[2]*fSkin[22]+3.31456303681194*coeff[0]*fSkin[22]+3.705794133009818*coeff[2]*fEdge[22]-3.31456303681194*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-12.2291206389324*coeff[2]*fSkin[23])+1.72229747539442*coeff[1]*fSkin[23]+10.9380580214794*coeff[0]*fSkin[23]-10.00564415912651*coeff[2]*fEdge[23]-1.72229747539442*coeff[1]*fEdge[23]+8.949320199392238*coeff[0]*fEdge[23]+6.418623720763661*coeff[2]*fSkin[13]-5.740991584648071*coeff[0]*fSkin[13]-6.418623720763661*coeff[2]*fEdge[13]+5.740991584648071*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 6.418623720763661*coeff[2]*fSkin[34]-0.9943689110435818*coeff[1]*fSkin[34]-5.740991584648071*coeff[0]*fSkin[34]+6.418623720763661*coeff[2]*fEdge[34]+0.9943689110435818*coeff[1]*fEdge[34]-5.740991584648071*coeff[0]*fEdge[34]-3.705794133009818*coeff[2]*fSkin[24]+3.31456303681194*coeff[0]*fSkin[24]+3.705794133009818*coeff[2]*fEdge[24]-3.31456303681194*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = (-34.53800131965151*coeff[2]*fSkin[25])+11.53939308514262*coeff[1]*fSkin[25]+14.38520357976382*coeff[0]*fSkin[25]+3.409330602369041*coeff[2]*fEdge[25]+0.516689242618324*coeff[1]*fEdge[25]-0.4640388251536773*coeff[0]*fEdge[25]+42.48333772639572*coeff[2]*fSkin[8]-14.89729241469946*coeff[1]*fSkin[8]-11.0400327997135*coeff[0]*fSkin[8]+11.48198316929615*coeff[2]*fEdge[8]-4.669300607592371*coeff[1]*fEdge[8]-3.337684334797105*coeff[0]*fEdge[8]-26.18504799081433*coeff[2]*fSkin[4]+9.756308055560769*coeff[1]*fSkin[4]+2.964635306407854*coeff[0]*fSkin[4]+10.27514541411701*coeff[2]*fEdge[4]-5.648388874272023*coeff[1]*fEdge[4]-2.964635306407854*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 6.418623720763661*coeff[2]*fSkin[36]-0.9943689110435818*coeff[1]*fSkin[36]-5.740991584648071*coeff[0]*fSkin[36]+6.418623720763661*coeff[2]*fEdge[36]+0.9943689110435818*coeff[1]*fEdge[36]-5.740991584648071*coeff[0]*fEdge[36]-3.705794133009818*coeff[2]*fSkin[26]+3.31456303681194*coeff[0]*fSkin[26]+3.705794133009818*coeff[2]*fEdge[26]-3.31456303681194*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 6.418623720763661*coeff[2]*fSkin[39]-0.9943689110435818*coeff[1]*fSkin[39]-5.740991584648071*coeff[0]*fSkin[39]+6.418623720763661*coeff[2]*fEdge[39]+0.9943689110435818*coeff[1]*fEdge[39]-5.740991584648071*coeff[0]*fEdge[39]-3.705794133009818*coeff[2]*fSkin[27]+3.31456303681194*coeff[0]*fSkin[27]+3.705794133009818*coeff[2]*fEdge[27]-3.31456303681194*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = (-12.2291206389324*coeff[2]*fSkin[28])+1.72229747539442*coeff[1]*fSkin[28]+10.9380580214794*coeff[0]*fSkin[28]-10.00564415912651*coeff[2]*fEdge[28]-1.72229747539442*coeff[1]*fEdge[28]+8.949320199392238*coeff[0]*fEdge[28]+6.418623720763661*coeff[2]*fSkin[14]-5.740991584648071*coeff[0]*fSkin[14]-6.418623720763661*coeff[2]*fEdge[14]+5.740991584648071*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 6.418623720763661*coeff[2]*fSkin[41]-0.9943689110435818*coeff[1]*fSkin[41]-5.740991584648071*coeff[0]*fSkin[41]+6.418623720763661*coeff[2]*fEdge[41]+0.9943689110435818*coeff[1]*fEdge[41]-5.740991584648071*coeff[0]*fEdge[41]-3.705794133009818*coeff[2]*fSkin[29]+3.31456303681194*coeff[0]*fSkin[29]+3.705794133009818*coeff[2]*fEdge[29]-3.31456303681194*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 6.418623720763661*coeff[2]*fSkin[42]-0.9943689110435818*coeff[1]*fSkin[42]-5.740991584648071*coeff[0]*fSkin[42]+6.418623720763661*coeff[2]*fEdge[42]+0.9943689110435818*coeff[1]*fEdge[42]-5.740991584648071*coeff[0]*fEdge[42]-3.705794133009818*coeff[2]*fSkin[30]+3.31456303681194*coeff[0]*fSkin[30]+3.705794133009818*coeff[2]*fEdge[30]-3.31456303681194*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 11.19493359006374*coeff[2]*fSkin[44]-3.112867071728246*coeff[1]*fSkin[44]-10.01305300439132*coeff[0]*fSkin[44]-7.176239480810093*coeff[2]*fEdge[44]-3.112867071728246*coeff[1]*fEdge[44]+6.418623720763666*coeff[0]*fEdge[44]-12.2291206389324*coeff[2]*fSkin[31]+1.72229747539442*coeff[1]*fSkin[31]+10.9380580214794*coeff[0]*fSkin[31]-10.00564415912651*coeff[2]*fEdge[31]-1.72229747539442*coeff[1]*fEdge[31]+8.949320199392238*coeff[0]*fEdge[31]+6.418623720763661*coeff[2]*fSkin[18]-5.74099158464807*coeff[0]*fSkin[18]-6.418623720763661*coeff[2]*fEdge[18]+5.74099158464807*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = (-34.53800131965151*coeff[2]*fSkin[32])+11.53939308514262*coeff[1]*fSkin[32]+14.38520357976382*coeff[0]*fSkin[32]+3.409330602369041*coeff[2]*fEdge[32]+0.516689242618324*coeff[1]*fEdge[32]-0.4640388251536773*coeff[0]*fEdge[32]+42.4833377263957*coeff[2]*fSkin[15]-14.89729241469947*coeff[1]*fSkin[15]-11.0400327997135*coeff[0]*fSkin[15]+11.48198316929614*coeff[2]*fEdge[15]-4.669300607592371*coeff[1]*fEdge[15]-3.337684334797107*coeff[0]*fEdge[15]-26.18504799081432*coeff[2]*fSkin[7]+9.756308055560766*coeff[1]*fSkin[7]+2.964635306407856*coeff[0]*fSkin[7]+10.27514541411701*coeff[2]*fEdge[7]-5.648388874272021*coeff[1]*fEdge[7]-2.964635306407856*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = (-12.2291206389324*coeff[2]*fSkin[33])+1.72229747539442*coeff[1]*fSkin[33]+10.9380580214794*coeff[0]*fSkin[33]-10.00564415912651*coeff[2]*fEdge[33]-1.72229747539442*coeff[1]*fEdge[33]+8.949320199392238*coeff[0]*fEdge[33]+6.418623720763661*coeff[2]*fSkin[22]-5.740991584648071*coeff[0]*fSkin[22]-6.418623720763661*coeff[2]*fEdge[22]+5.740991584648071*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = (-12.2291206389324*coeff[2]*fSkin[34])+1.72229747539442*coeff[1]*fSkin[34]+10.9380580214794*coeff[0]*fSkin[34]-10.00564415912651*coeff[2]*fEdge[34]-1.72229747539442*coeff[1]*fEdge[34]+8.949320199392238*coeff[0]*fEdge[34]+6.418623720763661*coeff[2]*fSkin[24]-5.740991584648071*coeff[0]*fSkin[24]-6.418623720763661*coeff[2]*fEdge[24]+5.740991584648071*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = (-34.53800131965151*coeff[2]*fSkin[35])+11.53939308514262*coeff[1]*fSkin[35]+14.38520357976382*coeff[0]*fSkin[35]+3.409330602369041*coeff[2]*fEdge[35]+0.516689242618324*coeff[1]*fEdge[35]-0.4640388251536773*coeff[0]*fEdge[35]+42.4833377263957*coeff[2]*fSkin[16]-14.89729241469947*coeff[1]*fSkin[16]-11.0400327997135*coeff[0]*fSkin[16]+11.48198316929614*coeff[2]*fEdge[16]-4.669300607592371*coeff[1]*fEdge[16]-3.337684334797107*coeff[0]*fEdge[16]-26.18504799081432*coeff[2]*fSkin[9]+9.756308055560766*coeff[1]*fSkin[9]+2.964635306407856*coeff[0]*fSkin[9]+10.27514541411701*coeff[2]*fEdge[9]-5.648388874272021*coeff[1]*fEdge[9]-2.964635306407856*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = (-12.2291206389324*coeff[2]*fSkin[36])+1.72229747539442*coeff[1]*fSkin[36]+10.9380580214794*coeff[0]*fSkin[36]-10.00564415912651*coeff[2]*fEdge[36]-1.72229747539442*coeff[1]*fEdge[36]+8.949320199392238*coeff[0]*fEdge[36]+6.418623720763661*coeff[2]*fSkin[26]-5.740991584648071*coeff[0]*fSkin[26]-6.418623720763661*coeff[2]*fEdge[26]+5.740991584648071*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = (-34.53800131965151*coeff[2]*fSkin[37])+11.53939308514262*coeff[1]*fSkin[37]+14.38520357976382*coeff[0]*fSkin[37]+3.409330602369041*coeff[2]*fEdge[37]+0.516689242618324*coeff[1]*fEdge[37]-0.4640388251536773*coeff[0]*fEdge[37]+42.4833377263957*coeff[2]*fSkin[17]-14.89729241469947*coeff[1]*fSkin[17]-11.0400327997135*coeff[0]*fSkin[17]+11.48198316929614*coeff[2]*fEdge[17]-4.669300607592371*coeff[1]*fEdge[17]-3.337684334797107*coeff[0]*fEdge[17]-26.18504799081432*coeff[2]*fSkin[10]+9.756308055560766*coeff[1]*fSkin[10]+2.964635306407856*coeff[0]*fSkin[10]+10.27514541411701*coeff[2]*fEdge[10]-5.648388874272021*coeff[1]*fEdge[10]-2.964635306407856*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 6.418623720763661*coeff[2]*fSkin[45]-0.9943689110435818*coeff[1]*fSkin[45]-5.740991584648071*coeff[0]*fSkin[45]+6.418623720763661*coeff[2]*fEdge[45]+0.9943689110435818*coeff[1]*fEdge[45]-5.740991584648071*coeff[0]*fEdge[45]-3.705794133009818*coeff[2]*fSkin[38]+3.31456303681194*coeff[0]*fSkin[38]+3.705794133009818*coeff[2]*fEdge[38]-3.31456303681194*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = (-12.2291206389324*coeff[2]*fSkin[39])+1.72229747539442*coeff[1]*fSkin[39]+10.9380580214794*coeff[0]*fSkin[39]-10.00564415912651*coeff[2]*fEdge[39]-1.72229747539442*coeff[1]*fEdge[39]+8.949320199392238*coeff[0]*fEdge[39]+6.418623720763661*coeff[2]*fSkin[27]-5.740991584648071*coeff[0]*fSkin[27]-6.418623720763661*coeff[2]*fEdge[27]+5.740991584648071*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 6.418623720763661*coeff[2]*fSkin[46]-0.9943689110435818*coeff[1]*fSkin[46]-5.740991584648071*coeff[0]*fSkin[46]+6.418623720763661*coeff[2]*fEdge[46]+0.9943689110435818*coeff[1]*fEdge[46]-5.740991584648071*coeff[0]*fEdge[46]-3.705794133009818*coeff[2]*fSkin[40]+3.31456303681194*coeff[0]*fSkin[40]+3.705794133009818*coeff[2]*fEdge[40]-3.31456303681194*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = (-12.2291206389324*coeff[2]*fSkin[41])+1.72229747539442*coeff[1]*fSkin[41]+10.9380580214794*coeff[0]*fSkin[41]-10.00564415912651*coeff[2]*fEdge[41]-1.72229747539442*coeff[1]*fEdge[41]+8.949320199392238*coeff[0]*fEdge[41]+6.418623720763661*coeff[2]*fSkin[29]-5.740991584648071*coeff[0]*fSkin[29]-6.418623720763661*coeff[2]*fEdge[29]+5.740991584648071*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = (-12.2291206389324*coeff[2]*fSkin[42])+1.72229747539442*coeff[1]*fSkin[42]+10.9380580214794*coeff[0]*fSkin[42]-10.00564415912651*coeff[2]*fEdge[42]-1.72229747539442*coeff[1]*fEdge[42]+8.949320199392238*coeff[0]*fEdge[42]+6.418623720763661*coeff[2]*fSkin[30]-5.740991584648071*coeff[0]*fSkin[30]-6.418623720763661*coeff[2]*fEdge[30]+5.740991584648071*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 6.418623720763661*coeff[2]*fSkin[47]-0.9943689110435818*coeff[1]*fSkin[47]-5.740991584648071*coeff[0]*fSkin[47]+6.418623720763661*coeff[2]*fEdge[47]+0.9943689110435818*coeff[1]*fEdge[47]-5.740991584648071*coeff[0]*fEdge[47]-3.705794133009818*coeff[2]*fSkin[43]+3.31456303681194*coeff[0]*fSkin[43]+3.705794133009818*coeff[2]*fEdge[43]-3.31456303681194*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = (-34.53800131965151*coeff[2]*fSkin[44])+11.53939308514262*coeff[1]*fSkin[44]+14.38520357976382*coeff[0]*fSkin[44]+3.409330602369041*coeff[2]*fEdge[44]+0.516689242618324*coeff[1]*fEdge[44]-0.4640388251536773*coeff[0]*fEdge[44]+42.48333772639572*coeff[2]*fSkin[31]-14.89729241469946*coeff[1]*fSkin[31]-11.0400327997135*coeff[0]*fSkin[31]+11.48198316929615*coeff[2]*fEdge[31]-4.669300607592371*coeff[1]*fEdge[31]-3.337684334797105*coeff[0]*fEdge[31]-26.18504799081433*coeff[2]*fSkin[18]+9.756308055560769*coeff[1]*fSkin[18]+2.964635306407854*coeff[0]*fSkin[18]+10.27514541411701*coeff[2]*fEdge[18]-5.648388874272023*coeff[1]*fEdge[18]-2.964635306407854*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = (-12.2291206389324*coeff[2]*fSkin[45])+1.72229747539442*coeff[1]*fSkin[45]+10.9380580214794*coeff[0]*fSkin[45]-10.00564415912651*coeff[2]*fEdge[45]-1.72229747539442*coeff[1]*fEdge[45]+8.949320199392238*coeff[0]*fEdge[45]+6.418623720763661*coeff[2]*fSkin[38]-5.740991584648071*coeff[0]*fSkin[38]-6.418623720763661*coeff[2]*fEdge[38]+5.740991584648071*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = (-12.2291206389324*coeff[2]*fSkin[46])+1.72229747539442*coeff[1]*fSkin[46]+10.9380580214794*coeff[0]*fSkin[46]-10.00564415912651*coeff[2]*fEdge[46]-1.72229747539442*coeff[1]*fEdge[46]+8.949320199392238*coeff[0]*fEdge[46]+6.418623720763661*coeff[2]*fSkin[40]-5.740991584648071*coeff[0]*fSkin[40]-6.418623720763661*coeff[2]*fEdge[40]+5.740991584648071*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = (-12.2291206389324*coeff[2]*fSkin[47])+1.72229747539442*coeff[1]*fSkin[47]+10.9380580214794*coeff[0]*fSkin[47]-10.00564415912651*coeff[2]*fEdge[47]-1.72229747539442*coeff[1]*fEdge[47]+8.949320199392238*coeff[0]*fEdge[47]+6.418623720763661*coeff[2]*fSkin[43]-5.740991584648071*coeff[0]*fSkin[43]-6.418623720763661*coeff[2]*fEdge[43]+5.740991584648071*coeff[0]*fEdge[43]; 

  boundSurf_incr[0] = (-3.59442928362765*coeff[1]*fSkin[11])-1.988737822087164*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-4.018694109253648*coeff[2]*fSkin[11])-6.225734143456494*coeff[1]*fSkin[11]+3.59442928362765*coeff[0]*fSkin[11]-2.223476479805891*fSkin[1]*coeff[2]-3.444594950788841*coeff[1]*fSkin[1]+1.988737822087164*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-3.594429283627651*coeff[1]*fSkin[19])-1.988737822087164*coeff[1]*fSkin[5]; 
  boundSurf_incr[3] = (-3.594429283627651*coeff[1]*fSkin[21])-1.988737822087164*coeff[1]*fSkin[6]; 
  boundSurf_incr[4] = (-3.594429283627651*coeff[1]*fSkin[25])-1.988737822087164*coeff[1]*fSkin[8]; 
  boundSurf_incr[5] = (-4.018694109253649*coeff[2]*fSkin[19])-6.225734143456493*coeff[1]*fSkin[19]+3.594429283627651*coeff[0]*fSkin[19]-2.223476479805891*coeff[2]*fSkin[5]-3.444594950788841*coeff[1]*fSkin[5]+1.988737822087164*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = (-4.018694109253649*coeff[2]*fSkin[21])-6.225734143456493*coeff[1]*fSkin[21]+3.594429283627651*coeff[0]*fSkin[21]-2.223476479805891*coeff[2]*fSkin[6]-3.444594950788841*coeff[1]*fSkin[6]+1.988737822087164*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = (-3.59442928362765*coeff[1]*fSkin[32])-1.988737822087164*coeff[1]*fSkin[15]; 
  boundSurf_incr[8] = (-4.018694109253649*coeff[2]*fSkin[25])-6.225734143456493*coeff[1]*fSkin[25]+3.594429283627651*coeff[0]*fSkin[25]-2.223476479805891*coeff[2]*fSkin[8]-3.444594950788841*coeff[1]*fSkin[8]+1.988737822087164*coeff[0]*fSkin[8]; 
  boundSurf_incr[9] = (-3.59442928362765*coeff[1]*fSkin[35])-1.988737822087164*coeff[1]*fSkin[16]; 
  boundSurf_incr[10] = (-3.59442928362765*coeff[1]*fSkin[37])-1.988737822087164*coeff[1]*fSkin[17]; 
  boundSurf_incr[11] = (-31.12867071728247*coeff[2]*fSkin[11])-12.05608232776095*coeff[1]*fSkin[11]+13.92116475461015*coeff[0]*fSkin[11]-31.00135455709956*fSkin[1]*coeff[2]-15.90990257669731*fSkin[0]*coeff[2]-10.2279918071071*coeff[1]*fSkin[1]+7.702348464916393*coeff[0]*fSkin[1]-4.107919181288745*fSkin[0]*coeff[1]; 
  boundSurf_incr[12] = -1.988737822087164*coeff[1]*fSkin[20]; 
  boundSurf_incr[13] = -1.988737822087164*coeff[1]*fSkin[23]; 
  boundSurf_incr[14] = -1.988737822087164*coeff[1]*fSkin[28]; 
  boundSurf_incr[15] = (-4.018694109253648*coeff[2]*fSkin[32])-6.225734143456494*coeff[1]*fSkin[32]+3.59442928362765*coeff[0]*fSkin[32]-2.223476479805891*coeff[2]*fSkin[15]-3.444594950788841*coeff[1]*fSkin[15]+1.988737822087164*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = (-4.018694109253648*coeff[2]*fSkin[35])-6.225734143456494*coeff[1]*fSkin[35]+3.59442928362765*coeff[0]*fSkin[35]-2.223476479805891*coeff[2]*fSkin[16]-3.444594950788841*coeff[1]*fSkin[16]+1.988737822087164*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = (-4.018694109253648*coeff[2]*fSkin[37])-6.225734143456494*coeff[1]*fSkin[37]+3.59442928362765*coeff[0]*fSkin[37]-2.223476479805891*coeff[2]*fSkin[17]-3.444594950788841*coeff[1]*fSkin[17]+1.988737822087164*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = (-3.594429283627651*coeff[1]*fSkin[44])-1.988737822087164*coeff[1]*fSkin[31]; 
  boundSurf_incr[19] = (-31.12867071728247*coeff[2]*fSkin[19])-12.05608232776095*coeff[1]*fSkin[19]+13.92116475461015*coeff[0]*fSkin[19]-31.00135455709958*coeff[2]*fSkin[5]-10.22799180710709*coeff[1]*fSkin[5]+7.702348464916396*coeff[0]*fSkin[5]-15.90990257669732*coeff[2]*fSkin[2]-4.107919181288745*coeff[1]*fSkin[2]; 
  boundSurf_incr[20] = (-2.223476479805891*coeff[2]*fSkin[20])-3.444594950788841*coeff[1]*fSkin[20]+1.988737822087164*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = (-31.12867071728247*coeff[2]*fSkin[21])-12.05608232776095*coeff[1]*fSkin[21]+13.92116475461015*coeff[0]*fSkin[21]-31.00135455709958*coeff[2]*fSkin[6]-10.22799180710709*coeff[1]*fSkin[6]+7.702348464916396*coeff[0]*fSkin[6]-15.90990257669732*coeff[2]*fSkin[3]-4.107919181288745*coeff[1]*fSkin[3]; 
  boundSurf_incr[22] = -1.988737822087164*coeff[1]*fSkin[33]; 
  boundSurf_incr[23] = (-2.223476479805891*coeff[2]*fSkin[23])-3.444594950788841*coeff[1]*fSkin[23]+1.988737822087164*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = -1.988737822087164*coeff[1]*fSkin[34]; 
  boundSurf_incr[25] = (-31.12867071728247*coeff[2]*fSkin[25])-12.05608232776095*coeff[1]*fSkin[25]+13.92116475461015*coeff[0]*fSkin[25]-31.00135455709958*coeff[2]*fSkin[8]-10.22799180710709*coeff[1]*fSkin[8]+7.702348464916396*coeff[0]*fSkin[8]-15.90990257669732*coeff[2]*fSkin[4]-4.107919181288745*coeff[1]*fSkin[4]; 
  boundSurf_incr[26] = -1.988737822087164*coeff[1]*fSkin[36]; 
  boundSurf_incr[27] = -1.988737822087164*coeff[1]*fSkin[39]; 
  boundSurf_incr[28] = (-2.223476479805891*coeff[2]*fSkin[28])-3.444594950788841*coeff[1]*fSkin[28]+1.988737822087164*coeff[0]*fSkin[28]; 
  boundSurf_incr[29] = -1.988737822087164*coeff[1]*fSkin[41]; 
  boundSurf_incr[30] = -1.988737822087164*coeff[1]*fSkin[42]; 
  boundSurf_incr[31] = (-4.018694109253649*coeff[2]*fSkin[44])-6.225734143456493*coeff[1]*fSkin[44]+3.594429283627651*coeff[0]*fSkin[44]-2.223476479805891*coeff[2]*fSkin[31]-3.444594950788841*coeff[1]*fSkin[31]+1.988737822087164*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = (-31.12867071728247*coeff[2]*fSkin[32])-12.05608232776095*coeff[1]*fSkin[32]+13.92116475461015*coeff[0]*fSkin[32]-31.00135455709956*coeff[2]*fSkin[15]-10.2279918071071*coeff[1]*fSkin[15]+7.702348464916393*coeff[0]*fSkin[15]-15.90990257669731*coeff[2]*fSkin[7]-4.107919181288745*coeff[1]*fSkin[7]; 
  boundSurf_incr[33] = (-2.223476479805891*coeff[2]*fSkin[33])-3.444594950788841*coeff[1]*fSkin[33]+1.988737822087164*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = (-2.223476479805891*coeff[2]*fSkin[34])-3.444594950788841*coeff[1]*fSkin[34]+1.988737822087164*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = (-31.12867071728247*coeff[2]*fSkin[35])-12.05608232776095*coeff[1]*fSkin[35]+13.92116475461015*coeff[0]*fSkin[35]-31.00135455709956*coeff[2]*fSkin[16]-10.2279918071071*coeff[1]*fSkin[16]+7.702348464916393*coeff[0]*fSkin[16]-15.90990257669731*coeff[2]*fSkin[9]-4.107919181288745*coeff[1]*fSkin[9]; 
  boundSurf_incr[36] = (-2.223476479805891*coeff[2]*fSkin[36])-3.444594950788841*coeff[1]*fSkin[36]+1.988737822087164*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = (-31.12867071728247*coeff[2]*fSkin[37])-12.05608232776095*coeff[1]*fSkin[37]+13.92116475461015*coeff[0]*fSkin[37]-31.00135455709956*coeff[2]*fSkin[17]-10.2279918071071*coeff[1]*fSkin[17]+7.702348464916393*coeff[0]*fSkin[17]-15.90990257669731*coeff[2]*fSkin[10]-4.107919181288745*coeff[1]*fSkin[10]; 
  boundSurf_incr[38] = -1.988737822087164*coeff[1]*fSkin[45]; 
  boundSurf_incr[39] = (-2.223476479805891*coeff[2]*fSkin[39])-3.444594950788841*coeff[1]*fSkin[39]+1.988737822087164*coeff[0]*fSkin[39]; 
  boundSurf_incr[40] = -1.988737822087164*coeff[1]*fSkin[46]; 
  boundSurf_incr[41] = (-2.223476479805891*coeff[2]*fSkin[41])-3.444594950788841*coeff[1]*fSkin[41]+1.988737822087164*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = (-2.223476479805891*coeff[2]*fSkin[42])-3.444594950788841*coeff[1]*fSkin[42]+1.988737822087164*coeff[0]*fSkin[42]; 
  boundSurf_incr[43] = -1.988737822087164*coeff[1]*fSkin[47]; 
  boundSurf_incr[44] = (-31.12867071728247*coeff[2]*fSkin[44])-12.05608232776095*coeff[1]*fSkin[44]+13.92116475461015*coeff[0]*fSkin[44]-31.00135455709958*coeff[2]*fSkin[31]-10.22799180710709*coeff[1]*fSkin[31]+7.702348464916396*coeff[0]*fSkin[31]-15.90990257669732*coeff[2]*fSkin[18]-4.107919181288745*coeff[1]*fSkin[18]; 
  boundSurf_incr[45] = (-2.223476479805891*coeff[2]*fSkin[45])-3.444594950788841*coeff[1]*fSkin[45]+1.988737822087164*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = (-2.223476479805891*coeff[2]*fSkin[46])-3.444594950788841*coeff[1]*fSkin[46]+1.988737822087164*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = (-2.223476479805891*coeff[2]*fSkin[47])-3.444594950788841*coeff[1]*fSkin[47]+1.988737822087164*coeff[0]*fSkin[47]; 

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
  out[20] += -1.0*(vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += -1.0*(vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += -1.0*(vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += -1.0*(vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += -1.0*(vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += -1.0*(vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += -1.0*(vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += -1.0*(vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += -1.0*(vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += -1.0*(vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += -1.0*(vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += -1.0*(vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += -1.0*(vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += -1.0*(vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += -1.0*(vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += -1.0*(vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += -1.0*(vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += -1.0*(vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += -1.0*(vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += -1.0*(vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += -1.0*(vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += -1.0*(vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += -1.0*(vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += -1.0*(vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += -1.0*(vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += -1.0*(vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += -1.0*(vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += -1.0*(vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 

  }

  return 0.;
}

