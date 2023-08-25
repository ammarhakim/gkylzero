#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order6_vlasov_boundary_surfx_1x3v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[11])+35.21807064562169*coeff[0]*fEdge[11]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-70.53065765632411*coeff[0]*fSkin[11])+51.46831774920949*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[19])+35.21807064562168*coeff[0]*fEdge[19]-34.09975027401226*coeff[0]*fSkin[5]-34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[21])+35.21807064562168*coeff[0]*fEdge[21]-34.09975027401226*coeff[0]*fSkin[6]-34.09975027401226*coeff[0]*fEdge[6]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-35.21807064562168*coeff[0]*fSkin[25])+35.21807064562168*coeff[0]*fEdge[25]-34.09975027401226*coeff[0]*fSkin[8]-34.09975027401226*coeff[0]*fEdge[8]-19.6875*coeff[0]*fSkin[4]+19.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-70.53065765632414*coeff[0]*fSkin[19])+51.4683177492095*coeff[0]*fEdge[19]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-70.53065765632414*coeff[0]*fSkin[21])+51.4683177492095*coeff[0]*fEdge[21]-61.5234375*coeff[0]*fSkin[6]-56.6015625*coeff[0]*fEdge[6]-34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-35.21807064562169*coeff[0]*fSkin[32])+35.21807064562169*coeff[0]*fEdge[32]-34.09975027401226*coeff[0]*fSkin[15]-34.09975027401226*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[7]+19.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-70.53065765632414*coeff[0]*fSkin[25])+51.4683177492095*coeff[0]*fEdge[25]-61.5234375*coeff[0]*fSkin[8]-56.6015625*coeff[0]*fEdge[8]-34.09975027401226*coeff[0]*fSkin[4]+34.09975027401226*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-35.21807064562169*coeff[0]*fSkin[35])+35.21807064562169*coeff[0]*fEdge[35]-34.09975027401226*coeff[0]*fSkin[16]-34.09975027401226*coeff[0]*fEdge[16]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-35.21807064562169*coeff[0]*fSkin[37])+35.21807064562169*coeff[0]*fEdge[37]-34.09975027401226*coeff[0]*fSkin[17]-34.09975027401226*coeff[0]*fEdge[17]-19.6875*coeff[0]*fSkin[10]+19.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]-31.316701275974*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = (-34.09975027401227*coeff[0]*fSkin[20])-34.09975027401227*coeff[0]*fEdge[20]-19.6875*coeff[0]*fSkin[12]+19.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = (-34.09975027401227*coeff[0]*fSkin[23])-34.09975027401227*coeff[0]*fEdge[23]-19.6875*coeff[0]*fSkin[13]+19.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = (-34.09975027401227*coeff[0]*fSkin[28])-34.09975027401227*coeff[0]*fEdge[28]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-70.53065765632411*coeff[0]*fSkin[32])+51.46831774920949*coeff[0]*fEdge[32]-61.5234375*coeff[0]*fSkin[15]-56.6015625*coeff[0]*fEdge[15]-34.09975027401226*coeff[0]*fSkin[7]+34.09975027401226*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = (-70.53065765632411*coeff[0]*fSkin[35])+51.46831774920949*coeff[0]*fEdge[35]-61.5234375*coeff[0]*fSkin[16]-56.6015625*coeff[0]*fEdge[16]-34.09975027401226*coeff[0]*fSkin[9]+34.09975027401226*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = (-70.53065765632411*coeff[0]*fSkin[37])+51.46831774920949*coeff[0]*fEdge[37]-61.5234375*coeff[0]*fSkin[17]-56.6015625*coeff[0]*fEdge[17]-34.09975027401226*coeff[0]*fSkin[10]+34.09975027401226*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = (-35.21807064562168*coeff[0]*fSkin[44])+35.21807064562168*coeff[0]*fEdge[44]-34.09975027401226*coeff[0]*fSkin[31]-34.09975027401226*coeff[0]*fEdge[31]-19.6875*coeff[0]*fSkin[18]+19.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-70.6640625*coeff[0]*fSkin[19])-3.1640625*coeff[0]*fEdge[19]-31.31670127597403*coeff[0]*fSkin[5]-12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = (-61.5234375*coeff[0]*fSkin[20])-56.6015625*coeff[0]*fEdge[20]-34.09975027401227*coeff[0]*fSkin[12]+34.09975027401227*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = (-70.6640625*coeff[0]*fSkin[21])-3.1640625*coeff[0]*fEdge[21]-31.31670127597403*coeff[0]*fSkin[6]-12.2543613688594*coeff[0]*fEdge[6]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = (-34.09975027401227*coeff[0]*fSkin[33])-34.09975027401227*coeff[0]*fEdge[33]-19.6875*coeff[0]*fSkin[22]+19.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-61.5234375*coeff[0]*fSkin[23])-56.6015625*coeff[0]*fEdge[23]-34.09975027401227*coeff[0]*fSkin[13]+34.09975027401227*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = (-34.09975027401227*coeff[0]*fSkin[34])-34.09975027401227*coeff[0]*fEdge[34]-19.6875*coeff[0]*fSkin[24]+19.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = (-70.6640625*coeff[0]*fSkin[25])-3.1640625*coeff[0]*fEdge[25]-31.31670127597403*coeff[0]*fSkin[8]-12.2543613688594*coeff[0]*fEdge[8]-12.57788237343632*coeff[0]*fSkin[4]+12.57788237343632*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = (-34.09975027401227*coeff[0]*fSkin[36])-34.09975027401227*coeff[0]*fEdge[36]-19.6875*coeff[0]*fSkin[26]+19.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = (-34.09975027401227*coeff[0]*fSkin[39])-34.09975027401227*coeff[0]*fEdge[39]-19.6875*coeff[0]*fSkin[27]+19.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = (-61.5234375*coeff[0]*fSkin[28])-56.6015625*coeff[0]*fEdge[28]-34.09975027401227*coeff[0]*fSkin[14]+34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = (-34.09975027401227*coeff[0]*fSkin[41])-34.09975027401227*coeff[0]*fEdge[41]-19.6875*coeff[0]*fSkin[29]+19.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = (-34.09975027401227*coeff[0]*fSkin[42])-34.09975027401227*coeff[0]*fEdge[42]-19.6875*coeff[0]*fSkin[30]+19.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = (-70.53065765632414*coeff[0]*fSkin[44])+51.4683177492095*coeff[0]*fEdge[44]-61.5234375*coeff[0]*fSkin[31]-56.6015625*coeff[0]*fEdge[31]-34.09975027401226*coeff[0]*fSkin[18]+34.09975027401226*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = (-70.6640625*coeff[0]*fSkin[32])-3.1640625*coeff[0]*fEdge[32]-31.316701275974*coeff[0]*fSkin[15]-12.2543613688594*coeff[0]*fEdge[15]-12.57788237343632*coeff[0]*fSkin[7]+12.57788237343632*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = (-61.5234375*coeff[0]*fSkin[33])-56.6015625*coeff[0]*fEdge[33]-34.09975027401227*coeff[0]*fSkin[22]+34.09975027401227*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = (-61.5234375*coeff[0]*fSkin[34])-56.6015625*coeff[0]*fEdge[34]-34.09975027401227*coeff[0]*fSkin[24]+34.09975027401227*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = (-70.6640625*coeff[0]*fSkin[35])-3.1640625*coeff[0]*fEdge[35]-31.316701275974*coeff[0]*fSkin[16]-12.2543613688594*coeff[0]*fEdge[16]-12.57788237343632*coeff[0]*fSkin[9]+12.57788237343632*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = (-61.5234375*coeff[0]*fSkin[36])-56.6015625*coeff[0]*fEdge[36]-34.09975027401227*coeff[0]*fSkin[26]+34.09975027401227*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = (-70.6640625*coeff[0]*fSkin[37])-3.1640625*coeff[0]*fEdge[37]-31.316701275974*coeff[0]*fSkin[17]-12.2543613688594*coeff[0]*fEdge[17]-12.57788237343632*coeff[0]*fSkin[10]+12.57788237343632*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = (-34.09975027401227*coeff[0]*fSkin[45])-34.09975027401227*coeff[0]*fEdge[45]-19.6875*coeff[0]*fSkin[38]+19.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = (-61.5234375*coeff[0]*fSkin[39])-56.6015625*coeff[0]*fEdge[39]-34.09975027401227*coeff[0]*fSkin[27]+34.09975027401227*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = (-34.09975027401227*coeff[0]*fSkin[46])-34.09975027401227*coeff[0]*fEdge[46]-19.6875*coeff[0]*fSkin[40]+19.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = (-61.5234375*coeff[0]*fSkin[41])-56.6015625*coeff[0]*fEdge[41]-34.09975027401227*coeff[0]*fSkin[29]+34.09975027401227*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = (-61.5234375*coeff[0]*fSkin[42])-56.6015625*coeff[0]*fEdge[42]-34.09975027401227*coeff[0]*fSkin[30]+34.09975027401227*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = (-34.09975027401227*coeff[0]*fSkin[47])-34.09975027401227*coeff[0]*fEdge[47]-19.6875*coeff[0]*fSkin[43]+19.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = (-70.6640625*coeff[0]*fSkin[44])-3.1640625*coeff[0]*fEdge[44]-31.31670127597403*coeff[0]*fSkin[31]-12.2543613688594*coeff[0]*fEdge[31]-12.57788237343632*coeff[0]*fSkin[18]+12.57788237343632*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = (-61.5234375*coeff[0]*fSkin[45])-56.6015625*coeff[0]*fEdge[45]-34.09975027401227*coeff[0]*fSkin[38]+34.09975027401227*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = (-61.5234375*coeff[0]*fSkin[46])-56.6015625*coeff[0]*fEdge[46]-34.09975027401227*coeff[0]*fSkin[40]+34.09975027401227*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = (-61.5234375*coeff[0]*fSkin[47])-56.6015625*coeff[0]*fEdge[47]-34.09975027401227*coeff[0]*fSkin[43]+34.09975027401227*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = 19.06233990711463*coeff[0]*fSkin[11]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 19.06233990711463*coeff[0]*fSkin[19]-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 19.06233990711463*coeff[0]*fSkin[21]-4.921875*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 19.06233990711463*coeff[0]*fSkin[25]-4.921875*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 19.06233990711463*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[11]; 
  boundSurf_incr[15] = 19.06233990711463*coeff[0]*fSkin[32]-4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 19.06233990711463*coeff[0]*fSkin[35]-4.921875*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 19.06233990711463*coeff[0]*fSkin[37]-4.921875*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = 19.06233990711463*coeff[0]*fSkin[5]-73.828125*coeff[0]*fSkin[19]; 
  boundSurf_incr[20] = -4.921875*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 19.06233990711463*coeff[0]*fSkin[6]-73.828125*coeff[0]*fSkin[21]; 
  boundSurf_incr[23] = -4.921875*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = 19.06233990711463*coeff[0]*fSkin[8]-73.828125*coeff[0]*fSkin[25]; 
  boundSurf_incr[28] = -4.921875*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = 19.06233990711463*coeff[0]*fSkin[44]-4.921875*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 19.06233990711463*coeff[0]*fSkin[15]-73.828125*coeff[0]*fSkin[32]; 
  boundSurf_incr[33] = -4.921875*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = -4.921875*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 19.06233990711463*coeff[0]*fSkin[16]-73.828125*coeff[0]*fSkin[35]; 
  boundSurf_incr[36] = -4.921875*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 19.06233990711463*coeff[0]*fSkin[17]-73.828125*coeff[0]*fSkin[37]; 
  boundSurf_incr[39] = -4.921875*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = -4.921875*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = -4.921875*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = 19.06233990711463*coeff[0]*fSkin[31]-73.828125*coeff[0]*fSkin[44]; 
  boundSurf_incr[45] = -4.921875*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = -4.921875*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = -4.921875*coeff[0]*fSkin[47]; 

  } else { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[11])+35.21807064562169*coeff[0]*fEdge[11]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[11]-51.46831774920949*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[19])+35.21807064562168*coeff[0]*fEdge[19]+34.09975027401226*coeff[0]*fSkin[5]+34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[21])+35.21807064562168*coeff[0]*fEdge[21]+34.09975027401226*coeff[0]*fSkin[6]+34.09975027401226*coeff[0]*fEdge[6]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-35.21807064562168*coeff[0]*fSkin[25])+35.21807064562168*coeff[0]*fEdge[25]+34.09975027401226*coeff[0]*fSkin[8]+34.09975027401226*coeff[0]*fEdge[8]-19.6875*coeff[0]*fSkin[4]+19.6875*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[0]*fSkin[19]-51.4683177492095*coeff[0]*fEdge[19]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 70.53065765632414*coeff[0]*fSkin[21]-51.4683177492095*coeff[0]*fEdge[21]-61.5234375*coeff[0]*fSkin[6]-56.6015625*coeff[0]*fEdge[6]+34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-35.21807064562169*coeff[0]*fSkin[32])+35.21807064562169*coeff[0]*fEdge[32]+34.09975027401226*coeff[0]*fSkin[15]+34.09975027401226*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[7]+19.6875*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 70.53065765632414*coeff[0]*fSkin[25]-51.4683177492095*coeff[0]*fEdge[25]-61.5234375*coeff[0]*fSkin[8]-56.6015625*coeff[0]*fEdge[8]+34.09975027401226*coeff[0]*fSkin[4]-34.09975027401226*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-35.21807064562169*coeff[0]*fSkin[35])+35.21807064562169*coeff[0]*fEdge[35]+34.09975027401226*coeff[0]*fSkin[16]+34.09975027401226*coeff[0]*fEdge[16]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-35.21807064562169*coeff[0]*fSkin[37])+35.21807064562169*coeff[0]*fEdge[37]+34.09975027401226*coeff[0]*fSkin[17]+34.09975027401226*coeff[0]*fEdge[17]-19.6875*coeff[0]*fSkin[10]+19.6875*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]+31.316701275974*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 34.09975027401227*coeff[0]*fSkin[20]+34.09975027401227*coeff[0]*fEdge[20]-19.6875*coeff[0]*fSkin[12]+19.6875*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 34.09975027401227*coeff[0]*fSkin[23]+34.09975027401227*coeff[0]*fEdge[23]-19.6875*coeff[0]*fSkin[13]+19.6875*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[0]*fSkin[28]+34.09975027401227*coeff[0]*fEdge[28]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 70.53065765632411*coeff[0]*fSkin[32]-51.46831774920949*coeff[0]*fEdge[32]-61.5234375*coeff[0]*fSkin[15]-56.6015625*coeff[0]*fEdge[15]+34.09975027401226*coeff[0]*fSkin[7]-34.09975027401226*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 70.53065765632411*coeff[0]*fSkin[35]-51.46831774920949*coeff[0]*fEdge[35]-61.5234375*coeff[0]*fSkin[16]-56.6015625*coeff[0]*fEdge[16]+34.09975027401226*coeff[0]*fSkin[9]-34.09975027401226*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 70.53065765632411*coeff[0]*fSkin[37]-51.46831774920949*coeff[0]*fEdge[37]-61.5234375*coeff[0]*fSkin[17]-56.6015625*coeff[0]*fEdge[17]+34.09975027401226*coeff[0]*fSkin[10]-34.09975027401226*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = (-35.21807064562168*coeff[0]*fSkin[44])+35.21807064562168*coeff[0]*fEdge[44]+34.09975027401226*coeff[0]*fSkin[31]+34.09975027401226*coeff[0]*fEdge[31]-19.6875*coeff[0]*fSkin[18]+19.6875*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-70.6640625*coeff[0]*fSkin[19])-3.1640625*coeff[0]*fEdge[19]+31.31670127597403*coeff[0]*fSkin[5]+12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = (-61.5234375*coeff[0]*fSkin[20])-56.6015625*coeff[0]*fEdge[20]+34.09975027401227*coeff[0]*fSkin[12]-34.09975027401227*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = (-70.6640625*coeff[0]*fSkin[21])-3.1640625*coeff[0]*fEdge[21]+31.31670127597403*coeff[0]*fSkin[6]+12.2543613688594*coeff[0]*fEdge[6]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 34.09975027401227*coeff[0]*fSkin[33]+34.09975027401227*coeff[0]*fEdge[33]-19.6875*coeff[0]*fSkin[22]+19.6875*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-61.5234375*coeff[0]*fSkin[23])-56.6015625*coeff[0]*fEdge[23]+34.09975027401227*coeff[0]*fSkin[13]-34.09975027401227*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 34.09975027401227*coeff[0]*fSkin[34]+34.09975027401227*coeff[0]*fEdge[34]-19.6875*coeff[0]*fSkin[24]+19.6875*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = (-70.6640625*coeff[0]*fSkin[25])-3.1640625*coeff[0]*fEdge[25]+31.31670127597403*coeff[0]*fSkin[8]+12.2543613688594*coeff[0]*fEdge[8]-12.57788237343632*coeff[0]*fSkin[4]+12.57788237343632*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 34.09975027401227*coeff[0]*fSkin[36]+34.09975027401227*coeff[0]*fEdge[36]-19.6875*coeff[0]*fSkin[26]+19.6875*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 34.09975027401227*coeff[0]*fSkin[39]+34.09975027401227*coeff[0]*fEdge[39]-19.6875*coeff[0]*fSkin[27]+19.6875*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = (-61.5234375*coeff[0]*fSkin[28])-56.6015625*coeff[0]*fEdge[28]+34.09975027401227*coeff[0]*fSkin[14]-34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 34.09975027401227*coeff[0]*fSkin[41]+34.09975027401227*coeff[0]*fEdge[41]-19.6875*coeff[0]*fSkin[29]+19.6875*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 34.09975027401227*coeff[0]*fSkin[42]+34.09975027401227*coeff[0]*fEdge[42]-19.6875*coeff[0]*fSkin[30]+19.6875*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 70.53065765632414*coeff[0]*fSkin[44]-51.4683177492095*coeff[0]*fEdge[44]-61.5234375*coeff[0]*fSkin[31]-56.6015625*coeff[0]*fEdge[31]+34.09975027401226*coeff[0]*fSkin[18]-34.09975027401226*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = (-70.6640625*coeff[0]*fSkin[32])-3.1640625*coeff[0]*fEdge[32]+31.316701275974*coeff[0]*fSkin[15]+12.2543613688594*coeff[0]*fEdge[15]-12.57788237343632*coeff[0]*fSkin[7]+12.57788237343632*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = (-61.5234375*coeff[0]*fSkin[33])-56.6015625*coeff[0]*fEdge[33]+34.09975027401227*coeff[0]*fSkin[22]-34.09975027401227*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = (-61.5234375*coeff[0]*fSkin[34])-56.6015625*coeff[0]*fEdge[34]+34.09975027401227*coeff[0]*fSkin[24]-34.09975027401227*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = (-70.6640625*coeff[0]*fSkin[35])-3.1640625*coeff[0]*fEdge[35]+31.316701275974*coeff[0]*fSkin[16]+12.2543613688594*coeff[0]*fEdge[16]-12.57788237343632*coeff[0]*fSkin[9]+12.57788237343632*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = (-61.5234375*coeff[0]*fSkin[36])-56.6015625*coeff[0]*fEdge[36]+34.09975027401227*coeff[0]*fSkin[26]-34.09975027401227*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = (-70.6640625*coeff[0]*fSkin[37])-3.1640625*coeff[0]*fEdge[37]+31.316701275974*coeff[0]*fSkin[17]+12.2543613688594*coeff[0]*fEdge[17]-12.57788237343632*coeff[0]*fSkin[10]+12.57788237343632*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 34.09975027401227*coeff[0]*fSkin[45]+34.09975027401227*coeff[0]*fEdge[45]-19.6875*coeff[0]*fSkin[38]+19.6875*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = (-61.5234375*coeff[0]*fSkin[39])-56.6015625*coeff[0]*fEdge[39]+34.09975027401227*coeff[0]*fSkin[27]-34.09975027401227*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 34.09975027401227*coeff[0]*fSkin[46]+34.09975027401227*coeff[0]*fEdge[46]-19.6875*coeff[0]*fSkin[40]+19.6875*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = (-61.5234375*coeff[0]*fSkin[41])-56.6015625*coeff[0]*fEdge[41]+34.09975027401227*coeff[0]*fSkin[29]-34.09975027401227*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = (-61.5234375*coeff[0]*fSkin[42])-56.6015625*coeff[0]*fEdge[42]+34.09975027401227*coeff[0]*fSkin[30]-34.09975027401227*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 34.09975027401227*coeff[0]*fSkin[47]+34.09975027401227*coeff[0]*fEdge[47]-19.6875*coeff[0]*fSkin[43]+19.6875*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = (-70.6640625*coeff[0]*fSkin[44])-3.1640625*coeff[0]*fEdge[44]+31.31670127597403*coeff[0]*fSkin[31]+12.2543613688594*coeff[0]*fEdge[31]-12.57788237343632*coeff[0]*fSkin[18]+12.57788237343632*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = (-61.5234375*coeff[0]*fSkin[45])-56.6015625*coeff[0]*fEdge[45]+34.09975027401227*coeff[0]*fSkin[38]-34.09975027401227*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = (-61.5234375*coeff[0]*fSkin[46])-56.6015625*coeff[0]*fEdge[46]+34.09975027401227*coeff[0]*fSkin[40]-34.09975027401227*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = (-61.5234375*coeff[0]*fSkin[47])-56.6015625*coeff[0]*fEdge[47]+34.09975027401227*coeff[0]*fSkin[43]-34.09975027401227*coeff[0]*fEdge[43]; 

  boundSurf_incr[1] = (-19.06233990711463*coeff[0]*fSkin[11])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = (-19.06233990711463*coeff[0]*fSkin[19])-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = (-19.06233990711463*coeff[0]*fSkin[21])-4.921875*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = (-19.06233990711463*coeff[0]*fSkin[25])-4.921875*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = (-73.828125*coeff[0]*fSkin[11])-19.06233990711463*coeff[0]*fSkin[1]; 
  boundSurf_incr[15] = (-19.06233990711463*coeff[0]*fSkin[32])-4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = (-19.06233990711463*coeff[0]*fSkin[35])-4.921875*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = (-19.06233990711463*coeff[0]*fSkin[37])-4.921875*coeff[0]*fSkin[17]; 
  boundSurf_incr[19] = (-73.828125*coeff[0]*fSkin[19])-19.06233990711463*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = -4.921875*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = (-73.828125*coeff[0]*fSkin[21])-19.06233990711463*coeff[0]*fSkin[6]; 
  boundSurf_incr[23] = -4.921875*coeff[0]*fSkin[23]; 
  boundSurf_incr[25] = (-73.828125*coeff[0]*fSkin[25])-19.06233990711463*coeff[0]*fSkin[8]; 
  boundSurf_incr[28] = -4.921875*coeff[0]*fSkin[28]; 
  boundSurf_incr[31] = (-19.06233990711463*coeff[0]*fSkin[44])-4.921875*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = (-73.828125*coeff[0]*fSkin[32])-19.06233990711463*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = -4.921875*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = -4.921875*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = (-73.828125*coeff[0]*fSkin[35])-19.06233990711463*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = -4.921875*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = (-73.828125*coeff[0]*fSkin[37])-19.06233990711463*coeff[0]*fSkin[17]; 
  boundSurf_incr[39] = -4.921875*coeff[0]*fSkin[39]; 
  boundSurf_incr[41] = -4.921875*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = -4.921875*coeff[0]*fSkin[42]; 
  boundSurf_incr[44] = (-73.828125*coeff[0]*fSkin[44])-19.06233990711463*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = -4.921875*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = -4.921875*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = -4.921875*coeff[0]*fSkin[47]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order6_vlasov_boundary_surfx_1x3v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[48] = {0.0}; 

  double edgeSurf_incr[48] = {0.0}; 
  double boundSurf_incr[48] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[11]+13.47910981360368*coeff[1]*fSkin[11]-24.90293657382598*coeff[0]*fSkin[11]-59.66213466261492*coeff[2]*fEdge[11]+13.47910981360368*coeff[1]*fEdge[11]+24.90293657382598*coeff[0]*fEdge[11]+65.46996195178933*fSkin[1]*coeff[2]+65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]-24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]-24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 110.8728999785158*coeff[2]*fSkin[11]+9.116253567204142*coeff[1]*fSkin[11]-49.87270631033365*coeff[0]*fSkin[11]-95.80279706881466*coeff[2]*fEdge[11]+37.57675250871956*coeff[1]*fEdge[11]+36.39359649672996*coeff[0]*fEdge[11]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]+65.46996195178933*fSkin[0]*coeff[2]-65.46996195178933*fEdge[0]*coeff[2]-11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]-23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]-24.11216465552189*coeff[0]*fSkin[0]+24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 59.66213466261492*coeff[2]*fSkin[19]+13.47910981360369*coeff[1]*fSkin[19]-24.90293657382598*coeff[0]*fSkin[19]-59.66213466261492*coeff[2]*fEdge[19]+13.47910981360369*coeff[1]*fEdge[19]+24.90293657382598*coeff[0]*fEdge[19]+65.46996195178933*coeff[2]*fSkin[5]+3.480291188652536*coeff[1]*fSkin[5]-24.11216465552189*coeff[0]*fSkin[5]+65.46996195178933*coeff[2]*fEdge[5]-3.480291188652536*coeff[1]*fEdge[5]-24.11216465552189*coeff[0]*fEdge[5]+37.79910015670014*coeff[2]*fSkin[2]-13.92116475461015*coeff[0]*fSkin[2]-37.79910015670014*coeff[2]*fEdge[2]+13.92116475461015*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 59.66213466261492*coeff[2]*fSkin[21]+13.47910981360369*coeff[1]*fSkin[21]-24.90293657382598*coeff[0]*fSkin[21]-59.66213466261492*coeff[2]*fEdge[21]+13.47910981360369*coeff[1]*fEdge[21]+24.90293657382598*coeff[0]*fEdge[21]+65.46996195178933*coeff[2]*fSkin[6]+3.480291188652536*coeff[1]*fSkin[6]-24.11216465552189*coeff[0]*fSkin[6]+65.46996195178933*coeff[2]*fEdge[6]-3.480291188652536*coeff[1]*fEdge[6]-24.11216465552189*coeff[0]*fEdge[6]+37.79910015670014*coeff[2]*fSkin[3]-13.92116475461015*coeff[0]*fSkin[3]-37.79910015670014*coeff[2]*fEdge[3]+13.92116475461015*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 59.66213466261492*coeff[2]*fSkin[25]+13.47910981360369*coeff[1]*fSkin[25]-24.90293657382598*coeff[0]*fSkin[25]-59.66213466261492*coeff[2]*fEdge[25]+13.47910981360369*coeff[1]*fEdge[25]+24.90293657382598*coeff[0]*fEdge[25]+65.46996195178933*coeff[2]*fSkin[8]+3.480291188652536*coeff[1]*fSkin[8]-24.11216465552189*coeff[0]*fSkin[8]+65.46996195178933*coeff[2]*fEdge[8]-3.480291188652536*coeff[1]*fEdge[8]-24.11216465552189*coeff[0]*fEdge[8]+37.79910015670014*coeff[2]*fSkin[4]-13.92116475461015*coeff[0]*fSkin[4]-37.79910015670014*coeff[2]*fEdge[4]+13.92116475461015*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 110.8728999785159*coeff[2]*fSkin[19]+9.116253567204131*coeff[1]*fSkin[19]-49.87270631033366*coeff[0]*fSkin[19]-95.80279706881471*coeff[2]*fEdge[19]+37.57675250871954*coeff[1]*fEdge[19]+36.39359649672998*coeff[0]*fEdge[19]+115.3428423899306*coeff[2]*fSkin[5]-11.19493359006374*coeff[1]*fSkin[5]-43.5036398581567*coeff[0]*fSkin[5]+111.4517585502703*coeff[2]*fEdge[5]-23.25101591782468*coeff[1]*fEdge[5]-40.02334866950417*coeff[0]*fEdge[5]+65.46996195178933*coeff[2]*fSkin[2]-9.94368911043582*coeff[1]*fSkin[2]-24.11216465552189*coeff[0]*fSkin[2]-65.46996195178933*coeff[2]*fEdge[2]+9.94368911043582*coeff[1]*fEdge[2]+24.11216465552189*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 110.8728999785159*coeff[2]*fSkin[21]+9.116253567204131*coeff[1]*fSkin[21]-49.87270631033366*coeff[0]*fSkin[21]-95.80279706881471*coeff[2]*fEdge[21]+37.57675250871954*coeff[1]*fEdge[21]+36.39359649672998*coeff[0]*fEdge[21]+115.3428423899306*coeff[2]*fSkin[6]-11.19493359006374*coeff[1]*fSkin[6]-43.5036398581567*coeff[0]*fSkin[6]+111.4517585502703*coeff[2]*fEdge[6]-23.25101591782468*coeff[1]*fEdge[6]-40.02334866950417*coeff[0]*fEdge[6]+65.46996195178933*coeff[2]*fSkin[3]-9.94368911043582*coeff[1]*fSkin[3]-24.11216465552189*coeff[0]*fSkin[3]-65.46996195178933*coeff[2]*fEdge[3]+9.94368911043582*coeff[1]*fEdge[3]+24.11216465552189*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 59.66213466261492*coeff[2]*fSkin[32]+13.47910981360368*coeff[1]*fSkin[32]-24.90293657382598*coeff[0]*fSkin[32]-59.66213466261492*coeff[2]*fEdge[32]+13.47910981360368*coeff[1]*fEdge[32]+24.90293657382598*coeff[0]*fEdge[32]+65.46996195178933*coeff[2]*fSkin[15]+3.480291188652536*coeff[1]*fSkin[15]-24.11216465552189*coeff[0]*fSkin[15]+65.46996195178933*coeff[2]*fEdge[15]-3.480291188652536*coeff[1]*fEdge[15]-24.11216465552189*coeff[0]*fEdge[15]+37.79910015670014*coeff[2]*fSkin[7]-13.92116475461015*coeff[0]*fSkin[7]-37.79910015670014*coeff[2]*fEdge[7]+13.92116475461015*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 110.8728999785159*coeff[2]*fSkin[25]+9.116253567204131*coeff[1]*fSkin[25]-49.87270631033366*coeff[0]*fSkin[25]-95.80279706881471*coeff[2]*fEdge[25]+37.57675250871954*coeff[1]*fEdge[25]+36.39359649672998*coeff[0]*fEdge[25]+115.3428423899306*coeff[2]*fSkin[8]-11.19493359006374*coeff[1]*fSkin[8]-43.5036398581567*coeff[0]*fSkin[8]+111.4517585502703*coeff[2]*fEdge[8]-23.25101591782468*coeff[1]*fEdge[8]-40.02334866950417*coeff[0]*fEdge[8]+65.46996195178933*coeff[2]*fSkin[4]-9.94368911043582*coeff[1]*fSkin[4]-24.11216465552189*coeff[0]*fSkin[4]-65.46996195178933*coeff[2]*fEdge[4]+9.94368911043582*coeff[1]*fEdge[4]+24.11216465552189*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 59.66213466261492*coeff[2]*fSkin[35]+13.47910981360368*coeff[1]*fSkin[35]-24.90293657382598*coeff[0]*fSkin[35]-59.66213466261492*coeff[2]*fEdge[35]+13.47910981360368*coeff[1]*fEdge[35]+24.90293657382598*coeff[0]*fEdge[35]+65.46996195178933*coeff[2]*fSkin[16]+3.480291188652536*coeff[1]*fSkin[16]-24.11216465552189*coeff[0]*fSkin[16]+65.46996195178933*coeff[2]*fEdge[16]-3.480291188652536*coeff[1]*fEdge[16]-24.11216465552189*coeff[0]*fEdge[16]+37.79910015670014*coeff[2]*fSkin[9]-13.92116475461015*coeff[0]*fSkin[9]-37.79910015670014*coeff[2]*fEdge[9]+13.92116475461015*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 59.66213466261492*coeff[2]*fSkin[37]+13.47910981360368*coeff[1]*fSkin[37]-24.90293657382598*coeff[0]*fSkin[37]-59.66213466261492*coeff[2]*fEdge[37]+13.47910981360368*coeff[1]*fEdge[37]+24.90293657382598*coeff[0]*fEdge[37]+65.46996195178933*coeff[2]*fSkin[17]+3.480291188652536*coeff[1]*fSkin[17]-24.11216465552189*coeff[0]*fSkin[17]+65.46996195178933*coeff[2]*fEdge[17]-3.480291188652536*coeff[1]*fEdge[17]-24.11216465552189*coeff[0]*fEdge[17]+37.79910015670014*coeff[2]*fSkin[10]-13.92116475461015*coeff[0]*fSkin[10]-37.79910015670014*coeff[2]*fEdge[10]+13.92116475461015*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 127.0160939089115*coeff[2]*fSkin[11]-24.97331339321913*coeff[1]*fSkin[11]-49.96703777993997*coeff[0]*fSkin[11]-68.64983631400692*coeff[2]*fEdge[11]+85.25372503202384*coeff[1]*fEdge[11]-2.237330049848051*coeff[0]*fEdge[11]+110.8728999785158*fSkin[1]*coeff[2]+95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]-22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]-8.665142023030945*coeff[0]*fEdge[1]-38.51174232458197*fSkin[0]*coeff[1]+38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = 65.46996195178934*coeff[2]*fSkin[20]+3.480291188652535*coeff[1]*fSkin[20]-24.1121646555219*coeff[0]*fSkin[20]+65.46996195178934*coeff[2]*fEdge[20]-3.480291188652535*coeff[1]*fEdge[20]-24.1121646555219*coeff[0]*fEdge[20]+37.79910015670014*coeff[2]*fSkin[12]-13.92116475461015*coeff[0]*fSkin[12]-37.79910015670014*coeff[2]*fEdge[12]+13.92116475461015*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = 65.46996195178934*coeff[2]*fSkin[23]+3.480291188652535*coeff[1]*fSkin[23]-24.1121646555219*coeff[0]*fSkin[23]+65.46996195178934*coeff[2]*fEdge[23]-3.480291188652535*coeff[1]*fEdge[23]-24.1121646555219*coeff[0]*fEdge[23]+37.79910015670014*coeff[2]*fSkin[13]-13.92116475461015*coeff[0]*fSkin[13]-37.79910015670014*coeff[2]*fEdge[13]+13.92116475461015*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = 65.46996195178934*coeff[2]*fSkin[28]+3.480291188652535*coeff[1]*fSkin[28]-24.1121646555219*coeff[0]*fSkin[28]+65.46996195178934*coeff[2]*fEdge[28]-3.480291188652535*coeff[1]*fEdge[28]-24.1121646555219*coeff[0]*fEdge[28]+37.79910015670014*coeff[2]*fSkin[14]-13.92116475461015*coeff[0]*fSkin[14]-37.79910015670014*coeff[2]*fEdge[14]+13.92116475461015*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 110.8728999785158*coeff[2]*fSkin[32]+9.116253567204142*coeff[1]*fSkin[32]-49.87270631033365*coeff[0]*fSkin[32]-95.80279706881466*coeff[2]*fEdge[32]+37.57675250871956*coeff[1]*fEdge[32]+36.39359649672996*coeff[0]*fEdge[32]+115.3428423899306*coeff[2]*fSkin[15]-11.19493359006374*coeff[1]*fSkin[15]-43.5036398581567*coeff[0]*fSkin[15]+111.4517585502703*coeff[2]*fEdge[15]-23.25101591782468*coeff[1]*fEdge[15]-40.02334866950417*coeff[0]*fEdge[15]+65.46996195178933*coeff[2]*fSkin[7]-9.94368911043582*coeff[1]*fSkin[7]-24.11216465552189*coeff[0]*fSkin[7]-65.46996195178933*coeff[2]*fEdge[7]+9.94368911043582*coeff[1]*fEdge[7]+24.11216465552189*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = 110.8728999785158*coeff[2]*fSkin[35]+9.116253567204142*coeff[1]*fSkin[35]-49.87270631033365*coeff[0]*fSkin[35]-95.80279706881466*coeff[2]*fEdge[35]+37.57675250871956*coeff[1]*fEdge[35]+36.39359649672996*coeff[0]*fEdge[35]+115.3428423899306*coeff[2]*fSkin[16]-11.19493359006374*coeff[1]*fSkin[16]-43.5036398581567*coeff[0]*fSkin[16]+111.4517585502703*coeff[2]*fEdge[16]-23.25101591782468*coeff[1]*fEdge[16]-40.02334866950417*coeff[0]*fEdge[16]+65.46996195178933*coeff[2]*fSkin[9]-9.94368911043582*coeff[1]*fSkin[9]-24.11216465552189*coeff[0]*fSkin[9]-65.46996195178933*coeff[2]*fEdge[9]+9.94368911043582*coeff[1]*fEdge[9]+24.11216465552189*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = 110.8728999785158*coeff[2]*fSkin[37]+9.116253567204142*coeff[1]*fSkin[37]-49.87270631033365*coeff[0]*fSkin[37]-95.80279706881466*coeff[2]*fEdge[37]+37.57675250871956*coeff[1]*fEdge[37]+36.39359649672996*coeff[0]*fEdge[37]+115.3428423899306*coeff[2]*fSkin[17]-11.19493359006374*coeff[1]*fSkin[17]-43.5036398581567*coeff[0]*fSkin[17]+111.4517585502703*coeff[2]*fEdge[17]-23.25101591782468*coeff[1]*fEdge[17]-40.02334866950417*coeff[0]*fEdge[17]+65.46996195178933*coeff[2]*fSkin[10]-9.94368911043582*coeff[1]*fSkin[10]-24.11216465552189*coeff[0]*fSkin[10]-65.46996195178933*coeff[2]*fEdge[10]+9.94368911043582*coeff[1]*fEdge[10]+24.11216465552189*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 59.66213466261492*coeff[2]*fSkin[44]+13.47910981360369*coeff[1]*fSkin[44]-24.90293657382598*coeff[0]*fSkin[44]-59.66213466261492*coeff[2]*fEdge[44]+13.47910981360369*coeff[1]*fEdge[44]+24.90293657382598*coeff[0]*fEdge[44]+65.46996195178933*coeff[2]*fSkin[31]+3.480291188652536*coeff[1]*fSkin[31]-24.11216465552189*coeff[0]*fSkin[31]+65.46996195178933*coeff[2]*fEdge[31]-3.480291188652536*coeff[1]*fEdge[31]-24.11216465552189*coeff[0]*fEdge[31]+37.79910015670014*coeff[2]*fSkin[18]-13.92116475461015*coeff[0]*fSkin[18]-37.79910015670014*coeff[2]*fEdge[18]+13.92116475461015*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 127.0160939089115*coeff[2]*fSkin[19]-24.97331339321913*coeff[1]*fSkin[19]-49.96703777993997*coeff[0]*fSkin[19]-68.64983631400692*coeff[2]*fEdge[19]+85.25372503202384*coeff[1]*fEdge[19]-2.237330049848051*coeff[0]*fEdge[19]+110.8728999785159*coeff[2]*fSkin[5]-58.92212671485608*coeff[1]*fSkin[5]-22.14425183663463*coeff[0]*fSkin[5]+95.80279706881468*coeff[2]*fEdge[5]-74.48646207349731*coeff[1]*fEdge[5]-8.665142023030945*coeff[0]*fEdge[5]+59.66213466261491*coeff[2]*fSkin[2]-38.51174232458198*coeff[1]*fSkin[2]-8.893905919223561*coeff[0]*fSkin[2]-59.66213466261491*coeff[2]*fEdge[2]+38.51174232458198*coeff[1]*fEdge[2]+8.893905919223561*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 115.3428423899306*coeff[2]*fSkin[20]-11.19493359006374*coeff[1]*fSkin[20]-43.5036398581567*coeff[0]*fSkin[20]+111.4517585502703*coeff[2]*fEdge[20]-23.25101591782468*coeff[1]*fEdge[20]-40.02334866950417*coeff[0]*fEdge[20]+65.46996195178934*coeff[2]*fSkin[12]-9.94368911043582*coeff[1]*fSkin[12]-24.1121646555219*coeff[0]*fSkin[12]-65.46996195178934*coeff[2]*fEdge[12]+9.94368911043582*coeff[1]*fEdge[12]+24.1121646555219*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 127.0160939089115*coeff[2]*fSkin[21]-24.97331339321913*coeff[1]*fSkin[21]-49.96703777993997*coeff[0]*fSkin[21]-68.64983631400692*coeff[2]*fEdge[21]+85.25372503202384*coeff[1]*fEdge[21]-2.237330049848051*coeff[0]*fEdge[21]+110.8728999785159*coeff[2]*fSkin[6]-58.92212671485608*coeff[1]*fSkin[6]-22.14425183663463*coeff[0]*fSkin[6]+95.80279706881468*coeff[2]*fEdge[6]-74.48646207349731*coeff[1]*fEdge[6]-8.665142023030945*coeff[0]*fEdge[6]+59.66213466261491*coeff[2]*fSkin[3]-38.51174232458198*coeff[1]*fSkin[3]-8.893905919223561*coeff[0]*fSkin[3]-59.66213466261491*coeff[2]*fEdge[3]+38.51174232458198*coeff[1]*fEdge[3]+8.893905919223561*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = 65.46996195178934*coeff[2]*fSkin[33]+3.480291188652535*coeff[1]*fSkin[33]-24.1121646555219*coeff[0]*fSkin[33]+65.46996195178934*coeff[2]*fEdge[33]-3.480291188652535*coeff[1]*fEdge[33]-24.1121646555219*coeff[0]*fEdge[33]+37.79910015670014*coeff[2]*fSkin[22]-13.92116475461015*coeff[0]*fSkin[22]-37.79910015670014*coeff[2]*fEdge[22]+13.92116475461015*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 115.3428423899306*coeff[2]*fSkin[23]-11.19493359006374*coeff[1]*fSkin[23]-43.5036398581567*coeff[0]*fSkin[23]+111.4517585502703*coeff[2]*fEdge[23]-23.25101591782468*coeff[1]*fEdge[23]-40.02334866950417*coeff[0]*fEdge[23]+65.46996195178934*coeff[2]*fSkin[13]-9.94368911043582*coeff[1]*fSkin[13]-24.1121646555219*coeff[0]*fSkin[13]-65.46996195178934*coeff[2]*fEdge[13]+9.94368911043582*coeff[1]*fEdge[13]+24.1121646555219*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = 65.46996195178934*coeff[2]*fSkin[34]+3.480291188652535*coeff[1]*fSkin[34]-24.1121646555219*coeff[0]*fSkin[34]+65.46996195178934*coeff[2]*fEdge[34]-3.480291188652535*coeff[1]*fEdge[34]-24.1121646555219*coeff[0]*fEdge[34]+37.79910015670014*coeff[2]*fSkin[24]-13.92116475461015*coeff[0]*fSkin[24]-37.79910015670014*coeff[2]*fEdge[24]+13.92116475461015*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 127.0160939089115*coeff[2]*fSkin[25]-24.97331339321913*coeff[1]*fSkin[25]-49.96703777993997*coeff[0]*fSkin[25]-68.64983631400692*coeff[2]*fEdge[25]+85.25372503202384*coeff[1]*fEdge[25]-2.237330049848051*coeff[0]*fEdge[25]+110.8728999785159*coeff[2]*fSkin[8]-58.92212671485608*coeff[1]*fSkin[8]-22.14425183663463*coeff[0]*fSkin[8]+95.80279706881468*coeff[2]*fEdge[8]-74.48646207349731*coeff[1]*fEdge[8]-8.665142023030945*coeff[0]*fEdge[8]+59.66213466261491*coeff[2]*fSkin[4]-38.51174232458198*coeff[1]*fSkin[4]-8.893905919223561*coeff[0]*fSkin[4]-59.66213466261491*coeff[2]*fEdge[4]+38.51174232458198*coeff[1]*fEdge[4]+8.893905919223561*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = 65.46996195178934*coeff[2]*fSkin[36]+3.480291188652535*coeff[1]*fSkin[36]-24.1121646555219*coeff[0]*fSkin[36]+65.46996195178934*coeff[2]*fEdge[36]-3.480291188652535*coeff[1]*fEdge[36]-24.1121646555219*coeff[0]*fEdge[36]+37.79910015670014*coeff[2]*fSkin[26]-13.92116475461015*coeff[0]*fSkin[26]-37.79910015670014*coeff[2]*fEdge[26]+13.92116475461015*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = 65.46996195178934*coeff[2]*fSkin[39]+3.480291188652535*coeff[1]*fSkin[39]-24.1121646555219*coeff[0]*fSkin[39]+65.46996195178934*coeff[2]*fEdge[39]-3.480291188652535*coeff[1]*fEdge[39]-24.1121646555219*coeff[0]*fEdge[39]+37.79910015670014*coeff[2]*fSkin[27]-13.92116475461015*coeff[0]*fSkin[27]-37.79910015670014*coeff[2]*fEdge[27]+13.92116475461015*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 115.3428423899306*coeff[2]*fSkin[28]-11.19493359006374*coeff[1]*fSkin[28]-43.5036398581567*coeff[0]*fSkin[28]+111.4517585502703*coeff[2]*fEdge[28]-23.25101591782468*coeff[1]*fEdge[28]-40.02334866950417*coeff[0]*fEdge[28]+65.46996195178934*coeff[2]*fSkin[14]-9.94368911043582*coeff[1]*fSkin[14]-24.1121646555219*coeff[0]*fSkin[14]-65.46996195178934*coeff[2]*fEdge[14]+9.94368911043582*coeff[1]*fEdge[14]+24.1121646555219*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = 65.46996195178934*coeff[2]*fSkin[41]+3.480291188652535*coeff[1]*fSkin[41]-24.1121646555219*coeff[0]*fSkin[41]+65.46996195178934*coeff[2]*fEdge[41]-3.480291188652535*coeff[1]*fEdge[41]-24.1121646555219*coeff[0]*fEdge[41]+37.79910015670014*coeff[2]*fSkin[29]-13.92116475461015*coeff[0]*fSkin[29]-37.79910015670014*coeff[2]*fEdge[29]+13.92116475461015*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = 65.46996195178934*coeff[2]*fSkin[42]+3.480291188652535*coeff[1]*fSkin[42]-24.1121646555219*coeff[0]*fSkin[42]+65.46996195178934*coeff[2]*fEdge[42]-3.480291188652535*coeff[1]*fEdge[42]-24.1121646555219*coeff[0]*fEdge[42]+37.79910015670014*coeff[2]*fSkin[30]-13.92116475461015*coeff[0]*fSkin[30]-37.79910015670014*coeff[2]*fEdge[30]+13.92116475461015*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = 110.8728999785159*coeff[2]*fSkin[44]+9.116253567204131*coeff[1]*fSkin[44]-49.87270631033366*coeff[0]*fSkin[44]-95.80279706881471*coeff[2]*fEdge[44]+37.57675250871954*coeff[1]*fEdge[44]+36.39359649672998*coeff[0]*fEdge[44]+115.3428423899306*coeff[2]*fSkin[31]-11.19493359006374*coeff[1]*fSkin[31]-43.5036398581567*coeff[0]*fSkin[31]+111.4517585502703*coeff[2]*fEdge[31]-23.25101591782468*coeff[1]*fEdge[31]-40.02334866950417*coeff[0]*fEdge[31]+65.46996195178933*coeff[2]*fSkin[18]-9.94368911043582*coeff[1]*fSkin[18]-24.11216465552189*coeff[0]*fSkin[18]-65.46996195178933*coeff[2]*fEdge[18]+9.94368911043582*coeff[1]*fEdge[18]+24.11216465552189*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 127.0160939089115*coeff[2]*fSkin[32]-24.97331339321913*coeff[1]*fSkin[32]-49.96703777993997*coeff[0]*fSkin[32]-68.64983631400692*coeff[2]*fEdge[32]+85.25372503202384*coeff[1]*fEdge[32]-2.237330049848051*coeff[0]*fEdge[32]+110.8728999785158*coeff[2]*fSkin[15]-58.92212671485613*coeff[1]*fSkin[15]-22.14425183663463*coeff[0]*fSkin[15]+95.80279706881466*coeff[2]*fEdge[15]-74.48646207349736*coeff[1]*fEdge[15]-8.665142023030945*coeff[0]*fEdge[15]+59.66213466261492*coeff[2]*fSkin[7]-38.51174232458197*coeff[1]*fSkin[7]-8.893905919223567*coeff[0]*fSkin[7]-59.66213466261492*coeff[2]*fEdge[7]+38.51174232458197*coeff[1]*fEdge[7]+8.893905919223567*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 115.3428423899306*coeff[2]*fSkin[33]-11.19493359006374*coeff[1]*fSkin[33]-43.5036398581567*coeff[0]*fSkin[33]+111.4517585502703*coeff[2]*fEdge[33]-23.25101591782468*coeff[1]*fEdge[33]-40.02334866950417*coeff[0]*fEdge[33]+65.46996195178934*coeff[2]*fSkin[22]-9.94368911043582*coeff[1]*fSkin[22]-24.1121646555219*coeff[0]*fSkin[22]-65.46996195178934*coeff[2]*fEdge[22]+9.94368911043582*coeff[1]*fEdge[22]+24.1121646555219*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 115.3428423899306*coeff[2]*fSkin[34]-11.19493359006374*coeff[1]*fSkin[34]-43.5036398581567*coeff[0]*fSkin[34]+111.4517585502703*coeff[2]*fEdge[34]-23.25101591782468*coeff[1]*fEdge[34]-40.02334866950417*coeff[0]*fEdge[34]+65.46996195178934*coeff[2]*fSkin[24]-9.94368911043582*coeff[1]*fSkin[24]-24.1121646555219*coeff[0]*fSkin[24]-65.46996195178934*coeff[2]*fEdge[24]+9.94368911043582*coeff[1]*fEdge[24]+24.1121646555219*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 127.0160939089115*coeff[2]*fSkin[35]-24.97331339321913*coeff[1]*fSkin[35]-49.96703777993997*coeff[0]*fSkin[35]-68.64983631400692*coeff[2]*fEdge[35]+85.25372503202384*coeff[1]*fEdge[35]-2.237330049848051*coeff[0]*fEdge[35]+110.8728999785158*coeff[2]*fSkin[16]-58.92212671485613*coeff[1]*fSkin[16]-22.14425183663463*coeff[0]*fSkin[16]+95.80279706881466*coeff[2]*fEdge[16]-74.48646207349736*coeff[1]*fEdge[16]-8.665142023030945*coeff[0]*fEdge[16]+59.66213466261492*coeff[2]*fSkin[9]-38.51174232458197*coeff[1]*fSkin[9]-8.893905919223567*coeff[0]*fSkin[9]-59.66213466261492*coeff[2]*fEdge[9]+38.51174232458197*coeff[1]*fEdge[9]+8.893905919223567*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 115.3428423899306*coeff[2]*fSkin[36]-11.19493359006374*coeff[1]*fSkin[36]-43.5036398581567*coeff[0]*fSkin[36]+111.4517585502703*coeff[2]*fEdge[36]-23.25101591782468*coeff[1]*fEdge[36]-40.02334866950417*coeff[0]*fEdge[36]+65.46996195178934*coeff[2]*fSkin[26]-9.94368911043582*coeff[1]*fSkin[26]-24.1121646555219*coeff[0]*fSkin[26]-65.46996195178934*coeff[2]*fEdge[26]+9.94368911043582*coeff[1]*fEdge[26]+24.1121646555219*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 127.0160939089115*coeff[2]*fSkin[37]-24.97331339321913*coeff[1]*fSkin[37]-49.96703777993997*coeff[0]*fSkin[37]-68.64983631400692*coeff[2]*fEdge[37]+85.25372503202384*coeff[1]*fEdge[37]-2.237330049848051*coeff[0]*fEdge[37]+110.8728999785158*coeff[2]*fSkin[17]-58.92212671485613*coeff[1]*fSkin[17]-22.14425183663463*coeff[0]*fSkin[17]+95.80279706881466*coeff[2]*fEdge[17]-74.48646207349736*coeff[1]*fEdge[17]-8.665142023030945*coeff[0]*fEdge[17]+59.66213466261492*coeff[2]*fSkin[10]-38.51174232458197*coeff[1]*fSkin[10]-8.893905919223567*coeff[0]*fSkin[10]-59.66213466261492*coeff[2]*fEdge[10]+38.51174232458197*coeff[1]*fEdge[10]+8.893905919223567*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = 65.46996195178934*coeff[2]*fSkin[45]+3.480291188652535*coeff[1]*fSkin[45]-24.1121646555219*coeff[0]*fSkin[45]+65.46996195178934*coeff[2]*fEdge[45]-3.480291188652535*coeff[1]*fEdge[45]-24.1121646555219*coeff[0]*fEdge[45]+37.79910015670014*coeff[2]*fSkin[38]-13.92116475461015*coeff[0]*fSkin[38]-37.79910015670014*coeff[2]*fEdge[38]+13.92116475461015*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 115.3428423899306*coeff[2]*fSkin[39]-11.19493359006374*coeff[1]*fSkin[39]-43.5036398581567*coeff[0]*fSkin[39]+111.4517585502703*coeff[2]*fEdge[39]-23.25101591782468*coeff[1]*fEdge[39]-40.02334866950417*coeff[0]*fEdge[39]+65.46996195178934*coeff[2]*fSkin[27]-9.94368911043582*coeff[1]*fSkin[27]-24.1121646555219*coeff[0]*fSkin[27]-65.46996195178934*coeff[2]*fEdge[27]+9.94368911043582*coeff[1]*fEdge[27]+24.1121646555219*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = 65.46996195178934*coeff[2]*fSkin[46]+3.480291188652535*coeff[1]*fSkin[46]-24.1121646555219*coeff[0]*fSkin[46]+65.46996195178934*coeff[2]*fEdge[46]-3.480291188652535*coeff[1]*fEdge[46]-24.1121646555219*coeff[0]*fEdge[46]+37.79910015670014*coeff[2]*fSkin[40]-13.92116475461015*coeff[0]*fSkin[40]-37.79910015670014*coeff[2]*fEdge[40]+13.92116475461015*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 115.3428423899306*coeff[2]*fSkin[41]-11.19493359006374*coeff[1]*fSkin[41]-43.5036398581567*coeff[0]*fSkin[41]+111.4517585502703*coeff[2]*fEdge[41]-23.25101591782468*coeff[1]*fEdge[41]-40.02334866950417*coeff[0]*fEdge[41]+65.46996195178934*coeff[2]*fSkin[29]-9.94368911043582*coeff[1]*fSkin[29]-24.1121646555219*coeff[0]*fSkin[29]-65.46996195178934*coeff[2]*fEdge[29]+9.94368911043582*coeff[1]*fEdge[29]+24.1121646555219*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 115.3428423899306*coeff[2]*fSkin[42]-11.19493359006374*coeff[1]*fSkin[42]-43.5036398581567*coeff[0]*fSkin[42]+111.4517585502703*coeff[2]*fEdge[42]-23.25101591782468*coeff[1]*fEdge[42]-40.02334866950417*coeff[0]*fEdge[42]+65.46996195178934*coeff[2]*fSkin[30]-9.94368911043582*coeff[1]*fSkin[30]-24.1121646555219*coeff[0]*fSkin[30]-65.46996195178934*coeff[2]*fEdge[30]+9.94368911043582*coeff[1]*fEdge[30]+24.1121646555219*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = 65.46996195178934*coeff[2]*fSkin[47]+3.480291188652535*coeff[1]*fSkin[47]-24.1121646555219*coeff[0]*fSkin[47]+65.46996195178934*coeff[2]*fEdge[47]-3.480291188652535*coeff[1]*fEdge[47]-24.1121646555219*coeff[0]*fEdge[47]+37.79910015670014*coeff[2]*fSkin[43]-13.92116475461015*coeff[0]*fSkin[43]-37.79910015670014*coeff[2]*fEdge[43]+13.92116475461015*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 127.0160939089115*coeff[2]*fSkin[44]-24.97331339321913*coeff[1]*fSkin[44]-49.96703777993997*coeff[0]*fSkin[44]-68.64983631400692*coeff[2]*fEdge[44]+85.25372503202384*coeff[1]*fEdge[44]-2.237330049848051*coeff[0]*fEdge[44]+110.8728999785159*coeff[2]*fSkin[31]-58.92212671485608*coeff[1]*fSkin[31]-22.14425183663463*coeff[0]*fSkin[31]+95.80279706881468*coeff[2]*fEdge[31]-74.48646207349731*coeff[1]*fEdge[31]-8.665142023030945*coeff[0]*fEdge[31]+59.66213466261491*coeff[2]*fSkin[18]-38.51174232458198*coeff[1]*fSkin[18]-8.893905919223561*coeff[0]*fSkin[18]-59.66213466261491*coeff[2]*fEdge[18]+38.51174232458198*coeff[1]*fEdge[18]+8.893905919223561*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 115.3428423899306*coeff[2]*fSkin[45]-11.19493359006374*coeff[1]*fSkin[45]-43.5036398581567*coeff[0]*fSkin[45]+111.4517585502703*coeff[2]*fEdge[45]-23.25101591782468*coeff[1]*fEdge[45]-40.02334866950417*coeff[0]*fEdge[45]+65.46996195178934*coeff[2]*fSkin[38]-9.94368911043582*coeff[1]*fSkin[38]-24.1121646555219*coeff[0]*fSkin[38]-65.46996195178934*coeff[2]*fEdge[38]+9.94368911043582*coeff[1]*fEdge[38]+24.1121646555219*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 115.3428423899306*coeff[2]*fSkin[46]-11.19493359006374*coeff[1]*fSkin[46]-43.5036398581567*coeff[0]*fSkin[46]+111.4517585502703*coeff[2]*fEdge[46]-23.25101591782468*coeff[1]*fEdge[46]-40.02334866950417*coeff[0]*fEdge[46]+65.46996195178934*coeff[2]*fSkin[40]-9.94368911043582*coeff[1]*fSkin[40]-24.1121646555219*coeff[0]*fSkin[40]-65.46996195178934*coeff[2]*fEdge[40]+9.94368911043582*coeff[1]*fEdge[40]+24.1121646555219*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 115.3428423899306*coeff[2]*fSkin[47]-11.19493359006374*coeff[1]*fSkin[47]-43.5036398581567*coeff[0]*fSkin[47]+111.4517585502703*coeff[2]*fEdge[47]-23.25101591782468*coeff[1]*fEdge[47]-40.02334866950417*coeff[0]*fEdge[47]+65.46996195178934*coeff[2]*fSkin[43]-9.94368911043582*coeff[1]*fSkin[43]-24.1121646555219*coeff[0]*fSkin[43]-65.46996195178934*coeff[2]*fEdge[43]+9.94368911043582*coeff[1]*fEdge[43]+24.1121646555219*coeff[0]*fEdge[43]; 

  boundSurf_incr[0] = 6.960582377305072*coeff[1]*fSkin[1]-26.95821962720737*coeff[1]*fSkin[11]; 
  boundSurf_incr[1] = (-15.07010290970117*coeff[2]*fSkin[11])+46.69300607592371*coeff[1]*fSkin[11]+13.47910981360368*coeff[0]*fSkin[11]+3.891083839660308*fSkin[1]*coeff[2]-12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 6.960582377305072*coeff[1]*fSkin[5]-26.95821962720738*coeff[1]*fSkin[19]; 
  boundSurf_incr[3] = 6.960582377305072*coeff[1]*fSkin[6]-26.95821962720738*coeff[1]*fSkin[21]; 
  boundSurf_incr[4] = 6.960582377305072*coeff[1]*fSkin[8]-26.95821962720738*coeff[1]*fSkin[25]; 
  boundSurf_incr[5] = (-15.07010290970118*coeff[2]*fSkin[19])+46.69300607592368*coeff[1]*fSkin[19]+13.47910981360369*coeff[0]*fSkin[19]+3.891083839660308*coeff[2]*fSkin[5]-12.05608232776094*coeff[1]*fSkin[5]-3.480291188652536*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = (-15.07010290970118*coeff[2]*fSkin[21])+46.69300607592368*coeff[1]*fSkin[21]+13.47910981360369*coeff[0]*fSkin[21]+3.891083839660308*coeff[2]*fSkin[6]-12.05608232776094*coeff[1]*fSkin[6]-3.480291188652536*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = 6.960582377305072*coeff[1]*fSkin[15]-26.95821962720737*coeff[1]*fSkin[32]; 
  boundSurf_incr[8] = (-15.07010290970118*coeff[2]*fSkin[25])+46.69300607592368*coeff[1]*fSkin[25]+13.47910981360369*coeff[0]*fSkin[25]+3.891083839660308*coeff[2]*fSkin[8]-12.05608232776094*coeff[1]*fSkin[8]-3.480291188652536*coeff[0]*fSkin[8]; 
  boundSurf_incr[9] = 6.960582377305072*coeff[1]*fSkin[16]-26.95821962720737*coeff[1]*fSkin[35]; 
  boundSurf_incr[10] = 6.960582377305072*coeff[1]*fSkin[17]-26.95821962720737*coeff[1]*fSkin[37]; 
  boundSurf_incr[11] = 58.36625759490461*coeff[2]*fSkin[11]-60.28041163880471*coeff[1]*fSkin[11]-52.20436782978803*coeff[0]*fSkin[11]-15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]+13.47910981360368*coeff[0]*fSkin[1]; 
  boundSurf_incr[12] = 6.960582377305072*coeff[1]*fSkin[20]; 
  boundSurf_incr[13] = 6.960582377305072*coeff[1]*fSkin[23]; 
  boundSurf_incr[14] = 6.960582377305072*coeff[1]*fSkin[28]; 
  boundSurf_incr[15] = (-15.07010290970117*coeff[2]*fSkin[32])+46.69300607592371*coeff[1]*fSkin[32]+13.47910981360368*coeff[0]*fSkin[32]+3.891083839660308*coeff[2]*fSkin[15]-12.05608232776094*coeff[1]*fSkin[15]-3.480291188652536*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = (-15.07010290970117*coeff[2]*fSkin[35])+46.69300607592371*coeff[1]*fSkin[35]+13.47910981360368*coeff[0]*fSkin[35]+3.891083839660308*coeff[2]*fSkin[16]-12.05608232776094*coeff[1]*fSkin[16]-3.480291188652536*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = (-15.07010290970117*coeff[2]*fSkin[37])+46.69300607592371*coeff[1]*fSkin[37]+13.47910981360368*coeff[0]*fSkin[37]+3.891083839660308*coeff[2]*fSkin[17]-12.05608232776094*coeff[1]*fSkin[17]-3.480291188652536*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = 6.960582377305072*coeff[1]*fSkin[31]-26.95821962720738*coeff[1]*fSkin[44]; 
  boundSurf_incr[19] = 58.36625759490461*coeff[2]*fSkin[19]-60.28041163880471*coeff[1]*fSkin[19]-52.20436782978803*coeff[0]*fSkin[19]-15.07010290970118*coeff[2]*fSkin[5]+15.56433535864123*coeff[1]*fSkin[5]+13.47910981360369*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 3.891083839660308*coeff[2]*fSkin[20]-12.05608232776094*coeff[1]*fSkin[20]-3.480291188652536*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 58.36625759490461*coeff[2]*fSkin[21]-60.28041163880471*coeff[1]*fSkin[21]-52.20436782978803*coeff[0]*fSkin[21]-15.07010290970118*coeff[2]*fSkin[6]+15.56433535864123*coeff[1]*fSkin[6]+13.47910981360369*coeff[0]*fSkin[6]; 
  boundSurf_incr[22] = 6.960582377305072*coeff[1]*fSkin[33]; 
  boundSurf_incr[23] = 3.891083839660308*coeff[2]*fSkin[23]-12.05608232776094*coeff[1]*fSkin[23]-3.480291188652536*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = 6.960582377305072*coeff[1]*fSkin[34]; 
  boundSurf_incr[25] = 58.36625759490461*coeff[2]*fSkin[25]-60.28041163880471*coeff[1]*fSkin[25]-52.20436782978803*coeff[0]*fSkin[25]-15.07010290970118*coeff[2]*fSkin[8]+15.56433535864123*coeff[1]*fSkin[8]+13.47910981360369*coeff[0]*fSkin[8]; 
  boundSurf_incr[26] = 6.960582377305072*coeff[1]*fSkin[36]; 
  boundSurf_incr[27] = 6.960582377305072*coeff[1]*fSkin[39]; 
  boundSurf_incr[28] = 3.891083839660308*coeff[2]*fSkin[28]-12.05608232776094*coeff[1]*fSkin[28]-3.480291188652536*coeff[0]*fSkin[28]; 
  boundSurf_incr[29] = 6.960582377305072*coeff[1]*fSkin[41]; 
  boundSurf_incr[30] = 6.960582377305072*coeff[1]*fSkin[42]; 
  boundSurf_incr[31] = (-15.07010290970118*coeff[2]*fSkin[44])+46.69300607592368*coeff[1]*fSkin[44]+13.47910981360369*coeff[0]*fSkin[44]+3.891083839660308*coeff[2]*fSkin[31]-12.05608232776094*coeff[1]*fSkin[31]-3.480291188652536*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 58.36625759490461*coeff[2]*fSkin[32]-60.28041163880471*coeff[1]*fSkin[32]-52.20436782978803*coeff[0]*fSkin[32]-15.07010290970117*coeff[2]*fSkin[15]+15.56433535864123*coeff[1]*fSkin[15]+13.47910981360368*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 3.891083839660308*coeff[2]*fSkin[33]-12.05608232776094*coeff[1]*fSkin[33]-3.480291188652536*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 3.891083839660308*coeff[2]*fSkin[34]-12.05608232776094*coeff[1]*fSkin[34]-3.480291188652536*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 58.36625759490461*coeff[2]*fSkin[35]-60.28041163880471*coeff[1]*fSkin[35]-52.20436782978803*coeff[0]*fSkin[35]-15.07010290970117*coeff[2]*fSkin[16]+15.56433535864123*coeff[1]*fSkin[16]+13.47910981360368*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 3.891083839660308*coeff[2]*fSkin[36]-12.05608232776094*coeff[1]*fSkin[36]-3.480291188652536*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 58.36625759490461*coeff[2]*fSkin[37]-60.28041163880471*coeff[1]*fSkin[37]-52.20436782978803*coeff[0]*fSkin[37]-15.07010290970117*coeff[2]*fSkin[17]+15.56433535864123*coeff[1]*fSkin[17]+13.47910981360368*coeff[0]*fSkin[17]; 
  boundSurf_incr[38] = 6.960582377305072*coeff[1]*fSkin[45]; 
  boundSurf_incr[39] = 3.891083839660308*coeff[2]*fSkin[39]-12.05608232776094*coeff[1]*fSkin[39]-3.480291188652536*coeff[0]*fSkin[39]; 
  boundSurf_incr[40] = 6.960582377305072*coeff[1]*fSkin[46]; 
  boundSurf_incr[41] = 3.891083839660308*coeff[2]*fSkin[41]-12.05608232776094*coeff[1]*fSkin[41]-3.480291188652536*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 3.891083839660308*coeff[2]*fSkin[42]-12.05608232776094*coeff[1]*fSkin[42]-3.480291188652536*coeff[0]*fSkin[42]; 
  boundSurf_incr[43] = 6.960582377305072*coeff[1]*fSkin[47]; 
  boundSurf_incr[44] = 58.36625759490461*coeff[2]*fSkin[44]-60.28041163880471*coeff[1]*fSkin[44]-52.20436782978803*coeff[0]*fSkin[44]-15.07010290970118*coeff[2]*fSkin[31]+15.56433535864123*coeff[1]*fSkin[31]+13.47910981360369*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 3.891083839660308*coeff[2]*fSkin[45]-12.05608232776094*coeff[1]*fSkin[45]-3.480291188652536*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 3.891083839660308*coeff[2]*fSkin[46]-12.05608232776094*coeff[1]*fSkin[46]-3.480291188652536*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 3.891083839660308*coeff[2]*fSkin[47]-12.05608232776094*coeff[1]*fSkin[47]-3.480291188652536*coeff[0]*fSkin[47]; 

  } else { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[11]-13.47910981360368*coeff[1]*fSkin[11]-24.90293657382598*coeff[0]*fSkin[11]-59.66213466261492*coeff[2]*fEdge[11]-13.47910981360368*coeff[1]*fEdge[11]+24.90293657382598*coeff[0]*fEdge[11]-65.46996195178933*fSkin[1]*coeff[2]-65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]+24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]+24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-110.8728999785158*coeff[2]*fSkin[11])+9.116253567204142*coeff[1]*fSkin[11]+49.87270631033365*coeff[0]*fSkin[11]+95.80279706881466*coeff[2]*fEdge[11]+37.57675250871956*coeff[1]*fEdge[11]-36.39359649672996*coeff[0]*fEdge[11]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]-65.46996195178933*fSkin[0]*coeff[2]+65.46996195178933*fEdge[0]*coeff[2]+11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]+23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]+24.11216465552189*coeff[0]*fSkin[0]-24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 59.66213466261492*coeff[2]*fSkin[19]-13.47910981360369*coeff[1]*fSkin[19]-24.90293657382598*coeff[0]*fSkin[19]-59.66213466261492*coeff[2]*fEdge[19]-13.47910981360369*coeff[1]*fEdge[19]+24.90293657382598*coeff[0]*fEdge[19]-65.46996195178933*coeff[2]*fSkin[5]+3.480291188652536*coeff[1]*fSkin[5]+24.11216465552189*coeff[0]*fSkin[5]-65.46996195178933*coeff[2]*fEdge[5]-3.480291188652536*coeff[1]*fEdge[5]+24.11216465552189*coeff[0]*fEdge[5]+37.79910015670014*coeff[2]*fSkin[2]-13.92116475461015*coeff[0]*fSkin[2]-37.79910015670014*coeff[2]*fEdge[2]+13.92116475461015*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 59.66213466261492*coeff[2]*fSkin[21]-13.47910981360369*coeff[1]*fSkin[21]-24.90293657382598*coeff[0]*fSkin[21]-59.66213466261492*coeff[2]*fEdge[21]-13.47910981360369*coeff[1]*fEdge[21]+24.90293657382598*coeff[0]*fEdge[21]-65.46996195178933*coeff[2]*fSkin[6]+3.480291188652536*coeff[1]*fSkin[6]+24.11216465552189*coeff[0]*fSkin[6]-65.46996195178933*coeff[2]*fEdge[6]-3.480291188652536*coeff[1]*fEdge[6]+24.11216465552189*coeff[0]*fEdge[6]+37.79910015670014*coeff[2]*fSkin[3]-13.92116475461015*coeff[0]*fSkin[3]-37.79910015670014*coeff[2]*fEdge[3]+13.92116475461015*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 59.66213466261492*coeff[2]*fSkin[25]-13.47910981360369*coeff[1]*fSkin[25]-24.90293657382598*coeff[0]*fSkin[25]-59.66213466261492*coeff[2]*fEdge[25]-13.47910981360369*coeff[1]*fEdge[25]+24.90293657382598*coeff[0]*fEdge[25]-65.46996195178933*coeff[2]*fSkin[8]+3.480291188652536*coeff[1]*fSkin[8]+24.11216465552189*coeff[0]*fSkin[8]-65.46996195178933*coeff[2]*fEdge[8]-3.480291188652536*coeff[1]*fEdge[8]+24.11216465552189*coeff[0]*fEdge[8]+37.79910015670014*coeff[2]*fSkin[4]-13.92116475461015*coeff[0]*fSkin[4]-37.79910015670014*coeff[2]*fEdge[4]+13.92116475461015*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-110.8728999785159*coeff[2]*fSkin[19])+9.116253567204131*coeff[1]*fSkin[19]+49.87270631033366*coeff[0]*fSkin[19]+95.80279706881471*coeff[2]*fEdge[19]+37.57675250871954*coeff[1]*fEdge[19]-36.39359649672998*coeff[0]*fEdge[19]+115.3428423899306*coeff[2]*fSkin[5]+11.19493359006374*coeff[1]*fSkin[5]-43.5036398581567*coeff[0]*fSkin[5]+111.4517585502703*coeff[2]*fEdge[5]+23.25101591782468*coeff[1]*fEdge[5]-40.02334866950417*coeff[0]*fEdge[5]-65.46996195178933*coeff[2]*fSkin[2]-9.94368911043582*coeff[1]*fSkin[2]+24.11216465552189*coeff[0]*fSkin[2]+65.46996195178933*coeff[2]*fEdge[2]+9.94368911043582*coeff[1]*fEdge[2]-24.11216465552189*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-110.8728999785159*coeff[2]*fSkin[21])+9.116253567204131*coeff[1]*fSkin[21]+49.87270631033366*coeff[0]*fSkin[21]+95.80279706881471*coeff[2]*fEdge[21]+37.57675250871954*coeff[1]*fEdge[21]-36.39359649672998*coeff[0]*fEdge[21]+115.3428423899306*coeff[2]*fSkin[6]+11.19493359006374*coeff[1]*fSkin[6]-43.5036398581567*coeff[0]*fSkin[6]+111.4517585502703*coeff[2]*fEdge[6]+23.25101591782468*coeff[1]*fEdge[6]-40.02334866950417*coeff[0]*fEdge[6]-65.46996195178933*coeff[2]*fSkin[3]-9.94368911043582*coeff[1]*fSkin[3]+24.11216465552189*coeff[0]*fSkin[3]+65.46996195178933*coeff[2]*fEdge[3]+9.94368911043582*coeff[1]*fEdge[3]-24.11216465552189*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 59.66213466261492*coeff[2]*fSkin[32]-13.47910981360368*coeff[1]*fSkin[32]-24.90293657382598*coeff[0]*fSkin[32]-59.66213466261492*coeff[2]*fEdge[32]-13.47910981360368*coeff[1]*fEdge[32]+24.90293657382598*coeff[0]*fEdge[32]-65.46996195178933*coeff[2]*fSkin[15]+3.480291188652536*coeff[1]*fSkin[15]+24.11216465552189*coeff[0]*fSkin[15]-65.46996195178933*coeff[2]*fEdge[15]-3.480291188652536*coeff[1]*fEdge[15]+24.11216465552189*coeff[0]*fEdge[15]+37.79910015670014*coeff[2]*fSkin[7]-13.92116475461015*coeff[0]*fSkin[7]-37.79910015670014*coeff[2]*fEdge[7]+13.92116475461015*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-110.8728999785159*coeff[2]*fSkin[25])+9.116253567204131*coeff[1]*fSkin[25]+49.87270631033366*coeff[0]*fSkin[25]+95.80279706881471*coeff[2]*fEdge[25]+37.57675250871954*coeff[1]*fEdge[25]-36.39359649672998*coeff[0]*fEdge[25]+115.3428423899306*coeff[2]*fSkin[8]+11.19493359006374*coeff[1]*fSkin[8]-43.5036398581567*coeff[0]*fSkin[8]+111.4517585502703*coeff[2]*fEdge[8]+23.25101591782468*coeff[1]*fEdge[8]-40.02334866950417*coeff[0]*fEdge[8]-65.46996195178933*coeff[2]*fSkin[4]-9.94368911043582*coeff[1]*fSkin[4]+24.11216465552189*coeff[0]*fSkin[4]+65.46996195178933*coeff[2]*fEdge[4]+9.94368911043582*coeff[1]*fEdge[4]-24.11216465552189*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 59.66213466261492*coeff[2]*fSkin[35]-13.47910981360368*coeff[1]*fSkin[35]-24.90293657382598*coeff[0]*fSkin[35]-59.66213466261492*coeff[2]*fEdge[35]-13.47910981360368*coeff[1]*fEdge[35]+24.90293657382598*coeff[0]*fEdge[35]-65.46996195178933*coeff[2]*fSkin[16]+3.480291188652536*coeff[1]*fSkin[16]+24.11216465552189*coeff[0]*fSkin[16]-65.46996195178933*coeff[2]*fEdge[16]-3.480291188652536*coeff[1]*fEdge[16]+24.11216465552189*coeff[0]*fEdge[16]+37.79910015670014*coeff[2]*fSkin[9]-13.92116475461015*coeff[0]*fSkin[9]-37.79910015670014*coeff[2]*fEdge[9]+13.92116475461015*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 59.66213466261492*coeff[2]*fSkin[37]-13.47910981360368*coeff[1]*fSkin[37]-24.90293657382598*coeff[0]*fSkin[37]-59.66213466261492*coeff[2]*fEdge[37]-13.47910981360368*coeff[1]*fEdge[37]+24.90293657382598*coeff[0]*fEdge[37]-65.46996195178933*coeff[2]*fSkin[17]+3.480291188652536*coeff[1]*fSkin[17]+24.11216465552189*coeff[0]*fSkin[17]-65.46996195178933*coeff[2]*fEdge[17]-3.480291188652536*coeff[1]*fEdge[17]+24.11216465552189*coeff[0]*fEdge[17]+37.79910015670014*coeff[2]*fSkin[10]-13.92116475461015*coeff[0]*fSkin[10]-37.79910015670014*coeff[2]*fEdge[10]+13.92116475461015*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 127.0160939089115*coeff[2]*fSkin[11]+24.97331339321913*coeff[1]*fSkin[11]-49.96703777993997*coeff[0]*fSkin[11]-68.64983631400692*coeff[2]*fEdge[11]-85.25372503202384*coeff[1]*fEdge[11]-2.237330049848051*coeff[0]*fEdge[11]-110.8728999785158*fSkin[1]*coeff[2]-95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]+22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]+8.665142023030945*coeff[0]*fEdge[1]+38.51174232458197*fSkin[0]*coeff[1]-38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 
  edgeSurf_incr[12] = (-65.46996195178934*coeff[2]*fSkin[20])+3.480291188652535*coeff[1]*fSkin[20]+24.1121646555219*coeff[0]*fSkin[20]-65.46996195178934*coeff[2]*fEdge[20]-3.480291188652535*coeff[1]*fEdge[20]+24.1121646555219*coeff[0]*fEdge[20]+37.79910015670014*coeff[2]*fSkin[12]-13.92116475461015*coeff[0]*fSkin[12]-37.79910015670014*coeff[2]*fEdge[12]+13.92116475461015*coeff[0]*fEdge[12]; 
  edgeSurf_incr[13] = (-65.46996195178934*coeff[2]*fSkin[23])+3.480291188652535*coeff[1]*fSkin[23]+24.1121646555219*coeff[0]*fSkin[23]-65.46996195178934*coeff[2]*fEdge[23]-3.480291188652535*coeff[1]*fEdge[23]+24.1121646555219*coeff[0]*fEdge[23]+37.79910015670014*coeff[2]*fSkin[13]-13.92116475461015*coeff[0]*fSkin[13]-37.79910015670014*coeff[2]*fEdge[13]+13.92116475461015*coeff[0]*fEdge[13]; 
  edgeSurf_incr[14] = (-65.46996195178934*coeff[2]*fSkin[28])+3.480291188652535*coeff[1]*fSkin[28]+24.1121646555219*coeff[0]*fSkin[28]-65.46996195178934*coeff[2]*fEdge[28]-3.480291188652535*coeff[1]*fEdge[28]+24.1121646555219*coeff[0]*fEdge[28]+37.79910015670014*coeff[2]*fSkin[14]-13.92116475461015*coeff[0]*fSkin[14]-37.79910015670014*coeff[2]*fEdge[14]+13.92116475461015*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-110.8728999785158*coeff[2]*fSkin[32])+9.116253567204142*coeff[1]*fSkin[32]+49.87270631033365*coeff[0]*fSkin[32]+95.80279706881466*coeff[2]*fEdge[32]+37.57675250871956*coeff[1]*fEdge[32]-36.39359649672996*coeff[0]*fEdge[32]+115.3428423899306*coeff[2]*fSkin[15]+11.19493359006374*coeff[1]*fSkin[15]-43.5036398581567*coeff[0]*fSkin[15]+111.4517585502703*coeff[2]*fEdge[15]+23.25101591782468*coeff[1]*fEdge[15]-40.02334866950417*coeff[0]*fEdge[15]-65.46996195178933*coeff[2]*fSkin[7]-9.94368911043582*coeff[1]*fSkin[7]+24.11216465552189*coeff[0]*fSkin[7]+65.46996195178933*coeff[2]*fEdge[7]+9.94368911043582*coeff[1]*fEdge[7]-24.11216465552189*coeff[0]*fEdge[7]; 
  edgeSurf_incr[16] = (-110.8728999785158*coeff[2]*fSkin[35])+9.116253567204142*coeff[1]*fSkin[35]+49.87270631033365*coeff[0]*fSkin[35]+95.80279706881466*coeff[2]*fEdge[35]+37.57675250871956*coeff[1]*fEdge[35]-36.39359649672996*coeff[0]*fEdge[35]+115.3428423899306*coeff[2]*fSkin[16]+11.19493359006374*coeff[1]*fSkin[16]-43.5036398581567*coeff[0]*fSkin[16]+111.4517585502703*coeff[2]*fEdge[16]+23.25101591782468*coeff[1]*fEdge[16]-40.02334866950417*coeff[0]*fEdge[16]-65.46996195178933*coeff[2]*fSkin[9]-9.94368911043582*coeff[1]*fSkin[9]+24.11216465552189*coeff[0]*fSkin[9]+65.46996195178933*coeff[2]*fEdge[9]+9.94368911043582*coeff[1]*fEdge[9]-24.11216465552189*coeff[0]*fEdge[9]; 
  edgeSurf_incr[17] = (-110.8728999785158*coeff[2]*fSkin[37])+9.116253567204142*coeff[1]*fSkin[37]+49.87270631033365*coeff[0]*fSkin[37]+95.80279706881466*coeff[2]*fEdge[37]+37.57675250871956*coeff[1]*fEdge[37]-36.39359649672996*coeff[0]*fEdge[37]+115.3428423899306*coeff[2]*fSkin[17]+11.19493359006374*coeff[1]*fSkin[17]-43.5036398581567*coeff[0]*fSkin[17]+111.4517585502703*coeff[2]*fEdge[17]+23.25101591782468*coeff[1]*fEdge[17]-40.02334866950417*coeff[0]*fEdge[17]-65.46996195178933*coeff[2]*fSkin[10]-9.94368911043582*coeff[1]*fSkin[10]+24.11216465552189*coeff[0]*fSkin[10]+65.46996195178933*coeff[2]*fEdge[10]+9.94368911043582*coeff[1]*fEdge[10]-24.11216465552189*coeff[0]*fEdge[10]; 
  edgeSurf_incr[18] = 59.66213466261492*coeff[2]*fSkin[44]-13.47910981360369*coeff[1]*fSkin[44]-24.90293657382598*coeff[0]*fSkin[44]-59.66213466261492*coeff[2]*fEdge[44]-13.47910981360369*coeff[1]*fEdge[44]+24.90293657382598*coeff[0]*fEdge[44]-65.46996195178933*coeff[2]*fSkin[31]+3.480291188652536*coeff[1]*fSkin[31]+24.11216465552189*coeff[0]*fSkin[31]-65.46996195178933*coeff[2]*fEdge[31]-3.480291188652536*coeff[1]*fEdge[31]+24.11216465552189*coeff[0]*fEdge[31]+37.79910015670014*coeff[2]*fSkin[18]-13.92116475461015*coeff[0]*fSkin[18]-37.79910015670014*coeff[2]*fEdge[18]+13.92116475461015*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 127.0160939089115*coeff[2]*fSkin[19]+24.97331339321913*coeff[1]*fSkin[19]-49.96703777993997*coeff[0]*fSkin[19]-68.64983631400692*coeff[2]*fEdge[19]-85.25372503202384*coeff[1]*fEdge[19]-2.237330049848051*coeff[0]*fEdge[19]-110.8728999785159*coeff[2]*fSkin[5]-58.92212671485608*coeff[1]*fSkin[5]+22.14425183663463*coeff[0]*fSkin[5]-95.80279706881468*coeff[2]*fEdge[5]-74.48646207349731*coeff[1]*fEdge[5]+8.665142023030945*coeff[0]*fEdge[5]+59.66213466261491*coeff[2]*fSkin[2]+38.51174232458198*coeff[1]*fSkin[2]-8.893905919223561*coeff[0]*fSkin[2]-59.66213466261491*coeff[2]*fEdge[2]-38.51174232458198*coeff[1]*fEdge[2]+8.893905919223561*coeff[0]*fEdge[2]; 
  edgeSurf_incr[20] = 115.3428423899306*coeff[2]*fSkin[20]+11.19493359006374*coeff[1]*fSkin[20]-43.5036398581567*coeff[0]*fSkin[20]+111.4517585502703*coeff[2]*fEdge[20]+23.25101591782468*coeff[1]*fEdge[20]-40.02334866950417*coeff[0]*fEdge[20]-65.46996195178934*coeff[2]*fSkin[12]-9.94368911043582*coeff[1]*fSkin[12]+24.1121646555219*coeff[0]*fSkin[12]+65.46996195178934*coeff[2]*fEdge[12]+9.94368911043582*coeff[1]*fEdge[12]-24.1121646555219*coeff[0]*fEdge[12]; 
  edgeSurf_incr[21] = 127.0160939089115*coeff[2]*fSkin[21]+24.97331339321913*coeff[1]*fSkin[21]-49.96703777993997*coeff[0]*fSkin[21]-68.64983631400692*coeff[2]*fEdge[21]-85.25372503202384*coeff[1]*fEdge[21]-2.237330049848051*coeff[0]*fEdge[21]-110.8728999785159*coeff[2]*fSkin[6]-58.92212671485608*coeff[1]*fSkin[6]+22.14425183663463*coeff[0]*fSkin[6]-95.80279706881468*coeff[2]*fEdge[6]-74.48646207349731*coeff[1]*fEdge[6]+8.665142023030945*coeff[0]*fEdge[6]+59.66213466261491*coeff[2]*fSkin[3]+38.51174232458198*coeff[1]*fSkin[3]-8.893905919223561*coeff[0]*fSkin[3]-59.66213466261491*coeff[2]*fEdge[3]-38.51174232458198*coeff[1]*fEdge[3]+8.893905919223561*coeff[0]*fEdge[3]; 
  edgeSurf_incr[22] = (-65.46996195178934*coeff[2]*fSkin[33])+3.480291188652535*coeff[1]*fSkin[33]+24.1121646555219*coeff[0]*fSkin[33]-65.46996195178934*coeff[2]*fEdge[33]-3.480291188652535*coeff[1]*fEdge[33]+24.1121646555219*coeff[0]*fEdge[33]+37.79910015670014*coeff[2]*fSkin[22]-13.92116475461015*coeff[0]*fSkin[22]-37.79910015670014*coeff[2]*fEdge[22]+13.92116475461015*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 115.3428423899306*coeff[2]*fSkin[23]+11.19493359006374*coeff[1]*fSkin[23]-43.5036398581567*coeff[0]*fSkin[23]+111.4517585502703*coeff[2]*fEdge[23]+23.25101591782468*coeff[1]*fEdge[23]-40.02334866950417*coeff[0]*fEdge[23]-65.46996195178934*coeff[2]*fSkin[13]-9.94368911043582*coeff[1]*fSkin[13]+24.1121646555219*coeff[0]*fSkin[13]+65.46996195178934*coeff[2]*fEdge[13]+9.94368911043582*coeff[1]*fEdge[13]-24.1121646555219*coeff[0]*fEdge[13]; 
  edgeSurf_incr[24] = (-65.46996195178934*coeff[2]*fSkin[34])+3.480291188652535*coeff[1]*fSkin[34]+24.1121646555219*coeff[0]*fSkin[34]-65.46996195178934*coeff[2]*fEdge[34]-3.480291188652535*coeff[1]*fEdge[34]+24.1121646555219*coeff[0]*fEdge[34]+37.79910015670014*coeff[2]*fSkin[24]-13.92116475461015*coeff[0]*fSkin[24]-37.79910015670014*coeff[2]*fEdge[24]+13.92116475461015*coeff[0]*fEdge[24]; 
  edgeSurf_incr[25] = 127.0160939089115*coeff[2]*fSkin[25]+24.97331339321913*coeff[1]*fSkin[25]-49.96703777993997*coeff[0]*fSkin[25]-68.64983631400692*coeff[2]*fEdge[25]-85.25372503202384*coeff[1]*fEdge[25]-2.237330049848051*coeff[0]*fEdge[25]-110.8728999785159*coeff[2]*fSkin[8]-58.92212671485608*coeff[1]*fSkin[8]+22.14425183663463*coeff[0]*fSkin[8]-95.80279706881468*coeff[2]*fEdge[8]-74.48646207349731*coeff[1]*fEdge[8]+8.665142023030945*coeff[0]*fEdge[8]+59.66213466261491*coeff[2]*fSkin[4]+38.51174232458198*coeff[1]*fSkin[4]-8.893905919223561*coeff[0]*fSkin[4]-59.66213466261491*coeff[2]*fEdge[4]-38.51174232458198*coeff[1]*fEdge[4]+8.893905919223561*coeff[0]*fEdge[4]; 
  edgeSurf_incr[26] = (-65.46996195178934*coeff[2]*fSkin[36])+3.480291188652535*coeff[1]*fSkin[36]+24.1121646555219*coeff[0]*fSkin[36]-65.46996195178934*coeff[2]*fEdge[36]-3.480291188652535*coeff[1]*fEdge[36]+24.1121646555219*coeff[0]*fEdge[36]+37.79910015670014*coeff[2]*fSkin[26]-13.92116475461015*coeff[0]*fSkin[26]-37.79910015670014*coeff[2]*fEdge[26]+13.92116475461015*coeff[0]*fEdge[26]; 
  edgeSurf_incr[27] = (-65.46996195178934*coeff[2]*fSkin[39])+3.480291188652535*coeff[1]*fSkin[39]+24.1121646555219*coeff[0]*fSkin[39]-65.46996195178934*coeff[2]*fEdge[39]-3.480291188652535*coeff[1]*fEdge[39]+24.1121646555219*coeff[0]*fEdge[39]+37.79910015670014*coeff[2]*fSkin[27]-13.92116475461015*coeff[0]*fSkin[27]-37.79910015670014*coeff[2]*fEdge[27]+13.92116475461015*coeff[0]*fEdge[27]; 
  edgeSurf_incr[28] = 115.3428423899306*coeff[2]*fSkin[28]+11.19493359006374*coeff[1]*fSkin[28]-43.5036398581567*coeff[0]*fSkin[28]+111.4517585502703*coeff[2]*fEdge[28]+23.25101591782468*coeff[1]*fEdge[28]-40.02334866950417*coeff[0]*fEdge[28]-65.46996195178934*coeff[2]*fSkin[14]-9.94368911043582*coeff[1]*fSkin[14]+24.1121646555219*coeff[0]*fSkin[14]+65.46996195178934*coeff[2]*fEdge[14]+9.94368911043582*coeff[1]*fEdge[14]-24.1121646555219*coeff[0]*fEdge[14]; 
  edgeSurf_incr[29] = (-65.46996195178934*coeff[2]*fSkin[41])+3.480291188652535*coeff[1]*fSkin[41]+24.1121646555219*coeff[0]*fSkin[41]-65.46996195178934*coeff[2]*fEdge[41]-3.480291188652535*coeff[1]*fEdge[41]+24.1121646555219*coeff[0]*fEdge[41]+37.79910015670014*coeff[2]*fSkin[29]-13.92116475461015*coeff[0]*fSkin[29]-37.79910015670014*coeff[2]*fEdge[29]+13.92116475461015*coeff[0]*fEdge[29]; 
  edgeSurf_incr[30] = (-65.46996195178934*coeff[2]*fSkin[42])+3.480291188652535*coeff[1]*fSkin[42]+24.1121646555219*coeff[0]*fSkin[42]-65.46996195178934*coeff[2]*fEdge[42]-3.480291188652535*coeff[1]*fEdge[42]+24.1121646555219*coeff[0]*fEdge[42]+37.79910015670014*coeff[2]*fSkin[30]-13.92116475461015*coeff[0]*fSkin[30]-37.79910015670014*coeff[2]*fEdge[30]+13.92116475461015*coeff[0]*fEdge[30]; 
  edgeSurf_incr[31] = (-110.8728999785159*coeff[2]*fSkin[44])+9.116253567204131*coeff[1]*fSkin[44]+49.87270631033366*coeff[0]*fSkin[44]+95.80279706881471*coeff[2]*fEdge[44]+37.57675250871954*coeff[1]*fEdge[44]-36.39359649672998*coeff[0]*fEdge[44]+115.3428423899306*coeff[2]*fSkin[31]+11.19493359006374*coeff[1]*fSkin[31]-43.5036398581567*coeff[0]*fSkin[31]+111.4517585502703*coeff[2]*fEdge[31]+23.25101591782468*coeff[1]*fEdge[31]-40.02334866950417*coeff[0]*fEdge[31]-65.46996195178933*coeff[2]*fSkin[18]-9.94368911043582*coeff[1]*fSkin[18]+24.11216465552189*coeff[0]*fSkin[18]+65.46996195178933*coeff[2]*fEdge[18]+9.94368911043582*coeff[1]*fEdge[18]-24.11216465552189*coeff[0]*fEdge[18]; 
  edgeSurf_incr[32] = 127.0160939089115*coeff[2]*fSkin[32]+24.97331339321913*coeff[1]*fSkin[32]-49.96703777993997*coeff[0]*fSkin[32]-68.64983631400692*coeff[2]*fEdge[32]-85.25372503202384*coeff[1]*fEdge[32]-2.237330049848051*coeff[0]*fEdge[32]-110.8728999785158*coeff[2]*fSkin[15]-58.92212671485613*coeff[1]*fSkin[15]+22.14425183663463*coeff[0]*fSkin[15]-95.80279706881466*coeff[2]*fEdge[15]-74.48646207349736*coeff[1]*fEdge[15]+8.665142023030945*coeff[0]*fEdge[15]+59.66213466261492*coeff[2]*fSkin[7]+38.51174232458197*coeff[1]*fSkin[7]-8.893905919223567*coeff[0]*fSkin[7]-59.66213466261492*coeff[2]*fEdge[7]-38.51174232458197*coeff[1]*fEdge[7]+8.893905919223567*coeff[0]*fEdge[7]; 
  edgeSurf_incr[33] = 115.3428423899306*coeff[2]*fSkin[33]+11.19493359006374*coeff[1]*fSkin[33]-43.5036398581567*coeff[0]*fSkin[33]+111.4517585502703*coeff[2]*fEdge[33]+23.25101591782468*coeff[1]*fEdge[33]-40.02334866950417*coeff[0]*fEdge[33]-65.46996195178934*coeff[2]*fSkin[22]-9.94368911043582*coeff[1]*fSkin[22]+24.1121646555219*coeff[0]*fSkin[22]+65.46996195178934*coeff[2]*fEdge[22]+9.94368911043582*coeff[1]*fEdge[22]-24.1121646555219*coeff[0]*fEdge[22]; 
  edgeSurf_incr[34] = 115.3428423899306*coeff[2]*fSkin[34]+11.19493359006374*coeff[1]*fSkin[34]-43.5036398581567*coeff[0]*fSkin[34]+111.4517585502703*coeff[2]*fEdge[34]+23.25101591782468*coeff[1]*fEdge[34]-40.02334866950417*coeff[0]*fEdge[34]-65.46996195178934*coeff[2]*fSkin[24]-9.94368911043582*coeff[1]*fSkin[24]+24.1121646555219*coeff[0]*fSkin[24]+65.46996195178934*coeff[2]*fEdge[24]+9.94368911043582*coeff[1]*fEdge[24]-24.1121646555219*coeff[0]*fEdge[24]; 
  edgeSurf_incr[35] = 127.0160939089115*coeff[2]*fSkin[35]+24.97331339321913*coeff[1]*fSkin[35]-49.96703777993997*coeff[0]*fSkin[35]-68.64983631400692*coeff[2]*fEdge[35]-85.25372503202384*coeff[1]*fEdge[35]-2.237330049848051*coeff[0]*fEdge[35]-110.8728999785158*coeff[2]*fSkin[16]-58.92212671485613*coeff[1]*fSkin[16]+22.14425183663463*coeff[0]*fSkin[16]-95.80279706881466*coeff[2]*fEdge[16]-74.48646207349736*coeff[1]*fEdge[16]+8.665142023030945*coeff[0]*fEdge[16]+59.66213466261492*coeff[2]*fSkin[9]+38.51174232458197*coeff[1]*fSkin[9]-8.893905919223567*coeff[0]*fSkin[9]-59.66213466261492*coeff[2]*fEdge[9]-38.51174232458197*coeff[1]*fEdge[9]+8.893905919223567*coeff[0]*fEdge[9]; 
  edgeSurf_incr[36] = 115.3428423899306*coeff[2]*fSkin[36]+11.19493359006374*coeff[1]*fSkin[36]-43.5036398581567*coeff[0]*fSkin[36]+111.4517585502703*coeff[2]*fEdge[36]+23.25101591782468*coeff[1]*fEdge[36]-40.02334866950417*coeff[0]*fEdge[36]-65.46996195178934*coeff[2]*fSkin[26]-9.94368911043582*coeff[1]*fSkin[26]+24.1121646555219*coeff[0]*fSkin[26]+65.46996195178934*coeff[2]*fEdge[26]+9.94368911043582*coeff[1]*fEdge[26]-24.1121646555219*coeff[0]*fEdge[26]; 
  edgeSurf_incr[37] = 127.0160939089115*coeff[2]*fSkin[37]+24.97331339321913*coeff[1]*fSkin[37]-49.96703777993997*coeff[0]*fSkin[37]-68.64983631400692*coeff[2]*fEdge[37]-85.25372503202384*coeff[1]*fEdge[37]-2.237330049848051*coeff[0]*fEdge[37]-110.8728999785158*coeff[2]*fSkin[17]-58.92212671485613*coeff[1]*fSkin[17]+22.14425183663463*coeff[0]*fSkin[17]-95.80279706881466*coeff[2]*fEdge[17]-74.48646207349736*coeff[1]*fEdge[17]+8.665142023030945*coeff[0]*fEdge[17]+59.66213466261492*coeff[2]*fSkin[10]+38.51174232458197*coeff[1]*fSkin[10]-8.893905919223567*coeff[0]*fSkin[10]-59.66213466261492*coeff[2]*fEdge[10]-38.51174232458197*coeff[1]*fEdge[10]+8.893905919223567*coeff[0]*fEdge[10]; 
  edgeSurf_incr[38] = (-65.46996195178934*coeff[2]*fSkin[45])+3.480291188652535*coeff[1]*fSkin[45]+24.1121646555219*coeff[0]*fSkin[45]-65.46996195178934*coeff[2]*fEdge[45]-3.480291188652535*coeff[1]*fEdge[45]+24.1121646555219*coeff[0]*fEdge[45]+37.79910015670014*coeff[2]*fSkin[38]-13.92116475461015*coeff[0]*fSkin[38]-37.79910015670014*coeff[2]*fEdge[38]+13.92116475461015*coeff[0]*fEdge[38]; 
  edgeSurf_incr[39] = 115.3428423899306*coeff[2]*fSkin[39]+11.19493359006374*coeff[1]*fSkin[39]-43.5036398581567*coeff[0]*fSkin[39]+111.4517585502703*coeff[2]*fEdge[39]+23.25101591782468*coeff[1]*fEdge[39]-40.02334866950417*coeff[0]*fEdge[39]-65.46996195178934*coeff[2]*fSkin[27]-9.94368911043582*coeff[1]*fSkin[27]+24.1121646555219*coeff[0]*fSkin[27]+65.46996195178934*coeff[2]*fEdge[27]+9.94368911043582*coeff[1]*fEdge[27]-24.1121646555219*coeff[0]*fEdge[27]; 
  edgeSurf_incr[40] = (-65.46996195178934*coeff[2]*fSkin[46])+3.480291188652535*coeff[1]*fSkin[46]+24.1121646555219*coeff[0]*fSkin[46]-65.46996195178934*coeff[2]*fEdge[46]-3.480291188652535*coeff[1]*fEdge[46]+24.1121646555219*coeff[0]*fEdge[46]+37.79910015670014*coeff[2]*fSkin[40]-13.92116475461015*coeff[0]*fSkin[40]-37.79910015670014*coeff[2]*fEdge[40]+13.92116475461015*coeff[0]*fEdge[40]; 
  edgeSurf_incr[41] = 115.3428423899306*coeff[2]*fSkin[41]+11.19493359006374*coeff[1]*fSkin[41]-43.5036398581567*coeff[0]*fSkin[41]+111.4517585502703*coeff[2]*fEdge[41]+23.25101591782468*coeff[1]*fEdge[41]-40.02334866950417*coeff[0]*fEdge[41]-65.46996195178934*coeff[2]*fSkin[29]-9.94368911043582*coeff[1]*fSkin[29]+24.1121646555219*coeff[0]*fSkin[29]+65.46996195178934*coeff[2]*fEdge[29]+9.94368911043582*coeff[1]*fEdge[29]-24.1121646555219*coeff[0]*fEdge[29]; 
  edgeSurf_incr[42] = 115.3428423899306*coeff[2]*fSkin[42]+11.19493359006374*coeff[1]*fSkin[42]-43.5036398581567*coeff[0]*fSkin[42]+111.4517585502703*coeff[2]*fEdge[42]+23.25101591782468*coeff[1]*fEdge[42]-40.02334866950417*coeff[0]*fEdge[42]-65.46996195178934*coeff[2]*fSkin[30]-9.94368911043582*coeff[1]*fSkin[30]+24.1121646555219*coeff[0]*fSkin[30]+65.46996195178934*coeff[2]*fEdge[30]+9.94368911043582*coeff[1]*fEdge[30]-24.1121646555219*coeff[0]*fEdge[30]; 
  edgeSurf_incr[43] = (-65.46996195178934*coeff[2]*fSkin[47])+3.480291188652535*coeff[1]*fSkin[47]+24.1121646555219*coeff[0]*fSkin[47]-65.46996195178934*coeff[2]*fEdge[47]-3.480291188652535*coeff[1]*fEdge[47]+24.1121646555219*coeff[0]*fEdge[47]+37.79910015670014*coeff[2]*fSkin[43]-13.92116475461015*coeff[0]*fSkin[43]-37.79910015670014*coeff[2]*fEdge[43]+13.92116475461015*coeff[0]*fEdge[43]; 
  edgeSurf_incr[44] = 127.0160939089115*coeff[2]*fSkin[44]+24.97331339321913*coeff[1]*fSkin[44]-49.96703777993997*coeff[0]*fSkin[44]-68.64983631400692*coeff[2]*fEdge[44]-85.25372503202384*coeff[1]*fEdge[44]-2.237330049848051*coeff[0]*fEdge[44]-110.8728999785159*coeff[2]*fSkin[31]-58.92212671485608*coeff[1]*fSkin[31]+22.14425183663463*coeff[0]*fSkin[31]-95.80279706881468*coeff[2]*fEdge[31]-74.48646207349731*coeff[1]*fEdge[31]+8.665142023030945*coeff[0]*fEdge[31]+59.66213466261491*coeff[2]*fSkin[18]+38.51174232458198*coeff[1]*fSkin[18]-8.893905919223561*coeff[0]*fSkin[18]-59.66213466261491*coeff[2]*fEdge[18]-38.51174232458198*coeff[1]*fEdge[18]+8.893905919223561*coeff[0]*fEdge[18]; 
  edgeSurf_incr[45] = 115.3428423899306*coeff[2]*fSkin[45]+11.19493359006374*coeff[1]*fSkin[45]-43.5036398581567*coeff[0]*fSkin[45]+111.4517585502703*coeff[2]*fEdge[45]+23.25101591782468*coeff[1]*fEdge[45]-40.02334866950417*coeff[0]*fEdge[45]-65.46996195178934*coeff[2]*fSkin[38]-9.94368911043582*coeff[1]*fSkin[38]+24.1121646555219*coeff[0]*fSkin[38]+65.46996195178934*coeff[2]*fEdge[38]+9.94368911043582*coeff[1]*fEdge[38]-24.1121646555219*coeff[0]*fEdge[38]; 
  edgeSurf_incr[46] = 115.3428423899306*coeff[2]*fSkin[46]+11.19493359006374*coeff[1]*fSkin[46]-43.5036398581567*coeff[0]*fSkin[46]+111.4517585502703*coeff[2]*fEdge[46]+23.25101591782468*coeff[1]*fEdge[46]-40.02334866950417*coeff[0]*fEdge[46]-65.46996195178934*coeff[2]*fSkin[40]-9.94368911043582*coeff[1]*fSkin[40]+24.1121646555219*coeff[0]*fSkin[40]+65.46996195178934*coeff[2]*fEdge[40]+9.94368911043582*coeff[1]*fEdge[40]-24.1121646555219*coeff[0]*fEdge[40]; 
  edgeSurf_incr[47] = 115.3428423899306*coeff[2]*fSkin[47]+11.19493359006374*coeff[1]*fSkin[47]-43.5036398581567*coeff[0]*fSkin[47]+111.4517585502703*coeff[2]*fEdge[47]+23.25101591782468*coeff[1]*fEdge[47]-40.02334866950417*coeff[0]*fEdge[47]-65.46996195178934*coeff[2]*fSkin[43]-9.94368911043582*coeff[1]*fSkin[43]+24.1121646555219*coeff[0]*fSkin[43]+65.46996195178934*coeff[2]*fEdge[43]+9.94368911043582*coeff[1]*fEdge[43]-24.1121646555219*coeff[0]*fEdge[43]; 

  boundSurf_incr[0] = 26.95821962720737*coeff[1]*fSkin[11]+6.960582377305072*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 15.07010290970117*coeff[2]*fSkin[11]+46.69300607592371*coeff[1]*fSkin[11]-13.47910981360368*coeff[0]*fSkin[11]+3.891083839660308*fSkin[1]*coeff[2]+12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 26.95821962720738*coeff[1]*fSkin[19]+6.960582377305072*coeff[1]*fSkin[5]; 
  boundSurf_incr[3] = 26.95821962720738*coeff[1]*fSkin[21]+6.960582377305072*coeff[1]*fSkin[6]; 
  boundSurf_incr[4] = 26.95821962720738*coeff[1]*fSkin[25]+6.960582377305072*coeff[1]*fSkin[8]; 
  boundSurf_incr[5] = 15.07010290970118*coeff[2]*fSkin[19]+46.69300607592368*coeff[1]*fSkin[19]-13.47910981360369*coeff[0]*fSkin[19]+3.891083839660308*coeff[2]*fSkin[5]+12.05608232776094*coeff[1]*fSkin[5]-3.480291188652536*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 15.07010290970118*coeff[2]*fSkin[21]+46.69300607592368*coeff[1]*fSkin[21]-13.47910981360369*coeff[0]*fSkin[21]+3.891083839660308*coeff[2]*fSkin[6]+12.05608232776094*coeff[1]*fSkin[6]-3.480291188652536*coeff[0]*fSkin[6]; 
  boundSurf_incr[7] = 26.95821962720737*coeff[1]*fSkin[32]+6.960582377305072*coeff[1]*fSkin[15]; 
  boundSurf_incr[8] = 15.07010290970118*coeff[2]*fSkin[25]+46.69300607592368*coeff[1]*fSkin[25]-13.47910981360369*coeff[0]*fSkin[25]+3.891083839660308*coeff[2]*fSkin[8]+12.05608232776094*coeff[1]*fSkin[8]-3.480291188652536*coeff[0]*fSkin[8]; 
  boundSurf_incr[9] = 26.95821962720737*coeff[1]*fSkin[35]+6.960582377305072*coeff[1]*fSkin[16]; 
  boundSurf_incr[10] = 26.95821962720737*coeff[1]*fSkin[37]+6.960582377305072*coeff[1]*fSkin[17]; 
  boundSurf_incr[11] = 58.36625759490461*coeff[2]*fSkin[11]+60.28041163880471*coeff[1]*fSkin[11]-52.20436782978803*coeff[0]*fSkin[11]+15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]-13.47910981360368*coeff[0]*fSkin[1]; 
  boundSurf_incr[12] = 6.960582377305072*coeff[1]*fSkin[20]; 
  boundSurf_incr[13] = 6.960582377305072*coeff[1]*fSkin[23]; 
  boundSurf_incr[14] = 6.960582377305072*coeff[1]*fSkin[28]; 
  boundSurf_incr[15] = 15.07010290970117*coeff[2]*fSkin[32]+46.69300607592371*coeff[1]*fSkin[32]-13.47910981360368*coeff[0]*fSkin[32]+3.891083839660308*coeff[2]*fSkin[15]+12.05608232776094*coeff[1]*fSkin[15]-3.480291188652536*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 15.07010290970117*coeff[2]*fSkin[35]+46.69300607592371*coeff[1]*fSkin[35]-13.47910981360368*coeff[0]*fSkin[35]+3.891083839660308*coeff[2]*fSkin[16]+12.05608232776094*coeff[1]*fSkin[16]-3.480291188652536*coeff[0]*fSkin[16]; 
  boundSurf_incr[17] = 15.07010290970117*coeff[2]*fSkin[37]+46.69300607592371*coeff[1]*fSkin[37]-13.47910981360368*coeff[0]*fSkin[37]+3.891083839660308*coeff[2]*fSkin[17]+12.05608232776094*coeff[1]*fSkin[17]-3.480291188652536*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = 26.95821962720738*coeff[1]*fSkin[44]+6.960582377305072*coeff[1]*fSkin[31]; 
  boundSurf_incr[19] = 58.36625759490461*coeff[2]*fSkin[19]+60.28041163880471*coeff[1]*fSkin[19]-52.20436782978803*coeff[0]*fSkin[19]+15.07010290970118*coeff[2]*fSkin[5]+15.56433535864123*coeff[1]*fSkin[5]-13.47910981360369*coeff[0]*fSkin[5]; 
  boundSurf_incr[20] = 3.891083839660308*coeff[2]*fSkin[20]+12.05608232776094*coeff[1]*fSkin[20]-3.480291188652536*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 58.36625759490461*coeff[2]*fSkin[21]+60.28041163880471*coeff[1]*fSkin[21]-52.20436782978803*coeff[0]*fSkin[21]+15.07010290970118*coeff[2]*fSkin[6]+15.56433535864123*coeff[1]*fSkin[6]-13.47910981360369*coeff[0]*fSkin[6]; 
  boundSurf_incr[22] = 6.960582377305072*coeff[1]*fSkin[33]; 
  boundSurf_incr[23] = 3.891083839660308*coeff[2]*fSkin[23]+12.05608232776094*coeff[1]*fSkin[23]-3.480291188652536*coeff[0]*fSkin[23]; 
  boundSurf_incr[24] = 6.960582377305072*coeff[1]*fSkin[34]; 
  boundSurf_incr[25] = 58.36625759490461*coeff[2]*fSkin[25]+60.28041163880471*coeff[1]*fSkin[25]-52.20436782978803*coeff[0]*fSkin[25]+15.07010290970118*coeff[2]*fSkin[8]+15.56433535864123*coeff[1]*fSkin[8]-13.47910981360369*coeff[0]*fSkin[8]; 
  boundSurf_incr[26] = 6.960582377305072*coeff[1]*fSkin[36]; 
  boundSurf_incr[27] = 6.960582377305072*coeff[1]*fSkin[39]; 
  boundSurf_incr[28] = 3.891083839660308*coeff[2]*fSkin[28]+12.05608232776094*coeff[1]*fSkin[28]-3.480291188652536*coeff[0]*fSkin[28]; 
  boundSurf_incr[29] = 6.960582377305072*coeff[1]*fSkin[41]; 
  boundSurf_incr[30] = 6.960582377305072*coeff[1]*fSkin[42]; 
  boundSurf_incr[31] = 15.07010290970118*coeff[2]*fSkin[44]+46.69300607592368*coeff[1]*fSkin[44]-13.47910981360369*coeff[0]*fSkin[44]+3.891083839660308*coeff[2]*fSkin[31]+12.05608232776094*coeff[1]*fSkin[31]-3.480291188652536*coeff[0]*fSkin[31]; 
  boundSurf_incr[32] = 58.36625759490461*coeff[2]*fSkin[32]+60.28041163880471*coeff[1]*fSkin[32]-52.20436782978803*coeff[0]*fSkin[32]+15.07010290970117*coeff[2]*fSkin[15]+15.56433535864123*coeff[1]*fSkin[15]-13.47910981360368*coeff[0]*fSkin[15]; 
  boundSurf_incr[33] = 3.891083839660308*coeff[2]*fSkin[33]+12.05608232776094*coeff[1]*fSkin[33]-3.480291188652536*coeff[0]*fSkin[33]; 
  boundSurf_incr[34] = 3.891083839660308*coeff[2]*fSkin[34]+12.05608232776094*coeff[1]*fSkin[34]-3.480291188652536*coeff[0]*fSkin[34]; 
  boundSurf_incr[35] = 58.36625759490461*coeff[2]*fSkin[35]+60.28041163880471*coeff[1]*fSkin[35]-52.20436782978803*coeff[0]*fSkin[35]+15.07010290970117*coeff[2]*fSkin[16]+15.56433535864123*coeff[1]*fSkin[16]-13.47910981360368*coeff[0]*fSkin[16]; 
  boundSurf_incr[36] = 3.891083839660308*coeff[2]*fSkin[36]+12.05608232776094*coeff[1]*fSkin[36]-3.480291188652536*coeff[0]*fSkin[36]; 
  boundSurf_incr[37] = 58.36625759490461*coeff[2]*fSkin[37]+60.28041163880471*coeff[1]*fSkin[37]-52.20436782978803*coeff[0]*fSkin[37]+15.07010290970117*coeff[2]*fSkin[17]+15.56433535864123*coeff[1]*fSkin[17]-13.47910981360368*coeff[0]*fSkin[17]; 
  boundSurf_incr[38] = 6.960582377305072*coeff[1]*fSkin[45]; 
  boundSurf_incr[39] = 3.891083839660308*coeff[2]*fSkin[39]+12.05608232776094*coeff[1]*fSkin[39]-3.480291188652536*coeff[0]*fSkin[39]; 
  boundSurf_incr[40] = 6.960582377305072*coeff[1]*fSkin[46]; 
  boundSurf_incr[41] = 3.891083839660308*coeff[2]*fSkin[41]+12.05608232776094*coeff[1]*fSkin[41]-3.480291188652536*coeff[0]*fSkin[41]; 
  boundSurf_incr[42] = 3.891083839660308*coeff[2]*fSkin[42]+12.05608232776094*coeff[1]*fSkin[42]-3.480291188652536*coeff[0]*fSkin[42]; 
  boundSurf_incr[43] = 6.960582377305072*coeff[1]*fSkin[47]; 
  boundSurf_incr[44] = 58.36625759490461*coeff[2]*fSkin[44]+60.28041163880471*coeff[1]*fSkin[44]-52.20436782978803*coeff[0]*fSkin[44]+15.07010290970118*coeff[2]*fSkin[31]+15.56433535864123*coeff[1]*fSkin[31]-13.47910981360369*coeff[0]*fSkin[31]; 
  boundSurf_incr[45] = 3.891083839660308*coeff[2]*fSkin[45]+12.05608232776094*coeff[1]*fSkin[45]-3.480291188652536*coeff[0]*fSkin[45]; 
  boundSurf_incr[46] = 3.891083839660308*coeff[2]*fSkin[46]+12.05608232776094*coeff[1]*fSkin[46]-3.480291188652536*coeff[0]*fSkin[46]; 
  boundSurf_incr[47] = 3.891083839660308*coeff[2]*fSkin[47]+12.05608232776094*coeff[1]*fSkin[47]-3.480291188652536*coeff[0]*fSkin[47]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += (vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += (vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += (vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += (vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += (vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 
  out[8] += (vol_incr[8]+edgeSurf_incr[8]+boundSurf_incr[8])*Jfac; 
  out[9] += (vol_incr[9]+edgeSurf_incr[9]+boundSurf_incr[9])*Jfac; 
  out[10] += (vol_incr[10]+edgeSurf_incr[10]+boundSurf_incr[10])*Jfac; 
  out[11] += (vol_incr[11]+edgeSurf_incr[11]+boundSurf_incr[11])*Jfac; 
  out[12] += (vol_incr[12]+edgeSurf_incr[12]+boundSurf_incr[12])*Jfac; 
  out[13] += (vol_incr[13]+edgeSurf_incr[13]+boundSurf_incr[13])*Jfac; 
  out[14] += (vol_incr[14]+edgeSurf_incr[14]+boundSurf_incr[14])*Jfac; 
  out[15] += (vol_incr[15]+edgeSurf_incr[15]+boundSurf_incr[15])*Jfac; 
  out[16] += (vol_incr[16]+edgeSurf_incr[16]+boundSurf_incr[16])*Jfac; 
  out[17] += (vol_incr[17]+edgeSurf_incr[17]+boundSurf_incr[17])*Jfac; 
  out[18] += (vol_incr[18]+edgeSurf_incr[18]+boundSurf_incr[18])*Jfac; 
  out[19] += (vol_incr[19]+edgeSurf_incr[19]+boundSurf_incr[19])*Jfac; 
  out[20] += (vol_incr[20]+edgeSurf_incr[20]+boundSurf_incr[20])*Jfac; 
  out[21] += (vol_incr[21]+edgeSurf_incr[21]+boundSurf_incr[21])*Jfac; 
  out[22] += (vol_incr[22]+edgeSurf_incr[22]+boundSurf_incr[22])*Jfac; 
  out[23] += (vol_incr[23]+edgeSurf_incr[23]+boundSurf_incr[23])*Jfac; 
  out[24] += (vol_incr[24]+edgeSurf_incr[24]+boundSurf_incr[24])*Jfac; 
  out[25] += (vol_incr[25]+edgeSurf_incr[25]+boundSurf_incr[25])*Jfac; 
  out[26] += (vol_incr[26]+edgeSurf_incr[26]+boundSurf_incr[26])*Jfac; 
  out[27] += (vol_incr[27]+edgeSurf_incr[27]+boundSurf_incr[27])*Jfac; 
  out[28] += (vol_incr[28]+edgeSurf_incr[28]+boundSurf_incr[28])*Jfac; 
  out[29] += (vol_incr[29]+edgeSurf_incr[29]+boundSurf_incr[29])*Jfac; 
  out[30] += (vol_incr[30]+edgeSurf_incr[30]+boundSurf_incr[30])*Jfac; 
  out[31] += (vol_incr[31]+edgeSurf_incr[31]+boundSurf_incr[31])*Jfac; 
  out[32] += (vol_incr[32]+edgeSurf_incr[32]+boundSurf_incr[32])*Jfac; 
  out[33] += (vol_incr[33]+edgeSurf_incr[33]+boundSurf_incr[33])*Jfac; 
  out[34] += (vol_incr[34]+edgeSurf_incr[34]+boundSurf_incr[34])*Jfac; 
  out[35] += (vol_incr[35]+edgeSurf_incr[35]+boundSurf_incr[35])*Jfac; 
  out[36] += (vol_incr[36]+edgeSurf_incr[36]+boundSurf_incr[36])*Jfac; 
  out[37] += (vol_incr[37]+edgeSurf_incr[37]+boundSurf_incr[37])*Jfac; 
  out[38] += (vol_incr[38]+edgeSurf_incr[38]+boundSurf_incr[38])*Jfac; 
  out[39] += (vol_incr[39]+edgeSurf_incr[39]+boundSurf_incr[39])*Jfac; 
  out[40] += (vol_incr[40]+edgeSurf_incr[40]+boundSurf_incr[40])*Jfac; 
  out[41] += (vol_incr[41]+edgeSurf_incr[41]+boundSurf_incr[41])*Jfac; 
  out[42] += (vol_incr[42]+edgeSurf_incr[42]+boundSurf_incr[42])*Jfac; 
  out[43] += (vol_incr[43]+edgeSurf_incr[43]+boundSurf_incr[43])*Jfac; 
  out[44] += (vol_incr[44]+edgeSurf_incr[44]+boundSurf_incr[44])*Jfac; 
  out[45] += (vol_incr[45]+edgeSurf_incr[45]+boundSurf_incr[45])*Jfac; 
  out[46] += (vol_incr[46]+edgeSurf_incr[46]+boundSurf_incr[46])*Jfac; 
  out[47] += (vol_incr[47]+edgeSurf_incr[47]+boundSurf_incr[47])*Jfac; 

  }

  return 0.;
}

