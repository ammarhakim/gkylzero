#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_surfx_1x2v_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[7])+35.21807064562169*coeff[0]*fEdge[7]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-70.53065765632411*coeff[0]*fSkin[7])+51.46831774920949*coeff[0]*fEdge[7]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[11])+35.21807064562168*coeff[0]*fEdge[11]-34.09975027401226*coeff[0]*fSkin[4]-34.09975027401226*coeff[0]*fEdge[4]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[13])+35.21807064562168*coeff[0]*fEdge[13]-34.09975027401226*coeff[0]*fSkin[5]-34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-70.53065765632414*coeff[0]*fSkin[11])+51.4683177492095*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[4]-56.6015625*coeff[0]*fEdge[4]-34.09975027401226*coeff[0]*fSkin[2]+34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-70.53065765632414*coeff[0]*fSkin[13])+51.4683177492095*coeff[0]*fEdge[13]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]-34.09975027401226*coeff[0]*fSkin[3]+34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-35.21807064562169*coeff[0]*fSkin[17])+35.21807064562169*coeff[0]*fEdge[17]-34.09975027401226*coeff[0]*fSkin[10]-34.09975027401226*coeff[0]*fEdge[10]-19.6875*coeff[0]*fSkin[6]+19.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-70.6640625*coeff[0]*fSkin[7])-3.1640625*coeff[0]*fEdge[7]-31.316701275974*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-34.09975027401227*coeff[0]*fSkin[12])-34.09975027401227*coeff[0]*fEdge[12]-19.6875*coeff[0]*fSkin[8]+19.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-34.09975027401227*coeff[0]*fSkin[15])-34.09975027401227*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-70.53065765632411*coeff[0]*fSkin[17])+51.46831774920949*coeff[0]*fEdge[17]-61.5234375*coeff[0]*fSkin[10]-56.6015625*coeff[0]*fEdge[10]-34.09975027401226*coeff[0]*fSkin[6]+34.09975027401226*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]-31.31670127597403*coeff[0]*fSkin[4]-12.2543613688594*coeff[0]*fEdge[4]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-61.5234375*coeff[0]*fSkin[12])-56.6015625*coeff[0]*fEdge[12]-34.09975027401227*coeff[0]*fSkin[8]+34.09975027401227*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-70.6640625*coeff[0]*fSkin[13])-3.1640625*coeff[0]*fEdge[13]-31.31670127597403*coeff[0]*fSkin[5]-12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-34.09975027401227*coeff[0]*fSkin[18])-34.09975027401227*coeff[0]*fEdge[18]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-61.5234375*coeff[0]*fSkin[15])-56.6015625*coeff[0]*fEdge[15]-34.09975027401227*coeff[0]*fSkin[9]+34.09975027401227*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-34.09975027401227*coeff[0]*fSkin[19])-34.09975027401227*coeff[0]*fEdge[19]-19.6875*coeff[0]*fSkin[16]+19.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-70.6640625*coeff[0]*fSkin[17])-3.1640625*coeff[0]*fEdge[17]-31.316701275974*coeff[0]*fSkin[10]-12.2543613688594*coeff[0]*fEdge[10]-12.57788237343632*coeff[0]*fSkin[6]+12.57788237343632*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-61.5234375*coeff[0]*fSkin[18])-56.6015625*coeff[0]*fEdge[18]-34.09975027401227*coeff[0]*fSkin[14]+34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-61.5234375*coeff[0]*fSkin[19])-56.6015625*coeff[0]*fEdge[19]-34.09975027401227*coeff[0]*fSkin[16]+34.09975027401227*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = 19.06233990711463*coeff[0]*fSkin[7]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 19.06233990711463*coeff[0]*fSkin[11]-4.921875*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 19.06233990711463*coeff[0]*fSkin[13]-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 19.06233990711463*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[7]; 
  boundSurf_incr[10] = 19.06233990711463*coeff[0]*fSkin[17]-4.921875*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 19.06233990711463*coeff[0]*fSkin[4]-73.828125*coeff[0]*fSkin[11]; 
  boundSurf_incr[12] = -4.921875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 19.06233990711463*coeff[0]*fSkin[5]-73.828125*coeff[0]*fSkin[13]; 
  boundSurf_incr[15] = -4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 19.06233990711463*coeff[0]*fSkin[10]-73.828125*coeff[0]*fSkin[17]; 
  boundSurf_incr[18] = -4.921875*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = -4.921875*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[7])+35.21807064562169*coeff[0]*fEdge[7]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[7]-51.46831774920949*coeff[0]*fEdge[7]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-35.21807064562168*coeff[0]*fSkin[11])+35.21807064562168*coeff[0]*fEdge[11]+34.09975027401226*coeff[0]*fSkin[4]+34.09975027401226*coeff[0]*fEdge[4]-19.6875*coeff[0]*fSkin[2]+19.6875*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-35.21807064562168*coeff[0]*fSkin[13])+35.21807064562168*coeff[0]*fEdge[13]+34.09975027401226*coeff[0]*fSkin[5]+34.09975027401226*coeff[0]*fEdge[5]-19.6875*coeff[0]*fSkin[3]+19.6875*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 70.53065765632414*coeff[0]*fSkin[11]-51.4683177492095*coeff[0]*fEdge[11]-61.5234375*coeff[0]*fSkin[4]-56.6015625*coeff[0]*fEdge[4]+34.09975027401226*coeff[0]*fSkin[2]-34.09975027401226*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 70.53065765632414*coeff[0]*fSkin[13]-51.4683177492095*coeff[0]*fEdge[13]-61.5234375*coeff[0]*fSkin[5]-56.6015625*coeff[0]*fEdge[5]+34.09975027401226*coeff[0]*fSkin[3]-34.09975027401226*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-35.21807064562169*coeff[0]*fSkin[17])+35.21807064562169*coeff[0]*fEdge[17]+34.09975027401226*coeff[0]*fSkin[10]+34.09975027401226*coeff[0]*fEdge[10]-19.6875*coeff[0]*fSkin[6]+19.6875*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = (-70.6640625*coeff[0]*fSkin[7])-3.1640625*coeff[0]*fEdge[7]+31.316701275974*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 34.09975027401227*coeff[0]*fSkin[12]+34.09975027401227*coeff[0]*fEdge[12]-19.6875*coeff[0]*fSkin[8]+19.6875*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 34.09975027401227*coeff[0]*fSkin[15]+34.09975027401227*coeff[0]*fEdge[15]-19.6875*coeff[0]*fSkin[9]+19.6875*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 70.53065765632411*coeff[0]*fSkin[17]-51.46831774920949*coeff[0]*fEdge[17]-61.5234375*coeff[0]*fSkin[10]-56.6015625*coeff[0]*fEdge[10]+34.09975027401226*coeff[0]*fSkin[6]-34.09975027401226*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = (-70.6640625*coeff[0]*fSkin[11])-3.1640625*coeff[0]*fEdge[11]+31.31670127597403*coeff[0]*fSkin[4]+12.2543613688594*coeff[0]*fEdge[4]-12.57788237343632*coeff[0]*fSkin[2]+12.57788237343632*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = (-61.5234375*coeff[0]*fSkin[12])-56.6015625*coeff[0]*fEdge[12]+34.09975027401227*coeff[0]*fSkin[8]-34.09975027401227*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = (-70.6640625*coeff[0]*fSkin[13])-3.1640625*coeff[0]*fEdge[13]+31.31670127597403*coeff[0]*fSkin[5]+12.2543613688594*coeff[0]*fEdge[5]-12.57788237343632*coeff[0]*fSkin[3]+12.57788237343632*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 34.09975027401227*coeff[0]*fSkin[18]+34.09975027401227*coeff[0]*fEdge[18]-19.6875*coeff[0]*fSkin[14]+19.6875*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-61.5234375*coeff[0]*fSkin[15])-56.6015625*coeff[0]*fEdge[15]+34.09975027401227*coeff[0]*fSkin[9]-34.09975027401227*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 34.09975027401227*coeff[0]*fSkin[19]+34.09975027401227*coeff[0]*fEdge[19]-19.6875*coeff[0]*fSkin[16]+19.6875*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-70.6640625*coeff[0]*fSkin[17])-3.1640625*coeff[0]*fEdge[17]+31.316701275974*coeff[0]*fSkin[10]+12.2543613688594*coeff[0]*fEdge[10]-12.57788237343632*coeff[0]*fSkin[6]+12.57788237343632*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = (-61.5234375*coeff[0]*fSkin[18])-56.6015625*coeff[0]*fEdge[18]+34.09975027401227*coeff[0]*fSkin[14]-34.09975027401227*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = (-61.5234375*coeff[0]*fSkin[19])-56.6015625*coeff[0]*fEdge[19]+34.09975027401227*coeff[0]*fSkin[16]-34.09975027401227*coeff[0]*fEdge[16]; 

  boundSurf_incr[1] = (-19.06233990711463*coeff[0]*fSkin[7])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = (-19.06233990711463*coeff[0]*fSkin[11])-4.921875*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = (-19.06233990711463*coeff[0]*fSkin[13])-4.921875*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = (-73.828125*coeff[0]*fSkin[7])-19.06233990711463*coeff[0]*fSkin[1]; 
  boundSurf_incr[10] = (-19.06233990711463*coeff[0]*fSkin[17])-4.921875*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = (-73.828125*coeff[0]*fSkin[11])-19.06233990711463*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = -4.921875*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = (-73.828125*coeff[0]*fSkin[13])-19.06233990711463*coeff[0]*fSkin[5]; 
  boundSurf_incr[15] = -4.921875*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = (-73.828125*coeff[0]*fSkin[17])-19.06233990711463*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = -4.921875*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = -4.921875*coeff[0]*fSkin[19]; 

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

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_gyrokinetic_order6_boundary_surfx_1x2v_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[20] = {0.0}; 

  double edgeSurf_incr[20] = {0.0}; 
  double boundSurf_incr[20] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[7]+13.47910981360368*coeff[1]*fSkin[7]-24.90293657382598*coeff[0]*fSkin[7]-59.66213466261492*coeff[2]*fEdge[7]+13.47910981360368*coeff[1]*fEdge[7]+24.90293657382598*coeff[0]*fEdge[7]+65.46996195178933*fSkin[1]*coeff[2]+65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]-24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]-24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 110.8728999785158*coeff[2]*fSkin[7]+9.116253567204142*coeff[1]*fSkin[7]-49.87270631033365*coeff[0]*fSkin[7]-95.80279706881466*coeff[2]*fEdge[7]+37.57675250871956*coeff[1]*fEdge[7]+36.39359649672996*coeff[0]*fEdge[7]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]+65.46996195178933*fSkin[0]*coeff[2]-65.46996195178933*fEdge[0]*coeff[2]-11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]-23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]-24.11216465552189*coeff[0]*fSkin[0]+24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 59.66213466261492*coeff[2]*fSkin[11]+13.47910981360369*coeff[1]*fSkin[11]-24.90293657382598*coeff[0]*fSkin[11]-59.66213466261492*coeff[2]*fEdge[11]+13.47910981360369*coeff[1]*fEdge[11]+24.90293657382598*coeff[0]*fEdge[11]+65.46996195178933*coeff[2]*fSkin[4]+3.480291188652536*coeff[1]*fSkin[4]-24.11216465552189*coeff[0]*fSkin[4]+65.46996195178933*coeff[2]*fEdge[4]-3.480291188652536*coeff[1]*fEdge[4]-24.11216465552189*coeff[0]*fEdge[4]+37.79910015670014*coeff[2]*fSkin[2]-13.92116475461015*coeff[0]*fSkin[2]-37.79910015670014*coeff[2]*fEdge[2]+13.92116475461015*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 59.66213466261492*coeff[2]*fSkin[13]+13.47910981360369*coeff[1]*fSkin[13]-24.90293657382598*coeff[0]*fSkin[13]-59.66213466261492*coeff[2]*fEdge[13]+13.47910981360369*coeff[1]*fEdge[13]+24.90293657382598*coeff[0]*fEdge[13]+65.46996195178933*coeff[2]*fSkin[5]+3.480291188652536*coeff[1]*fSkin[5]-24.11216465552189*coeff[0]*fSkin[5]+65.46996195178933*coeff[2]*fEdge[5]-3.480291188652536*coeff[1]*fEdge[5]-24.11216465552189*coeff[0]*fEdge[5]+37.79910015670014*coeff[2]*fSkin[3]-13.92116475461015*coeff[0]*fSkin[3]-37.79910015670014*coeff[2]*fEdge[3]+13.92116475461015*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 110.8728999785159*coeff[2]*fSkin[11]+9.116253567204131*coeff[1]*fSkin[11]-49.87270631033366*coeff[0]*fSkin[11]-95.80279706881471*coeff[2]*fEdge[11]+37.57675250871954*coeff[1]*fEdge[11]+36.39359649672998*coeff[0]*fEdge[11]+115.3428423899306*coeff[2]*fSkin[4]-11.19493359006374*coeff[1]*fSkin[4]-43.5036398581567*coeff[0]*fSkin[4]+111.4517585502703*coeff[2]*fEdge[4]-23.25101591782468*coeff[1]*fEdge[4]-40.02334866950417*coeff[0]*fEdge[4]+65.46996195178933*coeff[2]*fSkin[2]-9.94368911043582*coeff[1]*fSkin[2]-24.11216465552189*coeff[0]*fSkin[2]-65.46996195178933*coeff[2]*fEdge[2]+9.94368911043582*coeff[1]*fEdge[2]+24.11216465552189*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 110.8728999785159*coeff[2]*fSkin[13]+9.116253567204131*coeff[1]*fSkin[13]-49.87270631033366*coeff[0]*fSkin[13]-95.80279706881471*coeff[2]*fEdge[13]+37.57675250871954*coeff[1]*fEdge[13]+36.39359649672998*coeff[0]*fEdge[13]+115.3428423899306*coeff[2]*fSkin[5]-11.19493359006374*coeff[1]*fSkin[5]-43.5036398581567*coeff[0]*fSkin[5]+111.4517585502703*coeff[2]*fEdge[5]-23.25101591782468*coeff[1]*fEdge[5]-40.02334866950417*coeff[0]*fEdge[5]+65.46996195178933*coeff[2]*fSkin[3]-9.94368911043582*coeff[1]*fSkin[3]-24.11216465552189*coeff[0]*fSkin[3]-65.46996195178933*coeff[2]*fEdge[3]+9.94368911043582*coeff[1]*fEdge[3]+24.11216465552189*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 59.66213466261492*coeff[2]*fSkin[17]+13.47910981360368*coeff[1]*fSkin[17]-24.90293657382598*coeff[0]*fSkin[17]-59.66213466261492*coeff[2]*fEdge[17]+13.47910981360368*coeff[1]*fEdge[17]+24.90293657382598*coeff[0]*fEdge[17]+65.46996195178933*coeff[2]*fSkin[10]+3.480291188652536*coeff[1]*fSkin[10]-24.11216465552189*coeff[0]*fSkin[10]+65.46996195178933*coeff[2]*fEdge[10]-3.480291188652536*coeff[1]*fEdge[10]-24.11216465552189*coeff[0]*fEdge[10]+37.79910015670014*coeff[2]*fSkin[6]-13.92116475461015*coeff[0]*fSkin[6]-37.79910015670014*coeff[2]*fEdge[6]+13.92116475461015*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 127.0160939089115*coeff[2]*fSkin[7]-24.97331339321913*coeff[1]*fSkin[7]-49.96703777993997*coeff[0]*fSkin[7]-68.64983631400692*coeff[2]*fEdge[7]+85.25372503202384*coeff[1]*fEdge[7]-2.237330049848051*coeff[0]*fEdge[7]+110.8728999785158*fSkin[1]*coeff[2]+95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]-22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]-8.665142023030945*coeff[0]*fEdge[1]-38.51174232458197*fSkin[0]*coeff[1]+38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = 65.46996195178934*coeff[2]*fSkin[12]+3.480291188652535*coeff[1]*fSkin[12]-24.1121646555219*coeff[0]*fSkin[12]+65.46996195178934*coeff[2]*fEdge[12]-3.480291188652535*coeff[1]*fEdge[12]-24.1121646555219*coeff[0]*fEdge[12]+37.79910015670014*coeff[2]*fSkin[8]-13.92116475461015*coeff[0]*fSkin[8]-37.79910015670014*coeff[2]*fEdge[8]+13.92116475461015*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = 65.46996195178934*coeff[2]*fSkin[15]+3.480291188652535*coeff[1]*fSkin[15]-24.1121646555219*coeff[0]*fSkin[15]+65.46996195178934*coeff[2]*fEdge[15]-3.480291188652535*coeff[1]*fEdge[15]-24.1121646555219*coeff[0]*fEdge[15]+37.79910015670014*coeff[2]*fSkin[9]-13.92116475461015*coeff[0]*fSkin[9]-37.79910015670014*coeff[2]*fEdge[9]+13.92116475461015*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 110.8728999785158*coeff[2]*fSkin[17]+9.116253567204142*coeff[1]*fSkin[17]-49.87270631033365*coeff[0]*fSkin[17]-95.80279706881466*coeff[2]*fEdge[17]+37.57675250871956*coeff[1]*fEdge[17]+36.39359649672996*coeff[0]*fEdge[17]+115.3428423899306*coeff[2]*fSkin[10]-11.19493359006374*coeff[1]*fSkin[10]-43.5036398581567*coeff[0]*fSkin[10]+111.4517585502703*coeff[2]*fEdge[10]-23.25101591782468*coeff[1]*fEdge[10]-40.02334866950417*coeff[0]*fEdge[10]+65.46996195178933*coeff[2]*fSkin[6]-9.94368911043582*coeff[1]*fSkin[6]-24.11216465552189*coeff[0]*fSkin[6]-65.46996195178933*coeff[2]*fEdge[6]+9.94368911043582*coeff[1]*fEdge[6]+24.11216465552189*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 127.0160939089115*coeff[2]*fSkin[11]-24.97331339321913*coeff[1]*fSkin[11]-49.96703777993997*coeff[0]*fSkin[11]-68.64983631400692*coeff[2]*fEdge[11]+85.25372503202384*coeff[1]*fEdge[11]-2.237330049848051*coeff[0]*fEdge[11]+110.8728999785159*coeff[2]*fSkin[4]-58.92212671485608*coeff[1]*fSkin[4]-22.14425183663463*coeff[0]*fSkin[4]+95.80279706881468*coeff[2]*fEdge[4]-74.48646207349731*coeff[1]*fEdge[4]-8.665142023030945*coeff[0]*fEdge[4]+59.66213466261491*coeff[2]*fSkin[2]-38.51174232458198*coeff[1]*fSkin[2]-8.893905919223561*coeff[0]*fSkin[2]-59.66213466261491*coeff[2]*fEdge[2]+38.51174232458198*coeff[1]*fEdge[2]+8.893905919223561*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 115.3428423899306*coeff[2]*fSkin[12]-11.19493359006374*coeff[1]*fSkin[12]-43.5036398581567*coeff[0]*fSkin[12]+111.4517585502703*coeff[2]*fEdge[12]-23.25101591782468*coeff[1]*fEdge[12]-40.02334866950417*coeff[0]*fEdge[12]+65.46996195178934*coeff[2]*fSkin[8]-9.94368911043582*coeff[1]*fSkin[8]-24.1121646555219*coeff[0]*fSkin[8]-65.46996195178934*coeff[2]*fEdge[8]+9.94368911043582*coeff[1]*fEdge[8]+24.1121646555219*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 127.0160939089115*coeff[2]*fSkin[13]-24.97331339321913*coeff[1]*fSkin[13]-49.96703777993997*coeff[0]*fSkin[13]-68.64983631400692*coeff[2]*fEdge[13]+85.25372503202384*coeff[1]*fEdge[13]-2.237330049848051*coeff[0]*fEdge[13]+110.8728999785159*coeff[2]*fSkin[5]-58.92212671485608*coeff[1]*fSkin[5]-22.14425183663463*coeff[0]*fSkin[5]+95.80279706881468*coeff[2]*fEdge[5]-74.48646207349731*coeff[1]*fEdge[5]-8.665142023030945*coeff[0]*fEdge[5]+59.66213466261491*coeff[2]*fSkin[3]-38.51174232458198*coeff[1]*fSkin[3]-8.893905919223561*coeff[0]*fSkin[3]-59.66213466261491*coeff[2]*fEdge[3]+38.51174232458198*coeff[1]*fEdge[3]+8.893905919223561*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = 65.46996195178934*coeff[2]*fSkin[18]+3.480291188652535*coeff[1]*fSkin[18]-24.1121646555219*coeff[0]*fSkin[18]+65.46996195178934*coeff[2]*fEdge[18]-3.480291188652535*coeff[1]*fEdge[18]-24.1121646555219*coeff[0]*fEdge[18]+37.79910015670014*coeff[2]*fSkin[14]-13.92116475461015*coeff[0]*fSkin[14]-37.79910015670014*coeff[2]*fEdge[14]+13.92116475461015*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 115.3428423899306*coeff[2]*fSkin[15]-11.19493359006374*coeff[1]*fSkin[15]-43.5036398581567*coeff[0]*fSkin[15]+111.4517585502703*coeff[2]*fEdge[15]-23.25101591782468*coeff[1]*fEdge[15]-40.02334866950417*coeff[0]*fEdge[15]+65.46996195178934*coeff[2]*fSkin[9]-9.94368911043582*coeff[1]*fSkin[9]-24.1121646555219*coeff[0]*fSkin[9]-65.46996195178934*coeff[2]*fEdge[9]+9.94368911043582*coeff[1]*fEdge[9]+24.1121646555219*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = 65.46996195178934*coeff[2]*fSkin[19]+3.480291188652535*coeff[1]*fSkin[19]-24.1121646555219*coeff[0]*fSkin[19]+65.46996195178934*coeff[2]*fEdge[19]-3.480291188652535*coeff[1]*fEdge[19]-24.1121646555219*coeff[0]*fEdge[19]+37.79910015670014*coeff[2]*fSkin[16]-13.92116475461015*coeff[0]*fSkin[16]-37.79910015670014*coeff[2]*fEdge[16]+13.92116475461015*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 127.0160939089115*coeff[2]*fSkin[17]-24.97331339321913*coeff[1]*fSkin[17]-49.96703777993997*coeff[0]*fSkin[17]-68.64983631400692*coeff[2]*fEdge[17]+85.25372503202384*coeff[1]*fEdge[17]-2.237330049848051*coeff[0]*fEdge[17]+110.8728999785158*coeff[2]*fSkin[10]-58.92212671485613*coeff[1]*fSkin[10]-22.14425183663463*coeff[0]*fSkin[10]+95.80279706881466*coeff[2]*fEdge[10]-74.48646207349736*coeff[1]*fEdge[10]-8.665142023030945*coeff[0]*fEdge[10]+59.66213466261492*coeff[2]*fSkin[6]-38.51174232458197*coeff[1]*fSkin[6]-8.893905919223567*coeff[0]*fSkin[6]-59.66213466261492*coeff[2]*fEdge[6]+38.51174232458197*coeff[1]*fEdge[6]+8.893905919223567*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 115.3428423899306*coeff[2]*fSkin[18]-11.19493359006374*coeff[1]*fSkin[18]-43.5036398581567*coeff[0]*fSkin[18]+111.4517585502703*coeff[2]*fEdge[18]-23.25101591782468*coeff[1]*fEdge[18]-40.02334866950417*coeff[0]*fEdge[18]+65.46996195178934*coeff[2]*fSkin[14]-9.94368911043582*coeff[1]*fSkin[14]-24.1121646555219*coeff[0]*fSkin[14]-65.46996195178934*coeff[2]*fEdge[14]+9.94368911043582*coeff[1]*fEdge[14]+24.1121646555219*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 115.3428423899306*coeff[2]*fSkin[19]-11.19493359006374*coeff[1]*fSkin[19]-43.5036398581567*coeff[0]*fSkin[19]+111.4517585502703*coeff[2]*fEdge[19]-23.25101591782468*coeff[1]*fEdge[19]-40.02334866950417*coeff[0]*fEdge[19]+65.46996195178934*coeff[2]*fSkin[16]-9.94368911043582*coeff[1]*fSkin[16]-24.1121646555219*coeff[0]*fSkin[16]-65.46996195178934*coeff[2]*fEdge[16]+9.94368911043582*coeff[1]*fEdge[16]+24.1121646555219*coeff[0]*fEdge[16]; 

  boundSurf_incr[0] = 6.960582377305072*coeff[1]*fSkin[1]-26.95821962720737*coeff[1]*fSkin[7]; 
  boundSurf_incr[1] = (-15.07010290970117*coeff[2]*fSkin[7])+46.69300607592371*coeff[1]*fSkin[7]+13.47910981360368*coeff[0]*fSkin[7]+3.891083839660308*fSkin[1]*coeff[2]-12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 6.960582377305072*coeff[1]*fSkin[4]-26.95821962720738*coeff[1]*fSkin[11]; 
  boundSurf_incr[3] = 6.960582377305072*coeff[1]*fSkin[5]-26.95821962720738*coeff[1]*fSkin[13]; 
  boundSurf_incr[4] = (-15.07010290970118*coeff[2]*fSkin[11])+46.69300607592368*coeff[1]*fSkin[11]+13.47910981360369*coeff[0]*fSkin[11]+3.891083839660308*coeff[2]*fSkin[4]-12.05608232776094*coeff[1]*fSkin[4]-3.480291188652536*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = (-15.07010290970118*coeff[2]*fSkin[13])+46.69300607592368*coeff[1]*fSkin[13]+13.47910981360369*coeff[0]*fSkin[13]+3.891083839660308*coeff[2]*fSkin[5]-12.05608232776094*coeff[1]*fSkin[5]-3.480291188652536*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 6.960582377305072*coeff[1]*fSkin[10]-26.95821962720737*coeff[1]*fSkin[17]; 
  boundSurf_incr[7] = 58.36625759490461*coeff[2]*fSkin[7]-60.28041163880471*coeff[1]*fSkin[7]-52.20436782978803*coeff[0]*fSkin[7]-15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]+13.47910981360368*coeff[0]*fSkin[1]; 
  boundSurf_incr[8] = 6.960582377305072*coeff[1]*fSkin[12]; 
  boundSurf_incr[9] = 6.960582377305072*coeff[1]*fSkin[15]; 
  boundSurf_incr[10] = (-15.07010290970117*coeff[2]*fSkin[17])+46.69300607592371*coeff[1]*fSkin[17]+13.47910981360368*coeff[0]*fSkin[17]+3.891083839660308*coeff[2]*fSkin[10]-12.05608232776094*coeff[1]*fSkin[10]-3.480291188652536*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 58.36625759490461*coeff[2]*fSkin[11]-60.28041163880471*coeff[1]*fSkin[11]-52.20436782978803*coeff[0]*fSkin[11]-15.07010290970118*coeff[2]*fSkin[4]+15.56433535864123*coeff[1]*fSkin[4]+13.47910981360369*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 3.891083839660308*coeff[2]*fSkin[12]-12.05608232776094*coeff[1]*fSkin[12]-3.480291188652536*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 58.36625759490461*coeff[2]*fSkin[13]-60.28041163880471*coeff[1]*fSkin[13]-52.20436782978803*coeff[0]*fSkin[13]-15.07010290970118*coeff[2]*fSkin[5]+15.56433535864123*coeff[1]*fSkin[5]+13.47910981360369*coeff[0]*fSkin[5]; 
  boundSurf_incr[14] = 6.960582377305072*coeff[1]*fSkin[18]; 
  boundSurf_incr[15] = 3.891083839660308*coeff[2]*fSkin[15]-12.05608232776094*coeff[1]*fSkin[15]-3.480291188652536*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 6.960582377305072*coeff[1]*fSkin[19]; 
  boundSurf_incr[17] = 58.36625759490461*coeff[2]*fSkin[17]-60.28041163880471*coeff[1]*fSkin[17]-52.20436782978803*coeff[0]*fSkin[17]-15.07010290970117*coeff[2]*fSkin[10]+15.56433535864123*coeff[1]*fSkin[10]+13.47910981360368*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 3.891083839660308*coeff[2]*fSkin[18]-12.05608232776094*coeff[1]*fSkin[18]-3.480291188652536*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 3.891083839660308*coeff[2]*fSkin[19]-12.05608232776094*coeff[1]*fSkin[19]-3.480291188652536*coeff[0]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[7]-13.47910981360368*coeff[1]*fSkin[7]-24.90293657382598*coeff[0]*fSkin[7]-59.66213466261492*coeff[2]*fEdge[7]-13.47910981360368*coeff[1]*fEdge[7]+24.90293657382598*coeff[0]*fEdge[7]-65.46996195178933*fSkin[1]*coeff[2]-65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]+24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]+24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-110.8728999785158*coeff[2]*fSkin[7])+9.116253567204142*coeff[1]*fSkin[7]+49.87270631033365*coeff[0]*fSkin[7]+95.80279706881466*coeff[2]*fEdge[7]+37.57675250871956*coeff[1]*fEdge[7]-36.39359649672996*coeff[0]*fEdge[7]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]-65.46996195178933*fSkin[0]*coeff[2]+65.46996195178933*fEdge[0]*coeff[2]+11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]+23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]+24.11216465552189*coeff[0]*fSkin[0]-24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 59.66213466261492*coeff[2]*fSkin[11]-13.47910981360369*coeff[1]*fSkin[11]-24.90293657382598*coeff[0]*fSkin[11]-59.66213466261492*coeff[2]*fEdge[11]-13.47910981360369*coeff[1]*fEdge[11]+24.90293657382598*coeff[0]*fEdge[11]-65.46996195178933*coeff[2]*fSkin[4]+3.480291188652536*coeff[1]*fSkin[4]+24.11216465552189*coeff[0]*fSkin[4]-65.46996195178933*coeff[2]*fEdge[4]-3.480291188652536*coeff[1]*fEdge[4]+24.11216465552189*coeff[0]*fEdge[4]+37.79910015670014*coeff[2]*fSkin[2]-13.92116475461015*coeff[0]*fSkin[2]-37.79910015670014*coeff[2]*fEdge[2]+13.92116475461015*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 59.66213466261492*coeff[2]*fSkin[13]-13.47910981360369*coeff[1]*fSkin[13]-24.90293657382598*coeff[0]*fSkin[13]-59.66213466261492*coeff[2]*fEdge[13]-13.47910981360369*coeff[1]*fEdge[13]+24.90293657382598*coeff[0]*fEdge[13]-65.46996195178933*coeff[2]*fSkin[5]+3.480291188652536*coeff[1]*fSkin[5]+24.11216465552189*coeff[0]*fSkin[5]-65.46996195178933*coeff[2]*fEdge[5]-3.480291188652536*coeff[1]*fEdge[5]+24.11216465552189*coeff[0]*fEdge[5]+37.79910015670014*coeff[2]*fSkin[3]-13.92116475461015*coeff[0]*fSkin[3]-37.79910015670014*coeff[2]*fEdge[3]+13.92116475461015*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-110.8728999785159*coeff[2]*fSkin[11])+9.116253567204131*coeff[1]*fSkin[11]+49.87270631033366*coeff[0]*fSkin[11]+95.80279706881471*coeff[2]*fEdge[11]+37.57675250871954*coeff[1]*fEdge[11]-36.39359649672998*coeff[0]*fEdge[11]+115.3428423899306*coeff[2]*fSkin[4]+11.19493359006374*coeff[1]*fSkin[4]-43.5036398581567*coeff[0]*fSkin[4]+111.4517585502703*coeff[2]*fEdge[4]+23.25101591782468*coeff[1]*fEdge[4]-40.02334866950417*coeff[0]*fEdge[4]-65.46996195178933*coeff[2]*fSkin[2]-9.94368911043582*coeff[1]*fSkin[2]+24.11216465552189*coeff[0]*fSkin[2]+65.46996195178933*coeff[2]*fEdge[2]+9.94368911043582*coeff[1]*fEdge[2]-24.11216465552189*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = (-110.8728999785159*coeff[2]*fSkin[13])+9.116253567204131*coeff[1]*fSkin[13]+49.87270631033366*coeff[0]*fSkin[13]+95.80279706881471*coeff[2]*fEdge[13]+37.57675250871954*coeff[1]*fEdge[13]-36.39359649672998*coeff[0]*fEdge[13]+115.3428423899306*coeff[2]*fSkin[5]+11.19493359006374*coeff[1]*fSkin[5]-43.5036398581567*coeff[0]*fSkin[5]+111.4517585502703*coeff[2]*fEdge[5]+23.25101591782468*coeff[1]*fEdge[5]-40.02334866950417*coeff[0]*fEdge[5]-65.46996195178933*coeff[2]*fSkin[3]-9.94368911043582*coeff[1]*fSkin[3]+24.11216465552189*coeff[0]*fSkin[3]+65.46996195178933*coeff[2]*fEdge[3]+9.94368911043582*coeff[1]*fEdge[3]-24.11216465552189*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 59.66213466261492*coeff[2]*fSkin[17]-13.47910981360368*coeff[1]*fSkin[17]-24.90293657382598*coeff[0]*fSkin[17]-59.66213466261492*coeff[2]*fEdge[17]-13.47910981360368*coeff[1]*fEdge[17]+24.90293657382598*coeff[0]*fEdge[17]-65.46996195178933*coeff[2]*fSkin[10]+3.480291188652536*coeff[1]*fSkin[10]+24.11216465552189*coeff[0]*fSkin[10]-65.46996195178933*coeff[2]*fEdge[10]-3.480291188652536*coeff[1]*fEdge[10]+24.11216465552189*coeff[0]*fEdge[10]+37.79910015670014*coeff[2]*fSkin[6]-13.92116475461015*coeff[0]*fSkin[6]-37.79910015670014*coeff[2]*fEdge[6]+13.92116475461015*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 127.0160939089115*coeff[2]*fSkin[7]+24.97331339321913*coeff[1]*fSkin[7]-49.96703777993997*coeff[0]*fSkin[7]-68.64983631400692*coeff[2]*fEdge[7]-85.25372503202384*coeff[1]*fEdge[7]-2.237330049848051*coeff[0]*fEdge[7]-110.8728999785158*fSkin[1]*coeff[2]-95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]+22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]+8.665142023030945*coeff[0]*fEdge[1]+38.51174232458197*fSkin[0]*coeff[1]-38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 
  edgeSurf_incr[8] = (-65.46996195178934*coeff[2]*fSkin[12])+3.480291188652535*coeff[1]*fSkin[12]+24.1121646555219*coeff[0]*fSkin[12]-65.46996195178934*coeff[2]*fEdge[12]-3.480291188652535*coeff[1]*fEdge[12]+24.1121646555219*coeff[0]*fEdge[12]+37.79910015670014*coeff[2]*fSkin[8]-13.92116475461015*coeff[0]*fSkin[8]-37.79910015670014*coeff[2]*fEdge[8]+13.92116475461015*coeff[0]*fEdge[8]; 
  edgeSurf_incr[9] = (-65.46996195178934*coeff[2]*fSkin[15])+3.480291188652535*coeff[1]*fSkin[15]+24.1121646555219*coeff[0]*fSkin[15]-65.46996195178934*coeff[2]*fEdge[15]-3.480291188652535*coeff[1]*fEdge[15]+24.1121646555219*coeff[0]*fEdge[15]+37.79910015670014*coeff[2]*fSkin[9]-13.92116475461015*coeff[0]*fSkin[9]-37.79910015670014*coeff[2]*fEdge[9]+13.92116475461015*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-110.8728999785158*coeff[2]*fSkin[17])+9.116253567204142*coeff[1]*fSkin[17]+49.87270631033365*coeff[0]*fSkin[17]+95.80279706881466*coeff[2]*fEdge[17]+37.57675250871956*coeff[1]*fEdge[17]-36.39359649672996*coeff[0]*fEdge[17]+115.3428423899306*coeff[2]*fSkin[10]+11.19493359006374*coeff[1]*fSkin[10]-43.5036398581567*coeff[0]*fSkin[10]+111.4517585502703*coeff[2]*fEdge[10]+23.25101591782468*coeff[1]*fEdge[10]-40.02334866950417*coeff[0]*fEdge[10]-65.46996195178933*coeff[2]*fSkin[6]-9.94368911043582*coeff[1]*fSkin[6]+24.11216465552189*coeff[0]*fSkin[6]+65.46996195178933*coeff[2]*fEdge[6]+9.94368911043582*coeff[1]*fEdge[6]-24.11216465552189*coeff[0]*fEdge[6]; 
  edgeSurf_incr[11] = 127.0160939089115*coeff[2]*fSkin[11]+24.97331339321913*coeff[1]*fSkin[11]-49.96703777993997*coeff[0]*fSkin[11]-68.64983631400692*coeff[2]*fEdge[11]-85.25372503202384*coeff[1]*fEdge[11]-2.237330049848051*coeff[0]*fEdge[11]-110.8728999785159*coeff[2]*fSkin[4]-58.92212671485608*coeff[1]*fSkin[4]+22.14425183663463*coeff[0]*fSkin[4]-95.80279706881468*coeff[2]*fEdge[4]-74.48646207349731*coeff[1]*fEdge[4]+8.665142023030945*coeff[0]*fEdge[4]+59.66213466261491*coeff[2]*fSkin[2]+38.51174232458198*coeff[1]*fSkin[2]-8.893905919223561*coeff[0]*fSkin[2]-59.66213466261491*coeff[2]*fEdge[2]-38.51174232458198*coeff[1]*fEdge[2]+8.893905919223561*coeff[0]*fEdge[2]; 
  edgeSurf_incr[12] = 115.3428423899306*coeff[2]*fSkin[12]+11.19493359006374*coeff[1]*fSkin[12]-43.5036398581567*coeff[0]*fSkin[12]+111.4517585502703*coeff[2]*fEdge[12]+23.25101591782468*coeff[1]*fEdge[12]-40.02334866950417*coeff[0]*fEdge[12]-65.46996195178934*coeff[2]*fSkin[8]-9.94368911043582*coeff[1]*fSkin[8]+24.1121646555219*coeff[0]*fSkin[8]+65.46996195178934*coeff[2]*fEdge[8]+9.94368911043582*coeff[1]*fEdge[8]-24.1121646555219*coeff[0]*fEdge[8]; 
  edgeSurf_incr[13] = 127.0160939089115*coeff[2]*fSkin[13]+24.97331339321913*coeff[1]*fSkin[13]-49.96703777993997*coeff[0]*fSkin[13]-68.64983631400692*coeff[2]*fEdge[13]-85.25372503202384*coeff[1]*fEdge[13]-2.237330049848051*coeff[0]*fEdge[13]-110.8728999785159*coeff[2]*fSkin[5]-58.92212671485608*coeff[1]*fSkin[5]+22.14425183663463*coeff[0]*fSkin[5]-95.80279706881468*coeff[2]*fEdge[5]-74.48646207349731*coeff[1]*fEdge[5]+8.665142023030945*coeff[0]*fEdge[5]+59.66213466261491*coeff[2]*fSkin[3]+38.51174232458198*coeff[1]*fSkin[3]-8.893905919223561*coeff[0]*fSkin[3]-59.66213466261491*coeff[2]*fEdge[3]-38.51174232458198*coeff[1]*fEdge[3]+8.893905919223561*coeff[0]*fEdge[3]; 
  edgeSurf_incr[14] = (-65.46996195178934*coeff[2]*fSkin[18])+3.480291188652535*coeff[1]*fSkin[18]+24.1121646555219*coeff[0]*fSkin[18]-65.46996195178934*coeff[2]*fEdge[18]-3.480291188652535*coeff[1]*fEdge[18]+24.1121646555219*coeff[0]*fEdge[18]+37.79910015670014*coeff[2]*fSkin[14]-13.92116475461015*coeff[0]*fSkin[14]-37.79910015670014*coeff[2]*fEdge[14]+13.92116475461015*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 115.3428423899306*coeff[2]*fSkin[15]+11.19493359006374*coeff[1]*fSkin[15]-43.5036398581567*coeff[0]*fSkin[15]+111.4517585502703*coeff[2]*fEdge[15]+23.25101591782468*coeff[1]*fEdge[15]-40.02334866950417*coeff[0]*fEdge[15]-65.46996195178934*coeff[2]*fSkin[9]-9.94368911043582*coeff[1]*fSkin[9]+24.1121646555219*coeff[0]*fSkin[9]+65.46996195178934*coeff[2]*fEdge[9]+9.94368911043582*coeff[1]*fEdge[9]-24.1121646555219*coeff[0]*fEdge[9]; 
  edgeSurf_incr[16] = (-65.46996195178934*coeff[2]*fSkin[19])+3.480291188652535*coeff[1]*fSkin[19]+24.1121646555219*coeff[0]*fSkin[19]-65.46996195178934*coeff[2]*fEdge[19]-3.480291188652535*coeff[1]*fEdge[19]+24.1121646555219*coeff[0]*fEdge[19]+37.79910015670014*coeff[2]*fSkin[16]-13.92116475461015*coeff[0]*fSkin[16]-37.79910015670014*coeff[2]*fEdge[16]+13.92116475461015*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 127.0160939089115*coeff[2]*fSkin[17]+24.97331339321913*coeff[1]*fSkin[17]-49.96703777993997*coeff[0]*fSkin[17]-68.64983631400692*coeff[2]*fEdge[17]-85.25372503202384*coeff[1]*fEdge[17]-2.237330049848051*coeff[0]*fEdge[17]-110.8728999785158*coeff[2]*fSkin[10]-58.92212671485613*coeff[1]*fSkin[10]+22.14425183663463*coeff[0]*fSkin[10]-95.80279706881466*coeff[2]*fEdge[10]-74.48646207349736*coeff[1]*fEdge[10]+8.665142023030945*coeff[0]*fEdge[10]+59.66213466261492*coeff[2]*fSkin[6]+38.51174232458197*coeff[1]*fSkin[6]-8.893905919223567*coeff[0]*fSkin[6]-59.66213466261492*coeff[2]*fEdge[6]-38.51174232458197*coeff[1]*fEdge[6]+8.893905919223567*coeff[0]*fEdge[6]; 
  edgeSurf_incr[18] = 115.3428423899306*coeff[2]*fSkin[18]+11.19493359006374*coeff[1]*fSkin[18]-43.5036398581567*coeff[0]*fSkin[18]+111.4517585502703*coeff[2]*fEdge[18]+23.25101591782468*coeff[1]*fEdge[18]-40.02334866950417*coeff[0]*fEdge[18]-65.46996195178934*coeff[2]*fSkin[14]-9.94368911043582*coeff[1]*fSkin[14]+24.1121646555219*coeff[0]*fSkin[14]+65.46996195178934*coeff[2]*fEdge[14]+9.94368911043582*coeff[1]*fEdge[14]-24.1121646555219*coeff[0]*fEdge[14]; 
  edgeSurf_incr[19] = 115.3428423899306*coeff[2]*fSkin[19]+11.19493359006374*coeff[1]*fSkin[19]-43.5036398581567*coeff[0]*fSkin[19]+111.4517585502703*coeff[2]*fEdge[19]+23.25101591782468*coeff[1]*fEdge[19]-40.02334866950417*coeff[0]*fEdge[19]-65.46996195178934*coeff[2]*fSkin[16]-9.94368911043582*coeff[1]*fSkin[16]+24.1121646555219*coeff[0]*fSkin[16]+65.46996195178934*coeff[2]*fEdge[16]+9.94368911043582*coeff[1]*fEdge[16]-24.1121646555219*coeff[0]*fEdge[16]; 

  boundSurf_incr[0] = 26.95821962720737*coeff[1]*fSkin[7]+6.960582377305072*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 15.07010290970117*coeff[2]*fSkin[7]+46.69300607592371*coeff[1]*fSkin[7]-13.47910981360368*coeff[0]*fSkin[7]+3.891083839660308*fSkin[1]*coeff[2]+12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 26.95821962720738*coeff[1]*fSkin[11]+6.960582377305072*coeff[1]*fSkin[4]; 
  boundSurf_incr[3] = 26.95821962720738*coeff[1]*fSkin[13]+6.960582377305072*coeff[1]*fSkin[5]; 
  boundSurf_incr[4] = 15.07010290970118*coeff[2]*fSkin[11]+46.69300607592368*coeff[1]*fSkin[11]-13.47910981360369*coeff[0]*fSkin[11]+3.891083839660308*coeff[2]*fSkin[4]+12.05608232776094*coeff[1]*fSkin[4]-3.480291188652536*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 15.07010290970118*coeff[2]*fSkin[13]+46.69300607592368*coeff[1]*fSkin[13]-13.47910981360369*coeff[0]*fSkin[13]+3.891083839660308*coeff[2]*fSkin[5]+12.05608232776094*coeff[1]*fSkin[5]-3.480291188652536*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 26.95821962720737*coeff[1]*fSkin[17]+6.960582377305072*coeff[1]*fSkin[10]; 
  boundSurf_incr[7] = 58.36625759490461*coeff[2]*fSkin[7]+60.28041163880471*coeff[1]*fSkin[7]-52.20436782978803*coeff[0]*fSkin[7]+15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]-13.47910981360368*coeff[0]*fSkin[1]; 
  boundSurf_incr[8] = 6.960582377305072*coeff[1]*fSkin[12]; 
  boundSurf_incr[9] = 6.960582377305072*coeff[1]*fSkin[15]; 
  boundSurf_incr[10] = 15.07010290970117*coeff[2]*fSkin[17]+46.69300607592371*coeff[1]*fSkin[17]-13.47910981360368*coeff[0]*fSkin[17]+3.891083839660308*coeff[2]*fSkin[10]+12.05608232776094*coeff[1]*fSkin[10]-3.480291188652536*coeff[0]*fSkin[10]; 
  boundSurf_incr[11] = 58.36625759490461*coeff[2]*fSkin[11]+60.28041163880471*coeff[1]*fSkin[11]-52.20436782978803*coeff[0]*fSkin[11]+15.07010290970118*coeff[2]*fSkin[4]+15.56433535864123*coeff[1]*fSkin[4]-13.47910981360369*coeff[0]*fSkin[4]; 
  boundSurf_incr[12] = 3.891083839660308*coeff[2]*fSkin[12]+12.05608232776094*coeff[1]*fSkin[12]-3.480291188652536*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 58.36625759490461*coeff[2]*fSkin[13]+60.28041163880471*coeff[1]*fSkin[13]-52.20436782978803*coeff[0]*fSkin[13]+15.07010290970118*coeff[2]*fSkin[5]+15.56433535864123*coeff[1]*fSkin[5]-13.47910981360369*coeff[0]*fSkin[5]; 
  boundSurf_incr[14] = 6.960582377305072*coeff[1]*fSkin[18]; 
  boundSurf_incr[15] = 3.891083839660308*coeff[2]*fSkin[15]+12.05608232776094*coeff[1]*fSkin[15]-3.480291188652536*coeff[0]*fSkin[15]; 
  boundSurf_incr[16] = 6.960582377305072*coeff[1]*fSkin[19]; 
  boundSurf_incr[17] = 58.36625759490461*coeff[2]*fSkin[17]+60.28041163880471*coeff[1]*fSkin[17]-52.20436782978803*coeff[0]*fSkin[17]+15.07010290970117*coeff[2]*fSkin[10]+15.56433535864123*coeff[1]*fSkin[10]-13.47910981360368*coeff[0]*fSkin[10]; 
  boundSurf_incr[18] = 3.891083839660308*coeff[2]*fSkin[18]+12.05608232776094*coeff[1]*fSkin[18]-3.480291188652536*coeff[0]*fSkin[18]; 
  boundSurf_incr[19] = 3.891083839660308*coeff[2]*fSkin[19]+12.05608232776094*coeff[1]*fSkin[19]-3.480291188652536*coeff[0]*fSkin[19]; 

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

  }

  return 0.;
}

