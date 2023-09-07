#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_ser_p2_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[3] = {0.0}; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[2])+35.21807064562169*coeff[0]*fEdge[2]-34.09975027401226*coeff[0]*fSkin[1]-34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-70.53065765632411*coeff[0]*fSkin[2])+51.46831774920949*coeff[0]*fEdge[2]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]-34.09975027401226*coeff[0]*fSkin[0]+34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-70.6640625*coeff[0]*fSkin[2])-3.1640625*coeff[0]*fEdge[2]-31.316701275974*coeff[0]*fSkin[1]-12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = 19.06233990711463*coeff[0]*fSkin[2]-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 19.06233990711463*coeff[0]*fSkin[1]-73.828125*coeff[0]*fSkin[2]; 

  } else { 

  edgeSurf_incr[0] = (-35.21807064562169*coeff[0]*fSkin[2])+35.21807064562169*coeff[0]*fEdge[2]+34.09975027401226*coeff[0]*fSkin[1]+34.09975027401226*coeff[0]*fEdge[1]-19.6875*coeff[0]*fSkin[0]+19.6875*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 70.53065765632411*coeff[0]*fSkin[2]-51.46831774920949*coeff[0]*fEdge[2]-61.5234375*coeff[0]*fSkin[1]-56.6015625*coeff[0]*fEdge[1]+34.09975027401226*coeff[0]*fSkin[0]-34.09975027401226*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-70.6640625*coeff[0]*fSkin[2])-3.1640625*coeff[0]*fEdge[2]+31.316701275974*coeff[0]*fSkin[1]+12.2543613688594*coeff[0]*fEdge[1]-12.57788237343632*coeff[0]*fSkin[0]+12.57788237343632*coeff[0]*fEdge[0]; 

  boundSurf_incr[1] = (-19.06233990711463*coeff[0]*fSkin[2])-4.921875*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-73.828125*coeff[0]*fSkin[2])-19.06233990711463*coeff[0]*fSkin[1]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order6_boundary_surfx_1x_ser_p2_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],6.);

  double vol_incr[3] = {0.0}; 

  double edgeSurf_incr[3] = {0.0}; 
  double boundSurf_incr[3] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[2]+13.47910981360368*coeff[1]*fSkin[2]-24.90293657382598*coeff[0]*fSkin[2]-59.66213466261492*coeff[2]*fEdge[2]+13.47910981360368*coeff[1]*fEdge[2]+24.90293657382598*coeff[0]*fEdge[2]+65.46996195178933*fSkin[1]*coeff[2]+65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]-24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]-24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 110.8728999785158*coeff[2]*fSkin[2]+9.116253567204142*coeff[1]*fSkin[2]-49.87270631033365*coeff[0]*fSkin[2]-95.80279706881466*coeff[2]*fEdge[2]+37.57675250871956*coeff[1]*fEdge[2]+36.39359649672996*coeff[0]*fEdge[2]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]+65.46996195178933*fSkin[0]*coeff[2]-65.46996195178933*fEdge[0]*coeff[2]-11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]-23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]-24.11216465552189*coeff[0]*fSkin[0]+24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 127.0160939089115*coeff[2]*fSkin[2]-24.97331339321913*coeff[1]*fSkin[2]-49.96703777993997*coeff[0]*fSkin[2]-68.64983631400692*coeff[2]*fEdge[2]+85.25372503202384*coeff[1]*fEdge[2]-2.237330049848051*coeff[0]*fEdge[2]+110.8728999785158*fSkin[1]*coeff[2]+95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]-22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]-8.665142023030945*coeff[0]*fEdge[1]-38.51174232458197*fSkin[0]*coeff[1]+38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = 6.960582377305072*coeff[1]*fSkin[1]-26.95821962720737*coeff[1]*fSkin[2]; 
  boundSurf_incr[1] = (-15.07010290970117*coeff[2]*fSkin[2])+46.69300607592371*coeff[1]*fSkin[2]+13.47910981360368*coeff[0]*fSkin[2]+3.891083839660308*fSkin[1]*coeff[2]-12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 58.36625759490461*coeff[2]*fSkin[2]-60.28041163880471*coeff[1]*fSkin[2]-52.20436782978803*coeff[0]*fSkin[2]-15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]+13.47910981360368*coeff[0]*fSkin[1]; 

  } else { 

  edgeSurf_incr[0] = 59.66213466261492*coeff[2]*fSkin[2]-13.47910981360368*coeff[1]*fSkin[2]-24.90293657382598*coeff[0]*fSkin[2]-59.66213466261492*coeff[2]*fEdge[2]-13.47910981360368*coeff[1]*fEdge[2]+24.90293657382598*coeff[0]*fEdge[2]-65.46996195178933*fSkin[1]*coeff[2]-65.46996195178933*fEdge[1]*coeff[2]+37.79910015670014*fSkin[0]*coeff[2]-37.79910015670014*fEdge[0]*coeff[2]+3.480291188652536*coeff[1]*fSkin[1]+24.11216465552189*coeff[0]*fSkin[1]-3.480291188652536*coeff[1]*fEdge[1]+24.11216465552189*coeff[0]*fEdge[1]-13.92116475461015*coeff[0]*fSkin[0]+13.92116475461015*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-110.8728999785158*coeff[2]*fSkin[2])+9.116253567204142*coeff[1]*fSkin[2]+49.87270631033365*coeff[0]*fSkin[2]+95.80279706881466*coeff[2]*fEdge[2]+37.57675250871956*coeff[1]*fEdge[2]-36.39359649672996*coeff[0]*fEdge[2]+115.3428423899306*fSkin[1]*coeff[2]+111.4517585502703*fEdge[1]*coeff[2]-65.46996195178933*fSkin[0]*coeff[2]+65.46996195178933*fEdge[0]*coeff[2]+11.19493359006374*coeff[1]*fSkin[1]-43.5036398581567*coeff[0]*fSkin[1]+23.25101591782468*coeff[1]*fEdge[1]-40.02334866950417*coeff[0]*fEdge[1]-9.94368911043582*fSkin[0]*coeff[1]+9.94368911043582*fEdge[0]*coeff[1]+24.11216465552189*coeff[0]*fSkin[0]-24.11216465552189*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 127.0160939089115*coeff[2]*fSkin[2]+24.97331339321913*coeff[1]*fSkin[2]-49.96703777993997*coeff[0]*fSkin[2]-68.64983631400692*coeff[2]*fEdge[2]-85.25372503202384*coeff[1]*fEdge[2]-2.237330049848051*coeff[0]*fEdge[2]-110.8728999785158*fSkin[1]*coeff[2]-95.80279706881466*fEdge[1]*coeff[2]+59.66213466261492*fSkin[0]*coeff[2]-59.66213466261492*fEdge[0]*coeff[2]-58.92212671485613*coeff[1]*fSkin[1]+22.14425183663463*coeff[0]*fSkin[1]-74.48646207349736*coeff[1]*fEdge[1]+8.665142023030945*coeff[0]*fEdge[1]+38.51174232458197*fSkin[0]*coeff[1]-38.51174232458197*fEdge[0]*coeff[1]-8.893905919223567*coeff[0]*fSkin[0]+8.893905919223567*coeff[0]*fEdge[0]; 

  boundSurf_incr[0] = 26.95821962720737*coeff[1]*fSkin[2]+6.960582377305072*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 15.07010290970117*coeff[2]*fSkin[2]+46.69300607592371*coeff[1]*fSkin[2]-13.47910981360368*coeff[0]*fSkin[2]+3.891083839660308*fSkin[1]*coeff[2]+12.05608232776094*coeff[1]*fSkin[1]-3.480291188652536*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = 58.36625759490461*coeff[2]*fSkin[2]+60.28041163880471*coeff[1]*fSkin[2]-52.20436782978803*coeff[0]*fSkin[2]+15.07010290970117*fSkin[1]*coeff[2]+15.56433535864123*coeff[1]*fSkin[1]-13.47910981360368*coeff[0]*fSkin[1]; 

  out[0] += (vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += (vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += (vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 

  }

  return 0.;
}

