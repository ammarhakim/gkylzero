#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[0]*fSkin[4]+1.623797632095822*coeff[0]*fEdge[4]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 1.623797632095822*coeff[0]*fSkin[5]+1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*fSkin[4]+2.0625*coeff[0]*fEdge[4]+1.623797632095822*coeff[0]*fSkin[2]-1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]+1.623797632095822*coeff[0]*fSkin[3]-1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = 1.623797632095822*coeff[0]*fSkin[7]+1.623797632095822*coeff[0]*fEdge[7]+0.9375*coeff[0]*fSkin[6]-0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*fSkin[7]+2.0625*coeff[0]*fEdge[7]+1.623797632095822*coeff[0]*fSkin[6]-1.623797632095822*coeff[0]*fEdge[6]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[0]*fSkin[4])-1.623797632095822*coeff[0]*fEdge[4]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-1.623797632095822*coeff[0]*fSkin[5])-1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[0]*fSkin[4]+2.0625*coeff[0]*fEdge[4]-1.623797632095822*coeff[0]*fSkin[2]+1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]-1.623797632095822*coeff[0]*fSkin[3]+1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[6] = (-1.623797632095822*coeff[0]*fSkin[7])-1.623797632095822*coeff[0]*fEdge[7]+0.9375*coeff[0]*fSkin[6]-0.9375*coeff[0]*fEdge[6]; 
  edgeSurf_incr[7] = 3.5625*coeff[0]*fSkin[7]+2.0625*coeff[0]*fEdge[7]-1.623797632095822*coeff[0]*fSkin[6]+1.623797632095822*coeff[0]*fEdge[6]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[4] = 1.5*coeff[0]*fSkin[4]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[7] = 1.5*coeff[0]*fSkin[7]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfx_3x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.2651650429449552*coeff[7]*fSkin[7])+0.5740991584648069*coeff[6]*fSkin[7]+0.2651650429449552*coeff[7]*fEdge[7]+0.5740991584648069*coeff[6]*fEdge[7]+0.3314563036811939*coeff[6]*fSkin[6]-0.3314563036811939*coeff[6]*fEdge[6]-0.2651650429449552*coeff[5]*fSkin[5]+0.5740991584648069*coeff[3]*fSkin[5]+0.2651650429449552*coeff[5]*fEdge[5]+0.5740991584648069*coeff[3]*fEdge[5]-0.2651650429449552*coeff[4]*fSkin[4]+0.5740991584648069*coeff[2]*fSkin[4]+0.2651650429449552*coeff[4]*fEdge[4]+0.5740991584648069*coeff[2]*fEdge[4]+0.3314563036811939*coeff[3]*fSkin[3]-0.3314563036811939*coeff[3]*fEdge[3]+0.3314563036811939*coeff[2]*fSkin[2]-0.3314563036811939*coeff[2]*fEdge[2]-0.2651650429449552*coeff[1]*fSkin[1]+0.5740991584648069*coeff[0]*fSkin[1]+0.2651650429449552*coeff[1]*fEdge[1]+0.5740991584648069*coeff[0]*fEdge[1]+0.3314563036811939*coeff[0]*fSkin[0]-0.3314563036811939*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.4592793267718456*coeff[7]*fSkin[7])+1.259533953988537*coeff[6]*fSkin[7]+0.4592793267718456*coeff[7]*fEdge[7]+0.7292038680986265*coeff[6]*fEdge[7]+0.5740991584648069*coeff[6]*fSkin[6]-0.5740991584648069*coeff[6]*fEdge[6]-0.4592793267718456*coeff[5]*fSkin[5]+1.259533953988537*coeff[3]*fSkin[5]+0.4592793267718456*coeff[5]*fEdge[5]+0.7292038680986265*coeff[3]*fEdge[5]-0.4592793267718456*coeff[4]*fSkin[4]+1.259533953988537*coeff[2]*fSkin[4]+0.4592793267718456*coeff[4]*fEdge[4]+0.7292038680986265*coeff[2]*fEdge[4]+0.5740991584648069*coeff[3]*fSkin[3]-0.5740991584648069*coeff[3]*fEdge[3]+0.5740991584648069*coeff[2]*fSkin[2]-0.5740991584648069*coeff[2]*fEdge[2]-0.4592793267718456*coeff[1]*fSkin[1]+1.259533953988537*coeff[0]*fSkin[1]+0.4592793267718456*coeff[1]*fEdge[1]+0.7292038680986265*coeff[0]*fEdge[1]+0.5740991584648069*coeff[0]*fSkin[0]-0.5740991584648069*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.2651650429449552*coeff[5]*fSkin[7])+0.5740991584648069*coeff[3]*fSkin[7]+0.2651650429449552*coeff[5]*fEdge[7]+0.5740991584648069*coeff[3]*fEdge[7]-0.2651650429449552*fSkin[5]*coeff[7]+0.2651650429449552*fEdge[5]*coeff[7]+0.3314563036811939*coeff[3]*fSkin[6]-0.3314563036811939*coeff[3]*fEdge[6]+0.5740991584648069*fSkin[5]*coeff[6]+0.5740991584648069*fEdge[5]*coeff[6]+0.3314563036811939*fSkin[3]*coeff[6]-0.3314563036811939*fEdge[3]*coeff[6]-0.2651650429449552*coeff[1]*fSkin[4]+0.5740991584648069*coeff[0]*fSkin[4]+0.2651650429449552*coeff[1]*fEdge[4]+0.5740991584648069*coeff[0]*fEdge[4]-0.2651650429449552*fSkin[1]*coeff[4]+0.2651650429449552*fEdge[1]*coeff[4]+0.3314563036811939*coeff[0]*fSkin[2]-0.3314563036811939*coeff[0]*fEdge[2]+0.5740991584648069*fSkin[1]*coeff[2]+0.5740991584648069*fEdge[1]*coeff[2]+0.3314563036811939*fSkin[0]*coeff[2]-0.3314563036811939*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-0.2651650429449552*coeff[4]*fSkin[7])+0.5740991584648069*coeff[2]*fSkin[7]+0.2651650429449552*coeff[4]*fEdge[7]+0.5740991584648069*coeff[2]*fEdge[7]-0.2651650429449552*fSkin[4]*coeff[7]+0.2651650429449552*fEdge[4]*coeff[7]+0.3314563036811939*coeff[2]*fSkin[6]-0.3314563036811939*coeff[2]*fEdge[6]+0.5740991584648069*fSkin[4]*coeff[6]+0.5740991584648069*fEdge[4]*coeff[6]+0.3314563036811939*fSkin[2]*coeff[6]-0.3314563036811939*fEdge[2]*coeff[6]-0.2651650429449552*coeff[1]*fSkin[5]+0.5740991584648069*coeff[0]*fSkin[5]+0.2651650429449552*coeff[1]*fEdge[5]+0.5740991584648069*coeff[0]*fEdge[5]-0.2651650429449552*fSkin[1]*coeff[5]+0.2651650429449552*fEdge[1]*coeff[5]+0.3314563036811939*coeff[0]*fSkin[3]-0.3314563036811939*coeff[0]*fEdge[3]+0.5740991584648069*fSkin[1]*coeff[3]+0.5740991584648069*fEdge[1]*coeff[3]+0.3314563036811939*fSkin[0]*coeff[3]-0.3314563036811939*fEdge[0]*coeff[3]; 
  edgeSurf_incr[4] = (-0.4592793267718456*coeff[5]*fSkin[7])+1.259533953988537*coeff[3]*fSkin[7]+0.4592793267718456*coeff[5]*fEdge[7]+0.7292038680986265*coeff[3]*fEdge[7]-0.4592793267718456*fSkin[5]*coeff[7]+0.4592793267718456*fEdge[5]*coeff[7]+0.5740991584648069*coeff[3]*fSkin[6]-0.5740991584648069*coeff[3]*fEdge[6]+1.259533953988537*fSkin[5]*coeff[6]+0.7292038680986265*fEdge[5]*coeff[6]+0.5740991584648069*fSkin[3]*coeff[6]-0.5740991584648069*fEdge[3]*coeff[6]-0.4592793267718456*coeff[1]*fSkin[4]+1.259533953988537*coeff[0]*fSkin[4]+0.4592793267718456*coeff[1]*fEdge[4]+0.7292038680986265*coeff[0]*fEdge[4]-0.4592793267718456*fSkin[1]*coeff[4]+0.4592793267718456*fEdge[1]*coeff[4]+0.5740991584648069*coeff[0]*fSkin[2]-0.5740991584648069*coeff[0]*fEdge[2]+1.259533953988537*fSkin[1]*coeff[2]+0.7292038680986265*fEdge[1]*coeff[2]+0.5740991584648069*fSkin[0]*coeff[2]-0.5740991584648069*fEdge[0]*coeff[2]; 
  edgeSurf_incr[5] = (-0.4592793267718456*coeff[4]*fSkin[7])+1.259533953988537*coeff[2]*fSkin[7]+0.4592793267718456*coeff[4]*fEdge[7]+0.7292038680986265*coeff[2]*fEdge[7]-0.4592793267718456*fSkin[4]*coeff[7]+0.4592793267718456*fEdge[4]*coeff[7]+0.5740991584648069*coeff[2]*fSkin[6]-0.5740991584648069*coeff[2]*fEdge[6]+1.259533953988537*fSkin[4]*coeff[6]+0.7292038680986265*fEdge[4]*coeff[6]+0.5740991584648069*fSkin[2]*coeff[6]-0.5740991584648069*fEdge[2]*coeff[6]-0.4592793267718456*coeff[1]*fSkin[5]+1.259533953988537*coeff[0]*fSkin[5]+0.4592793267718456*coeff[1]*fEdge[5]+0.7292038680986265*coeff[0]*fEdge[5]-0.4592793267718456*fSkin[1]*coeff[5]+0.4592793267718456*fEdge[1]*coeff[5]+0.5740991584648069*coeff[0]*fSkin[3]-0.5740991584648069*coeff[0]*fEdge[3]+1.259533953988537*fSkin[1]*coeff[3]+0.7292038680986265*fEdge[1]*coeff[3]+0.5740991584648069*fSkin[0]*coeff[3]-0.5740991584648069*fEdge[0]*coeff[3]; 
  edgeSurf_incr[6] = (-0.2651650429449552*coeff[1]*fSkin[7])+0.5740991584648069*coeff[0]*fSkin[7]+0.2651650429449552*coeff[1]*fEdge[7]+0.5740991584648069*coeff[0]*fEdge[7]-0.2651650429449552*fSkin[1]*coeff[7]+0.2651650429449552*fEdge[1]*coeff[7]+0.3314563036811939*coeff[0]*fSkin[6]-0.3314563036811939*coeff[0]*fEdge[6]+0.5740991584648069*fSkin[1]*coeff[6]+0.5740991584648069*fEdge[1]*coeff[6]+0.3314563036811939*fSkin[0]*coeff[6]-0.3314563036811939*fEdge[0]*coeff[6]-0.2651650429449552*coeff[4]*fSkin[5]+0.5740991584648069*coeff[2]*fSkin[5]+0.2651650429449552*coeff[4]*fEdge[5]+0.5740991584648069*coeff[2]*fEdge[5]-0.2651650429449552*fSkin[4]*coeff[5]+0.2651650429449552*fEdge[4]*coeff[5]+0.5740991584648069*coeff[3]*fSkin[4]+0.5740991584648069*coeff[3]*fEdge[4]+0.3314563036811939*coeff[2]*fSkin[3]-0.3314563036811939*coeff[2]*fEdge[3]+0.3314563036811939*fSkin[2]*coeff[3]-0.3314563036811939*fEdge[2]*coeff[3]; 
  edgeSurf_incr[7] = (-0.4592793267718456*coeff[1]*fSkin[7])+1.259533953988537*coeff[0]*fSkin[7]+0.4592793267718456*coeff[1]*fEdge[7]+0.7292038680986265*coeff[0]*fEdge[7]-0.4592793267718456*fSkin[1]*coeff[7]+0.4592793267718456*fEdge[1]*coeff[7]+0.5740991584648069*coeff[0]*fSkin[6]-0.5740991584648069*coeff[0]*fEdge[6]+1.259533953988537*fSkin[1]*coeff[6]+0.7292038680986265*fEdge[1]*coeff[6]+0.5740991584648069*fSkin[0]*coeff[6]-0.5740991584648069*fEdge[0]*coeff[6]-0.4592793267718456*coeff[4]*fSkin[5]+1.259533953988537*coeff[2]*fSkin[5]+0.4592793267718456*coeff[4]*fEdge[5]+0.7292038680986265*coeff[2]*fEdge[5]-0.4592793267718456*fSkin[4]*coeff[5]+0.4592793267718456*fEdge[4]*coeff[5]+1.259533953988537*coeff[3]*fSkin[4]+0.7292038680986265*coeff[3]*fEdge[4]+0.5740991584648069*coeff[2]*fSkin[3]-0.5740991584648069*coeff[2]*fEdge[3]+0.5740991584648069*fSkin[2]*coeff[3]-0.5740991584648069*fEdge[2]*coeff[3]; 

  boundSurf_incr[0] = (-0.5303300858899105*coeff[7]*fSkin[7])-0.5303300858899105*coeff[5]*fSkin[5]-0.5303300858899105*coeff[4]*fSkin[4]-0.5303300858899105*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = 0.9185586535436913*coeff[7]*fSkin[7]+0.5303300858899105*coeff[6]*fSkin[7]+0.9185586535436913*coeff[5]*fSkin[5]+0.5303300858899105*coeff[3]*fSkin[5]+0.9185586535436913*coeff[4]*fSkin[4]+0.5303300858899105*coeff[2]*fSkin[4]+0.9185586535436913*coeff[1]*fSkin[1]+0.5303300858899105*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-0.5303300858899105*coeff[5]*fSkin[7])-0.5303300858899105*fSkin[5]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[4]-0.5303300858899105*fSkin[1]*coeff[4]; 
  boundSurf_incr[3] = (-0.5303300858899105*coeff[4]*fSkin[7])-0.5303300858899105*fSkin[4]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[5]-0.5303300858899105*fSkin[1]*coeff[5]; 
  boundSurf_incr[4] = 0.9185586535436913*coeff[5]*fSkin[7]+0.5303300858899105*coeff[3]*fSkin[7]+0.9185586535436913*fSkin[5]*coeff[7]+0.5303300858899105*fSkin[5]*coeff[6]+0.9185586535436913*coeff[1]*fSkin[4]+0.5303300858899105*coeff[0]*fSkin[4]+0.9185586535436913*fSkin[1]*coeff[4]+0.5303300858899105*fSkin[1]*coeff[2]; 
  boundSurf_incr[5] = 0.9185586535436913*coeff[4]*fSkin[7]+0.5303300858899105*coeff[2]*fSkin[7]+0.9185586535436913*fSkin[4]*coeff[7]+0.5303300858899105*fSkin[4]*coeff[6]+0.9185586535436913*coeff[1]*fSkin[5]+0.5303300858899105*coeff[0]*fSkin[5]+0.9185586535436913*fSkin[1]*coeff[5]+0.5303300858899105*fSkin[1]*coeff[3]; 
  boundSurf_incr[6] = (-0.5303300858899105*coeff[1]*fSkin[7])-0.5303300858899105*fSkin[1]*coeff[7]-0.5303300858899105*coeff[4]*fSkin[5]-0.5303300858899105*fSkin[4]*coeff[5]; 
  boundSurf_incr[7] = 0.9185586535436913*coeff[1]*fSkin[7]+0.5303300858899105*coeff[0]*fSkin[7]+0.9185586535436913*fSkin[1]*coeff[7]+0.5303300858899105*fSkin[1]*coeff[6]+0.9185586535436913*coeff[4]*fSkin[5]+0.5303300858899105*coeff[2]*fSkin[5]+0.9185586535436913*fSkin[4]*coeff[5]+0.5303300858899105*coeff[3]*fSkin[4]; 

  } else { 

  edgeSurf_incr[0] = (-0.2651650429449552*coeff[7]*fSkin[7])-0.5740991584648069*coeff[6]*fSkin[7]+0.2651650429449552*coeff[7]*fEdge[7]-0.5740991584648069*coeff[6]*fEdge[7]+0.3314563036811939*coeff[6]*fSkin[6]-0.3314563036811939*coeff[6]*fEdge[6]-0.2651650429449552*coeff[5]*fSkin[5]-0.5740991584648069*coeff[3]*fSkin[5]+0.2651650429449552*coeff[5]*fEdge[5]-0.5740991584648069*coeff[3]*fEdge[5]-0.2651650429449552*coeff[4]*fSkin[4]-0.5740991584648069*coeff[2]*fSkin[4]+0.2651650429449552*coeff[4]*fEdge[4]-0.5740991584648069*coeff[2]*fEdge[4]+0.3314563036811939*coeff[3]*fSkin[3]-0.3314563036811939*coeff[3]*fEdge[3]+0.3314563036811939*coeff[2]*fSkin[2]-0.3314563036811939*coeff[2]*fEdge[2]-0.2651650429449552*coeff[1]*fSkin[1]-0.5740991584648069*coeff[0]*fSkin[1]+0.2651650429449552*coeff[1]*fEdge[1]-0.5740991584648069*coeff[0]*fEdge[1]+0.3314563036811939*coeff[0]*fSkin[0]-0.3314563036811939*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.4592793267718456*coeff[7]*fSkin[7]+1.259533953988537*coeff[6]*fSkin[7]-0.4592793267718456*coeff[7]*fEdge[7]+0.7292038680986265*coeff[6]*fEdge[7]-0.5740991584648069*coeff[6]*fSkin[6]+0.5740991584648069*coeff[6]*fEdge[6]+0.4592793267718456*coeff[5]*fSkin[5]+1.259533953988537*coeff[3]*fSkin[5]-0.4592793267718456*coeff[5]*fEdge[5]+0.7292038680986265*coeff[3]*fEdge[5]+0.4592793267718456*coeff[4]*fSkin[4]+1.259533953988537*coeff[2]*fSkin[4]-0.4592793267718456*coeff[4]*fEdge[4]+0.7292038680986265*coeff[2]*fEdge[4]-0.5740991584648069*coeff[3]*fSkin[3]+0.5740991584648069*coeff[3]*fEdge[3]-0.5740991584648069*coeff[2]*fSkin[2]+0.5740991584648069*coeff[2]*fEdge[2]+0.4592793267718456*coeff[1]*fSkin[1]+1.259533953988537*coeff[0]*fSkin[1]-0.4592793267718456*coeff[1]*fEdge[1]+0.7292038680986265*coeff[0]*fEdge[1]-0.5740991584648069*coeff[0]*fSkin[0]+0.5740991584648069*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.2651650429449552*coeff[5]*fSkin[7])-0.5740991584648069*coeff[3]*fSkin[7]+0.2651650429449552*coeff[5]*fEdge[7]-0.5740991584648069*coeff[3]*fEdge[7]-0.2651650429449552*fSkin[5]*coeff[7]+0.2651650429449552*fEdge[5]*coeff[7]+0.3314563036811939*coeff[3]*fSkin[6]-0.3314563036811939*coeff[3]*fEdge[6]-0.5740991584648069*fSkin[5]*coeff[6]-0.5740991584648069*fEdge[5]*coeff[6]+0.3314563036811939*fSkin[3]*coeff[6]-0.3314563036811939*fEdge[3]*coeff[6]-0.2651650429449552*coeff[1]*fSkin[4]-0.5740991584648069*coeff[0]*fSkin[4]+0.2651650429449552*coeff[1]*fEdge[4]-0.5740991584648069*coeff[0]*fEdge[4]-0.2651650429449552*fSkin[1]*coeff[4]+0.2651650429449552*fEdge[1]*coeff[4]+0.3314563036811939*coeff[0]*fSkin[2]-0.3314563036811939*coeff[0]*fEdge[2]-0.5740991584648069*fSkin[1]*coeff[2]-0.5740991584648069*fEdge[1]*coeff[2]+0.3314563036811939*fSkin[0]*coeff[2]-0.3314563036811939*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-0.2651650429449552*coeff[4]*fSkin[7])-0.5740991584648069*coeff[2]*fSkin[7]+0.2651650429449552*coeff[4]*fEdge[7]-0.5740991584648069*coeff[2]*fEdge[7]-0.2651650429449552*fSkin[4]*coeff[7]+0.2651650429449552*fEdge[4]*coeff[7]+0.3314563036811939*coeff[2]*fSkin[6]-0.3314563036811939*coeff[2]*fEdge[6]-0.5740991584648069*fSkin[4]*coeff[6]-0.5740991584648069*fEdge[4]*coeff[6]+0.3314563036811939*fSkin[2]*coeff[6]-0.3314563036811939*fEdge[2]*coeff[6]-0.2651650429449552*coeff[1]*fSkin[5]-0.5740991584648069*coeff[0]*fSkin[5]+0.2651650429449552*coeff[1]*fEdge[5]-0.5740991584648069*coeff[0]*fEdge[5]-0.2651650429449552*fSkin[1]*coeff[5]+0.2651650429449552*fEdge[1]*coeff[5]+0.3314563036811939*coeff[0]*fSkin[3]-0.3314563036811939*coeff[0]*fEdge[3]-0.5740991584648069*fSkin[1]*coeff[3]-0.5740991584648069*fEdge[1]*coeff[3]+0.3314563036811939*fSkin[0]*coeff[3]-0.3314563036811939*fEdge[0]*coeff[3]; 
  edgeSurf_incr[4] = 0.4592793267718456*coeff[5]*fSkin[7]+1.259533953988537*coeff[3]*fSkin[7]-0.4592793267718456*coeff[5]*fEdge[7]+0.7292038680986265*coeff[3]*fEdge[7]+0.4592793267718456*fSkin[5]*coeff[7]-0.4592793267718456*fEdge[5]*coeff[7]-0.5740991584648069*coeff[3]*fSkin[6]+0.5740991584648069*coeff[3]*fEdge[6]+1.259533953988537*fSkin[5]*coeff[6]+0.7292038680986265*fEdge[5]*coeff[6]-0.5740991584648069*fSkin[3]*coeff[6]+0.5740991584648069*fEdge[3]*coeff[6]+0.4592793267718456*coeff[1]*fSkin[4]+1.259533953988537*coeff[0]*fSkin[4]-0.4592793267718456*coeff[1]*fEdge[4]+0.7292038680986265*coeff[0]*fEdge[4]+0.4592793267718456*fSkin[1]*coeff[4]-0.4592793267718456*fEdge[1]*coeff[4]-0.5740991584648069*coeff[0]*fSkin[2]+0.5740991584648069*coeff[0]*fEdge[2]+1.259533953988537*fSkin[1]*coeff[2]+0.7292038680986265*fEdge[1]*coeff[2]-0.5740991584648069*fSkin[0]*coeff[2]+0.5740991584648069*fEdge[0]*coeff[2]; 
  edgeSurf_incr[5] = 0.4592793267718456*coeff[4]*fSkin[7]+1.259533953988537*coeff[2]*fSkin[7]-0.4592793267718456*coeff[4]*fEdge[7]+0.7292038680986265*coeff[2]*fEdge[7]+0.4592793267718456*fSkin[4]*coeff[7]-0.4592793267718456*fEdge[4]*coeff[7]-0.5740991584648069*coeff[2]*fSkin[6]+0.5740991584648069*coeff[2]*fEdge[6]+1.259533953988537*fSkin[4]*coeff[6]+0.7292038680986265*fEdge[4]*coeff[6]-0.5740991584648069*fSkin[2]*coeff[6]+0.5740991584648069*fEdge[2]*coeff[6]+0.4592793267718456*coeff[1]*fSkin[5]+1.259533953988537*coeff[0]*fSkin[5]-0.4592793267718456*coeff[1]*fEdge[5]+0.7292038680986265*coeff[0]*fEdge[5]+0.4592793267718456*fSkin[1]*coeff[5]-0.4592793267718456*fEdge[1]*coeff[5]-0.5740991584648069*coeff[0]*fSkin[3]+0.5740991584648069*coeff[0]*fEdge[3]+1.259533953988537*fSkin[1]*coeff[3]+0.7292038680986265*fEdge[1]*coeff[3]-0.5740991584648069*fSkin[0]*coeff[3]+0.5740991584648069*fEdge[0]*coeff[3]; 
  edgeSurf_incr[6] = (-0.2651650429449552*coeff[1]*fSkin[7])-0.5740991584648069*coeff[0]*fSkin[7]+0.2651650429449552*coeff[1]*fEdge[7]-0.5740991584648069*coeff[0]*fEdge[7]-0.2651650429449552*fSkin[1]*coeff[7]+0.2651650429449552*fEdge[1]*coeff[7]+0.3314563036811939*coeff[0]*fSkin[6]-0.3314563036811939*coeff[0]*fEdge[6]-0.5740991584648069*fSkin[1]*coeff[6]-0.5740991584648069*fEdge[1]*coeff[6]+0.3314563036811939*fSkin[0]*coeff[6]-0.3314563036811939*fEdge[0]*coeff[6]-0.2651650429449552*coeff[4]*fSkin[5]-0.5740991584648069*coeff[2]*fSkin[5]+0.2651650429449552*coeff[4]*fEdge[5]-0.5740991584648069*coeff[2]*fEdge[5]-0.2651650429449552*fSkin[4]*coeff[5]+0.2651650429449552*fEdge[4]*coeff[5]-0.5740991584648069*coeff[3]*fSkin[4]-0.5740991584648069*coeff[3]*fEdge[4]+0.3314563036811939*coeff[2]*fSkin[3]-0.3314563036811939*coeff[2]*fEdge[3]+0.3314563036811939*fSkin[2]*coeff[3]-0.3314563036811939*fEdge[2]*coeff[3]; 
  edgeSurf_incr[7] = 0.4592793267718456*coeff[1]*fSkin[7]+1.259533953988537*coeff[0]*fSkin[7]-0.4592793267718456*coeff[1]*fEdge[7]+0.7292038680986265*coeff[0]*fEdge[7]+0.4592793267718456*fSkin[1]*coeff[7]-0.4592793267718456*fEdge[1]*coeff[7]-0.5740991584648069*coeff[0]*fSkin[6]+0.5740991584648069*coeff[0]*fEdge[6]+1.259533953988537*fSkin[1]*coeff[6]+0.7292038680986265*fEdge[1]*coeff[6]-0.5740991584648069*fSkin[0]*coeff[6]+0.5740991584648069*fEdge[0]*coeff[6]+0.4592793267718456*coeff[4]*fSkin[5]+1.259533953988537*coeff[2]*fSkin[5]-0.4592793267718456*coeff[4]*fEdge[5]+0.7292038680986265*coeff[2]*fEdge[5]+0.4592793267718456*fSkin[4]*coeff[5]-0.4592793267718456*fEdge[4]*coeff[5]+1.259533953988537*coeff[3]*fSkin[4]+0.7292038680986265*coeff[3]*fEdge[4]-0.5740991584648069*coeff[2]*fSkin[3]+0.5740991584648069*coeff[2]*fEdge[3]-0.5740991584648069*fSkin[2]*coeff[3]+0.5740991584648069*fEdge[2]*coeff[3]; 

  boundSurf_incr[0] = (-0.5303300858899105*coeff[7]*fSkin[7])-0.5303300858899105*coeff[5]*fSkin[5]-0.5303300858899105*coeff[4]*fSkin[4]-0.5303300858899105*coeff[1]*fSkin[1]; 
  boundSurf_incr[1] = (-0.9185586535436913*coeff[7]*fSkin[7])+0.5303300858899105*coeff[6]*fSkin[7]-0.9185586535436913*coeff[5]*fSkin[5]+0.5303300858899105*coeff[3]*fSkin[5]-0.9185586535436913*coeff[4]*fSkin[4]+0.5303300858899105*coeff[2]*fSkin[4]-0.9185586535436913*coeff[1]*fSkin[1]+0.5303300858899105*coeff[0]*fSkin[1]; 
  boundSurf_incr[2] = (-0.5303300858899105*coeff[5]*fSkin[7])-0.5303300858899105*fSkin[5]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[4]-0.5303300858899105*fSkin[1]*coeff[4]; 
  boundSurf_incr[3] = (-0.5303300858899105*coeff[4]*fSkin[7])-0.5303300858899105*fSkin[4]*coeff[7]-0.5303300858899105*coeff[1]*fSkin[5]-0.5303300858899105*fSkin[1]*coeff[5]; 
  boundSurf_incr[4] = (-0.9185586535436913*coeff[5]*fSkin[7])+0.5303300858899105*coeff[3]*fSkin[7]-0.9185586535436913*fSkin[5]*coeff[7]+0.5303300858899105*fSkin[5]*coeff[6]-0.9185586535436913*coeff[1]*fSkin[4]+0.5303300858899105*coeff[0]*fSkin[4]-0.9185586535436913*fSkin[1]*coeff[4]+0.5303300858899105*fSkin[1]*coeff[2]; 
  boundSurf_incr[5] = (-0.9185586535436913*coeff[4]*fSkin[7])+0.5303300858899105*coeff[2]*fSkin[7]-0.9185586535436913*fSkin[4]*coeff[7]+0.5303300858899105*fSkin[4]*coeff[6]-0.9185586535436913*coeff[1]*fSkin[5]+0.5303300858899105*coeff[0]*fSkin[5]-0.9185586535436913*fSkin[1]*coeff[5]+0.5303300858899105*fSkin[1]*coeff[3]; 
  boundSurf_incr[6] = (-0.5303300858899105*coeff[1]*fSkin[7])-0.5303300858899105*fSkin[1]*coeff[7]-0.5303300858899105*coeff[4]*fSkin[5]-0.5303300858899105*fSkin[4]*coeff[5]; 
  boundSurf_incr[7] = (-0.9185586535436913*coeff[1]*fSkin[7])+0.5303300858899105*coeff[0]*fSkin[7]-0.9185586535436913*fSkin[1]*coeff[7]+0.5303300858899105*fSkin[1]*coeff[6]-0.9185586535436913*coeff[4]*fSkin[5]+0.5303300858899105*coeff[2]*fSkin[5]-0.9185586535436913*fSkin[4]*coeff[5]+0.5303300858899105*coeff[3]*fSkin[4]; 

  out[0] += -1.0*(vol_incr[0]+edgeSurf_incr[0]+boundSurf_incr[0])*Jfac; 
  out[1] += -1.0*(vol_incr[1]+edgeSurf_incr[1]+boundSurf_incr[1])*Jfac; 
  out[2] += -1.0*(vol_incr[2]+edgeSurf_incr[2]+boundSurf_incr[2])*Jfac; 
  out[3] += -1.0*(vol_incr[3]+edgeSurf_incr[3]+boundSurf_incr[3])*Jfac; 
  out[4] += -1.0*(vol_incr[4]+edgeSurf_incr[4]+boundSurf_incr[4])*Jfac; 
  out[5] += -1.0*(vol_incr[5]+edgeSurf_incr[5]+boundSurf_incr[5])*Jfac; 
  out[6] += -1.0*(vol_incr[6]+edgeSurf_incr[6]+boundSurf_incr[6])*Jfac; 
  out[7] += -1.0*(vol_incr[7]+edgeSurf_incr[7]+boundSurf_incr[7])*Jfac; 

  }

  return 0.;
}

