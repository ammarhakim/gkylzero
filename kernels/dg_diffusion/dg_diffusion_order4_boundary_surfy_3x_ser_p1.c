#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order4_boundary_surfy_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[1]*fSkin[2]+1.623797632095822*coeff[1]*fEdge[2]+0.9375*fSkin[0]*coeff[1]-0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 1.623797632095822*coeff[1]*fSkin[4]+1.623797632095822*coeff[1]*fEdge[4]+0.9375*coeff[1]*fSkin[1]-0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*fSkin[2]+2.0625*coeff[1]*fEdge[2]+1.623797632095822*fSkin[0]*coeff[1]-1.623797632095822*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 1.623797632095822*coeff[1]*fSkin[6]+1.623797632095822*coeff[1]*fEdge[6]+0.9375*coeff[1]*fSkin[3]-0.9375*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[1]*fSkin[4]+2.0625*coeff[1]*fEdge[4]+1.623797632095822*coeff[1]*fSkin[1]-1.623797632095822*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = 1.623797632095822*coeff[1]*fSkin[7]+1.623797632095822*coeff[1]*fEdge[7]+0.9375*coeff[1]*fSkin[5]-0.9375*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = 3.5625*coeff[1]*fSkin[6]+2.0625*coeff[1]*fEdge[6]+1.623797632095822*coeff[1]*fSkin[3]-1.623797632095822*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = 3.5625*coeff[1]*fSkin[7]+2.0625*coeff[1]*fEdge[7]+1.623797632095822*coeff[1]*fSkin[5]-1.623797632095822*coeff[1]*fEdge[5]; 

  boundSurf_incr[2] = 1.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 1.5*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 1.5*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = 1.5*coeff[1]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[1]*fSkin[2])-1.623797632095822*coeff[1]*fEdge[2]+0.9375*fSkin[0]*coeff[1]-0.9375*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = (-1.623797632095822*coeff[1]*fSkin[4])-1.623797632095822*coeff[1]*fEdge[4]+0.9375*coeff[1]*fSkin[1]-0.9375*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = 3.5625*coeff[1]*fSkin[2]+2.0625*coeff[1]*fEdge[2]-1.623797632095822*fSkin[0]*coeff[1]+1.623797632095822*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = (-1.623797632095822*coeff[1]*fSkin[6])-1.623797632095822*coeff[1]*fEdge[6]+0.9375*coeff[1]*fSkin[3]-0.9375*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 3.5625*coeff[1]*fSkin[4]+2.0625*coeff[1]*fEdge[4]-1.623797632095822*coeff[1]*fSkin[1]+1.623797632095822*coeff[1]*fEdge[1]; 
  edgeSurf_incr[5] = (-1.623797632095822*coeff[1]*fSkin[7])-1.623797632095822*coeff[1]*fEdge[7]+0.9375*coeff[1]*fSkin[5]-0.9375*coeff[1]*fEdge[5]; 
  edgeSurf_incr[6] = 3.5625*coeff[1]*fSkin[6]+2.0625*coeff[1]*fEdge[6]-1.623797632095822*coeff[1]*fSkin[3]+1.623797632095822*coeff[1]*fEdge[3]; 
  edgeSurf_incr[7] = 3.5625*coeff[1]*fSkin[7]+2.0625*coeff[1]*fEdge[7]-1.623797632095822*coeff[1]*fSkin[5]+1.623797632095822*coeff[1]*fEdge[5]; 

  boundSurf_incr[2] = 1.5*coeff[1]*fSkin[2]; 
  boundSurf_incr[4] = 1.5*coeff[1]*fSkin[4]; 
  boundSurf_incr[6] = 1.5*coeff[1]*fSkin[6]; 
  boundSurf_incr[7] = 1.5*coeff[1]*fSkin[7]; 

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

GKYL_CU_DH double dg_diffusion_order4_boundary_surfy_3x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.2651650429449552*fSkin[7]*coeff[15])+0.2651650429449552*fEdge[7]*coeff[15]-0.2651650429449552*fSkin[6]*coeff[14]+0.2651650429449552*fEdge[6]*coeff[14]+0.5740991584648069*fSkin[7]*coeff[13]+0.5740991584648069*fEdge[7]*coeff[13]+0.3314563036811939*fSkin[5]*coeff[13]-0.3314563036811939*fEdge[5]*coeff[13]-0.2651650429449552*fSkin[4]*coeff[12]+0.2651650429449552*fEdge[4]*coeff[12]+0.5740991584648069*fSkin[6]*coeff[11]+0.5740991584648069*fEdge[6]*coeff[11]+0.3314563036811939*fSkin[3]*coeff[11]-0.3314563036811939*fEdge[3]*coeff[11]-0.2651650429449552*fSkin[2]*coeff[10]+0.2651650429449552*fEdge[2]*coeff[10]+0.5740991584648069*fSkin[4]*coeff[9]+0.5740991584648069*fEdge[4]*coeff[9]+0.3314563036811939*fSkin[1]*coeff[9]-0.3314563036811939*fEdge[1]*coeff[9]+0.5740991584648069*fSkin[2]*coeff[8]+0.5740991584648069*fEdge[2]*coeff[8]+0.3314563036811939*fSkin[0]*coeff[8]-0.3314563036811939*fEdge[0]*coeff[8]; 
  edgeSurf_incr[1] = (-0.2651650429449552*fSkin[6]*coeff[15])+0.2651650429449552*fEdge[6]*coeff[15]-0.2651650429449552*fSkin[7]*coeff[14]+0.2651650429449552*fEdge[7]*coeff[14]+0.5740991584648069*fSkin[6]*coeff[13]+0.5740991584648069*fEdge[6]*coeff[13]+0.3314563036811939*fSkin[3]*coeff[13]-0.3314563036811939*fEdge[3]*coeff[13]-0.2651650429449552*fSkin[2]*coeff[12]+0.2651650429449552*fEdge[2]*coeff[12]+0.5740991584648069*fSkin[7]*coeff[11]+0.5740991584648069*fEdge[7]*coeff[11]+0.3314563036811939*fSkin[5]*coeff[11]-0.3314563036811939*fEdge[5]*coeff[11]-0.2651650429449552*fSkin[4]*coeff[10]+0.2651650429449552*fEdge[4]*coeff[10]+0.5740991584648069*fSkin[2]*coeff[9]+0.5740991584648069*fEdge[2]*coeff[9]+0.3314563036811939*fSkin[0]*coeff[9]-0.3314563036811939*fEdge[0]*coeff[9]+0.5740991584648069*fSkin[4]*coeff[8]+0.5740991584648069*fEdge[4]*coeff[8]+0.3314563036811939*fSkin[1]*coeff[8]-0.3314563036811939*fEdge[1]*coeff[8]; 
  edgeSurf_incr[2] = (-0.4592793267718456*fSkin[7]*coeff[15])+0.4592793267718456*fEdge[7]*coeff[15]-0.4592793267718456*fSkin[6]*coeff[14]+0.4592793267718456*fEdge[6]*coeff[14]+1.259533953988537*fSkin[7]*coeff[13]+0.7292038680986265*fEdge[7]*coeff[13]+0.5740991584648069*fSkin[5]*coeff[13]-0.5740991584648069*fEdge[5]*coeff[13]-0.4592793267718456*fSkin[4]*coeff[12]+0.4592793267718456*fEdge[4]*coeff[12]+1.259533953988537*fSkin[6]*coeff[11]+0.7292038680986265*fEdge[6]*coeff[11]+0.5740991584648069*fSkin[3]*coeff[11]-0.5740991584648069*fEdge[3]*coeff[11]-0.4592793267718456*fSkin[2]*coeff[10]+0.4592793267718456*fEdge[2]*coeff[10]+1.259533953988537*fSkin[4]*coeff[9]+0.7292038680986265*fEdge[4]*coeff[9]+0.5740991584648069*fSkin[1]*coeff[9]-0.5740991584648069*fEdge[1]*coeff[9]+1.259533953988537*fSkin[2]*coeff[8]+0.7292038680986265*fEdge[2]*coeff[8]+0.5740991584648069*fSkin[0]*coeff[8]-0.5740991584648069*fEdge[0]*coeff[8]; 
  edgeSurf_incr[3] = (-0.2651650429449552*fSkin[4]*coeff[15])+0.2651650429449552*fEdge[4]*coeff[15]-0.2651650429449552*fSkin[2]*coeff[14]+0.2651650429449552*fEdge[2]*coeff[14]+0.5740991584648069*fSkin[4]*coeff[13]+0.5740991584648069*fEdge[4]*coeff[13]+0.3314563036811939*fSkin[1]*coeff[13]-0.3314563036811939*fEdge[1]*coeff[13]-0.2651650429449552*fSkin[7]*coeff[12]+0.2651650429449552*fEdge[7]*coeff[12]+0.5740991584648069*fSkin[2]*coeff[11]+0.5740991584648069*fEdge[2]*coeff[11]+0.3314563036811939*fSkin[0]*coeff[11]-0.3314563036811939*fEdge[0]*coeff[11]-0.2651650429449552*fSkin[6]*coeff[10]+0.2651650429449552*fEdge[6]*coeff[10]+0.5740991584648069*fSkin[7]*coeff[9]+0.5740991584648069*fEdge[7]*coeff[9]+0.3314563036811939*fSkin[5]*coeff[9]-0.3314563036811939*fEdge[5]*coeff[9]+0.5740991584648069*fSkin[6]*coeff[8]+0.5740991584648069*fEdge[6]*coeff[8]+0.3314563036811939*fSkin[3]*coeff[8]-0.3314563036811939*fEdge[3]*coeff[8]; 
  edgeSurf_incr[4] = (-0.4592793267718456*fSkin[6]*coeff[15])+0.4592793267718456*fEdge[6]*coeff[15]-0.4592793267718456*fSkin[7]*coeff[14]+0.4592793267718456*fEdge[7]*coeff[14]+1.259533953988537*fSkin[6]*coeff[13]+0.7292038680986265*fEdge[6]*coeff[13]+0.5740991584648069*fSkin[3]*coeff[13]-0.5740991584648069*fEdge[3]*coeff[13]-0.4592793267718456*fSkin[2]*coeff[12]+0.4592793267718456*fEdge[2]*coeff[12]+1.259533953988537*fSkin[7]*coeff[11]+0.7292038680986265*fEdge[7]*coeff[11]+0.5740991584648069*fSkin[5]*coeff[11]-0.5740991584648069*fEdge[5]*coeff[11]-0.4592793267718456*fSkin[4]*coeff[10]+0.4592793267718456*fEdge[4]*coeff[10]+1.259533953988537*fSkin[2]*coeff[9]+0.7292038680986265*fEdge[2]*coeff[9]+0.5740991584648069*fSkin[0]*coeff[9]-0.5740991584648069*fEdge[0]*coeff[9]+1.259533953988537*fSkin[4]*coeff[8]+0.7292038680986265*fEdge[4]*coeff[8]+0.5740991584648069*fSkin[1]*coeff[8]-0.5740991584648069*fEdge[1]*coeff[8]; 
  edgeSurf_incr[5] = (-0.2651650429449552*fSkin[2]*coeff[15])+0.2651650429449552*fEdge[2]*coeff[15]-0.2651650429449552*fSkin[4]*coeff[14]+0.2651650429449552*fEdge[4]*coeff[14]+0.5740991584648069*fSkin[2]*coeff[13]+0.5740991584648069*fEdge[2]*coeff[13]+0.3314563036811939*fSkin[0]*coeff[13]-0.3314563036811939*fEdge[0]*coeff[13]-0.2651650429449552*fSkin[6]*coeff[12]+0.2651650429449552*fEdge[6]*coeff[12]+0.5740991584648069*fSkin[4]*coeff[11]+0.5740991584648069*fEdge[4]*coeff[11]+0.3314563036811939*fSkin[1]*coeff[11]-0.3314563036811939*fEdge[1]*coeff[11]-0.2651650429449552*fSkin[7]*coeff[10]+0.2651650429449552*fEdge[7]*coeff[10]+0.5740991584648069*fSkin[6]*coeff[9]+0.5740991584648069*fEdge[6]*coeff[9]+0.3314563036811939*fSkin[3]*coeff[9]-0.3314563036811939*fEdge[3]*coeff[9]+0.5740991584648069*fSkin[7]*coeff[8]+0.5740991584648069*fEdge[7]*coeff[8]+0.3314563036811939*fSkin[5]*coeff[8]-0.3314563036811939*fEdge[5]*coeff[8]; 
  edgeSurf_incr[6] = (-0.4592793267718456*fSkin[4]*coeff[15])+0.4592793267718456*fEdge[4]*coeff[15]-0.4592793267718456*fSkin[2]*coeff[14]+0.4592793267718456*fEdge[2]*coeff[14]+1.259533953988537*fSkin[4]*coeff[13]+0.7292038680986265*fEdge[4]*coeff[13]+0.5740991584648069*fSkin[1]*coeff[13]-0.5740991584648069*fEdge[1]*coeff[13]-0.4592793267718456*fSkin[7]*coeff[12]+0.4592793267718456*fEdge[7]*coeff[12]+1.259533953988537*fSkin[2]*coeff[11]+0.7292038680986265*fEdge[2]*coeff[11]+0.5740991584648069*fSkin[0]*coeff[11]-0.5740991584648069*fEdge[0]*coeff[11]-0.4592793267718456*fSkin[6]*coeff[10]+0.4592793267718456*fEdge[6]*coeff[10]+1.259533953988537*fSkin[7]*coeff[9]+0.7292038680986265*fEdge[7]*coeff[9]+0.5740991584648069*fSkin[5]*coeff[9]-0.5740991584648069*fEdge[5]*coeff[9]+1.259533953988537*fSkin[6]*coeff[8]+0.7292038680986265*fEdge[6]*coeff[8]+0.5740991584648069*fSkin[3]*coeff[8]-0.5740991584648069*fEdge[3]*coeff[8]; 
  edgeSurf_incr[7] = (-0.4592793267718456*fSkin[2]*coeff[15])+0.4592793267718456*fEdge[2]*coeff[15]-0.4592793267718456*fSkin[4]*coeff[14]+0.4592793267718456*fEdge[4]*coeff[14]+1.259533953988537*fSkin[2]*coeff[13]+0.7292038680986265*fEdge[2]*coeff[13]+0.5740991584648069*fSkin[0]*coeff[13]-0.5740991584648069*fEdge[0]*coeff[13]-0.4592793267718456*fSkin[6]*coeff[12]+0.4592793267718456*fEdge[6]*coeff[12]+1.259533953988537*fSkin[4]*coeff[11]+0.7292038680986265*fEdge[4]*coeff[11]+0.5740991584648069*fSkin[1]*coeff[11]-0.5740991584648069*fEdge[1]*coeff[11]-0.4592793267718456*fSkin[7]*coeff[10]+0.4592793267718456*fEdge[7]*coeff[10]+1.259533953988537*fSkin[6]*coeff[9]+0.7292038680986265*fEdge[6]*coeff[9]+0.5740991584648069*fSkin[3]*coeff[9]-0.5740991584648069*fEdge[3]*coeff[9]+1.259533953988537*fSkin[7]*coeff[8]+0.7292038680986265*fEdge[7]*coeff[8]+0.5740991584648069*fSkin[5]*coeff[8]-0.5740991584648069*fEdge[5]*coeff[8]; 

  boundSurf_incr[0] = (-0.5303300858899105*fSkin[7]*coeff[15])-0.5303300858899105*fSkin[6]*coeff[14]-0.5303300858899105*fSkin[4]*coeff[12]-0.5303300858899105*fSkin[2]*coeff[10]; 
  boundSurf_incr[1] = (-0.5303300858899105*fSkin[6]*coeff[15])-0.5303300858899105*fSkin[7]*coeff[14]-0.5303300858899105*fSkin[2]*coeff[12]-0.5303300858899105*fSkin[4]*coeff[10]; 
  boundSurf_incr[2] = 0.9185586535436913*fSkin[7]*coeff[15]+0.9185586535436913*fSkin[6]*coeff[14]+0.5303300858899105*fSkin[7]*coeff[13]+0.9185586535436913*fSkin[4]*coeff[12]+0.5303300858899105*fSkin[6]*coeff[11]+0.9185586535436913*fSkin[2]*coeff[10]+0.5303300858899105*fSkin[4]*coeff[9]+0.5303300858899105*fSkin[2]*coeff[8]; 
  boundSurf_incr[3] = (-0.5303300858899105*fSkin[4]*coeff[15])-0.5303300858899105*fSkin[2]*coeff[14]-0.5303300858899105*fSkin[7]*coeff[12]-0.5303300858899105*fSkin[6]*coeff[10]; 
  boundSurf_incr[4] = 0.9185586535436913*fSkin[6]*coeff[15]+0.9185586535436913*fSkin[7]*coeff[14]+0.5303300858899105*fSkin[6]*coeff[13]+0.9185586535436913*fSkin[2]*coeff[12]+0.5303300858899105*fSkin[7]*coeff[11]+0.9185586535436913*fSkin[4]*coeff[10]+0.5303300858899105*fSkin[2]*coeff[9]+0.5303300858899105*fSkin[4]*coeff[8]; 
  boundSurf_incr[5] = (-0.5303300858899105*fSkin[2]*coeff[15])-0.5303300858899105*fSkin[4]*coeff[14]-0.5303300858899105*fSkin[6]*coeff[12]-0.5303300858899105*fSkin[7]*coeff[10]; 
  boundSurf_incr[6] = 0.9185586535436913*fSkin[4]*coeff[15]+0.9185586535436913*fSkin[2]*coeff[14]+0.5303300858899105*fSkin[4]*coeff[13]+0.9185586535436913*fSkin[7]*coeff[12]+0.5303300858899105*fSkin[2]*coeff[11]+0.9185586535436913*fSkin[6]*coeff[10]+0.5303300858899105*fSkin[7]*coeff[9]+0.5303300858899105*fSkin[6]*coeff[8]; 
  boundSurf_incr[7] = 0.9185586535436913*fSkin[2]*coeff[15]+0.9185586535436913*fSkin[4]*coeff[14]+0.5303300858899105*fSkin[2]*coeff[13]+0.9185586535436913*fSkin[6]*coeff[12]+0.5303300858899105*fSkin[4]*coeff[11]+0.9185586535436913*fSkin[7]*coeff[10]+0.5303300858899105*fSkin[6]*coeff[9]+0.5303300858899105*fSkin[7]*coeff[8]; 

  } else { 

  edgeSurf_incr[0] = (-0.2651650429449552*fSkin[7]*coeff[15])+0.2651650429449552*fEdge[7]*coeff[15]-0.2651650429449552*fSkin[6]*coeff[14]+0.2651650429449552*fEdge[6]*coeff[14]-0.5740991584648069*fSkin[7]*coeff[13]-0.5740991584648069*fEdge[7]*coeff[13]+0.3314563036811939*fSkin[5]*coeff[13]-0.3314563036811939*fEdge[5]*coeff[13]-0.2651650429449552*fSkin[4]*coeff[12]+0.2651650429449552*fEdge[4]*coeff[12]-0.5740991584648069*fSkin[6]*coeff[11]-0.5740991584648069*fEdge[6]*coeff[11]+0.3314563036811939*fSkin[3]*coeff[11]-0.3314563036811939*fEdge[3]*coeff[11]-0.2651650429449552*fSkin[2]*coeff[10]+0.2651650429449552*fEdge[2]*coeff[10]-0.5740991584648069*fSkin[4]*coeff[9]-0.5740991584648069*fEdge[4]*coeff[9]+0.3314563036811939*fSkin[1]*coeff[9]-0.3314563036811939*fEdge[1]*coeff[9]-0.5740991584648069*fSkin[2]*coeff[8]-0.5740991584648069*fEdge[2]*coeff[8]+0.3314563036811939*fSkin[0]*coeff[8]-0.3314563036811939*fEdge[0]*coeff[8]; 
  edgeSurf_incr[1] = (-0.2651650429449552*fSkin[6]*coeff[15])+0.2651650429449552*fEdge[6]*coeff[15]-0.2651650429449552*fSkin[7]*coeff[14]+0.2651650429449552*fEdge[7]*coeff[14]-0.5740991584648069*fSkin[6]*coeff[13]-0.5740991584648069*fEdge[6]*coeff[13]+0.3314563036811939*fSkin[3]*coeff[13]-0.3314563036811939*fEdge[3]*coeff[13]-0.2651650429449552*fSkin[2]*coeff[12]+0.2651650429449552*fEdge[2]*coeff[12]-0.5740991584648069*fSkin[7]*coeff[11]-0.5740991584648069*fEdge[7]*coeff[11]+0.3314563036811939*fSkin[5]*coeff[11]-0.3314563036811939*fEdge[5]*coeff[11]-0.2651650429449552*fSkin[4]*coeff[10]+0.2651650429449552*fEdge[4]*coeff[10]-0.5740991584648069*fSkin[2]*coeff[9]-0.5740991584648069*fEdge[2]*coeff[9]+0.3314563036811939*fSkin[0]*coeff[9]-0.3314563036811939*fEdge[0]*coeff[9]-0.5740991584648069*fSkin[4]*coeff[8]-0.5740991584648069*fEdge[4]*coeff[8]+0.3314563036811939*fSkin[1]*coeff[8]-0.3314563036811939*fEdge[1]*coeff[8]; 
  edgeSurf_incr[2] = 0.4592793267718456*fSkin[7]*coeff[15]-0.4592793267718456*fEdge[7]*coeff[15]+0.4592793267718456*fSkin[6]*coeff[14]-0.4592793267718456*fEdge[6]*coeff[14]+1.259533953988537*fSkin[7]*coeff[13]+0.7292038680986265*fEdge[7]*coeff[13]-0.5740991584648069*fSkin[5]*coeff[13]+0.5740991584648069*fEdge[5]*coeff[13]+0.4592793267718456*fSkin[4]*coeff[12]-0.4592793267718456*fEdge[4]*coeff[12]+1.259533953988537*fSkin[6]*coeff[11]+0.7292038680986265*fEdge[6]*coeff[11]-0.5740991584648069*fSkin[3]*coeff[11]+0.5740991584648069*fEdge[3]*coeff[11]+0.4592793267718456*fSkin[2]*coeff[10]-0.4592793267718456*fEdge[2]*coeff[10]+1.259533953988537*fSkin[4]*coeff[9]+0.7292038680986265*fEdge[4]*coeff[9]-0.5740991584648069*fSkin[1]*coeff[9]+0.5740991584648069*fEdge[1]*coeff[9]+1.259533953988537*fSkin[2]*coeff[8]+0.7292038680986265*fEdge[2]*coeff[8]-0.5740991584648069*fSkin[0]*coeff[8]+0.5740991584648069*fEdge[0]*coeff[8]; 
  edgeSurf_incr[3] = (-0.2651650429449552*fSkin[4]*coeff[15])+0.2651650429449552*fEdge[4]*coeff[15]-0.2651650429449552*fSkin[2]*coeff[14]+0.2651650429449552*fEdge[2]*coeff[14]-0.5740991584648069*fSkin[4]*coeff[13]-0.5740991584648069*fEdge[4]*coeff[13]+0.3314563036811939*fSkin[1]*coeff[13]-0.3314563036811939*fEdge[1]*coeff[13]-0.2651650429449552*fSkin[7]*coeff[12]+0.2651650429449552*fEdge[7]*coeff[12]-0.5740991584648069*fSkin[2]*coeff[11]-0.5740991584648069*fEdge[2]*coeff[11]+0.3314563036811939*fSkin[0]*coeff[11]-0.3314563036811939*fEdge[0]*coeff[11]-0.2651650429449552*fSkin[6]*coeff[10]+0.2651650429449552*fEdge[6]*coeff[10]-0.5740991584648069*fSkin[7]*coeff[9]-0.5740991584648069*fEdge[7]*coeff[9]+0.3314563036811939*fSkin[5]*coeff[9]-0.3314563036811939*fEdge[5]*coeff[9]-0.5740991584648069*fSkin[6]*coeff[8]-0.5740991584648069*fEdge[6]*coeff[8]+0.3314563036811939*fSkin[3]*coeff[8]-0.3314563036811939*fEdge[3]*coeff[8]; 
  edgeSurf_incr[4] = 0.4592793267718456*fSkin[6]*coeff[15]-0.4592793267718456*fEdge[6]*coeff[15]+0.4592793267718456*fSkin[7]*coeff[14]-0.4592793267718456*fEdge[7]*coeff[14]+1.259533953988537*fSkin[6]*coeff[13]+0.7292038680986265*fEdge[6]*coeff[13]-0.5740991584648069*fSkin[3]*coeff[13]+0.5740991584648069*fEdge[3]*coeff[13]+0.4592793267718456*fSkin[2]*coeff[12]-0.4592793267718456*fEdge[2]*coeff[12]+1.259533953988537*fSkin[7]*coeff[11]+0.7292038680986265*fEdge[7]*coeff[11]-0.5740991584648069*fSkin[5]*coeff[11]+0.5740991584648069*fEdge[5]*coeff[11]+0.4592793267718456*fSkin[4]*coeff[10]-0.4592793267718456*fEdge[4]*coeff[10]+1.259533953988537*fSkin[2]*coeff[9]+0.7292038680986265*fEdge[2]*coeff[9]-0.5740991584648069*fSkin[0]*coeff[9]+0.5740991584648069*fEdge[0]*coeff[9]+1.259533953988537*fSkin[4]*coeff[8]+0.7292038680986265*fEdge[4]*coeff[8]-0.5740991584648069*fSkin[1]*coeff[8]+0.5740991584648069*fEdge[1]*coeff[8]; 
  edgeSurf_incr[5] = (-0.2651650429449552*fSkin[2]*coeff[15])+0.2651650429449552*fEdge[2]*coeff[15]-0.2651650429449552*fSkin[4]*coeff[14]+0.2651650429449552*fEdge[4]*coeff[14]-0.5740991584648069*fSkin[2]*coeff[13]-0.5740991584648069*fEdge[2]*coeff[13]+0.3314563036811939*fSkin[0]*coeff[13]-0.3314563036811939*fEdge[0]*coeff[13]-0.2651650429449552*fSkin[6]*coeff[12]+0.2651650429449552*fEdge[6]*coeff[12]-0.5740991584648069*fSkin[4]*coeff[11]-0.5740991584648069*fEdge[4]*coeff[11]+0.3314563036811939*fSkin[1]*coeff[11]-0.3314563036811939*fEdge[1]*coeff[11]-0.2651650429449552*fSkin[7]*coeff[10]+0.2651650429449552*fEdge[7]*coeff[10]-0.5740991584648069*fSkin[6]*coeff[9]-0.5740991584648069*fEdge[6]*coeff[9]+0.3314563036811939*fSkin[3]*coeff[9]-0.3314563036811939*fEdge[3]*coeff[9]-0.5740991584648069*fSkin[7]*coeff[8]-0.5740991584648069*fEdge[7]*coeff[8]+0.3314563036811939*fSkin[5]*coeff[8]-0.3314563036811939*fEdge[5]*coeff[8]; 
  edgeSurf_incr[6] = 0.4592793267718456*fSkin[4]*coeff[15]-0.4592793267718456*fEdge[4]*coeff[15]+0.4592793267718456*fSkin[2]*coeff[14]-0.4592793267718456*fEdge[2]*coeff[14]+1.259533953988537*fSkin[4]*coeff[13]+0.7292038680986265*fEdge[4]*coeff[13]-0.5740991584648069*fSkin[1]*coeff[13]+0.5740991584648069*fEdge[1]*coeff[13]+0.4592793267718456*fSkin[7]*coeff[12]-0.4592793267718456*fEdge[7]*coeff[12]+1.259533953988537*fSkin[2]*coeff[11]+0.7292038680986265*fEdge[2]*coeff[11]-0.5740991584648069*fSkin[0]*coeff[11]+0.5740991584648069*fEdge[0]*coeff[11]+0.4592793267718456*fSkin[6]*coeff[10]-0.4592793267718456*fEdge[6]*coeff[10]+1.259533953988537*fSkin[7]*coeff[9]+0.7292038680986265*fEdge[7]*coeff[9]-0.5740991584648069*fSkin[5]*coeff[9]+0.5740991584648069*fEdge[5]*coeff[9]+1.259533953988537*fSkin[6]*coeff[8]+0.7292038680986265*fEdge[6]*coeff[8]-0.5740991584648069*fSkin[3]*coeff[8]+0.5740991584648069*fEdge[3]*coeff[8]; 
  edgeSurf_incr[7] = 0.4592793267718456*fSkin[2]*coeff[15]-0.4592793267718456*fEdge[2]*coeff[15]+0.4592793267718456*fSkin[4]*coeff[14]-0.4592793267718456*fEdge[4]*coeff[14]+1.259533953988537*fSkin[2]*coeff[13]+0.7292038680986265*fEdge[2]*coeff[13]-0.5740991584648069*fSkin[0]*coeff[13]+0.5740991584648069*fEdge[0]*coeff[13]+0.4592793267718456*fSkin[6]*coeff[12]-0.4592793267718456*fEdge[6]*coeff[12]+1.259533953988537*fSkin[4]*coeff[11]+0.7292038680986265*fEdge[4]*coeff[11]-0.5740991584648069*fSkin[1]*coeff[11]+0.5740991584648069*fEdge[1]*coeff[11]+0.4592793267718456*fSkin[7]*coeff[10]-0.4592793267718456*fEdge[7]*coeff[10]+1.259533953988537*fSkin[6]*coeff[9]+0.7292038680986265*fEdge[6]*coeff[9]-0.5740991584648069*fSkin[3]*coeff[9]+0.5740991584648069*fEdge[3]*coeff[9]+1.259533953988537*fSkin[7]*coeff[8]+0.7292038680986265*fEdge[7]*coeff[8]-0.5740991584648069*fSkin[5]*coeff[8]+0.5740991584648069*fEdge[5]*coeff[8]; 

  boundSurf_incr[0] = (-0.5303300858899105*fSkin[7]*coeff[15])-0.5303300858899105*fSkin[6]*coeff[14]-0.5303300858899105*fSkin[4]*coeff[12]-0.5303300858899105*fSkin[2]*coeff[10]; 
  boundSurf_incr[1] = (-0.5303300858899105*fSkin[6]*coeff[15])-0.5303300858899105*fSkin[7]*coeff[14]-0.5303300858899105*fSkin[2]*coeff[12]-0.5303300858899105*fSkin[4]*coeff[10]; 
  boundSurf_incr[2] = (-0.9185586535436913*fSkin[7]*coeff[15])-0.9185586535436913*fSkin[6]*coeff[14]+0.5303300858899105*fSkin[7]*coeff[13]-0.9185586535436913*fSkin[4]*coeff[12]+0.5303300858899105*fSkin[6]*coeff[11]-0.9185586535436913*fSkin[2]*coeff[10]+0.5303300858899105*fSkin[4]*coeff[9]+0.5303300858899105*fSkin[2]*coeff[8]; 
  boundSurf_incr[3] = (-0.5303300858899105*fSkin[4]*coeff[15])-0.5303300858899105*fSkin[2]*coeff[14]-0.5303300858899105*fSkin[7]*coeff[12]-0.5303300858899105*fSkin[6]*coeff[10]; 
  boundSurf_incr[4] = (-0.9185586535436913*fSkin[6]*coeff[15])-0.9185586535436913*fSkin[7]*coeff[14]+0.5303300858899105*fSkin[6]*coeff[13]-0.9185586535436913*fSkin[2]*coeff[12]+0.5303300858899105*fSkin[7]*coeff[11]-0.9185586535436913*fSkin[4]*coeff[10]+0.5303300858899105*fSkin[2]*coeff[9]+0.5303300858899105*fSkin[4]*coeff[8]; 
  boundSurf_incr[5] = (-0.5303300858899105*fSkin[2]*coeff[15])-0.5303300858899105*fSkin[4]*coeff[14]-0.5303300858899105*fSkin[6]*coeff[12]-0.5303300858899105*fSkin[7]*coeff[10]; 
  boundSurf_incr[6] = (-0.9185586535436913*fSkin[4]*coeff[15])-0.9185586535436913*fSkin[2]*coeff[14]+0.5303300858899105*fSkin[4]*coeff[13]-0.9185586535436913*fSkin[7]*coeff[12]+0.5303300858899105*fSkin[2]*coeff[11]-0.9185586535436913*fSkin[6]*coeff[10]+0.5303300858899105*fSkin[7]*coeff[9]+0.5303300858899105*fSkin[6]*coeff[8]; 
  boundSurf_incr[7] = (-0.9185586535436913*fSkin[2]*coeff[15])-0.9185586535436913*fSkin[4]*coeff[14]+0.5303300858899105*fSkin[2]*coeff[13]-0.9185586535436913*fSkin[6]*coeff[12]+0.5303300858899105*fSkin[4]*coeff[11]-0.9185586535436913*fSkin[7]*coeff[10]+0.5303300858899105*fSkin[6]*coeff[9]+0.5303300858899105*fSkin[7]*coeff[8]; 

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

