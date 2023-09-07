#include <gkyl_dg_diffusion_fluid_kernels.h>

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[2]*fSkin[3]+1.623797632095822*coeff[2]*fEdge[3]+0.9375*fSkin[0]*coeff[2]-0.9375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = 1.623797632095822*coeff[2]*fSkin[5]+1.623797632095822*coeff[2]*fEdge[5]+0.9375*fSkin[1]*coeff[2]-0.9375*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[2]*fSkin[6]+1.623797632095822*coeff[2]*fEdge[6]+0.9375*coeff[2]*fSkin[2]-0.9375*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[2]*fSkin[3]+2.0625*coeff[2]*fEdge[3]+1.623797632095822*fSkin[0]*coeff[2]-1.623797632095822*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = 1.623797632095822*coeff[2]*fSkin[7]+1.623797632095822*coeff[2]*fEdge[7]+0.9375*coeff[2]*fSkin[4]-0.9375*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[2]*fSkin[5]+2.0625*coeff[2]*fEdge[5]+1.623797632095822*fSkin[1]*coeff[2]-1.623797632095822*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[2]*fSkin[6]+2.0625*coeff[2]*fEdge[6]+1.623797632095822*coeff[2]*fSkin[2]-1.623797632095822*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = 3.5625*coeff[2]*fSkin[7]+2.0625*coeff[2]*fEdge[7]+1.623797632095822*coeff[2]*fSkin[4]-1.623797632095822*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = 1.5*coeff[2]*fSkin[3]; 
  boundSurf_incr[5] = 1.5*coeff[2]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[2]*fSkin[6]; 
  boundSurf_incr[7] = 1.5*coeff[2]*fSkin[7]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[2]*fSkin[3])-1.623797632095822*coeff[2]*fEdge[3]+0.9375*fSkin[0]*coeff[2]-0.9375*fEdge[0]*coeff[2]; 
  edgeSurf_incr[1] = (-1.623797632095822*coeff[2]*fSkin[5])-1.623797632095822*coeff[2]*fEdge[5]+0.9375*fSkin[1]*coeff[2]-0.9375*fEdge[1]*coeff[2]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[2]*fSkin[6])-1.623797632095822*coeff[2]*fEdge[6]+0.9375*coeff[2]*fSkin[2]-0.9375*coeff[2]*fEdge[2]; 
  edgeSurf_incr[3] = 3.5625*coeff[2]*fSkin[3]+2.0625*coeff[2]*fEdge[3]-1.623797632095822*fSkin[0]*coeff[2]+1.623797632095822*fEdge[0]*coeff[2]; 
  edgeSurf_incr[4] = (-1.623797632095822*coeff[2]*fSkin[7])-1.623797632095822*coeff[2]*fEdge[7]+0.9375*coeff[2]*fSkin[4]-0.9375*coeff[2]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[2]*fSkin[5]+2.0625*coeff[2]*fEdge[5]-1.623797632095822*fSkin[1]*coeff[2]+1.623797632095822*fEdge[1]*coeff[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[2]*fSkin[6]+2.0625*coeff[2]*fEdge[6]-1.623797632095822*coeff[2]*fSkin[2]+1.623797632095822*coeff[2]*fEdge[2]; 
  edgeSurf_incr[7] = 3.5625*coeff[2]*fSkin[7]+2.0625*coeff[2]*fEdge[7]-1.623797632095822*coeff[2]*fSkin[4]+1.623797632095822*coeff[2]*fEdge[4]; 

  boundSurf_incr[3] = 1.5*coeff[2]*fSkin[3]; 
  boundSurf_incr[5] = 1.5*coeff[2]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[2]*fSkin[6]; 
  boundSurf_incr[7] = 1.5*coeff[2]*fSkin[7]; 

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

GKYL_CU_DH double dg_diffusion_fluid_order4_boundary_surfz_3x_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[2],4.);

  double vol_incr[8] = {0.0}; 

  double edgeSurf_incr[8] = {0.0}; 
  double boundSurf_incr[8] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.2651650429449552*fSkin[7]*coeff[23])+0.2651650429449552*fEdge[7]*coeff[23]-0.2651650429449552*fSkin[6]*coeff[22]+0.2651650429449552*fEdge[6]*coeff[22]-0.2651650429449552*fSkin[5]*coeff[21]+0.2651650429449552*fEdge[5]*coeff[21]+0.5740991584648069*fSkin[7]*coeff[20]+0.5740991584648069*fEdge[7]*coeff[20]+0.3314563036811939*fSkin[4]*coeff[20]-0.3314563036811939*fEdge[4]*coeff[20]-0.2651650429449552*fSkin[3]*coeff[19]+0.2651650429449552*fEdge[3]*coeff[19]+0.5740991584648069*fSkin[6]*coeff[18]+0.5740991584648069*fEdge[6]*coeff[18]+0.3314563036811939*fSkin[2]*coeff[18]-0.3314563036811939*fEdge[2]*coeff[18]+0.5740991584648069*fSkin[5]*coeff[17]+0.5740991584648069*fEdge[5]*coeff[17]+0.3314563036811939*fSkin[1]*coeff[17]-0.3314563036811939*fEdge[1]*coeff[17]+0.5740991584648069*fSkin[3]*coeff[16]+0.5740991584648069*fEdge[3]*coeff[16]+0.3314563036811939*fSkin[0]*coeff[16]-0.3314563036811939*fEdge[0]*coeff[16]; 
  edgeSurf_incr[1] = (-0.2651650429449552*fSkin[6]*coeff[23])+0.2651650429449552*fEdge[6]*coeff[23]-0.2651650429449552*fSkin[7]*coeff[22]+0.2651650429449552*fEdge[7]*coeff[22]-0.2651650429449552*fSkin[3]*coeff[21]+0.2651650429449552*fEdge[3]*coeff[21]+0.5740991584648069*fSkin[6]*coeff[20]+0.5740991584648069*fEdge[6]*coeff[20]+0.3314563036811939*fSkin[2]*coeff[20]-0.3314563036811939*fEdge[2]*coeff[20]-0.2651650429449552*fSkin[5]*coeff[19]+0.2651650429449552*fEdge[5]*coeff[19]+0.5740991584648069*fSkin[7]*coeff[18]+0.5740991584648069*fEdge[7]*coeff[18]+0.3314563036811939*fSkin[4]*coeff[18]-0.3314563036811939*fEdge[4]*coeff[18]+0.5740991584648069*fSkin[3]*coeff[17]+0.5740991584648069*fEdge[3]*coeff[17]+0.3314563036811939*fSkin[0]*coeff[17]-0.3314563036811939*fEdge[0]*coeff[17]+0.5740991584648069*fSkin[5]*coeff[16]+0.5740991584648069*fEdge[5]*coeff[16]+0.3314563036811939*fSkin[1]*coeff[16]-0.3314563036811939*fEdge[1]*coeff[16]; 
  edgeSurf_incr[2] = (-0.2651650429449552*fSkin[5]*coeff[23])+0.2651650429449552*fEdge[5]*coeff[23]-0.2651650429449552*fSkin[3]*coeff[22]+0.2651650429449552*fEdge[3]*coeff[22]-0.2651650429449552*fSkin[7]*coeff[21]+0.2651650429449552*fEdge[7]*coeff[21]+0.5740991584648069*fSkin[5]*coeff[20]+0.5740991584648069*fEdge[5]*coeff[20]+0.3314563036811939*fSkin[1]*coeff[20]-0.3314563036811939*fEdge[1]*coeff[20]-0.2651650429449552*fSkin[6]*coeff[19]+0.2651650429449552*fEdge[6]*coeff[19]+0.5740991584648069*fSkin[3]*coeff[18]+0.5740991584648069*fEdge[3]*coeff[18]+0.3314563036811939*fSkin[0]*coeff[18]-0.3314563036811939*fEdge[0]*coeff[18]+0.5740991584648069*fSkin[7]*coeff[17]+0.5740991584648069*fEdge[7]*coeff[17]+0.3314563036811939*fSkin[4]*coeff[17]-0.3314563036811939*fEdge[4]*coeff[17]+0.5740991584648069*fSkin[6]*coeff[16]+0.5740991584648069*fEdge[6]*coeff[16]+0.3314563036811939*fSkin[2]*coeff[16]-0.3314563036811939*fEdge[2]*coeff[16]; 
  edgeSurf_incr[3] = (-0.4592793267718456*fSkin[7]*coeff[23])+0.4592793267718456*fEdge[7]*coeff[23]-0.4592793267718456*fSkin[6]*coeff[22]+0.4592793267718456*fEdge[6]*coeff[22]-0.4592793267718456*fSkin[5]*coeff[21]+0.4592793267718456*fEdge[5]*coeff[21]+1.259533953988537*fSkin[7]*coeff[20]+0.7292038680986265*fEdge[7]*coeff[20]+0.5740991584648069*fSkin[4]*coeff[20]-0.5740991584648069*fEdge[4]*coeff[20]-0.4592793267718456*fSkin[3]*coeff[19]+0.4592793267718456*fEdge[3]*coeff[19]+1.259533953988537*fSkin[6]*coeff[18]+0.7292038680986265*fEdge[6]*coeff[18]+0.5740991584648069*fSkin[2]*coeff[18]-0.5740991584648069*fEdge[2]*coeff[18]+1.259533953988537*fSkin[5]*coeff[17]+0.7292038680986265*fEdge[5]*coeff[17]+0.5740991584648069*fSkin[1]*coeff[17]-0.5740991584648069*fEdge[1]*coeff[17]+1.259533953988537*fSkin[3]*coeff[16]+0.7292038680986265*fEdge[3]*coeff[16]+0.5740991584648069*fSkin[0]*coeff[16]-0.5740991584648069*fEdge[0]*coeff[16]; 
  edgeSurf_incr[4] = (-0.2651650429449552*fSkin[3]*coeff[23])+0.2651650429449552*fEdge[3]*coeff[23]-0.2651650429449552*fSkin[5]*coeff[22]+0.2651650429449552*fEdge[5]*coeff[22]-0.2651650429449552*fSkin[6]*coeff[21]+0.2651650429449552*fEdge[6]*coeff[21]+0.5740991584648069*fSkin[3]*coeff[20]+0.5740991584648069*fEdge[3]*coeff[20]+0.3314563036811939*fSkin[0]*coeff[20]-0.3314563036811939*fEdge[0]*coeff[20]-0.2651650429449552*fSkin[7]*coeff[19]+0.2651650429449552*fEdge[7]*coeff[19]+0.5740991584648069*fSkin[5]*coeff[18]+0.5740991584648069*fEdge[5]*coeff[18]+0.3314563036811939*fSkin[1]*coeff[18]-0.3314563036811939*fEdge[1]*coeff[18]+0.5740991584648069*fSkin[6]*coeff[17]+0.5740991584648069*fEdge[6]*coeff[17]+0.3314563036811939*fSkin[2]*coeff[17]-0.3314563036811939*fEdge[2]*coeff[17]+0.5740991584648069*fSkin[7]*coeff[16]+0.5740991584648069*fEdge[7]*coeff[16]+0.3314563036811939*fSkin[4]*coeff[16]-0.3314563036811939*fEdge[4]*coeff[16]; 
  edgeSurf_incr[5] = (-0.4592793267718456*fSkin[6]*coeff[23])+0.4592793267718456*fEdge[6]*coeff[23]-0.4592793267718456*fSkin[7]*coeff[22]+0.4592793267718456*fEdge[7]*coeff[22]-0.4592793267718456*fSkin[3]*coeff[21]+0.4592793267718456*fEdge[3]*coeff[21]+1.259533953988537*fSkin[6]*coeff[20]+0.7292038680986265*fEdge[6]*coeff[20]+0.5740991584648069*fSkin[2]*coeff[20]-0.5740991584648069*fEdge[2]*coeff[20]-0.4592793267718456*fSkin[5]*coeff[19]+0.4592793267718456*fEdge[5]*coeff[19]+1.259533953988537*fSkin[7]*coeff[18]+0.7292038680986265*fEdge[7]*coeff[18]+0.5740991584648069*fSkin[4]*coeff[18]-0.5740991584648069*fEdge[4]*coeff[18]+1.259533953988537*fSkin[3]*coeff[17]+0.7292038680986265*fEdge[3]*coeff[17]+0.5740991584648069*fSkin[0]*coeff[17]-0.5740991584648069*fEdge[0]*coeff[17]+1.259533953988537*fSkin[5]*coeff[16]+0.7292038680986265*fEdge[5]*coeff[16]+0.5740991584648069*fSkin[1]*coeff[16]-0.5740991584648069*fEdge[1]*coeff[16]; 
  edgeSurf_incr[6] = (-0.4592793267718456*fSkin[5]*coeff[23])+0.4592793267718456*fEdge[5]*coeff[23]-0.4592793267718456*fSkin[3]*coeff[22]+0.4592793267718456*fEdge[3]*coeff[22]-0.4592793267718456*fSkin[7]*coeff[21]+0.4592793267718456*fEdge[7]*coeff[21]+1.259533953988537*fSkin[5]*coeff[20]+0.7292038680986265*fEdge[5]*coeff[20]+0.5740991584648069*fSkin[1]*coeff[20]-0.5740991584648069*fEdge[1]*coeff[20]-0.4592793267718456*fSkin[6]*coeff[19]+0.4592793267718456*fEdge[6]*coeff[19]+1.259533953988537*fSkin[3]*coeff[18]+0.7292038680986265*fEdge[3]*coeff[18]+0.5740991584648069*fSkin[0]*coeff[18]-0.5740991584648069*fEdge[0]*coeff[18]+1.259533953988537*fSkin[7]*coeff[17]+0.7292038680986265*fEdge[7]*coeff[17]+0.5740991584648069*fSkin[4]*coeff[17]-0.5740991584648069*fEdge[4]*coeff[17]+1.259533953988537*fSkin[6]*coeff[16]+0.7292038680986265*fEdge[6]*coeff[16]+0.5740991584648069*fSkin[2]*coeff[16]-0.5740991584648069*fEdge[2]*coeff[16]; 
  edgeSurf_incr[7] = (-0.4592793267718456*fSkin[3]*coeff[23])+0.4592793267718456*fEdge[3]*coeff[23]-0.4592793267718456*fSkin[5]*coeff[22]+0.4592793267718456*fEdge[5]*coeff[22]-0.4592793267718456*fSkin[6]*coeff[21]+0.4592793267718456*fEdge[6]*coeff[21]+1.259533953988537*fSkin[3]*coeff[20]+0.7292038680986265*fEdge[3]*coeff[20]+0.5740991584648069*fSkin[0]*coeff[20]-0.5740991584648069*fEdge[0]*coeff[20]-0.4592793267718456*fSkin[7]*coeff[19]+0.4592793267718456*fEdge[7]*coeff[19]+1.259533953988537*fSkin[5]*coeff[18]+0.7292038680986265*fEdge[5]*coeff[18]+0.5740991584648069*fSkin[1]*coeff[18]-0.5740991584648069*fEdge[1]*coeff[18]+1.259533953988537*fSkin[6]*coeff[17]+0.7292038680986265*fEdge[6]*coeff[17]+0.5740991584648069*fSkin[2]*coeff[17]-0.5740991584648069*fEdge[2]*coeff[17]+1.259533953988537*fSkin[7]*coeff[16]+0.7292038680986265*fEdge[7]*coeff[16]+0.5740991584648069*fSkin[4]*coeff[16]-0.5740991584648069*fEdge[4]*coeff[16]; 

  boundSurf_incr[0] = (-0.5303300858899105*fSkin[7]*coeff[23])-0.5303300858899105*fSkin[6]*coeff[22]-0.5303300858899105*fSkin[5]*coeff[21]-0.5303300858899105*fSkin[3]*coeff[19]; 
  boundSurf_incr[1] = (-0.5303300858899105*fSkin[6]*coeff[23])-0.5303300858899105*fSkin[7]*coeff[22]-0.5303300858899105*fSkin[3]*coeff[21]-0.5303300858899105*fSkin[5]*coeff[19]; 
  boundSurf_incr[2] = (-0.5303300858899105*fSkin[5]*coeff[23])-0.5303300858899105*fSkin[3]*coeff[22]-0.5303300858899105*fSkin[7]*coeff[21]-0.5303300858899105*fSkin[6]*coeff[19]; 
  boundSurf_incr[3] = 0.9185586535436913*fSkin[7]*coeff[23]+0.9185586535436913*fSkin[6]*coeff[22]+0.9185586535436913*fSkin[5]*coeff[21]+0.5303300858899105*fSkin[7]*coeff[20]+0.9185586535436913*fSkin[3]*coeff[19]+0.5303300858899105*fSkin[6]*coeff[18]+0.5303300858899105*fSkin[5]*coeff[17]+0.5303300858899105*fSkin[3]*coeff[16]; 
  boundSurf_incr[4] = (-0.5303300858899105*fSkin[3]*coeff[23])-0.5303300858899105*fSkin[5]*coeff[22]-0.5303300858899105*fSkin[6]*coeff[21]-0.5303300858899105*fSkin[7]*coeff[19]; 
  boundSurf_incr[5] = 0.9185586535436913*fSkin[6]*coeff[23]+0.9185586535436913*fSkin[7]*coeff[22]+0.9185586535436913*fSkin[3]*coeff[21]+0.5303300858899105*fSkin[6]*coeff[20]+0.9185586535436913*fSkin[5]*coeff[19]+0.5303300858899105*fSkin[7]*coeff[18]+0.5303300858899105*fSkin[3]*coeff[17]+0.5303300858899105*fSkin[5]*coeff[16]; 
  boundSurf_incr[6] = 0.9185586535436913*fSkin[5]*coeff[23]+0.9185586535436913*fSkin[3]*coeff[22]+0.9185586535436913*fSkin[7]*coeff[21]+0.5303300858899105*fSkin[5]*coeff[20]+0.9185586535436913*fSkin[6]*coeff[19]+0.5303300858899105*fSkin[3]*coeff[18]+0.5303300858899105*fSkin[7]*coeff[17]+0.5303300858899105*fSkin[6]*coeff[16]; 
  boundSurf_incr[7] = 0.9185586535436913*fSkin[3]*coeff[23]+0.9185586535436913*fSkin[5]*coeff[22]+0.9185586535436913*fSkin[6]*coeff[21]+0.5303300858899105*fSkin[3]*coeff[20]+0.9185586535436913*fSkin[7]*coeff[19]+0.5303300858899105*fSkin[5]*coeff[18]+0.5303300858899105*fSkin[6]*coeff[17]+0.5303300858899105*fSkin[7]*coeff[16]; 

  } else { 

  edgeSurf_incr[0] = (-0.2651650429449552*fSkin[7]*coeff[23])+0.2651650429449552*fEdge[7]*coeff[23]-0.2651650429449552*fSkin[6]*coeff[22]+0.2651650429449552*fEdge[6]*coeff[22]-0.2651650429449552*fSkin[5]*coeff[21]+0.2651650429449552*fEdge[5]*coeff[21]-0.5740991584648069*fSkin[7]*coeff[20]-0.5740991584648069*fEdge[7]*coeff[20]+0.3314563036811939*fSkin[4]*coeff[20]-0.3314563036811939*fEdge[4]*coeff[20]-0.2651650429449552*fSkin[3]*coeff[19]+0.2651650429449552*fEdge[3]*coeff[19]-0.5740991584648069*fSkin[6]*coeff[18]-0.5740991584648069*fEdge[6]*coeff[18]+0.3314563036811939*fSkin[2]*coeff[18]-0.3314563036811939*fEdge[2]*coeff[18]-0.5740991584648069*fSkin[5]*coeff[17]-0.5740991584648069*fEdge[5]*coeff[17]+0.3314563036811939*fSkin[1]*coeff[17]-0.3314563036811939*fEdge[1]*coeff[17]-0.5740991584648069*fSkin[3]*coeff[16]-0.5740991584648069*fEdge[3]*coeff[16]+0.3314563036811939*fSkin[0]*coeff[16]-0.3314563036811939*fEdge[0]*coeff[16]; 
  edgeSurf_incr[1] = (-0.2651650429449552*fSkin[6]*coeff[23])+0.2651650429449552*fEdge[6]*coeff[23]-0.2651650429449552*fSkin[7]*coeff[22]+0.2651650429449552*fEdge[7]*coeff[22]-0.2651650429449552*fSkin[3]*coeff[21]+0.2651650429449552*fEdge[3]*coeff[21]-0.5740991584648069*fSkin[6]*coeff[20]-0.5740991584648069*fEdge[6]*coeff[20]+0.3314563036811939*fSkin[2]*coeff[20]-0.3314563036811939*fEdge[2]*coeff[20]-0.2651650429449552*fSkin[5]*coeff[19]+0.2651650429449552*fEdge[5]*coeff[19]-0.5740991584648069*fSkin[7]*coeff[18]-0.5740991584648069*fEdge[7]*coeff[18]+0.3314563036811939*fSkin[4]*coeff[18]-0.3314563036811939*fEdge[4]*coeff[18]-0.5740991584648069*fSkin[3]*coeff[17]-0.5740991584648069*fEdge[3]*coeff[17]+0.3314563036811939*fSkin[0]*coeff[17]-0.3314563036811939*fEdge[0]*coeff[17]-0.5740991584648069*fSkin[5]*coeff[16]-0.5740991584648069*fEdge[5]*coeff[16]+0.3314563036811939*fSkin[1]*coeff[16]-0.3314563036811939*fEdge[1]*coeff[16]; 
  edgeSurf_incr[2] = (-0.2651650429449552*fSkin[5]*coeff[23])+0.2651650429449552*fEdge[5]*coeff[23]-0.2651650429449552*fSkin[3]*coeff[22]+0.2651650429449552*fEdge[3]*coeff[22]-0.2651650429449552*fSkin[7]*coeff[21]+0.2651650429449552*fEdge[7]*coeff[21]-0.5740991584648069*fSkin[5]*coeff[20]-0.5740991584648069*fEdge[5]*coeff[20]+0.3314563036811939*fSkin[1]*coeff[20]-0.3314563036811939*fEdge[1]*coeff[20]-0.2651650429449552*fSkin[6]*coeff[19]+0.2651650429449552*fEdge[6]*coeff[19]-0.5740991584648069*fSkin[3]*coeff[18]-0.5740991584648069*fEdge[3]*coeff[18]+0.3314563036811939*fSkin[0]*coeff[18]-0.3314563036811939*fEdge[0]*coeff[18]-0.5740991584648069*fSkin[7]*coeff[17]-0.5740991584648069*fEdge[7]*coeff[17]+0.3314563036811939*fSkin[4]*coeff[17]-0.3314563036811939*fEdge[4]*coeff[17]-0.5740991584648069*fSkin[6]*coeff[16]-0.5740991584648069*fEdge[6]*coeff[16]+0.3314563036811939*fSkin[2]*coeff[16]-0.3314563036811939*fEdge[2]*coeff[16]; 
  edgeSurf_incr[3] = 0.4592793267718456*fSkin[7]*coeff[23]-0.4592793267718456*fEdge[7]*coeff[23]+0.4592793267718456*fSkin[6]*coeff[22]-0.4592793267718456*fEdge[6]*coeff[22]+0.4592793267718456*fSkin[5]*coeff[21]-0.4592793267718456*fEdge[5]*coeff[21]+1.259533953988537*fSkin[7]*coeff[20]+0.7292038680986265*fEdge[7]*coeff[20]-0.5740991584648069*fSkin[4]*coeff[20]+0.5740991584648069*fEdge[4]*coeff[20]+0.4592793267718456*fSkin[3]*coeff[19]-0.4592793267718456*fEdge[3]*coeff[19]+1.259533953988537*fSkin[6]*coeff[18]+0.7292038680986265*fEdge[6]*coeff[18]-0.5740991584648069*fSkin[2]*coeff[18]+0.5740991584648069*fEdge[2]*coeff[18]+1.259533953988537*fSkin[5]*coeff[17]+0.7292038680986265*fEdge[5]*coeff[17]-0.5740991584648069*fSkin[1]*coeff[17]+0.5740991584648069*fEdge[1]*coeff[17]+1.259533953988537*fSkin[3]*coeff[16]+0.7292038680986265*fEdge[3]*coeff[16]-0.5740991584648069*fSkin[0]*coeff[16]+0.5740991584648069*fEdge[0]*coeff[16]; 
  edgeSurf_incr[4] = (-0.2651650429449552*fSkin[3]*coeff[23])+0.2651650429449552*fEdge[3]*coeff[23]-0.2651650429449552*fSkin[5]*coeff[22]+0.2651650429449552*fEdge[5]*coeff[22]-0.2651650429449552*fSkin[6]*coeff[21]+0.2651650429449552*fEdge[6]*coeff[21]-0.5740991584648069*fSkin[3]*coeff[20]-0.5740991584648069*fEdge[3]*coeff[20]+0.3314563036811939*fSkin[0]*coeff[20]-0.3314563036811939*fEdge[0]*coeff[20]-0.2651650429449552*fSkin[7]*coeff[19]+0.2651650429449552*fEdge[7]*coeff[19]-0.5740991584648069*fSkin[5]*coeff[18]-0.5740991584648069*fEdge[5]*coeff[18]+0.3314563036811939*fSkin[1]*coeff[18]-0.3314563036811939*fEdge[1]*coeff[18]-0.5740991584648069*fSkin[6]*coeff[17]-0.5740991584648069*fEdge[6]*coeff[17]+0.3314563036811939*fSkin[2]*coeff[17]-0.3314563036811939*fEdge[2]*coeff[17]-0.5740991584648069*fSkin[7]*coeff[16]-0.5740991584648069*fEdge[7]*coeff[16]+0.3314563036811939*fSkin[4]*coeff[16]-0.3314563036811939*fEdge[4]*coeff[16]; 
  edgeSurf_incr[5] = 0.4592793267718456*fSkin[6]*coeff[23]-0.4592793267718456*fEdge[6]*coeff[23]+0.4592793267718456*fSkin[7]*coeff[22]-0.4592793267718456*fEdge[7]*coeff[22]+0.4592793267718456*fSkin[3]*coeff[21]-0.4592793267718456*fEdge[3]*coeff[21]+1.259533953988537*fSkin[6]*coeff[20]+0.7292038680986265*fEdge[6]*coeff[20]-0.5740991584648069*fSkin[2]*coeff[20]+0.5740991584648069*fEdge[2]*coeff[20]+0.4592793267718456*fSkin[5]*coeff[19]-0.4592793267718456*fEdge[5]*coeff[19]+1.259533953988537*fSkin[7]*coeff[18]+0.7292038680986265*fEdge[7]*coeff[18]-0.5740991584648069*fSkin[4]*coeff[18]+0.5740991584648069*fEdge[4]*coeff[18]+1.259533953988537*fSkin[3]*coeff[17]+0.7292038680986265*fEdge[3]*coeff[17]-0.5740991584648069*fSkin[0]*coeff[17]+0.5740991584648069*fEdge[0]*coeff[17]+1.259533953988537*fSkin[5]*coeff[16]+0.7292038680986265*fEdge[5]*coeff[16]-0.5740991584648069*fSkin[1]*coeff[16]+0.5740991584648069*fEdge[1]*coeff[16]; 
  edgeSurf_incr[6] = 0.4592793267718456*fSkin[5]*coeff[23]-0.4592793267718456*fEdge[5]*coeff[23]+0.4592793267718456*fSkin[3]*coeff[22]-0.4592793267718456*fEdge[3]*coeff[22]+0.4592793267718456*fSkin[7]*coeff[21]-0.4592793267718456*fEdge[7]*coeff[21]+1.259533953988537*fSkin[5]*coeff[20]+0.7292038680986265*fEdge[5]*coeff[20]-0.5740991584648069*fSkin[1]*coeff[20]+0.5740991584648069*fEdge[1]*coeff[20]+0.4592793267718456*fSkin[6]*coeff[19]-0.4592793267718456*fEdge[6]*coeff[19]+1.259533953988537*fSkin[3]*coeff[18]+0.7292038680986265*fEdge[3]*coeff[18]-0.5740991584648069*fSkin[0]*coeff[18]+0.5740991584648069*fEdge[0]*coeff[18]+1.259533953988537*fSkin[7]*coeff[17]+0.7292038680986265*fEdge[7]*coeff[17]-0.5740991584648069*fSkin[4]*coeff[17]+0.5740991584648069*fEdge[4]*coeff[17]+1.259533953988537*fSkin[6]*coeff[16]+0.7292038680986265*fEdge[6]*coeff[16]-0.5740991584648069*fSkin[2]*coeff[16]+0.5740991584648069*fEdge[2]*coeff[16]; 
  edgeSurf_incr[7] = 0.4592793267718456*fSkin[3]*coeff[23]-0.4592793267718456*fEdge[3]*coeff[23]+0.4592793267718456*fSkin[5]*coeff[22]-0.4592793267718456*fEdge[5]*coeff[22]+0.4592793267718456*fSkin[6]*coeff[21]-0.4592793267718456*fEdge[6]*coeff[21]+1.259533953988537*fSkin[3]*coeff[20]+0.7292038680986265*fEdge[3]*coeff[20]-0.5740991584648069*fSkin[0]*coeff[20]+0.5740991584648069*fEdge[0]*coeff[20]+0.4592793267718456*fSkin[7]*coeff[19]-0.4592793267718456*fEdge[7]*coeff[19]+1.259533953988537*fSkin[5]*coeff[18]+0.7292038680986265*fEdge[5]*coeff[18]-0.5740991584648069*fSkin[1]*coeff[18]+0.5740991584648069*fEdge[1]*coeff[18]+1.259533953988537*fSkin[6]*coeff[17]+0.7292038680986265*fEdge[6]*coeff[17]-0.5740991584648069*fSkin[2]*coeff[17]+0.5740991584648069*fEdge[2]*coeff[17]+1.259533953988537*fSkin[7]*coeff[16]+0.7292038680986265*fEdge[7]*coeff[16]-0.5740991584648069*fSkin[4]*coeff[16]+0.5740991584648069*fEdge[4]*coeff[16]; 

  boundSurf_incr[0] = (-0.5303300858899105*fSkin[7]*coeff[23])-0.5303300858899105*fSkin[6]*coeff[22]-0.5303300858899105*fSkin[5]*coeff[21]-0.5303300858899105*fSkin[3]*coeff[19]; 
  boundSurf_incr[1] = (-0.5303300858899105*fSkin[6]*coeff[23])-0.5303300858899105*fSkin[7]*coeff[22]-0.5303300858899105*fSkin[3]*coeff[21]-0.5303300858899105*fSkin[5]*coeff[19]; 
  boundSurf_incr[2] = (-0.5303300858899105*fSkin[5]*coeff[23])-0.5303300858899105*fSkin[3]*coeff[22]-0.5303300858899105*fSkin[7]*coeff[21]-0.5303300858899105*fSkin[6]*coeff[19]; 
  boundSurf_incr[3] = (-0.9185586535436913*fSkin[7]*coeff[23])-0.9185586535436913*fSkin[6]*coeff[22]-0.9185586535436913*fSkin[5]*coeff[21]+0.5303300858899105*fSkin[7]*coeff[20]-0.9185586535436913*fSkin[3]*coeff[19]+0.5303300858899105*fSkin[6]*coeff[18]+0.5303300858899105*fSkin[5]*coeff[17]+0.5303300858899105*fSkin[3]*coeff[16]; 
  boundSurf_incr[4] = (-0.5303300858899105*fSkin[3]*coeff[23])-0.5303300858899105*fSkin[5]*coeff[22]-0.5303300858899105*fSkin[6]*coeff[21]-0.5303300858899105*fSkin[7]*coeff[19]; 
  boundSurf_incr[5] = (-0.9185586535436913*fSkin[6]*coeff[23])-0.9185586535436913*fSkin[7]*coeff[22]-0.9185586535436913*fSkin[3]*coeff[21]+0.5303300858899105*fSkin[6]*coeff[20]-0.9185586535436913*fSkin[5]*coeff[19]+0.5303300858899105*fSkin[7]*coeff[18]+0.5303300858899105*fSkin[3]*coeff[17]+0.5303300858899105*fSkin[5]*coeff[16]; 
  boundSurf_incr[6] = (-0.9185586535436913*fSkin[5]*coeff[23])-0.9185586535436913*fSkin[3]*coeff[22]-0.9185586535436913*fSkin[7]*coeff[21]+0.5303300858899105*fSkin[5]*coeff[20]-0.9185586535436913*fSkin[6]*coeff[19]+0.5303300858899105*fSkin[3]*coeff[18]+0.5303300858899105*fSkin[7]*coeff[17]+0.5303300858899105*fSkin[6]*coeff[16]; 
  boundSurf_incr[7] = (-0.9185586535436913*fSkin[3]*coeff[23])-0.9185586535436913*fSkin[5]*coeff[22]-0.9185586535436913*fSkin[6]*coeff[21]+0.5303300858899105*fSkin[3]*coeff[20]-0.9185586535436913*fSkin[7]*coeff[19]+0.5303300858899105*fSkin[5]*coeff[18]+0.5303300858899105*fSkin[6]*coeff[17]+0.5303300858899105*fSkin[7]*coeff[16]; 

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

