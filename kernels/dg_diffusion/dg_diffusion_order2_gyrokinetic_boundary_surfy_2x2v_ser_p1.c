#include <gkyl_dg_diffusion_kernels.h>

GKYL_CU_DH double dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[24] = {0.0}; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5412658773652741*coeff[1]*fSkin[2])-0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = (-0.5412658773652741*coeff[1]*fSkin[5])-0.5412658773652741*coeff[1]*fEdge[5]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = (-1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]-1.407291281149712*fSkin[0]*coeff[1]+0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = (-0.5412658773652741*coeff[1]*fSkin[7])-0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.5412658773652741*coeff[1]*fSkin[9])-0.5412658773652741*coeff[1]*fEdge[9]-0.5625*coeff[1]*fSkin[4]+0.5625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[1]*fSkin[5])-0.4375*coeff[1]*fEdge[5]-1.407291281149712*coeff[1]*fSkin[1]+0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = (-0.5412658773652741*coeff[1]*fSkin[11])-0.5412658773652741*coeff[1]*fEdge[11]-0.5625*coeff[1]*fSkin[6]+0.5625*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]-1.407291281149712*coeff[1]*fSkin[3]+0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = (-0.5412658773652741*coeff[1]*fSkin[12])-0.5412658773652741*coeff[1]*fEdge[12]-0.5625*coeff[1]*fSkin[8]+0.5625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = (-1.4375*coeff[1]*fSkin[9])-0.4375*coeff[1]*fEdge[9]-1.407291281149712*coeff[1]*fSkin[4]+0.5412658773652739*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = (-0.5412658773652741*coeff[1]*fSkin[14])-0.5412658773652741*coeff[1]*fEdge[14]-0.5625*coeff[1]*fSkin[10]+0.5625*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[1]*fSkin[11])-0.4375*coeff[1]*fEdge[11]-1.407291281149712*coeff[1]*fSkin[6]+0.5412658773652739*coeff[1]*fEdge[6]; 
  edgeSurf_incr[12] = (-1.4375*coeff[1]*fSkin[12])-0.4375*coeff[1]*fEdge[12]-1.407291281149712*coeff[1]*fSkin[8]+0.5412658773652739*coeff[1]*fEdge[8]; 
  edgeSurf_incr[13] = (-0.5412658773652741*coeff[1]*fSkin[15])-0.5412658773652741*coeff[1]*fEdge[15]-0.5625*coeff[1]*fSkin[13]+0.5625*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = (-1.4375*coeff[1]*fSkin[14])-0.4375*coeff[1]*fEdge[14]-1.407291281149712*coeff[1]*fSkin[10]+0.5412658773652739*coeff[1]*fEdge[10]; 
  edgeSurf_incr[15] = (-1.4375*coeff[1]*fSkin[15])-0.4375*coeff[1]*fEdge[15]-1.407291281149712*coeff[1]*fSkin[13]+0.5412658773652739*coeff[1]*fEdge[13]; 
  edgeSurf_incr[16] = (-0.5412658773652742*coeff[1]*fSkin[18])-0.5412658773652742*coeff[1]*fEdge[18]-0.5625*coeff[1]*fSkin[16]+0.5625*coeff[1]*fEdge[16]; 
  edgeSurf_incr[17] = (-0.5412658773652742*coeff[1]*fSkin[20])-0.5412658773652742*coeff[1]*fEdge[20]-0.5625*coeff[1]*fSkin[17]+0.5625*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = (-1.4375*coeff[1]*fSkin[18])-0.4375*coeff[1]*fEdge[18]-1.407291281149713*coeff[1]*fSkin[16]+0.5412658773652742*coeff[1]*fEdge[16]; 
  edgeSurf_incr[19] = (-0.5412658773652742*coeff[1]*fSkin[22])-0.5412658773652742*coeff[1]*fEdge[22]-0.5625*coeff[1]*fSkin[19]+0.5625*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = (-1.4375*coeff[1]*fSkin[20])-0.4375*coeff[1]*fEdge[20]-1.407291281149713*coeff[1]*fSkin[17]+0.5412658773652742*coeff[1]*fEdge[17]; 
  edgeSurf_incr[21] = (-0.5412658773652742*coeff[1]*fSkin[23])-0.5412658773652742*coeff[1]*fEdge[23]-0.5625*coeff[1]*fSkin[21]+0.5625*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = (-1.4375*coeff[1]*fSkin[22])-0.4375*coeff[1]*fEdge[22]-1.407291281149713*coeff[1]*fSkin[19]+0.5412658773652742*coeff[1]*fEdge[19]; 
  edgeSurf_incr[23] = (-1.4375*coeff[1]*fSkin[23])-0.4375*coeff[1]*fEdge[23]-1.407291281149713*coeff[1]*fSkin[21]+0.5412658773652742*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[0]*coeff[1]-1.0*coeff[1]*fSkin[2]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[1]*fSkin[1]-1.0*coeff[1]*fSkin[5]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[1]*fSkin[3]-1.0*coeff[1]*fSkin[7]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[1]*fSkin[4]-1.0*coeff[1]*fSkin[9]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[1]*fSkin[6]-1.0*coeff[1]*fSkin[11]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[1]*fSkin[8]-1.0*coeff[1]*fSkin[12]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[1]*fSkin[10]-1.0*coeff[1]*fSkin[14]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[1]*fSkin[13]-1.0*coeff[1]*fSkin[15]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[1]*fSkin[16]-1.0*coeff[1]*fSkin[18]; 
  boundSurf_incr[20] = 0.8660254037844387*coeff[1]*fSkin[17]-1.0*coeff[1]*fSkin[20]; 
  boundSurf_incr[22] = 0.8660254037844387*coeff[1]*fSkin[19]-1.0*coeff[1]*fSkin[22]; 
  boundSurf_incr[23] = 0.8660254037844387*coeff[1]*fSkin[21]-1.0*coeff[1]*fSkin[23]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[1]*fSkin[2]+0.5412658773652741*coeff[1]*fEdge[2]-0.5625*fSkin[0]*coeff[1]+0.5625*fEdge[0]*coeff[1]; 
  edgeSurf_incr[1] = 0.5412658773652741*coeff[1]*fSkin[5]+0.5412658773652741*coeff[1]*fEdge[5]-0.5625*coeff[1]*fSkin[1]+0.5625*coeff[1]*fEdge[1]; 
  edgeSurf_incr[2] = (-1.4375*coeff[1]*fSkin[2])-0.4375*coeff[1]*fEdge[2]+1.407291281149712*fSkin[0]*coeff[1]-0.5412658773652739*fEdge[0]*coeff[1]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[1]*fSkin[7]+0.5412658773652741*coeff[1]*fEdge[7]-0.5625*coeff[1]*fSkin[3]+0.5625*coeff[1]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[1]*fSkin[9]+0.5412658773652741*coeff[1]*fEdge[9]-0.5625*coeff[1]*fSkin[4]+0.5625*coeff[1]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[1]*fSkin[5])-0.4375*coeff[1]*fEdge[5]+1.407291281149712*coeff[1]*fSkin[1]-0.5412658773652739*coeff[1]*fEdge[1]; 
  edgeSurf_incr[6] = 0.5412658773652741*coeff[1]*fSkin[11]+0.5412658773652741*coeff[1]*fEdge[11]-0.5625*coeff[1]*fSkin[6]+0.5625*coeff[1]*fEdge[6]; 
  edgeSurf_incr[7] = (-1.4375*coeff[1]*fSkin[7])-0.4375*coeff[1]*fEdge[7]+1.407291281149712*coeff[1]*fSkin[3]-0.5412658773652739*coeff[1]*fEdge[3]; 
  edgeSurf_incr[8] = 0.5412658773652741*coeff[1]*fSkin[12]+0.5412658773652741*coeff[1]*fEdge[12]-0.5625*coeff[1]*fSkin[8]+0.5625*coeff[1]*fEdge[8]; 
  edgeSurf_incr[9] = (-1.4375*coeff[1]*fSkin[9])-0.4375*coeff[1]*fEdge[9]+1.407291281149712*coeff[1]*fSkin[4]-0.5412658773652739*coeff[1]*fEdge[4]; 
  edgeSurf_incr[10] = 0.5412658773652741*coeff[1]*fSkin[14]+0.5412658773652741*coeff[1]*fEdge[14]-0.5625*coeff[1]*fSkin[10]+0.5625*coeff[1]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[1]*fSkin[11])-0.4375*coeff[1]*fEdge[11]+1.407291281149712*coeff[1]*fSkin[6]-0.5412658773652739*coeff[1]*fEdge[6]; 
  edgeSurf_incr[12] = (-1.4375*coeff[1]*fSkin[12])-0.4375*coeff[1]*fEdge[12]+1.407291281149712*coeff[1]*fSkin[8]-0.5412658773652739*coeff[1]*fEdge[8]; 
  edgeSurf_incr[13] = 0.5412658773652741*coeff[1]*fSkin[15]+0.5412658773652741*coeff[1]*fEdge[15]-0.5625*coeff[1]*fSkin[13]+0.5625*coeff[1]*fEdge[13]; 
  edgeSurf_incr[14] = (-1.4375*coeff[1]*fSkin[14])-0.4375*coeff[1]*fEdge[14]+1.407291281149712*coeff[1]*fSkin[10]-0.5412658773652739*coeff[1]*fEdge[10]; 
  edgeSurf_incr[15] = (-1.4375*coeff[1]*fSkin[15])-0.4375*coeff[1]*fEdge[15]+1.407291281149712*coeff[1]*fSkin[13]-0.5412658773652739*coeff[1]*fEdge[13]; 
  edgeSurf_incr[16] = 0.5412658773652742*coeff[1]*fSkin[18]+0.5412658773652742*coeff[1]*fEdge[18]-0.5625*coeff[1]*fSkin[16]+0.5625*coeff[1]*fEdge[16]; 
  edgeSurf_incr[17] = 0.5412658773652742*coeff[1]*fSkin[20]+0.5412658773652742*coeff[1]*fEdge[20]-0.5625*coeff[1]*fSkin[17]+0.5625*coeff[1]*fEdge[17]; 
  edgeSurf_incr[18] = (-1.4375*coeff[1]*fSkin[18])-0.4375*coeff[1]*fEdge[18]+1.407291281149713*coeff[1]*fSkin[16]-0.5412658773652742*coeff[1]*fEdge[16]; 
  edgeSurf_incr[19] = 0.5412658773652742*coeff[1]*fSkin[22]+0.5412658773652742*coeff[1]*fEdge[22]-0.5625*coeff[1]*fSkin[19]+0.5625*coeff[1]*fEdge[19]; 
  edgeSurf_incr[20] = (-1.4375*coeff[1]*fSkin[20])-0.4375*coeff[1]*fEdge[20]+1.407291281149713*coeff[1]*fSkin[17]-0.5412658773652742*coeff[1]*fEdge[17]; 
  edgeSurf_incr[21] = 0.5412658773652742*coeff[1]*fSkin[23]+0.5412658773652742*coeff[1]*fEdge[23]-0.5625*coeff[1]*fSkin[21]+0.5625*coeff[1]*fEdge[21]; 
  edgeSurf_incr[22] = (-1.4375*coeff[1]*fSkin[22])-0.4375*coeff[1]*fEdge[22]+1.407291281149713*coeff[1]*fSkin[19]-0.5412658773652742*coeff[1]*fEdge[19]; 
  edgeSurf_incr[23] = (-1.4375*coeff[1]*fSkin[23])-0.4375*coeff[1]*fEdge[23]+1.407291281149713*coeff[1]*fSkin[21]-0.5412658773652742*coeff[1]*fEdge[21]; 

  boundSurf_incr[2] = (-1.0*coeff[1]*fSkin[2])-0.8660254037844386*fSkin[0]*coeff[1]; 
  boundSurf_incr[5] = (-1.0*coeff[1]*fSkin[5])-0.8660254037844386*coeff[1]*fSkin[1]; 
  boundSurf_incr[7] = (-1.0*coeff[1]*fSkin[7])-0.8660254037844386*coeff[1]*fSkin[3]; 
  boundSurf_incr[9] = (-1.0*coeff[1]*fSkin[9])-0.8660254037844386*coeff[1]*fSkin[4]; 
  boundSurf_incr[11] = (-1.0*coeff[1]*fSkin[11])-0.8660254037844386*coeff[1]*fSkin[6]; 
  boundSurf_incr[12] = (-1.0*coeff[1]*fSkin[12])-0.8660254037844386*coeff[1]*fSkin[8]; 
  boundSurf_incr[14] = (-1.0*coeff[1]*fSkin[14])-0.8660254037844386*coeff[1]*fSkin[10]; 
  boundSurf_incr[15] = (-1.0*coeff[1]*fSkin[15])-0.8660254037844386*coeff[1]*fSkin[13]; 
  boundSurf_incr[18] = (-1.0*coeff[1]*fSkin[18])-0.8660254037844387*coeff[1]*fSkin[16]; 
  boundSurf_incr[20] = (-1.0*coeff[1]*fSkin[20])-0.8660254037844387*coeff[1]*fSkin[17]; 
  boundSurf_incr[22] = (-1.0*coeff[1]*fSkin[22])-0.8660254037844387*coeff[1]*fSkin[19]; 
  boundSurf_incr[23] = (-1.0*coeff[1]*fSkin[23])-0.8660254037844387*coeff[1]*fSkin[21]; 

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

  }

  return 0.;
}

GKYL_CU_DH double dg_diffusion_order2_gyrokinetic_boundary_surfy_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[1],2.);

  double vol_incr[24] = {0.0}; 
  vol_incr[2] = 1.5*fSkin[1]*coeff[7]+1.5*fSkin[0]*coeff[6]; 
  vol_incr[5] = 1.5*fSkin[0]*coeff[7]+1.5*fSkin[1]*coeff[6]; 
  vol_incr[7] = 1.5*fSkin[6]*coeff[7]+1.5*fSkin[3]*coeff[6]; 
  vol_incr[9] = 1.5*coeff[7]*fSkin[8]+1.5*fSkin[4]*coeff[6]; 
  vol_incr[11] = 1.5*fSkin[3]*coeff[7]+1.5*coeff[6]*fSkin[6]; 
  vol_incr[12] = 1.5*coeff[6]*fSkin[8]+1.5*fSkin[4]*coeff[7]; 
  vol_incr[14] = 1.5*coeff[7]*fSkin[13]+1.5*coeff[6]*fSkin[10]; 
  vol_incr[15] = 1.5*coeff[6]*fSkin[13]+1.5*coeff[7]*fSkin[10]; 
  vol_incr[18] = 1.5*coeff[7]*fSkin[17]+1.5*coeff[6]*fSkin[16]; 
  vol_incr[20] = 1.5*coeff[6]*fSkin[17]+1.5*coeff[7]*fSkin[16]; 
  vol_incr[22] = 1.5*coeff[7]*fSkin[21]+1.5*coeff[6]*fSkin[19]; 
  vol_incr[23] = 1.5*coeff[6]*fSkin[21]+1.5*coeff[7]*fSkin[19]; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.270632938682637*coeff[5]*fSkin[5])-0.270632938682637*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]-0.270632938682637*fSkin[2]*coeff[4]-0.270632938682637*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = (-0.270632938682637*coeff[4]*fSkin[5])-0.270632938682637*coeff[4]*fEdge[5]-0.270632938682637*fSkin[2]*coeff[5]-0.270632938682637*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = (-0.4330127018922193*fSkin[5]*coeff[7])+0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]-0.4330127018922193*fSkin[2]*coeff[6]+0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]-0.7036456405748562*fSkin[1]*coeff[5]+0.2706329386826369*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]-0.7036456405748562*fSkin[0]*coeff[4]+0.2706329386826369*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = (-0.270632938682637*coeff[5]*fSkin[11])-0.270632938682637*coeff[5]*fEdge[11]-0.270632938682637*coeff[4]*fSkin[7]-0.270632938682637*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = (-0.270632938682637*coeff[5]*fSkin[12])-0.270632938682637*coeff[5]*fEdge[12]-0.270632938682637*coeff[4]*fSkin[9]-0.270632938682637*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = (-0.4330127018922193*fSkin[2]*coeff[7])+0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]-0.4330127018922193*fSkin[5]*coeff[6]+0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]-0.7036456405748562*fSkin[0]*coeff[5]+0.2706329386826369*fEdge[0]*coeff[5]-0.7036456405748562*fSkin[1]*coeff[4]+0.2706329386826369*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = (-0.270632938682637*coeff[4]*fSkin[11])-0.270632938682637*coeff[4]*fEdge[11]-0.270632938682637*coeff[5]*fSkin[7]-0.270632938682637*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = (-0.4330127018922193*coeff[7]*fSkin[11])-0.71875*coeff[5]*fSkin[11]+0.4330127018922193*coeff[7]*fEdge[11]-0.21875*coeff[5]*fEdge[11]-0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]+0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]-0.7036456405748562*coeff[5]*fSkin[6]+0.2706329386826369*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]-0.7036456405748562*fSkin[3]*coeff[4]+0.2706329386826369*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = (-0.270632938682637*coeff[4]*fSkin[12])-0.270632938682637*coeff[4]*fEdge[12]-0.270632938682637*coeff[5]*fSkin[9]-0.270632938682637*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = (-0.4330127018922193*coeff[7]*fSkin[12])-0.71875*coeff[5]*fSkin[12]+0.4330127018922193*coeff[7]*fEdge[12]-0.21875*coeff[5]*fEdge[12]-0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]+0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]-0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]+0.2706329386826369*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]-0.7036456405748562*coeff[4]*fSkin[4]+0.2706329386826369*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = (-0.270632938682637*coeff[5]*fSkin[15])-0.270632938682637*coeff[5]*fEdge[15]-0.270632938682637*coeff[4]*fSkin[14]-0.270632938682637*coeff[4]*fEdge[14]-0.28125*coeff[5]*fSkin[13]+0.28125*coeff[5]*fEdge[13]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = (-0.4330127018922193*coeff[6]*fSkin[11])-0.71875*coeff[4]*fSkin[11]+0.4330127018922193*coeff[6]*fEdge[11]-0.21875*coeff[4]*fEdge[11]-0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]+0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]-0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]+0.2706329386826369*coeff[4]*fEdge[6]-0.7036456405748562*fSkin[3]*coeff[5]+0.2706329386826369*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = (-0.4330127018922193*coeff[6]*fSkin[12])-0.71875*coeff[4]*fSkin[12]+0.4330127018922193*coeff[6]*fEdge[12]-0.21875*coeff[4]*fEdge[12]-0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]+0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]-0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]+0.2706329386826369*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]-0.7036456405748562*fSkin[4]*coeff[5]+0.2706329386826369*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = (-0.270632938682637*coeff[4]*fSkin[15])-0.270632938682637*coeff[4]*fEdge[15]-0.270632938682637*coeff[5]*fSkin[14]-0.270632938682637*coeff[5]*fEdge[14]-0.28125*coeff[4]*fSkin[13]+0.28125*coeff[4]*fEdge[13]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = (-0.4330127018922193*coeff[7]*fSkin[15])-0.71875*coeff[5]*fSkin[15]+0.4330127018922193*coeff[7]*fEdge[15]-0.21875*coeff[5]*fEdge[15]-0.4330127018922193*coeff[6]*fSkin[14]-0.71875*coeff[4]*fSkin[14]+0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]-0.375*coeff[7]*fSkin[13]-0.7036456405748562*coeff[5]*fSkin[13]-0.375*coeff[7]*fEdge[13]+0.2706329386826369*coeff[5]*fEdge[13]-0.375*coeff[6]*fSkin[10]-0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]+0.2706329386826369*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = (-0.4330127018922193*coeff[6]*fSkin[15])-0.71875*coeff[4]*fSkin[15]+0.4330127018922193*coeff[6]*fEdge[15]-0.21875*coeff[4]*fEdge[15]-0.4330127018922193*coeff[7]*fSkin[14]-0.71875*coeff[5]*fSkin[14]+0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]-0.375*coeff[6]*fSkin[13]-0.7036456405748562*coeff[4]*fSkin[13]-0.375*coeff[6]*fEdge[13]+0.2706329386826369*coeff[4]*fEdge[13]-0.375*coeff[7]*fSkin[10]-0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]+0.2706329386826369*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = (-0.270632938682637*coeff[5]*fSkin[20])-0.270632938682637*coeff[5]*fEdge[20]-0.2706329386826371*coeff[4]*fSkin[18]-0.2706329386826371*coeff[4]*fEdge[18]-0.2812499999999999*coeff[5]*fSkin[17]+0.2812499999999999*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[16]+0.28125*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = (-0.2706329386826371*coeff[4]*fSkin[20])-0.2706329386826371*coeff[4]*fEdge[20]-0.270632938682637*coeff[5]*fSkin[18]-0.270632938682637*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.2812499999999999*coeff[5]*fSkin[16]+0.2812499999999999*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = (-0.4330127018922194*coeff[7]*fSkin[20])-0.7187500000000001*coeff[5]*fSkin[20]+0.4330127018922194*coeff[7]*fEdge[20]-0.21875*coeff[5]*fEdge[20]-0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]+0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]-0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]+0.2706329386826369*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[16]-0.7036456405748563*coeff[4]*fSkin[16]-0.375*coeff[6]*fEdge[16]+0.2706329386826371*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = (-0.270632938682637*coeff[5]*fSkin[23])-0.270632938682637*coeff[5]*fEdge[23]-0.2706329386826371*coeff[4]*fSkin[22]-0.2706329386826371*coeff[4]*fEdge[22]-0.2812499999999999*coeff[5]*fSkin[21]+0.2812499999999999*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = (-0.4330127018922193*coeff[6]*fSkin[20])-0.71875*coeff[4]*fSkin[20]+0.4330127018922193*coeff[6]*fEdge[20]-0.21875*coeff[4]*fEdge[20]-0.4330127018922194*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]+0.4330127018922194*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]-0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]+0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[16]-0.7036456405748562*coeff[5]*fSkin[16]-0.375*coeff[7]*fEdge[16]+0.2706329386826369*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = (-0.2706329386826371*coeff[4]*fSkin[23])-0.2706329386826371*coeff[4]*fEdge[23]-0.270632938682637*coeff[5]*fSkin[22]-0.270632938682637*coeff[5]*fEdge[22]-0.28125*coeff[4]*fSkin[21]+0.28125*coeff[4]*fEdge[21]-0.2812499999999999*coeff[5]*fSkin[19]+0.2812499999999999*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = (-0.4330127018922194*coeff[7]*fSkin[23])-0.7187500000000001*coeff[5]*fSkin[23]+0.4330127018922194*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]-0.4330127018922193*coeff[6]*fSkin[22]-0.71875*coeff[4]*fSkin[22]+0.4330127018922193*coeff[6]*fEdge[22]-0.21875*coeff[4]*fEdge[22]-0.375*coeff[7]*fSkin[21]-0.7036456405748562*coeff[5]*fSkin[21]-0.375*coeff[7]*fEdge[21]+0.2706329386826369*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[19]-0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]+0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = (-0.4330127018922193*coeff[6]*fSkin[23])-0.71875*coeff[4]*fSkin[23]+0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]-0.4330127018922194*coeff[7]*fSkin[22]-0.7187500000000001*coeff[5]*fSkin[22]+0.4330127018922194*coeff[7]*fEdge[22]-0.21875*coeff[5]*fEdge[22]-0.375*coeff[6]*fSkin[21]-0.7036456405748563*coeff[4]*fSkin[21]-0.375*coeff[6]*fEdge[21]+0.2706329386826371*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[19]-0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]+0.2706329386826369*coeff[5]*fEdge[19]; 

  boundSurf_incr[2] = 0.8660254037844386*fSkin[5]*coeff[7]-0.75*fSkin[1]*coeff[7]+0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]+0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]+0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = 0.8660254037844386*fSkin[2]*coeff[7]-0.75*fSkin[0]*coeff[7]+0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]+0.4330127018922193*fSkin[0]*coeff[5]+0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = 0.8660254037844386*coeff[7]*fSkin[11]-0.5*coeff[5]*fSkin[11]+0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]+0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]+0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = 0.8660254037844386*coeff[7]*fSkin[12]-0.5*coeff[5]*fSkin[12]+0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]+0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]+0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[6]*fSkin[11]-0.5*coeff[4]*fSkin[11]+0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]+0.4330127018922193*coeff[4]*fSkin[6]+0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[6]*fSkin[12]-0.5*coeff[4]*fSkin[12]+0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]+0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]+0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[14] = 0.8660254037844386*coeff[7]*fSkin[15]-0.5*coeff[5]*fSkin[15]+0.8660254037844386*coeff[6]*fSkin[14]-0.5*coeff[4]*fSkin[14]-0.75*coeff[7]*fSkin[13]+0.4330127018922193*coeff[5]*fSkin[13]-0.75*coeff[6]*fSkin[10]+0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[6]*fSkin[15]-0.5*coeff[4]*fSkin[15]+0.8660254037844386*coeff[7]*fSkin[14]-0.5*coeff[5]*fSkin[14]-0.75*coeff[6]*fSkin[13]+0.4330127018922193*coeff[4]*fSkin[13]-0.75*coeff[7]*fSkin[10]+0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[18] = 0.8660254037844387*coeff[7]*fSkin[20]-0.5000000000000001*coeff[5]*fSkin[20]+0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]+0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[16]+0.4330127018922194*coeff[4]*fSkin[16]; 
  boundSurf_incr[20] = 0.8660254037844386*coeff[6]*fSkin[20]-0.5*coeff[4]*fSkin[20]+0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]+0.4330127018922194*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[16]+0.4330127018922193*coeff[5]*fSkin[16]; 
  boundSurf_incr[22] = 0.8660254037844387*coeff[7]*fSkin[23]-0.5000000000000001*coeff[5]*fSkin[23]+0.8660254037844386*coeff[6]*fSkin[22]-0.5*coeff[4]*fSkin[22]-0.75*coeff[7]*fSkin[21]+0.4330127018922193*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[19]+0.4330127018922194*coeff[4]*fSkin[19]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[6]*fSkin[23]-0.5*coeff[4]*fSkin[23]+0.8660254037844387*coeff[7]*fSkin[22]-0.5000000000000001*coeff[5]*fSkin[22]-0.75*coeff[6]*fSkin[21]+0.4330127018922194*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[19]+0.4330127018922193*coeff[5]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 0.270632938682637*coeff[5]*fSkin[5]+0.270632938682637*coeff[5]*fEdge[5]-0.28125*fSkin[1]*coeff[5]+0.28125*fEdge[1]*coeff[5]+0.270632938682637*fSkin[2]*coeff[4]+0.270632938682637*fEdge[2]*coeff[4]-0.28125*fSkin[0]*coeff[4]+0.28125*fEdge[0]*coeff[4]; 
  edgeSurf_incr[1] = 0.270632938682637*coeff[4]*fSkin[5]+0.270632938682637*coeff[4]*fEdge[5]+0.270632938682637*fSkin[2]*coeff[5]+0.270632938682637*fEdge[2]*coeff[5]-0.28125*fSkin[0]*coeff[5]+0.28125*fEdge[0]*coeff[5]-0.28125*fSkin[1]*coeff[4]+0.28125*fEdge[1]*coeff[4]; 
  edgeSurf_incr[2] = 0.4330127018922193*fSkin[5]*coeff[7]-0.4330127018922193*fEdge[5]*coeff[7]-0.375*fSkin[1]*coeff[7]-0.375*fEdge[1]*coeff[7]+0.4330127018922193*fSkin[2]*coeff[6]-0.4330127018922193*fEdge[2]*coeff[6]-0.375*fSkin[0]*coeff[6]-0.375*fEdge[0]*coeff[6]-0.71875*coeff[5]*fSkin[5]-0.21875*coeff[5]*fEdge[5]+0.7036456405748562*fSkin[1]*coeff[5]-0.2706329386826369*fEdge[1]*coeff[5]-0.71875*fSkin[2]*coeff[4]-0.21875*fEdge[2]*coeff[4]+0.7036456405748562*fSkin[0]*coeff[4]-0.2706329386826369*fEdge[0]*coeff[4]; 
  edgeSurf_incr[3] = 0.270632938682637*coeff[5]*fSkin[11]+0.270632938682637*coeff[5]*fEdge[11]+0.270632938682637*coeff[4]*fSkin[7]+0.270632938682637*coeff[4]*fEdge[7]-0.28125*coeff[5]*fSkin[6]+0.28125*coeff[5]*fEdge[6]-0.28125*fSkin[3]*coeff[4]+0.28125*fEdge[3]*coeff[4]; 
  edgeSurf_incr[4] = 0.270632938682637*coeff[5]*fSkin[12]+0.270632938682637*coeff[5]*fEdge[12]+0.270632938682637*coeff[4]*fSkin[9]+0.270632938682637*coeff[4]*fEdge[9]-0.28125*coeff[5]*fSkin[8]+0.28125*coeff[5]*fEdge[8]-0.28125*coeff[4]*fSkin[4]+0.28125*coeff[4]*fEdge[4]; 
  edgeSurf_incr[5] = 0.4330127018922193*fSkin[2]*coeff[7]-0.4330127018922193*fEdge[2]*coeff[7]-0.375*fSkin[0]*coeff[7]-0.375*fEdge[0]*coeff[7]+0.4330127018922193*fSkin[5]*coeff[6]-0.4330127018922193*fEdge[5]*coeff[6]-0.375*fSkin[1]*coeff[6]-0.375*fEdge[1]*coeff[6]-0.71875*coeff[4]*fSkin[5]-0.21875*coeff[4]*fEdge[5]-0.71875*fSkin[2]*coeff[5]-0.21875*fEdge[2]*coeff[5]+0.7036456405748562*fSkin[0]*coeff[5]-0.2706329386826369*fEdge[0]*coeff[5]+0.7036456405748562*fSkin[1]*coeff[4]-0.2706329386826369*fEdge[1]*coeff[4]; 
  edgeSurf_incr[6] = 0.270632938682637*coeff[4]*fSkin[11]+0.270632938682637*coeff[4]*fEdge[11]+0.270632938682637*coeff[5]*fSkin[7]+0.270632938682637*coeff[5]*fEdge[7]-0.28125*coeff[4]*fSkin[6]+0.28125*coeff[4]*fEdge[6]-0.28125*fSkin[3]*coeff[5]+0.28125*fEdge[3]*coeff[5]; 
  edgeSurf_incr[7] = 0.4330127018922193*coeff[7]*fSkin[11]-0.71875*coeff[5]*fSkin[11]-0.4330127018922193*coeff[7]*fEdge[11]-0.21875*coeff[5]*fEdge[11]+0.4330127018922193*coeff[6]*fSkin[7]-0.71875*coeff[4]*fSkin[7]-0.4330127018922193*coeff[6]*fEdge[7]-0.21875*coeff[4]*fEdge[7]-0.375*fSkin[6]*coeff[7]-0.375*fEdge[6]*coeff[7]+0.7036456405748562*coeff[5]*fSkin[6]-0.2706329386826369*coeff[5]*fEdge[6]-0.375*fSkin[3]*coeff[6]-0.375*fEdge[3]*coeff[6]+0.7036456405748562*fSkin[3]*coeff[4]-0.2706329386826369*fEdge[3]*coeff[4]; 
  edgeSurf_incr[8] = 0.270632938682637*coeff[4]*fSkin[12]+0.270632938682637*coeff[4]*fEdge[12]+0.270632938682637*coeff[5]*fSkin[9]+0.270632938682637*coeff[5]*fEdge[9]-0.28125*coeff[4]*fSkin[8]+0.28125*coeff[4]*fEdge[8]-0.28125*fSkin[4]*coeff[5]+0.28125*fEdge[4]*coeff[5]; 
  edgeSurf_incr[9] = 0.4330127018922193*coeff[7]*fSkin[12]-0.71875*coeff[5]*fSkin[12]-0.4330127018922193*coeff[7]*fEdge[12]-0.21875*coeff[5]*fEdge[12]+0.4330127018922193*coeff[6]*fSkin[9]-0.71875*coeff[4]*fSkin[9]-0.4330127018922193*coeff[6]*fEdge[9]-0.21875*coeff[4]*fEdge[9]-0.375*coeff[7]*fSkin[8]+0.7036456405748562*coeff[5]*fSkin[8]-0.375*coeff[7]*fEdge[8]-0.2706329386826369*coeff[5]*fEdge[8]-0.375*fSkin[4]*coeff[6]-0.375*fEdge[4]*coeff[6]+0.7036456405748562*coeff[4]*fSkin[4]-0.2706329386826369*coeff[4]*fEdge[4]; 
  edgeSurf_incr[10] = 0.270632938682637*coeff[5]*fSkin[15]+0.270632938682637*coeff[5]*fEdge[15]+0.270632938682637*coeff[4]*fSkin[14]+0.270632938682637*coeff[4]*fEdge[14]-0.28125*coeff[5]*fSkin[13]+0.28125*coeff[5]*fEdge[13]-0.28125*coeff[4]*fSkin[10]+0.28125*coeff[4]*fEdge[10]; 
  edgeSurf_incr[11] = 0.4330127018922193*coeff[6]*fSkin[11]-0.71875*coeff[4]*fSkin[11]-0.4330127018922193*coeff[6]*fEdge[11]-0.21875*coeff[4]*fEdge[11]+0.4330127018922193*coeff[7]*fSkin[7]-0.71875*coeff[5]*fSkin[7]-0.4330127018922193*coeff[7]*fEdge[7]-0.21875*coeff[5]*fEdge[7]-0.375*fSkin[3]*coeff[7]-0.375*fEdge[3]*coeff[7]-0.375*coeff[6]*fSkin[6]+0.7036456405748562*coeff[4]*fSkin[6]-0.375*coeff[6]*fEdge[6]-0.2706329386826369*coeff[4]*fEdge[6]+0.7036456405748562*fSkin[3]*coeff[5]-0.2706329386826369*fEdge[3]*coeff[5]; 
  edgeSurf_incr[12] = 0.4330127018922193*coeff[6]*fSkin[12]-0.71875*coeff[4]*fSkin[12]-0.4330127018922193*coeff[6]*fEdge[12]-0.21875*coeff[4]*fEdge[12]+0.4330127018922193*coeff[7]*fSkin[9]-0.71875*coeff[5]*fSkin[9]-0.4330127018922193*coeff[7]*fEdge[9]-0.21875*coeff[5]*fEdge[9]-0.375*coeff[6]*fSkin[8]+0.7036456405748562*coeff[4]*fSkin[8]-0.375*coeff[6]*fEdge[8]-0.2706329386826369*coeff[4]*fEdge[8]-0.375*fSkin[4]*coeff[7]-0.375*fEdge[4]*coeff[7]+0.7036456405748562*fSkin[4]*coeff[5]-0.2706329386826369*fEdge[4]*coeff[5]; 
  edgeSurf_incr[13] = 0.270632938682637*coeff[4]*fSkin[15]+0.270632938682637*coeff[4]*fEdge[15]+0.270632938682637*coeff[5]*fSkin[14]+0.270632938682637*coeff[5]*fEdge[14]-0.28125*coeff[4]*fSkin[13]+0.28125*coeff[4]*fEdge[13]-0.28125*coeff[5]*fSkin[10]+0.28125*coeff[5]*fEdge[10]; 
  edgeSurf_incr[14] = 0.4330127018922193*coeff[7]*fSkin[15]-0.71875*coeff[5]*fSkin[15]-0.4330127018922193*coeff[7]*fEdge[15]-0.21875*coeff[5]*fEdge[15]+0.4330127018922193*coeff[6]*fSkin[14]-0.71875*coeff[4]*fSkin[14]-0.4330127018922193*coeff[6]*fEdge[14]-0.21875*coeff[4]*fEdge[14]-0.375*coeff[7]*fSkin[13]+0.7036456405748562*coeff[5]*fSkin[13]-0.375*coeff[7]*fEdge[13]-0.2706329386826369*coeff[5]*fEdge[13]-0.375*coeff[6]*fSkin[10]+0.7036456405748562*coeff[4]*fSkin[10]-0.375*coeff[6]*fEdge[10]-0.2706329386826369*coeff[4]*fEdge[10]; 
  edgeSurf_incr[15] = 0.4330127018922193*coeff[6]*fSkin[15]-0.71875*coeff[4]*fSkin[15]-0.4330127018922193*coeff[6]*fEdge[15]-0.21875*coeff[4]*fEdge[15]+0.4330127018922193*coeff[7]*fSkin[14]-0.71875*coeff[5]*fSkin[14]-0.4330127018922193*coeff[7]*fEdge[14]-0.21875*coeff[5]*fEdge[14]-0.375*coeff[6]*fSkin[13]+0.7036456405748562*coeff[4]*fSkin[13]-0.375*coeff[6]*fEdge[13]-0.2706329386826369*coeff[4]*fEdge[13]-0.375*coeff[7]*fSkin[10]+0.7036456405748562*coeff[5]*fSkin[10]-0.375*coeff[7]*fEdge[10]-0.2706329386826369*coeff[5]*fEdge[10]; 
  edgeSurf_incr[16] = 0.270632938682637*coeff[5]*fSkin[20]+0.270632938682637*coeff[5]*fEdge[20]+0.2706329386826371*coeff[4]*fSkin[18]+0.2706329386826371*coeff[4]*fEdge[18]-0.2812499999999999*coeff[5]*fSkin[17]+0.2812499999999999*coeff[5]*fEdge[17]-0.28125*coeff[4]*fSkin[16]+0.28125*coeff[4]*fEdge[16]; 
  edgeSurf_incr[17] = 0.2706329386826371*coeff[4]*fSkin[20]+0.2706329386826371*coeff[4]*fEdge[20]+0.270632938682637*coeff[5]*fSkin[18]+0.270632938682637*coeff[5]*fEdge[18]-0.28125*coeff[4]*fSkin[17]+0.28125*coeff[4]*fEdge[17]-0.2812499999999999*coeff[5]*fSkin[16]+0.2812499999999999*coeff[5]*fEdge[16]; 
  edgeSurf_incr[18] = 0.4330127018922194*coeff[7]*fSkin[20]-0.7187500000000001*coeff[5]*fSkin[20]-0.4330127018922194*coeff[7]*fEdge[20]-0.21875*coeff[5]*fEdge[20]+0.4330127018922193*coeff[6]*fSkin[18]-0.71875*coeff[4]*fSkin[18]-0.4330127018922193*coeff[6]*fEdge[18]-0.21875*coeff[4]*fEdge[18]-0.375*coeff[7]*fSkin[17]+0.7036456405748562*coeff[5]*fSkin[17]-0.375*coeff[7]*fEdge[17]-0.2706329386826369*coeff[5]*fEdge[17]-0.375*coeff[6]*fSkin[16]+0.7036456405748563*coeff[4]*fSkin[16]-0.375*coeff[6]*fEdge[16]-0.2706329386826371*coeff[4]*fEdge[16]; 
  edgeSurf_incr[19] = 0.270632938682637*coeff[5]*fSkin[23]+0.270632938682637*coeff[5]*fEdge[23]+0.2706329386826371*coeff[4]*fSkin[22]+0.2706329386826371*coeff[4]*fEdge[22]-0.2812499999999999*coeff[5]*fSkin[21]+0.2812499999999999*coeff[5]*fEdge[21]-0.28125*coeff[4]*fSkin[19]+0.28125*coeff[4]*fEdge[19]; 
  edgeSurf_incr[20] = 0.4330127018922193*coeff[6]*fSkin[20]-0.71875*coeff[4]*fSkin[20]-0.4330127018922193*coeff[6]*fEdge[20]-0.21875*coeff[4]*fEdge[20]+0.4330127018922194*coeff[7]*fSkin[18]-0.7187500000000001*coeff[5]*fSkin[18]-0.4330127018922194*coeff[7]*fEdge[18]-0.21875*coeff[5]*fEdge[18]-0.375*coeff[6]*fSkin[17]+0.7036456405748563*coeff[4]*fSkin[17]-0.375*coeff[6]*fEdge[17]-0.2706329386826371*coeff[4]*fEdge[17]-0.375*coeff[7]*fSkin[16]+0.7036456405748562*coeff[5]*fSkin[16]-0.375*coeff[7]*fEdge[16]-0.2706329386826369*coeff[5]*fEdge[16]; 
  edgeSurf_incr[21] = 0.2706329386826371*coeff[4]*fSkin[23]+0.2706329386826371*coeff[4]*fEdge[23]+0.270632938682637*coeff[5]*fSkin[22]+0.270632938682637*coeff[5]*fEdge[22]-0.28125*coeff[4]*fSkin[21]+0.28125*coeff[4]*fEdge[21]-0.2812499999999999*coeff[5]*fSkin[19]+0.2812499999999999*coeff[5]*fEdge[19]; 
  edgeSurf_incr[22] = 0.4330127018922194*coeff[7]*fSkin[23]-0.7187500000000001*coeff[5]*fSkin[23]-0.4330127018922194*coeff[7]*fEdge[23]-0.21875*coeff[5]*fEdge[23]+0.4330127018922193*coeff[6]*fSkin[22]-0.71875*coeff[4]*fSkin[22]-0.4330127018922193*coeff[6]*fEdge[22]-0.21875*coeff[4]*fEdge[22]-0.375*coeff[7]*fSkin[21]+0.7036456405748562*coeff[5]*fSkin[21]-0.375*coeff[7]*fEdge[21]-0.2706329386826369*coeff[5]*fEdge[21]-0.375*coeff[6]*fSkin[19]+0.7036456405748563*coeff[4]*fSkin[19]-0.375*coeff[6]*fEdge[19]-0.2706329386826371*coeff[4]*fEdge[19]; 
  edgeSurf_incr[23] = 0.4330127018922193*coeff[6]*fSkin[23]-0.71875*coeff[4]*fSkin[23]-0.4330127018922193*coeff[6]*fEdge[23]-0.21875*coeff[4]*fEdge[23]+0.4330127018922194*coeff[7]*fSkin[22]-0.7187500000000001*coeff[5]*fSkin[22]-0.4330127018922194*coeff[7]*fEdge[22]-0.21875*coeff[5]*fEdge[22]-0.375*coeff[6]*fSkin[21]+0.7036456405748563*coeff[4]*fSkin[21]-0.375*coeff[6]*fEdge[21]-0.2706329386826371*coeff[4]*fEdge[21]-0.375*coeff[7]*fSkin[19]+0.7036456405748562*coeff[5]*fSkin[19]-0.375*coeff[7]*fEdge[19]-0.2706329386826369*coeff[5]*fEdge[19]; 

  boundSurf_incr[2] = (-0.8660254037844386*fSkin[5]*coeff[7])-0.75*fSkin[1]*coeff[7]-0.8660254037844386*fSkin[2]*coeff[6]-0.75*fSkin[0]*coeff[6]-0.5*coeff[5]*fSkin[5]-0.4330127018922193*fSkin[1]*coeff[5]-0.5*fSkin[2]*coeff[4]-0.4330127018922193*fSkin[0]*coeff[4]; 
  boundSurf_incr[5] = (-0.8660254037844386*fSkin[2]*coeff[7])-0.75*fSkin[0]*coeff[7]-0.8660254037844386*fSkin[5]*coeff[6]-0.75*fSkin[1]*coeff[6]-0.5*coeff[4]*fSkin[5]-0.5*fSkin[2]*coeff[5]-0.4330127018922193*fSkin[0]*coeff[5]-0.4330127018922193*fSkin[1]*coeff[4]; 
  boundSurf_incr[7] = (-0.8660254037844386*coeff[7]*fSkin[11])-0.5*coeff[5]*fSkin[11]-0.8660254037844386*coeff[6]*fSkin[7]-0.5*coeff[4]*fSkin[7]-0.75*fSkin[6]*coeff[7]-0.4330127018922193*coeff[5]*fSkin[6]-0.75*fSkin[3]*coeff[6]-0.4330127018922193*fSkin[3]*coeff[4]; 
  boundSurf_incr[9] = (-0.8660254037844386*coeff[7]*fSkin[12])-0.5*coeff[5]*fSkin[12]-0.8660254037844386*coeff[6]*fSkin[9]-0.5*coeff[4]*fSkin[9]-0.75*coeff[7]*fSkin[8]-0.4330127018922193*coeff[5]*fSkin[8]-0.75*fSkin[4]*coeff[6]-0.4330127018922193*coeff[4]*fSkin[4]; 
  boundSurf_incr[11] = (-0.8660254037844386*coeff[6]*fSkin[11])-0.5*coeff[4]*fSkin[11]-0.8660254037844386*coeff[7]*fSkin[7]-0.5*coeff[5]*fSkin[7]-0.75*fSkin[3]*coeff[7]-0.75*coeff[6]*fSkin[6]-0.4330127018922193*coeff[4]*fSkin[6]-0.4330127018922193*fSkin[3]*coeff[5]; 
  boundSurf_incr[12] = (-0.8660254037844386*coeff[6]*fSkin[12])-0.5*coeff[4]*fSkin[12]-0.8660254037844386*coeff[7]*fSkin[9]-0.5*coeff[5]*fSkin[9]-0.75*coeff[6]*fSkin[8]-0.4330127018922193*coeff[4]*fSkin[8]-0.75*fSkin[4]*coeff[7]-0.4330127018922193*fSkin[4]*coeff[5]; 
  boundSurf_incr[14] = (-0.8660254037844386*coeff[7]*fSkin[15])-0.5*coeff[5]*fSkin[15]-0.8660254037844386*coeff[6]*fSkin[14]-0.5*coeff[4]*fSkin[14]-0.75*coeff[7]*fSkin[13]-0.4330127018922193*coeff[5]*fSkin[13]-0.75*coeff[6]*fSkin[10]-0.4330127018922193*coeff[4]*fSkin[10]; 
  boundSurf_incr[15] = (-0.8660254037844386*coeff[6]*fSkin[15])-0.5*coeff[4]*fSkin[15]-0.8660254037844386*coeff[7]*fSkin[14]-0.5*coeff[5]*fSkin[14]-0.75*coeff[6]*fSkin[13]-0.4330127018922193*coeff[4]*fSkin[13]-0.75*coeff[7]*fSkin[10]-0.4330127018922193*coeff[5]*fSkin[10]; 
  boundSurf_incr[18] = (-0.8660254037844387*coeff[7]*fSkin[20])-0.5000000000000001*coeff[5]*fSkin[20]-0.8660254037844386*coeff[6]*fSkin[18]-0.5*coeff[4]*fSkin[18]-0.75*coeff[7]*fSkin[17]-0.4330127018922193*coeff[5]*fSkin[17]-0.75*coeff[6]*fSkin[16]-0.4330127018922194*coeff[4]*fSkin[16]; 
  boundSurf_incr[20] = (-0.8660254037844386*coeff[6]*fSkin[20])-0.5*coeff[4]*fSkin[20]-0.8660254037844387*coeff[7]*fSkin[18]-0.5000000000000001*coeff[5]*fSkin[18]-0.75*coeff[6]*fSkin[17]-0.4330127018922194*coeff[4]*fSkin[17]-0.75*coeff[7]*fSkin[16]-0.4330127018922193*coeff[5]*fSkin[16]; 
  boundSurf_incr[22] = (-0.8660254037844387*coeff[7]*fSkin[23])-0.5000000000000001*coeff[5]*fSkin[23]-0.8660254037844386*coeff[6]*fSkin[22]-0.5*coeff[4]*fSkin[22]-0.75*coeff[7]*fSkin[21]-0.4330127018922193*coeff[5]*fSkin[21]-0.75*coeff[6]*fSkin[19]-0.4330127018922194*coeff[4]*fSkin[19]; 
  boundSurf_incr[23] = (-0.8660254037844386*coeff[6]*fSkin[23])-0.5*coeff[4]*fSkin[23]-0.8660254037844387*coeff[7]*fSkin[22]-0.5000000000000001*coeff[5]*fSkin[22]-0.75*coeff[6]*fSkin[21]-0.4330127018922194*coeff[4]*fSkin[21]-0.75*coeff[7]*fSkin[19]-0.4330127018922193*coeff[5]*fSkin[19]; 

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

  }

  return 0.;
}

