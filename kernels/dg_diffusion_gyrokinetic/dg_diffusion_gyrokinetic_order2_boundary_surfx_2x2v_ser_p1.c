#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[24] = {0.0}; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.5412658773652741*coeff[0]*fSkin[1])-0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]-1.407291281149712*coeff[0]*fSkin[0]+0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.5412658773652741*coeff[0]*fSkin[5])-0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-0.5412658773652741*coeff[0]*fSkin[6])-0.5412658773652741*coeff[0]*fEdge[6]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.5412658773652741*coeff[0]*fSkin[8])-0.5412658773652741*coeff[0]*fEdge[8]-0.5625*coeff[0]*fSkin[4]+0.5625*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]-1.407291281149712*coeff[0]*fSkin[2]+0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-1.4375*coeff[0]*fSkin[6])-0.4375*coeff[0]*fEdge[6]-1.407291281149712*coeff[0]*fSkin[3]+0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-0.5412658773652741*coeff[0]*fSkin[11])-0.5412658773652741*coeff[0]*fEdge[11]-0.5625*coeff[0]*fSkin[7]+0.5625*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-1.4375*coeff[0]*fSkin[8])-0.4375*coeff[0]*fEdge[8]-1.407291281149712*coeff[0]*fSkin[4]+0.5412658773652739*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-0.5412658773652741*coeff[0]*fSkin[12])-0.5412658773652741*coeff[0]*fEdge[12]-0.5625*coeff[0]*fSkin[9]+0.5625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-0.5412658773652741*coeff[0]*fSkin[13])-0.5412658773652741*coeff[0]*fEdge[13]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[0]*fSkin[11])-0.4375*coeff[0]*fEdge[11]-1.407291281149712*coeff[0]*fSkin[7]+0.5412658773652739*coeff[0]*fEdge[7]; 
  edgeSurf_incr[12] = (-1.4375*coeff[0]*fSkin[12])-0.4375*coeff[0]*fEdge[12]-1.407291281149712*coeff[0]*fSkin[9]+0.5412658773652739*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = (-1.4375*coeff[0]*fSkin[13])-0.4375*coeff[0]*fEdge[13]-1.407291281149712*coeff[0]*fSkin[10]+0.5412658773652739*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = (-0.5412658773652741*coeff[0]*fSkin[15])-0.5412658773652741*coeff[0]*fEdge[15]-0.5625*coeff[0]*fSkin[14]+0.5625*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-1.4375*coeff[0]*fSkin[15])-0.4375*coeff[0]*fEdge[15]-1.407291281149712*coeff[0]*fSkin[14]+0.5412658773652739*coeff[0]*fEdge[14]; 
  edgeSurf_incr[16] = (-0.5412658773652742*coeff[0]*fSkin[17])-0.5412658773652742*coeff[0]*fEdge[17]-0.5625*coeff[0]*fSkin[16]+0.5625*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-1.4375*coeff[0]*fSkin[17])-0.4375*coeff[0]*fEdge[17]-1.407291281149713*coeff[0]*fSkin[16]+0.5412658773652742*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = (-0.5412658773652742*coeff[0]*fSkin[20])-0.5412658773652742*coeff[0]*fEdge[20]-0.5625*coeff[0]*fSkin[18]+0.5625*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-0.5412658773652742*coeff[0]*fSkin[21])-0.5412658773652742*coeff[0]*fEdge[21]-0.5625*coeff[0]*fSkin[19]+0.5625*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = (-1.4375*coeff[0]*fSkin[20])-0.4375*coeff[0]*fEdge[20]-1.407291281149713*coeff[0]*fSkin[18]+0.5412658773652742*coeff[0]*fEdge[18]; 
  edgeSurf_incr[21] = (-1.4375*coeff[0]*fSkin[21])-0.4375*coeff[0]*fEdge[21]-1.407291281149713*coeff[0]*fSkin[19]+0.5412658773652742*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = (-0.5412658773652742*coeff[0]*fSkin[23])-0.5412658773652742*coeff[0]*fEdge[23]-0.5625*coeff[0]*fSkin[22]+0.5625*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-1.4375*coeff[0]*fSkin[23])-0.4375*coeff[0]*fEdge[23]-1.407291281149713*coeff[0]*fSkin[22]+0.5412658773652742*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[0]*fSkin[0]-1.0*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[0]*fSkin[2]-1.0*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[0]*fSkin[3]-1.0*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 0.8660254037844386*coeff[0]*fSkin[4]-1.0*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[0]*fSkin[7]-1.0*coeff[0]*fSkin[11]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[0]*fSkin[9]-1.0*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 0.8660254037844386*coeff[0]*fSkin[10]-1.0*coeff[0]*fSkin[13]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[0]*fSkin[14]-1.0*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 0.8660254037844387*coeff[0]*fSkin[16]-1.0*coeff[0]*fSkin[17]; 
  boundSurf_incr[20] = 0.8660254037844387*coeff[0]*fSkin[18]-1.0*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 0.8660254037844387*coeff[0]*fSkin[19]-1.0*coeff[0]*fSkin[21]; 
  boundSurf_incr[23] = 0.8660254037844387*coeff[0]*fSkin[22]-1.0*coeff[0]*fSkin[23]; 

  } else { 

  edgeSurf_incr[0] = 0.5412658773652741*coeff[0]*fSkin[1]+0.5412658773652741*coeff[0]*fEdge[1]-0.5625*coeff[0]*fSkin[0]+0.5625*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-1.4375*coeff[0]*fSkin[1])-0.4375*coeff[0]*fEdge[1]+1.407291281149712*coeff[0]*fSkin[0]-0.5412658773652739*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.5412658773652741*coeff[0]*fSkin[5]+0.5412658773652741*coeff[0]*fEdge[5]-0.5625*coeff[0]*fSkin[2]+0.5625*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 0.5412658773652741*coeff[0]*fSkin[6]+0.5412658773652741*coeff[0]*fEdge[6]-0.5625*coeff[0]*fSkin[3]+0.5625*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.5412658773652741*coeff[0]*fSkin[8]+0.5412658773652741*coeff[0]*fEdge[8]-0.5625*coeff[0]*fSkin[4]+0.5625*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-1.4375*coeff[0]*fSkin[5])-0.4375*coeff[0]*fEdge[5]+1.407291281149712*coeff[0]*fSkin[2]-0.5412658773652739*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = (-1.4375*coeff[0]*fSkin[6])-0.4375*coeff[0]*fEdge[6]+1.407291281149712*coeff[0]*fSkin[3]-0.5412658773652739*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 0.5412658773652741*coeff[0]*fSkin[11]+0.5412658773652741*coeff[0]*fEdge[11]-0.5625*coeff[0]*fSkin[7]+0.5625*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = (-1.4375*coeff[0]*fSkin[8])-0.4375*coeff[0]*fEdge[8]+1.407291281149712*coeff[0]*fSkin[4]-0.5412658773652739*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 0.5412658773652741*coeff[0]*fSkin[12]+0.5412658773652741*coeff[0]*fEdge[12]-0.5625*coeff[0]*fSkin[9]+0.5625*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 0.5412658773652741*coeff[0]*fSkin[13]+0.5412658773652741*coeff[0]*fEdge[13]-0.5625*coeff[0]*fSkin[10]+0.5625*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-1.4375*coeff[0]*fSkin[11])-0.4375*coeff[0]*fEdge[11]+1.407291281149712*coeff[0]*fSkin[7]-0.5412658773652739*coeff[0]*fEdge[7]; 
  edgeSurf_incr[12] = (-1.4375*coeff[0]*fSkin[12])-0.4375*coeff[0]*fEdge[12]+1.407291281149712*coeff[0]*fSkin[9]-0.5412658773652739*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = (-1.4375*coeff[0]*fSkin[13])-0.4375*coeff[0]*fEdge[13]+1.407291281149712*coeff[0]*fSkin[10]-0.5412658773652739*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = 0.5412658773652741*coeff[0]*fSkin[15]+0.5412658773652741*coeff[0]*fEdge[15]-0.5625*coeff[0]*fSkin[14]+0.5625*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = (-1.4375*coeff[0]*fSkin[15])-0.4375*coeff[0]*fEdge[15]+1.407291281149712*coeff[0]*fSkin[14]-0.5412658773652739*coeff[0]*fEdge[14]; 
  edgeSurf_incr[16] = 0.5412658773652742*coeff[0]*fSkin[17]+0.5412658773652742*coeff[0]*fEdge[17]-0.5625*coeff[0]*fSkin[16]+0.5625*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-1.4375*coeff[0]*fSkin[17])-0.4375*coeff[0]*fEdge[17]+1.407291281149713*coeff[0]*fSkin[16]-0.5412658773652742*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = 0.5412658773652742*coeff[0]*fSkin[20]+0.5412658773652742*coeff[0]*fEdge[20]-0.5625*coeff[0]*fSkin[18]+0.5625*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 0.5412658773652742*coeff[0]*fSkin[21]+0.5412658773652742*coeff[0]*fEdge[21]-0.5625*coeff[0]*fSkin[19]+0.5625*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = (-1.4375*coeff[0]*fSkin[20])-0.4375*coeff[0]*fEdge[20]+1.407291281149713*coeff[0]*fSkin[18]-0.5412658773652742*coeff[0]*fEdge[18]; 
  edgeSurf_incr[21] = (-1.4375*coeff[0]*fSkin[21])-0.4375*coeff[0]*fEdge[21]+1.407291281149713*coeff[0]*fSkin[19]-0.5412658773652742*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = 0.5412658773652742*coeff[0]*fSkin[23]+0.5412658773652742*coeff[0]*fEdge[23]-0.5625*coeff[0]*fSkin[22]+0.5625*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = (-1.4375*coeff[0]*fSkin[23])-0.4375*coeff[0]*fEdge[23]+1.407291281149713*coeff[0]*fSkin[22]-0.5412658773652742*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = (-1.0*coeff[0]*fSkin[1])-0.8660254037844386*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = (-1.0*coeff[0]*fSkin[5])-0.8660254037844386*coeff[0]*fSkin[2]; 
  boundSurf_incr[6] = (-1.0*coeff[0]*fSkin[6])-0.8660254037844386*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = (-1.0*coeff[0]*fSkin[8])-0.8660254037844386*coeff[0]*fSkin[4]; 
  boundSurf_incr[11] = (-1.0*coeff[0]*fSkin[11])-0.8660254037844386*coeff[0]*fSkin[7]; 
  boundSurf_incr[12] = (-1.0*coeff[0]*fSkin[12])-0.8660254037844386*coeff[0]*fSkin[9]; 
  boundSurf_incr[13] = (-1.0*coeff[0]*fSkin[13])-0.8660254037844386*coeff[0]*fSkin[10]; 
  boundSurf_incr[15] = (-1.0*coeff[0]*fSkin[15])-0.8660254037844386*coeff[0]*fSkin[14]; 
  boundSurf_incr[17] = (-1.0*coeff[0]*fSkin[17])-0.8660254037844387*coeff[0]*fSkin[16]; 
  boundSurf_incr[20] = (-1.0*coeff[0]*fSkin[20])-0.8660254037844387*coeff[0]*fSkin[18]; 
  boundSurf_incr[21] = (-1.0*coeff[0]*fSkin[21])-0.8660254037844387*coeff[0]*fSkin[19]; 
  boundSurf_incr[23] = (-1.0*coeff[0]*fSkin[23])-0.8660254037844387*coeff[0]*fSkin[22]; 

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

GKYL_CU_DH double dg_diffusion_gyrokinetic_order2_boundary_surfx_2x2v_ser_p1_varcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],2.);

  double vol_incr[24] = {0.0}; 
  vol_incr[1] = 1.5*fSkin[2]*coeff[3]+1.5*fSkin[0]*coeff[1]; 
  vol_incr[5] = 1.5*fSkin[0]*coeff[3]+1.5*coeff[1]*fSkin[2]; 
  vol_incr[6] = 1.5*coeff[3]*fSkin[7]+1.5*coeff[1]*fSkin[3]; 
  vol_incr[8] = 1.5*coeff[3]*fSkin[9]+1.5*coeff[1]*fSkin[4]; 
  vol_incr[11] = 1.5*coeff[1]*fSkin[7]+1.5*coeff[3]*fSkin[3]; 
  vol_incr[12] = 1.5*coeff[1]*fSkin[9]+1.5*coeff[3]*fSkin[4]; 
  vol_incr[13] = 1.5*coeff[3]*fSkin[14]+1.5*coeff[1]*fSkin[10]; 
  vol_incr[15] = 1.5*coeff[1]*fSkin[14]+1.5*coeff[3]*fSkin[10]; 
  vol_incr[17] = 1.5*coeff[3]*fSkin[18]+1.5*coeff[1]*fSkin[16]; 
  vol_incr[20] = 1.5*coeff[1]*fSkin[18]+1.5*coeff[3]*fSkin[16]; 
  vol_incr[21] = 1.5*coeff[3]*fSkin[22]+1.5*coeff[1]*fSkin[19]; 
  vol_incr[23] = 1.5*coeff[1]*fSkin[22]+1.5*coeff[3]*fSkin[19]; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = (-0.270632938682637*coeff[2]*fSkin[5])-0.270632938682637*coeff[2]*fEdge[5]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]-0.270632938682637*coeff[0]*fSkin[1]-0.270632938682637*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = (-0.4330127018922193*coeff[3]*fSkin[5])-0.71875*coeff[2]*fSkin[5]+0.4330127018922193*coeff[3]*fEdge[5]-0.21875*coeff[2]*fEdge[5]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]-0.7036456405748562*coeff[2]*fSkin[2]+0.2706329386826369*coeff[2]*fEdge[2]-0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]+0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]-0.7036456405748562*coeff[0]*fSkin[0]+0.2706329386826369*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-0.270632938682637*coeff[0]*fSkin[5])-0.270632938682637*coeff[0]*fEdge[5]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]-0.270632938682637*fSkin[1]*coeff[2]-0.270632938682637*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = (-0.270632938682637*coeff[2]*fSkin[11])-0.270632938682637*coeff[2]*fEdge[11]-0.28125*coeff[2]*fSkin[7]+0.28125*coeff[2]*fEdge[7]-0.270632938682637*coeff[0]*fSkin[6]-0.270632938682637*coeff[0]*fEdge[6]-0.28125*coeff[0]*fSkin[3]+0.28125*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-0.270632938682637*coeff[2]*fSkin[12])-0.270632938682637*coeff[2]*fEdge[12]-0.28125*coeff[2]*fSkin[9]+0.28125*coeff[2]*fEdge[9]-0.270632938682637*coeff[0]*fSkin[8]-0.270632938682637*coeff[0]*fEdge[8]-0.28125*coeff[0]*fSkin[4]+0.28125*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = (-0.4330127018922193*coeff[1]*fSkin[5])-0.71875*coeff[0]*fSkin[5]+0.4330127018922193*coeff[1]*fEdge[5]-0.21875*coeff[0]*fEdge[5]-0.4330127018922193*fSkin[1]*coeff[3]+0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]-0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]+0.2706329386826369*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]-0.7036456405748562*fSkin[0]*coeff[2]+0.2706329386826369*fEdge[0]*coeff[2]; 
  edgeSurf_incr[6] = (-0.4330127018922193*coeff[3]*fSkin[11])-0.71875*coeff[2]*fSkin[11]+0.4330127018922193*coeff[3]*fEdge[11]-0.21875*coeff[2]*fEdge[11]-0.375*coeff[3]*fSkin[7]-0.7036456405748562*coeff[2]*fSkin[7]-0.375*coeff[3]*fEdge[7]+0.2706329386826369*coeff[2]*fEdge[7]-0.4330127018922193*coeff[1]*fSkin[6]-0.71875*coeff[0]*fSkin[6]+0.4330127018922193*coeff[1]*fEdge[6]-0.21875*coeff[0]*fEdge[6]-0.375*coeff[1]*fSkin[3]-0.7036456405748562*coeff[0]*fSkin[3]-0.375*coeff[1]*fEdge[3]+0.2706329386826369*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-0.270632938682637*coeff[0]*fSkin[11])-0.270632938682637*coeff[0]*fEdge[11]-0.28125*coeff[0]*fSkin[7]+0.28125*coeff[0]*fEdge[7]-0.270632938682637*coeff[2]*fSkin[6]-0.270632938682637*coeff[2]*fEdge[6]-0.28125*coeff[2]*fSkin[3]+0.28125*coeff[2]*fEdge[3]; 
  edgeSurf_incr[8] = (-0.4330127018922193*coeff[3]*fSkin[12])-0.71875*coeff[2]*fSkin[12]+0.4330127018922193*coeff[3]*fEdge[12]-0.21875*coeff[2]*fEdge[12]-0.375*coeff[3]*fSkin[9]-0.7036456405748562*coeff[2]*fSkin[9]-0.375*coeff[3]*fEdge[9]+0.2706329386826369*coeff[2]*fEdge[9]-0.4330127018922193*coeff[1]*fSkin[8]-0.71875*coeff[0]*fSkin[8]+0.4330127018922193*coeff[1]*fEdge[8]-0.21875*coeff[0]*fEdge[8]-0.375*coeff[1]*fSkin[4]-0.7036456405748562*coeff[0]*fSkin[4]-0.375*coeff[1]*fEdge[4]+0.2706329386826369*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-0.270632938682637*coeff[0]*fSkin[12])-0.270632938682637*coeff[0]*fEdge[12]-0.28125*coeff[0]*fSkin[9]+0.28125*coeff[0]*fEdge[9]-0.270632938682637*coeff[2]*fSkin[8]-0.270632938682637*coeff[2]*fEdge[8]-0.28125*coeff[2]*fSkin[4]+0.28125*coeff[2]*fEdge[4]; 
  edgeSurf_incr[10] = (-0.270632938682637*coeff[2]*fSkin[15])-0.270632938682637*coeff[2]*fEdge[15]-0.28125*coeff[2]*fSkin[14]+0.28125*coeff[2]*fEdge[14]-0.270632938682637*coeff[0]*fSkin[13]-0.270632938682637*coeff[0]*fEdge[13]-0.28125*coeff[0]*fSkin[10]+0.28125*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = (-0.4330127018922193*coeff[1]*fSkin[11])-0.71875*coeff[0]*fSkin[11]+0.4330127018922193*coeff[1]*fEdge[11]-0.21875*coeff[0]*fEdge[11]-0.375*coeff[1]*fSkin[7]-0.7036456405748562*coeff[0]*fSkin[7]-0.375*coeff[1]*fEdge[7]+0.2706329386826369*coeff[0]*fEdge[7]-0.4330127018922193*coeff[3]*fSkin[6]-0.71875*coeff[2]*fSkin[6]+0.4330127018922193*coeff[3]*fEdge[6]-0.21875*coeff[2]*fEdge[6]-0.375*coeff[3]*fSkin[3]-0.7036456405748562*coeff[2]*fSkin[3]-0.375*coeff[3]*fEdge[3]+0.2706329386826369*coeff[2]*fEdge[3]; 
  edgeSurf_incr[12] = (-0.4330127018922193*coeff[1]*fSkin[12])-0.71875*coeff[0]*fSkin[12]+0.4330127018922193*coeff[1]*fEdge[12]-0.21875*coeff[0]*fEdge[12]-0.375*coeff[1]*fSkin[9]-0.7036456405748562*coeff[0]*fSkin[9]-0.375*coeff[1]*fEdge[9]+0.2706329386826369*coeff[0]*fEdge[9]-0.4330127018922193*coeff[3]*fSkin[8]-0.71875*coeff[2]*fSkin[8]+0.4330127018922193*coeff[3]*fEdge[8]-0.21875*coeff[2]*fEdge[8]-0.375*coeff[3]*fSkin[4]-0.7036456405748562*coeff[2]*fSkin[4]-0.375*coeff[3]*fEdge[4]+0.2706329386826369*coeff[2]*fEdge[4]; 
  edgeSurf_incr[13] = (-0.4330127018922193*coeff[3]*fSkin[15])-0.71875*coeff[2]*fSkin[15]+0.4330127018922193*coeff[3]*fEdge[15]-0.21875*coeff[2]*fEdge[15]-0.375*coeff[3]*fSkin[14]-0.7036456405748562*coeff[2]*fSkin[14]-0.375*coeff[3]*fEdge[14]+0.2706329386826369*coeff[2]*fEdge[14]-0.4330127018922193*coeff[1]*fSkin[13]-0.71875*coeff[0]*fSkin[13]+0.4330127018922193*coeff[1]*fEdge[13]-0.21875*coeff[0]*fEdge[13]-0.375*coeff[1]*fSkin[10]-0.7036456405748562*coeff[0]*fSkin[10]-0.375*coeff[1]*fEdge[10]+0.2706329386826369*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = (-0.270632938682637*coeff[0]*fSkin[15])-0.270632938682637*coeff[0]*fEdge[15]-0.28125*coeff[0]*fSkin[14]+0.28125*coeff[0]*fEdge[14]-0.270632938682637*coeff[2]*fSkin[13]-0.270632938682637*coeff[2]*fEdge[13]-0.28125*coeff[2]*fSkin[10]+0.28125*coeff[2]*fEdge[10]; 
  edgeSurf_incr[15] = (-0.4330127018922193*coeff[1]*fSkin[15])-0.71875*coeff[0]*fSkin[15]+0.4330127018922193*coeff[1]*fEdge[15]-0.21875*coeff[0]*fEdge[15]-0.375*coeff[1]*fSkin[14]-0.7036456405748562*coeff[0]*fSkin[14]-0.375*coeff[1]*fEdge[14]+0.2706329386826369*coeff[0]*fEdge[14]-0.4330127018922193*coeff[3]*fSkin[13]-0.71875*coeff[2]*fSkin[13]+0.4330127018922193*coeff[3]*fEdge[13]-0.21875*coeff[2]*fEdge[13]-0.375*coeff[3]*fSkin[10]-0.7036456405748562*coeff[2]*fSkin[10]-0.375*coeff[3]*fEdge[10]+0.2706329386826369*coeff[2]*fEdge[10]; 
  edgeSurf_incr[16] = (-0.270632938682637*coeff[2]*fSkin[20])-0.270632938682637*coeff[2]*fEdge[20]-0.2812499999999999*coeff[2]*fSkin[18]+0.2812499999999999*coeff[2]*fEdge[18]-0.2706329386826371*coeff[0]*fSkin[17]-0.2706329386826371*coeff[0]*fEdge[17]-0.28125*coeff[0]*fSkin[16]+0.28125*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = (-0.4330127018922194*coeff[3]*fSkin[20])-0.7187500000000001*coeff[2]*fSkin[20]+0.4330127018922194*coeff[3]*fEdge[20]-0.21875*coeff[2]*fEdge[20]-0.375*coeff[3]*fSkin[18]-0.7036456405748562*coeff[2]*fSkin[18]-0.375*coeff[3]*fEdge[18]+0.2706329386826369*coeff[2]*fEdge[18]-0.4330127018922193*coeff[1]*fSkin[17]-0.71875*coeff[0]*fSkin[17]+0.4330127018922193*coeff[1]*fEdge[17]-0.21875*coeff[0]*fEdge[17]-0.375*coeff[1]*fSkin[16]-0.7036456405748563*coeff[0]*fSkin[16]-0.375*coeff[1]*fEdge[16]+0.2706329386826371*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = (-0.2706329386826371*coeff[0]*fSkin[20])-0.2706329386826371*coeff[0]*fEdge[20]-0.28125*coeff[0]*fSkin[18]+0.28125*coeff[0]*fEdge[18]-0.270632938682637*coeff[2]*fSkin[17]-0.270632938682637*coeff[2]*fEdge[17]-0.2812499999999999*coeff[2]*fSkin[16]+0.2812499999999999*coeff[2]*fEdge[16]; 
  edgeSurf_incr[19] = (-0.270632938682637*coeff[2]*fSkin[23])-0.270632938682637*coeff[2]*fEdge[23]-0.2812499999999999*coeff[2]*fSkin[22]+0.2812499999999999*coeff[2]*fEdge[22]-0.2706329386826371*coeff[0]*fSkin[21]-0.2706329386826371*coeff[0]*fEdge[21]-0.28125*coeff[0]*fSkin[19]+0.28125*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = (-0.4330127018922193*coeff[1]*fSkin[20])-0.71875*coeff[0]*fSkin[20]+0.4330127018922193*coeff[1]*fEdge[20]-0.21875*coeff[0]*fEdge[20]-0.375*coeff[1]*fSkin[18]-0.7036456405748563*coeff[0]*fSkin[18]-0.375*coeff[1]*fEdge[18]+0.2706329386826371*coeff[0]*fEdge[18]-0.4330127018922194*coeff[3]*fSkin[17]-0.7187500000000001*coeff[2]*fSkin[17]+0.4330127018922194*coeff[3]*fEdge[17]-0.21875*coeff[2]*fEdge[17]-0.375*coeff[3]*fSkin[16]-0.7036456405748562*coeff[2]*fSkin[16]-0.375*coeff[3]*fEdge[16]+0.2706329386826369*coeff[2]*fEdge[16]; 
  edgeSurf_incr[21] = (-0.4330127018922194*coeff[3]*fSkin[23])-0.7187500000000001*coeff[2]*fSkin[23]+0.4330127018922194*coeff[3]*fEdge[23]-0.21875*coeff[2]*fEdge[23]-0.375*coeff[3]*fSkin[22]-0.7036456405748562*coeff[2]*fSkin[22]-0.375*coeff[3]*fEdge[22]+0.2706329386826369*coeff[2]*fEdge[22]-0.4330127018922193*coeff[1]*fSkin[21]-0.71875*coeff[0]*fSkin[21]+0.4330127018922193*coeff[1]*fEdge[21]-0.21875*coeff[0]*fEdge[21]-0.375*coeff[1]*fSkin[19]-0.7036456405748563*coeff[0]*fSkin[19]-0.375*coeff[1]*fEdge[19]+0.2706329386826371*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = (-0.2706329386826371*coeff[0]*fSkin[23])-0.2706329386826371*coeff[0]*fEdge[23]-0.28125*coeff[0]*fSkin[22]+0.28125*coeff[0]*fEdge[22]-0.270632938682637*coeff[2]*fSkin[21]-0.270632938682637*coeff[2]*fEdge[21]-0.2812499999999999*coeff[2]*fSkin[19]+0.2812499999999999*coeff[2]*fEdge[19]; 
  edgeSurf_incr[23] = (-0.4330127018922193*coeff[1]*fSkin[23])-0.71875*coeff[0]*fSkin[23]+0.4330127018922193*coeff[1]*fEdge[23]-0.21875*coeff[0]*fEdge[23]-0.375*coeff[1]*fSkin[22]-0.7036456405748563*coeff[0]*fSkin[22]-0.375*coeff[1]*fEdge[22]+0.2706329386826371*coeff[0]*fEdge[22]-0.4330127018922194*coeff[3]*fSkin[21]-0.7187500000000001*coeff[2]*fSkin[21]+0.4330127018922194*coeff[3]*fEdge[21]-0.21875*coeff[2]*fEdge[21]-0.375*coeff[3]*fSkin[19]-0.7036456405748562*coeff[2]*fSkin[19]-0.375*coeff[3]*fEdge[19]+0.2706329386826369*coeff[2]*fEdge[19]; 

  boundSurf_incr[1] = 0.8660254037844386*coeff[3]*fSkin[5]-0.5*coeff[2]*fSkin[5]-0.75*fSkin[2]*coeff[3]+0.4330127018922193*coeff[2]*fSkin[2]+0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]+0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = 0.8660254037844386*coeff[1]*fSkin[5]-0.5*coeff[0]*fSkin[5]+0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]+0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]+0.4330127018922193*fSkin[0]*coeff[2]; 
  boundSurf_incr[6] = 0.8660254037844386*coeff[3]*fSkin[11]-0.5*coeff[2]*fSkin[11]-0.75*coeff[3]*fSkin[7]+0.4330127018922193*coeff[2]*fSkin[7]+0.8660254037844386*coeff[1]*fSkin[6]-0.5*coeff[0]*fSkin[6]-0.75*coeff[1]*fSkin[3]+0.4330127018922193*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = 0.8660254037844386*coeff[3]*fSkin[12]-0.5*coeff[2]*fSkin[12]-0.75*coeff[3]*fSkin[9]+0.4330127018922193*coeff[2]*fSkin[9]+0.8660254037844386*coeff[1]*fSkin[8]-0.5*coeff[0]*fSkin[8]-0.75*coeff[1]*fSkin[4]+0.4330127018922193*coeff[0]*fSkin[4]; 
  boundSurf_incr[11] = 0.8660254037844386*coeff[1]*fSkin[11]-0.5*coeff[0]*fSkin[11]-0.75*coeff[1]*fSkin[7]+0.4330127018922193*coeff[0]*fSkin[7]+0.8660254037844386*coeff[3]*fSkin[6]-0.5*coeff[2]*fSkin[6]-0.75*coeff[3]*fSkin[3]+0.4330127018922193*coeff[2]*fSkin[3]; 
  boundSurf_incr[12] = 0.8660254037844386*coeff[1]*fSkin[12]-0.5*coeff[0]*fSkin[12]-0.75*coeff[1]*fSkin[9]+0.4330127018922193*coeff[0]*fSkin[9]+0.8660254037844386*coeff[3]*fSkin[8]-0.5*coeff[2]*fSkin[8]-0.75*coeff[3]*fSkin[4]+0.4330127018922193*coeff[2]*fSkin[4]; 
  boundSurf_incr[13] = 0.8660254037844386*coeff[3]*fSkin[15]-0.5*coeff[2]*fSkin[15]-0.75*coeff[3]*fSkin[14]+0.4330127018922193*coeff[2]*fSkin[14]+0.8660254037844386*coeff[1]*fSkin[13]-0.5*coeff[0]*fSkin[13]-0.75*coeff[1]*fSkin[10]+0.4330127018922193*coeff[0]*fSkin[10]; 
  boundSurf_incr[15] = 0.8660254037844386*coeff[1]*fSkin[15]-0.5*coeff[0]*fSkin[15]-0.75*coeff[1]*fSkin[14]+0.4330127018922193*coeff[0]*fSkin[14]+0.8660254037844386*coeff[3]*fSkin[13]-0.5*coeff[2]*fSkin[13]-0.75*coeff[3]*fSkin[10]+0.4330127018922193*coeff[2]*fSkin[10]; 
  boundSurf_incr[17] = 0.8660254037844387*coeff[3]*fSkin[20]-0.5000000000000001*coeff[2]*fSkin[20]-0.75*coeff[3]*fSkin[18]+0.4330127018922193*coeff[2]*fSkin[18]+0.8660254037844386*coeff[1]*fSkin[17]-0.5*coeff[0]*fSkin[17]-0.75*coeff[1]*fSkin[16]+0.4330127018922194*coeff[0]*fSkin[16]; 
  boundSurf_incr[20] = 0.8660254037844386*coeff[1]*fSkin[20]-0.5*coeff[0]*fSkin[20]-0.75*coeff[1]*fSkin[18]+0.4330127018922194*coeff[0]*fSkin[18]+0.8660254037844387*coeff[3]*fSkin[17]-0.5000000000000001*coeff[2]*fSkin[17]-0.75*coeff[3]*fSkin[16]+0.4330127018922193*coeff[2]*fSkin[16]; 
  boundSurf_incr[21] = 0.8660254037844387*coeff[3]*fSkin[23]-0.5000000000000001*coeff[2]*fSkin[23]-0.75*coeff[3]*fSkin[22]+0.4330127018922193*coeff[2]*fSkin[22]+0.8660254037844386*coeff[1]*fSkin[21]-0.5*coeff[0]*fSkin[21]-0.75*coeff[1]*fSkin[19]+0.4330127018922194*coeff[0]*fSkin[19]; 
  boundSurf_incr[23] = 0.8660254037844386*coeff[1]*fSkin[23]-0.5*coeff[0]*fSkin[23]-0.75*coeff[1]*fSkin[22]+0.4330127018922194*coeff[0]*fSkin[22]+0.8660254037844387*coeff[3]*fSkin[21]-0.5000000000000001*coeff[2]*fSkin[21]-0.75*coeff[3]*fSkin[19]+0.4330127018922193*coeff[2]*fSkin[19]; 

  } else { 

  edgeSurf_incr[0] = 0.270632938682637*coeff[2]*fSkin[5]+0.270632938682637*coeff[2]*fEdge[5]-0.28125*coeff[2]*fSkin[2]+0.28125*coeff[2]*fEdge[2]+0.270632938682637*coeff[0]*fSkin[1]+0.270632938682637*coeff[0]*fEdge[1]-0.28125*coeff[0]*fSkin[0]+0.28125*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 0.4330127018922193*coeff[3]*fSkin[5]-0.71875*coeff[2]*fSkin[5]-0.4330127018922193*coeff[3]*fEdge[5]-0.21875*coeff[2]*fEdge[5]-0.375*fSkin[2]*coeff[3]-0.375*fEdge[2]*coeff[3]+0.7036456405748562*coeff[2]*fSkin[2]-0.2706329386826369*coeff[2]*fEdge[2]+0.4330127018922193*coeff[1]*fSkin[1]-0.71875*coeff[0]*fSkin[1]-0.4330127018922193*coeff[1]*fEdge[1]-0.21875*coeff[0]*fEdge[1]-0.375*fSkin[0]*coeff[1]-0.375*fEdge[0]*coeff[1]+0.7036456405748562*coeff[0]*fSkin[0]-0.2706329386826369*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 0.270632938682637*coeff[0]*fSkin[5]+0.270632938682637*coeff[0]*fEdge[5]-0.28125*coeff[0]*fSkin[2]+0.28125*coeff[0]*fEdge[2]+0.270632938682637*fSkin[1]*coeff[2]+0.270632938682637*fEdge[1]*coeff[2]-0.28125*fSkin[0]*coeff[2]+0.28125*fEdge[0]*coeff[2]; 
  edgeSurf_incr[3] = 0.270632938682637*coeff[2]*fSkin[11]+0.270632938682637*coeff[2]*fEdge[11]-0.28125*coeff[2]*fSkin[7]+0.28125*coeff[2]*fEdge[7]+0.270632938682637*coeff[0]*fSkin[6]+0.270632938682637*coeff[0]*fEdge[6]-0.28125*coeff[0]*fSkin[3]+0.28125*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 0.270632938682637*coeff[2]*fSkin[12]+0.270632938682637*coeff[2]*fEdge[12]-0.28125*coeff[2]*fSkin[9]+0.28125*coeff[2]*fEdge[9]+0.270632938682637*coeff[0]*fSkin[8]+0.270632938682637*coeff[0]*fEdge[8]-0.28125*coeff[0]*fSkin[4]+0.28125*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 0.4330127018922193*coeff[1]*fSkin[5]-0.71875*coeff[0]*fSkin[5]-0.4330127018922193*coeff[1]*fEdge[5]-0.21875*coeff[0]*fEdge[5]+0.4330127018922193*fSkin[1]*coeff[3]-0.4330127018922193*fEdge[1]*coeff[3]-0.375*fSkin[0]*coeff[3]-0.375*fEdge[0]*coeff[3]-0.375*coeff[1]*fSkin[2]+0.7036456405748562*coeff[0]*fSkin[2]-0.375*coeff[1]*fEdge[2]-0.2706329386826369*coeff[0]*fEdge[2]-0.71875*fSkin[1]*coeff[2]-0.21875*fEdge[1]*coeff[2]+0.7036456405748562*fSkin[0]*coeff[2]-0.2706329386826369*fEdge[0]*coeff[2]; 
  edgeSurf_incr[6] = 0.4330127018922193*coeff[3]*fSkin[11]-0.71875*coeff[2]*fSkin[11]-0.4330127018922193*coeff[3]*fEdge[11]-0.21875*coeff[2]*fEdge[11]-0.375*coeff[3]*fSkin[7]+0.7036456405748562*coeff[2]*fSkin[7]-0.375*coeff[3]*fEdge[7]-0.2706329386826369*coeff[2]*fEdge[7]+0.4330127018922193*coeff[1]*fSkin[6]-0.71875*coeff[0]*fSkin[6]-0.4330127018922193*coeff[1]*fEdge[6]-0.21875*coeff[0]*fEdge[6]-0.375*coeff[1]*fSkin[3]+0.7036456405748562*coeff[0]*fSkin[3]-0.375*coeff[1]*fEdge[3]-0.2706329386826369*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 0.270632938682637*coeff[0]*fSkin[11]+0.270632938682637*coeff[0]*fEdge[11]-0.28125*coeff[0]*fSkin[7]+0.28125*coeff[0]*fEdge[7]+0.270632938682637*coeff[2]*fSkin[6]+0.270632938682637*coeff[2]*fEdge[6]-0.28125*coeff[2]*fSkin[3]+0.28125*coeff[2]*fEdge[3]; 
  edgeSurf_incr[8] = 0.4330127018922193*coeff[3]*fSkin[12]-0.71875*coeff[2]*fSkin[12]-0.4330127018922193*coeff[3]*fEdge[12]-0.21875*coeff[2]*fEdge[12]-0.375*coeff[3]*fSkin[9]+0.7036456405748562*coeff[2]*fSkin[9]-0.375*coeff[3]*fEdge[9]-0.2706329386826369*coeff[2]*fEdge[9]+0.4330127018922193*coeff[1]*fSkin[8]-0.71875*coeff[0]*fSkin[8]-0.4330127018922193*coeff[1]*fEdge[8]-0.21875*coeff[0]*fEdge[8]-0.375*coeff[1]*fSkin[4]+0.7036456405748562*coeff[0]*fSkin[4]-0.375*coeff[1]*fEdge[4]-0.2706329386826369*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 0.270632938682637*coeff[0]*fSkin[12]+0.270632938682637*coeff[0]*fEdge[12]-0.28125*coeff[0]*fSkin[9]+0.28125*coeff[0]*fEdge[9]+0.270632938682637*coeff[2]*fSkin[8]+0.270632938682637*coeff[2]*fEdge[8]-0.28125*coeff[2]*fSkin[4]+0.28125*coeff[2]*fEdge[4]; 
  edgeSurf_incr[10] = 0.270632938682637*coeff[2]*fSkin[15]+0.270632938682637*coeff[2]*fEdge[15]-0.28125*coeff[2]*fSkin[14]+0.28125*coeff[2]*fEdge[14]+0.270632938682637*coeff[0]*fSkin[13]+0.270632938682637*coeff[0]*fEdge[13]-0.28125*coeff[0]*fSkin[10]+0.28125*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 0.4330127018922193*coeff[1]*fSkin[11]-0.71875*coeff[0]*fSkin[11]-0.4330127018922193*coeff[1]*fEdge[11]-0.21875*coeff[0]*fEdge[11]-0.375*coeff[1]*fSkin[7]+0.7036456405748562*coeff[0]*fSkin[7]-0.375*coeff[1]*fEdge[7]-0.2706329386826369*coeff[0]*fEdge[7]+0.4330127018922193*coeff[3]*fSkin[6]-0.71875*coeff[2]*fSkin[6]-0.4330127018922193*coeff[3]*fEdge[6]-0.21875*coeff[2]*fEdge[6]-0.375*coeff[3]*fSkin[3]+0.7036456405748562*coeff[2]*fSkin[3]-0.375*coeff[3]*fEdge[3]-0.2706329386826369*coeff[2]*fEdge[3]; 
  edgeSurf_incr[12] = 0.4330127018922193*coeff[1]*fSkin[12]-0.71875*coeff[0]*fSkin[12]-0.4330127018922193*coeff[1]*fEdge[12]-0.21875*coeff[0]*fEdge[12]-0.375*coeff[1]*fSkin[9]+0.7036456405748562*coeff[0]*fSkin[9]-0.375*coeff[1]*fEdge[9]-0.2706329386826369*coeff[0]*fEdge[9]+0.4330127018922193*coeff[3]*fSkin[8]-0.71875*coeff[2]*fSkin[8]-0.4330127018922193*coeff[3]*fEdge[8]-0.21875*coeff[2]*fEdge[8]-0.375*coeff[3]*fSkin[4]+0.7036456405748562*coeff[2]*fSkin[4]-0.375*coeff[3]*fEdge[4]-0.2706329386826369*coeff[2]*fEdge[4]; 
  edgeSurf_incr[13] = 0.4330127018922193*coeff[3]*fSkin[15]-0.71875*coeff[2]*fSkin[15]-0.4330127018922193*coeff[3]*fEdge[15]-0.21875*coeff[2]*fEdge[15]-0.375*coeff[3]*fSkin[14]+0.7036456405748562*coeff[2]*fSkin[14]-0.375*coeff[3]*fEdge[14]-0.2706329386826369*coeff[2]*fEdge[14]+0.4330127018922193*coeff[1]*fSkin[13]-0.71875*coeff[0]*fSkin[13]-0.4330127018922193*coeff[1]*fEdge[13]-0.21875*coeff[0]*fEdge[13]-0.375*coeff[1]*fSkin[10]+0.7036456405748562*coeff[0]*fSkin[10]-0.375*coeff[1]*fEdge[10]-0.2706329386826369*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = 0.270632938682637*coeff[0]*fSkin[15]+0.270632938682637*coeff[0]*fEdge[15]-0.28125*coeff[0]*fSkin[14]+0.28125*coeff[0]*fEdge[14]+0.270632938682637*coeff[2]*fSkin[13]+0.270632938682637*coeff[2]*fEdge[13]-0.28125*coeff[2]*fSkin[10]+0.28125*coeff[2]*fEdge[10]; 
  edgeSurf_incr[15] = 0.4330127018922193*coeff[1]*fSkin[15]-0.71875*coeff[0]*fSkin[15]-0.4330127018922193*coeff[1]*fEdge[15]-0.21875*coeff[0]*fEdge[15]-0.375*coeff[1]*fSkin[14]+0.7036456405748562*coeff[0]*fSkin[14]-0.375*coeff[1]*fEdge[14]-0.2706329386826369*coeff[0]*fEdge[14]+0.4330127018922193*coeff[3]*fSkin[13]-0.71875*coeff[2]*fSkin[13]-0.4330127018922193*coeff[3]*fEdge[13]-0.21875*coeff[2]*fEdge[13]-0.375*coeff[3]*fSkin[10]+0.7036456405748562*coeff[2]*fSkin[10]-0.375*coeff[3]*fEdge[10]-0.2706329386826369*coeff[2]*fEdge[10]; 
  edgeSurf_incr[16] = 0.270632938682637*coeff[2]*fSkin[20]+0.270632938682637*coeff[2]*fEdge[20]-0.2812499999999999*coeff[2]*fSkin[18]+0.2812499999999999*coeff[2]*fEdge[18]+0.2706329386826371*coeff[0]*fSkin[17]+0.2706329386826371*coeff[0]*fEdge[17]-0.28125*coeff[0]*fSkin[16]+0.28125*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 0.4330127018922194*coeff[3]*fSkin[20]-0.7187500000000001*coeff[2]*fSkin[20]-0.4330127018922194*coeff[3]*fEdge[20]-0.21875*coeff[2]*fEdge[20]-0.375*coeff[3]*fSkin[18]+0.7036456405748562*coeff[2]*fSkin[18]-0.375*coeff[3]*fEdge[18]-0.2706329386826369*coeff[2]*fEdge[18]+0.4330127018922193*coeff[1]*fSkin[17]-0.71875*coeff[0]*fSkin[17]-0.4330127018922193*coeff[1]*fEdge[17]-0.21875*coeff[0]*fEdge[17]-0.375*coeff[1]*fSkin[16]+0.7036456405748563*coeff[0]*fSkin[16]-0.375*coeff[1]*fEdge[16]-0.2706329386826371*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = 0.2706329386826371*coeff[0]*fSkin[20]+0.2706329386826371*coeff[0]*fEdge[20]-0.28125*coeff[0]*fSkin[18]+0.28125*coeff[0]*fEdge[18]+0.270632938682637*coeff[2]*fSkin[17]+0.270632938682637*coeff[2]*fEdge[17]-0.2812499999999999*coeff[2]*fSkin[16]+0.2812499999999999*coeff[2]*fEdge[16]; 
  edgeSurf_incr[19] = 0.270632938682637*coeff[2]*fSkin[23]+0.270632938682637*coeff[2]*fEdge[23]-0.2812499999999999*coeff[2]*fSkin[22]+0.2812499999999999*coeff[2]*fEdge[22]+0.2706329386826371*coeff[0]*fSkin[21]+0.2706329386826371*coeff[0]*fEdge[21]-0.28125*coeff[0]*fSkin[19]+0.28125*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 0.4330127018922193*coeff[1]*fSkin[20]-0.71875*coeff[0]*fSkin[20]-0.4330127018922193*coeff[1]*fEdge[20]-0.21875*coeff[0]*fEdge[20]-0.375*coeff[1]*fSkin[18]+0.7036456405748563*coeff[0]*fSkin[18]-0.375*coeff[1]*fEdge[18]-0.2706329386826371*coeff[0]*fEdge[18]+0.4330127018922194*coeff[3]*fSkin[17]-0.7187500000000001*coeff[2]*fSkin[17]-0.4330127018922194*coeff[3]*fEdge[17]-0.21875*coeff[2]*fEdge[17]-0.375*coeff[3]*fSkin[16]+0.7036456405748562*coeff[2]*fSkin[16]-0.375*coeff[3]*fEdge[16]-0.2706329386826369*coeff[2]*fEdge[16]; 
  edgeSurf_incr[21] = 0.4330127018922194*coeff[3]*fSkin[23]-0.7187500000000001*coeff[2]*fSkin[23]-0.4330127018922194*coeff[3]*fEdge[23]-0.21875*coeff[2]*fEdge[23]-0.375*coeff[3]*fSkin[22]+0.7036456405748562*coeff[2]*fSkin[22]-0.375*coeff[3]*fEdge[22]-0.2706329386826369*coeff[2]*fEdge[22]+0.4330127018922193*coeff[1]*fSkin[21]-0.71875*coeff[0]*fSkin[21]-0.4330127018922193*coeff[1]*fEdge[21]-0.21875*coeff[0]*fEdge[21]-0.375*coeff[1]*fSkin[19]+0.7036456405748563*coeff[0]*fSkin[19]-0.375*coeff[1]*fEdge[19]-0.2706329386826371*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = 0.2706329386826371*coeff[0]*fSkin[23]+0.2706329386826371*coeff[0]*fEdge[23]-0.28125*coeff[0]*fSkin[22]+0.28125*coeff[0]*fEdge[22]+0.270632938682637*coeff[2]*fSkin[21]+0.270632938682637*coeff[2]*fEdge[21]-0.2812499999999999*coeff[2]*fSkin[19]+0.2812499999999999*coeff[2]*fEdge[19]; 
  edgeSurf_incr[23] = 0.4330127018922193*coeff[1]*fSkin[23]-0.71875*coeff[0]*fSkin[23]-0.4330127018922193*coeff[1]*fEdge[23]-0.21875*coeff[0]*fEdge[23]-0.375*coeff[1]*fSkin[22]+0.7036456405748563*coeff[0]*fSkin[22]-0.375*coeff[1]*fEdge[22]-0.2706329386826371*coeff[0]*fEdge[22]+0.4330127018922194*coeff[3]*fSkin[21]-0.7187500000000001*coeff[2]*fSkin[21]-0.4330127018922194*coeff[3]*fEdge[21]-0.21875*coeff[2]*fEdge[21]-0.375*coeff[3]*fSkin[19]+0.7036456405748562*coeff[2]*fSkin[19]-0.375*coeff[3]*fEdge[19]-0.2706329386826369*coeff[2]*fEdge[19]; 

  boundSurf_incr[1] = (-0.8660254037844386*coeff[3]*fSkin[5])-0.5*coeff[2]*fSkin[5]-0.75*fSkin[2]*coeff[3]-0.4330127018922193*coeff[2]*fSkin[2]-0.8660254037844386*coeff[1]*fSkin[1]-0.5*coeff[0]*fSkin[1]-0.75*fSkin[0]*coeff[1]-0.4330127018922193*coeff[0]*fSkin[0]; 
  boundSurf_incr[5] = (-0.8660254037844386*coeff[1]*fSkin[5])-0.5*coeff[0]*fSkin[5]-0.8660254037844386*fSkin[1]*coeff[3]-0.75*fSkin[0]*coeff[3]-0.75*coeff[1]*fSkin[2]-0.4330127018922193*coeff[0]*fSkin[2]-0.5*fSkin[1]*coeff[2]-0.4330127018922193*fSkin[0]*coeff[2]; 
  boundSurf_incr[6] = (-0.8660254037844386*coeff[3]*fSkin[11])-0.5*coeff[2]*fSkin[11]-0.75*coeff[3]*fSkin[7]-0.4330127018922193*coeff[2]*fSkin[7]-0.8660254037844386*coeff[1]*fSkin[6]-0.5*coeff[0]*fSkin[6]-0.75*coeff[1]*fSkin[3]-0.4330127018922193*coeff[0]*fSkin[3]; 
  boundSurf_incr[8] = (-0.8660254037844386*coeff[3]*fSkin[12])-0.5*coeff[2]*fSkin[12]-0.75*coeff[3]*fSkin[9]-0.4330127018922193*coeff[2]*fSkin[9]-0.8660254037844386*coeff[1]*fSkin[8]-0.5*coeff[0]*fSkin[8]-0.75*coeff[1]*fSkin[4]-0.4330127018922193*coeff[0]*fSkin[4]; 
  boundSurf_incr[11] = (-0.8660254037844386*coeff[1]*fSkin[11])-0.5*coeff[0]*fSkin[11]-0.75*coeff[1]*fSkin[7]-0.4330127018922193*coeff[0]*fSkin[7]-0.8660254037844386*coeff[3]*fSkin[6]-0.5*coeff[2]*fSkin[6]-0.75*coeff[3]*fSkin[3]-0.4330127018922193*coeff[2]*fSkin[3]; 
  boundSurf_incr[12] = (-0.8660254037844386*coeff[1]*fSkin[12])-0.5*coeff[0]*fSkin[12]-0.75*coeff[1]*fSkin[9]-0.4330127018922193*coeff[0]*fSkin[9]-0.8660254037844386*coeff[3]*fSkin[8]-0.5*coeff[2]*fSkin[8]-0.75*coeff[3]*fSkin[4]-0.4330127018922193*coeff[2]*fSkin[4]; 
  boundSurf_incr[13] = (-0.8660254037844386*coeff[3]*fSkin[15])-0.5*coeff[2]*fSkin[15]-0.75*coeff[3]*fSkin[14]-0.4330127018922193*coeff[2]*fSkin[14]-0.8660254037844386*coeff[1]*fSkin[13]-0.5*coeff[0]*fSkin[13]-0.75*coeff[1]*fSkin[10]-0.4330127018922193*coeff[0]*fSkin[10]; 
  boundSurf_incr[15] = (-0.8660254037844386*coeff[1]*fSkin[15])-0.5*coeff[0]*fSkin[15]-0.75*coeff[1]*fSkin[14]-0.4330127018922193*coeff[0]*fSkin[14]-0.8660254037844386*coeff[3]*fSkin[13]-0.5*coeff[2]*fSkin[13]-0.75*coeff[3]*fSkin[10]-0.4330127018922193*coeff[2]*fSkin[10]; 
  boundSurf_incr[17] = (-0.8660254037844387*coeff[3]*fSkin[20])-0.5000000000000001*coeff[2]*fSkin[20]-0.75*coeff[3]*fSkin[18]-0.4330127018922193*coeff[2]*fSkin[18]-0.8660254037844386*coeff[1]*fSkin[17]-0.5*coeff[0]*fSkin[17]-0.75*coeff[1]*fSkin[16]-0.4330127018922194*coeff[0]*fSkin[16]; 
  boundSurf_incr[20] = (-0.8660254037844386*coeff[1]*fSkin[20])-0.5*coeff[0]*fSkin[20]-0.75*coeff[1]*fSkin[18]-0.4330127018922194*coeff[0]*fSkin[18]-0.8660254037844387*coeff[3]*fSkin[17]-0.5000000000000001*coeff[2]*fSkin[17]-0.75*coeff[3]*fSkin[16]-0.4330127018922193*coeff[2]*fSkin[16]; 
  boundSurf_incr[21] = (-0.8660254037844387*coeff[3]*fSkin[23])-0.5000000000000001*coeff[2]*fSkin[23]-0.75*coeff[3]*fSkin[22]-0.4330127018922193*coeff[2]*fSkin[22]-0.8660254037844386*coeff[1]*fSkin[21]-0.5*coeff[0]*fSkin[21]-0.75*coeff[1]*fSkin[19]-0.4330127018922194*coeff[0]*fSkin[19]; 
  boundSurf_incr[23] = (-0.8660254037844386*coeff[1]*fSkin[23])-0.5*coeff[0]*fSkin[23]-0.75*coeff[1]*fSkin[22]-0.4330127018922194*coeff[0]*fSkin[22]-0.8660254037844387*coeff[3]*fSkin[21]-0.5000000000000001*coeff[2]*fSkin[21]-0.75*coeff[3]*fSkin[19]-0.4330127018922193*coeff[2]*fSkin[19]; 

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

