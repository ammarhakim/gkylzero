#include <gkyl_dg_diffusion_gyrokinetic_kernels.h>

GKYL_CU_DH double dg_diffusion_gyrokinetic_order4_boundary_surfx_2x2v_ser_p1_constcoeff(const double *w, const double *dx, const double *coeff, int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{
  // w[NDIM]: Cell-center coordinate.
  // dxv[NDIM]: Cell length.
  // coeff: Diffusion coefficient.
  // edge: -1 for lower boundary, +1 for upper boundary.
  // fSkin/Edge: scalar field in skind and egde cells.
  // out: Incremented output.

  const double Jfac = pow(2./dx[0],4.);

  double vol_incr[24] = {0.0}; 

  double edgeSurf_incr[24] = {0.0}; 
  double boundSurf_incr[24] = {0.0}; 

  if (edge == -1) { 

  edgeSurf_incr[0] = 1.623797632095822*coeff[0]*fSkin[1]+1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]+1.623797632095822*coeff[0]*fSkin[0]-1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = 1.623797632095822*coeff[0]*fSkin[5]+1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = 1.623797632095822*coeff[0]*fSkin[6]+1.623797632095822*coeff[0]*fEdge[6]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = 1.623797632095822*coeff[0]*fSkin[8]+1.623797632095822*coeff[0]*fEdge[8]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]+1.623797632095822*coeff[0]*fSkin[2]-1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*fSkin[6]+2.0625*coeff[0]*fEdge[6]+1.623797632095822*coeff[0]*fSkin[3]-1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = 1.623797632095822*coeff[0]*fSkin[11]+1.623797632095822*coeff[0]*fEdge[11]+0.9375*coeff[0]*fSkin[7]-0.9375*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*fSkin[8]+2.0625*coeff[0]*fEdge[8]+1.623797632095822*coeff[0]*fSkin[4]-1.623797632095822*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = 1.623797632095822*coeff[0]*fSkin[12]+1.623797632095822*coeff[0]*fEdge[12]+0.9375*coeff[0]*fSkin[9]-0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = 1.623797632095822*coeff[0]*fSkin[13]+1.623797632095822*coeff[0]*fEdge[13]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*fSkin[11]+2.0625*coeff[0]*fEdge[11]+1.623797632095822*coeff[0]*fSkin[7]-1.623797632095822*coeff[0]*fEdge[7]; 
  edgeSurf_incr[12] = 3.5625*coeff[0]*fSkin[12]+2.0625*coeff[0]*fEdge[12]+1.623797632095822*coeff[0]*fSkin[9]-1.623797632095822*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*fSkin[13]+2.0625*coeff[0]*fEdge[13]+1.623797632095822*coeff[0]*fSkin[10]-1.623797632095822*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = 1.623797632095822*coeff[0]*fSkin[15]+1.623797632095822*coeff[0]*fEdge[15]+0.9375*coeff[0]*fSkin[14]-0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*fSkin[15]+2.0625*coeff[0]*fEdge[15]+1.623797632095822*coeff[0]*fSkin[14]-1.623797632095822*coeff[0]*fEdge[14]; 
  edgeSurf_incr[16] = 1.623797632095823*coeff[0]*fSkin[17]+1.623797632095823*coeff[0]*fEdge[17]+0.9375*coeff[0]*fSkin[16]-0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*fSkin[17]+2.0625*coeff[0]*fEdge[17]+1.623797632095823*coeff[0]*fSkin[16]-1.623797632095823*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = 1.623797632095823*coeff[0]*fSkin[20]+1.623797632095823*coeff[0]*fEdge[20]+0.9375*coeff[0]*fSkin[18]-0.9375*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = 1.623797632095823*coeff[0]*fSkin[21]+1.623797632095823*coeff[0]*fEdge[21]+0.9375*coeff[0]*fSkin[19]-0.9375*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*fSkin[20]+2.0625*coeff[0]*fEdge[20]+1.623797632095823*coeff[0]*fSkin[18]-1.623797632095823*coeff[0]*fEdge[18]; 
  edgeSurf_incr[21] = 3.5625*coeff[0]*fSkin[21]+2.0625*coeff[0]*fEdge[21]+1.623797632095823*coeff[0]*fSkin[19]-1.623797632095823*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = 1.623797632095823*coeff[0]*fSkin[23]+1.623797632095823*coeff[0]*fEdge[23]+0.9375*coeff[0]*fSkin[22]-0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*fSkin[23]+2.0625*coeff[0]*fEdge[23]+1.623797632095823*coeff[0]*fSkin[22]-1.623797632095823*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[0]*fSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 1.5*coeff[0]*fSkin[13]; 
  boundSurf_incr[15] = 1.5*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*fSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 1.5*coeff[0]*fSkin[21]; 
  boundSurf_incr[23] = 1.5*coeff[0]*fSkin[23]; 

  } else { 

  edgeSurf_incr[0] = (-1.623797632095822*coeff[0]*fSkin[1])-1.623797632095822*coeff[0]*fEdge[1]+0.9375*coeff[0]*fSkin[0]-0.9375*coeff[0]*fEdge[0]; 
  edgeSurf_incr[1] = 3.5625*coeff[0]*fSkin[1]+2.0625*coeff[0]*fEdge[1]-1.623797632095822*coeff[0]*fSkin[0]+1.623797632095822*coeff[0]*fEdge[0]; 
  edgeSurf_incr[2] = (-1.623797632095822*coeff[0]*fSkin[5])-1.623797632095822*coeff[0]*fEdge[5]+0.9375*coeff[0]*fSkin[2]-0.9375*coeff[0]*fEdge[2]; 
  edgeSurf_incr[3] = (-1.623797632095822*coeff[0]*fSkin[6])-1.623797632095822*coeff[0]*fEdge[6]+0.9375*coeff[0]*fSkin[3]-0.9375*coeff[0]*fEdge[3]; 
  edgeSurf_incr[4] = (-1.623797632095822*coeff[0]*fSkin[8])-1.623797632095822*coeff[0]*fEdge[8]+0.9375*coeff[0]*fSkin[4]-0.9375*coeff[0]*fEdge[4]; 
  edgeSurf_incr[5] = 3.5625*coeff[0]*fSkin[5]+2.0625*coeff[0]*fEdge[5]-1.623797632095822*coeff[0]*fSkin[2]+1.623797632095822*coeff[0]*fEdge[2]; 
  edgeSurf_incr[6] = 3.5625*coeff[0]*fSkin[6]+2.0625*coeff[0]*fEdge[6]-1.623797632095822*coeff[0]*fSkin[3]+1.623797632095822*coeff[0]*fEdge[3]; 
  edgeSurf_incr[7] = (-1.623797632095822*coeff[0]*fSkin[11])-1.623797632095822*coeff[0]*fEdge[11]+0.9375*coeff[0]*fSkin[7]-0.9375*coeff[0]*fEdge[7]; 
  edgeSurf_incr[8] = 3.5625*coeff[0]*fSkin[8]+2.0625*coeff[0]*fEdge[8]-1.623797632095822*coeff[0]*fSkin[4]+1.623797632095822*coeff[0]*fEdge[4]; 
  edgeSurf_incr[9] = (-1.623797632095822*coeff[0]*fSkin[12])-1.623797632095822*coeff[0]*fEdge[12]+0.9375*coeff[0]*fSkin[9]-0.9375*coeff[0]*fEdge[9]; 
  edgeSurf_incr[10] = (-1.623797632095822*coeff[0]*fSkin[13])-1.623797632095822*coeff[0]*fEdge[13]+0.9375*coeff[0]*fSkin[10]-0.9375*coeff[0]*fEdge[10]; 
  edgeSurf_incr[11] = 3.5625*coeff[0]*fSkin[11]+2.0625*coeff[0]*fEdge[11]-1.623797632095822*coeff[0]*fSkin[7]+1.623797632095822*coeff[0]*fEdge[7]; 
  edgeSurf_incr[12] = 3.5625*coeff[0]*fSkin[12]+2.0625*coeff[0]*fEdge[12]-1.623797632095822*coeff[0]*fSkin[9]+1.623797632095822*coeff[0]*fEdge[9]; 
  edgeSurf_incr[13] = 3.5625*coeff[0]*fSkin[13]+2.0625*coeff[0]*fEdge[13]-1.623797632095822*coeff[0]*fSkin[10]+1.623797632095822*coeff[0]*fEdge[10]; 
  edgeSurf_incr[14] = (-1.623797632095822*coeff[0]*fSkin[15])-1.623797632095822*coeff[0]*fEdge[15]+0.9375*coeff[0]*fSkin[14]-0.9375*coeff[0]*fEdge[14]; 
  edgeSurf_incr[15] = 3.5625*coeff[0]*fSkin[15]+2.0625*coeff[0]*fEdge[15]-1.623797632095822*coeff[0]*fSkin[14]+1.623797632095822*coeff[0]*fEdge[14]; 
  edgeSurf_incr[16] = (-1.623797632095823*coeff[0]*fSkin[17])-1.623797632095823*coeff[0]*fEdge[17]+0.9375*coeff[0]*fSkin[16]-0.9375*coeff[0]*fEdge[16]; 
  edgeSurf_incr[17] = 3.5625*coeff[0]*fSkin[17]+2.0625*coeff[0]*fEdge[17]-1.623797632095823*coeff[0]*fSkin[16]+1.623797632095823*coeff[0]*fEdge[16]; 
  edgeSurf_incr[18] = (-1.623797632095823*coeff[0]*fSkin[20])-1.623797632095823*coeff[0]*fEdge[20]+0.9375*coeff[0]*fSkin[18]-0.9375*coeff[0]*fEdge[18]; 
  edgeSurf_incr[19] = (-1.623797632095823*coeff[0]*fSkin[21])-1.623797632095823*coeff[0]*fEdge[21]+0.9375*coeff[0]*fSkin[19]-0.9375*coeff[0]*fEdge[19]; 
  edgeSurf_incr[20] = 3.5625*coeff[0]*fSkin[20]+2.0625*coeff[0]*fEdge[20]-1.623797632095823*coeff[0]*fSkin[18]+1.623797632095823*coeff[0]*fEdge[18]; 
  edgeSurf_incr[21] = 3.5625*coeff[0]*fSkin[21]+2.0625*coeff[0]*fEdge[21]-1.623797632095823*coeff[0]*fSkin[19]+1.623797632095823*coeff[0]*fEdge[19]; 
  edgeSurf_incr[22] = (-1.623797632095823*coeff[0]*fSkin[23])-1.623797632095823*coeff[0]*fEdge[23]+0.9375*coeff[0]*fSkin[22]-0.9375*coeff[0]*fEdge[22]; 
  edgeSurf_incr[23] = 3.5625*coeff[0]*fSkin[23]+2.0625*coeff[0]*fEdge[23]-1.623797632095823*coeff[0]*fSkin[22]+1.623797632095823*coeff[0]*fEdge[22]; 

  boundSurf_incr[1] = 1.5*coeff[0]*fSkin[1]; 
  boundSurf_incr[5] = 1.5*coeff[0]*fSkin[5]; 
  boundSurf_incr[6] = 1.5*coeff[0]*fSkin[6]; 
  boundSurf_incr[8] = 1.5*coeff[0]*fSkin[8]; 
  boundSurf_incr[11] = 1.5*coeff[0]*fSkin[11]; 
  boundSurf_incr[12] = 1.5*coeff[0]*fSkin[12]; 
  boundSurf_incr[13] = 1.5*coeff[0]*fSkin[13]; 
  boundSurf_incr[15] = 1.5*coeff[0]*fSkin[15]; 
  boundSurf_incr[17] = 1.5*coeff[0]*fSkin[17]; 
  boundSurf_incr[20] = 1.5*coeff[0]*fSkin[20]; 
  boundSurf_incr[21] = 1.5*coeff[0]*fSkin[21]; 
  boundSurf_incr[23] = 1.5*coeff[0]*fSkin[23]; 

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

  }

  return 0.;
}

