#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_1x1v_ser_p2(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*0.25; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
 
  double J_inv_f[8] = {0.0}; 
  J_inv_f[0] = 0.7071067811865475*(jacobgeo_inv[2]*f[4]+f[1]*jacobgeo_inv[1]+f[0]*jacobgeo_inv[0]); 
  J_inv_f[1] = 0.6324555320336759*(jacobgeo_inv[1]*f[4]+f[1]*jacobgeo_inv[2])+0.7071067811865475*(f[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*f[1]); 
  J_inv_f[2] = 0.7071067811865475*(jacobgeo_inv[2]*f[6]+jacobgeo_inv[1]*f[3]+jacobgeo_inv[0]*f[2]); 
  J_inv_f[3] = 0.632455532033676*jacobgeo_inv[1]*f[6]+0.6324555320336759*jacobgeo_inv[2]*f[3]+0.7071067811865475*(jacobgeo_inv[0]*f[3]+jacobgeo_inv[1]*f[2]); 
  J_inv_f[4] = 0.4517539514526256*jacobgeo_inv[2]*f[4]+0.7071067811865475*(jacobgeo_inv[0]*f[4]+f[0]*jacobgeo_inv[2])+0.6324555320336759*f[1]*jacobgeo_inv[1]; 
  J_inv_f[5] = 0.7071067811865475*(jacobgeo_inv[1]*f[7]+jacobgeo_inv[0]*f[5]); 
  J_inv_f[6] = (0.4517539514526256*jacobgeo_inv[2]+0.7071067811865475*jacobgeo_inv[0])*f[6]+0.632455532033676*jacobgeo_inv[1]*f[3]+0.7071067811865475*f[2]*jacobgeo_inv[2]; 
  J_inv_f[7] = 0.6324555320336759*jacobgeo_inv[2]*f[7]+0.7071067811865475*(jacobgeo_inv[0]*f[7]+jacobgeo_inv[1]*f[5]); 
  out[0] += 2.0*J_inv_f[0]*volFact; 
  out[1] += volFact*(2.0*J_inv_f[0]*wx1+0.5773502691896258*J_inv_f[2]*dv1); 
  out[2] += volFact*(2.0*J_inv_f[0]*wx1_sq+1.154700538379252*J_inv_f[2]*dv1*wx1+0.149071198499986*J_inv_f[5]*dv1_sq+0.1666666666666667*J_inv_f[0]*dv1_sq); 
} 
