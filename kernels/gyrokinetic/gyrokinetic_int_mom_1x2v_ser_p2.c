#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_1x2v_ser_p2(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[0]*dxv[1]*dxv[2]/0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
 
  double J_inv_f[20] = {0.0}; 
  J_inv_f[0] = 0.7071067811865475*(jacobgeo_inv[2]*f[7]+f[1]*jacobgeo_inv[1]+f[0]*jacobgeo_inv[0]); 
  J_inv_f[1] = 0.6324555320336759*(jacobgeo_inv[1]*f[7]+f[1]*jacobgeo_inv[2])+0.7071067811865475*(f[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*f[1]); 
  J_inv_f[2] = 0.7071067811865475*(jacobgeo_inv[2]*f[11]+jacobgeo_inv[1]*f[4]+jacobgeo_inv[0]*f[2]); 
  J_inv_f[3] = 0.7071067811865475*(jacobgeo_inv[2]*f[13]+jacobgeo_inv[1]*f[5]+jacobgeo_inv[0]*f[3]); 
  J_inv_f[4] = 0.632455532033676*jacobgeo_inv[1]*f[11]+0.6324555320336759*jacobgeo_inv[2]*f[4]+0.7071067811865475*(jacobgeo_inv[0]*f[4]+jacobgeo_inv[1]*f[2]); 
  J_inv_f[5] = 0.632455532033676*jacobgeo_inv[1]*f[13]+0.6324555320336759*jacobgeo_inv[2]*f[5]+0.7071067811865475*(jacobgeo_inv[0]*f[5]+jacobgeo_inv[1]*f[3]); 
  J_inv_f[6] = 0.7071067811865475*(jacobgeo_inv[2]*f[17]+jacobgeo_inv[1]*f[10]+jacobgeo_inv[0]*f[6]); 
  J_inv_f[7] = 0.4517539514526256*jacobgeo_inv[2]*f[7]+0.7071067811865475*(jacobgeo_inv[0]*f[7]+f[0]*jacobgeo_inv[2])+0.6324555320336759*f[1]*jacobgeo_inv[1]; 
  J_inv_f[8] = 0.7071067811865475*(jacobgeo_inv[1]*f[12]+jacobgeo_inv[0]*f[8]); 
  J_inv_f[9] = 0.7071067811865475*(jacobgeo_inv[1]*f[15]+jacobgeo_inv[0]*f[9]); 
  J_inv_f[10] = 0.6324555320336759*(jacobgeo_inv[1]*f[17]+jacobgeo_inv[2]*f[10])+0.7071067811865475*(jacobgeo_inv[0]*f[10]+jacobgeo_inv[1]*f[6]); 
  J_inv_f[11] = (0.4517539514526256*jacobgeo_inv[2]+0.7071067811865475*jacobgeo_inv[0])*f[11]+0.632455532033676*jacobgeo_inv[1]*f[4]+0.7071067811865475*f[2]*jacobgeo_inv[2]; 
  J_inv_f[12] = 0.6324555320336759*jacobgeo_inv[2]*f[12]+0.7071067811865475*(jacobgeo_inv[0]*f[12]+jacobgeo_inv[1]*f[8]); 
  J_inv_f[13] = (0.4517539514526256*jacobgeo_inv[2]+0.7071067811865475*jacobgeo_inv[0])*f[13]+0.632455532033676*jacobgeo_inv[1]*f[5]+0.7071067811865475*jacobgeo_inv[2]*f[3]; 
  J_inv_f[14] = 0.7071067811865475*(jacobgeo_inv[1]*f[18]+jacobgeo_inv[0]*f[14]); 
  J_inv_f[15] = 0.6324555320336759*jacobgeo_inv[2]*f[15]+0.7071067811865475*(jacobgeo_inv[0]*f[15]+jacobgeo_inv[1]*f[9]); 
  J_inv_f[16] = 0.7071067811865475*(jacobgeo_inv[1]*f[19]+jacobgeo_inv[0]*f[16]); 
  J_inv_f[17] = (0.4517539514526256*jacobgeo_inv[2]+0.7071067811865475*jacobgeo_inv[0])*f[17]+0.6324555320336759*jacobgeo_inv[1]*f[10]+0.7071067811865475*jacobgeo_inv[2]*f[6]; 
  J_inv_f[18] = 0.6324555320336759*jacobgeo_inv[2]*f[18]+0.7071067811865475*(jacobgeo_inv[0]*f[18]+jacobgeo_inv[1]*f[14]); 
  J_inv_f[19] = 0.6324555320336759*jacobgeo_inv[2]*f[19]+0.7071067811865475*(jacobgeo_inv[0]*f[19]+jacobgeo_inv[1]*f[16]); 
  double tmp[3]; 
  tmp[0] = (4.0*J_inv_f[0]*wx2)/m_+(1.154700538379252*J_inv_f[3]*dv2)/m_; 
  tmp[1] = (4.0*J_inv_f[1]*wx2)/m_+(1.154700538379252*J_inv_f[5]*dv2)/m_; 
  tmp[2] = (4.0*J_inv_f[7]*wx2)/m_+(1.154700538379251*J_inv_f[13]*dv2)/m_; 
 
  out[0] += 2.828427124746191*J_inv_f[0]*volFact; 
  out[1] += volFact*(2.828427124746191*J_inv_f[0]*wx1+0.8164965809277261*J_inv_f[2]*dv1); 
  out[2] += volFact*(2.828427124746191*J_inv_f[0]*wx1_sq+1.632993161855453*J_inv_f[2]*dv1*wx1+0.210818510677892*J_inv_f[8]*dv1_sq+0.2357022603955158*J_inv_f[0]*dv1_sq); 
  out[3] += (bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
