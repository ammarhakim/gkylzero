#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *jacobgeo_inv, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[0]*dxv[1]*dxv[2]*dxv[3]/0.0625; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
 
  double J_inv_f[24] = {0.0}; 
  J_inv_f[0] = 0.5*(jacobgeo_inv[3]*f[5]+f[2]*jacobgeo_inv[2]+f[1]*jacobgeo_inv[1]+f[0]*jacobgeo_inv[0]); 
  J_inv_f[1] = 0.5*(jacobgeo_inv[2]*f[5]+f[2]*jacobgeo_inv[3]+f[0]*jacobgeo_inv[1]+jacobgeo_inv[0]*f[1]); 
  J_inv_f[2] = 0.5*(jacobgeo_inv[1]*f[5]+f[1]*jacobgeo_inv[3]+f[0]*jacobgeo_inv[2]+jacobgeo_inv[0]*f[2]); 
  J_inv_f[3] = 0.5*(jacobgeo_inv[3]*f[11]+jacobgeo_inv[2]*f[7]+jacobgeo_inv[1]*f[6]+jacobgeo_inv[0]*f[3]); 
  J_inv_f[4] = 0.5*(jacobgeo_inv[3]*f[12]+jacobgeo_inv[2]*f[9]+jacobgeo_inv[1]*f[8]+jacobgeo_inv[0]*f[4]); 
  J_inv_f[5] = 0.5*(jacobgeo_inv[0]*f[5]+f[0]*jacobgeo_inv[3]+f[1]*jacobgeo_inv[2]+jacobgeo_inv[1]*f[2]); 
  J_inv_f[6] = 0.5*(jacobgeo_inv[2]*f[11]+jacobgeo_inv[3]*f[7]+jacobgeo_inv[0]*f[6]+jacobgeo_inv[1]*f[3]); 
  J_inv_f[7] = 0.5*(jacobgeo_inv[1]*f[11]+jacobgeo_inv[0]*f[7]+jacobgeo_inv[3]*f[6]+jacobgeo_inv[2]*f[3]); 
  J_inv_f[8] = 0.5*(jacobgeo_inv[2]*f[12]+jacobgeo_inv[3]*f[9]+jacobgeo_inv[0]*f[8]+jacobgeo_inv[1]*f[4]); 
  J_inv_f[9] = 0.5*(jacobgeo_inv[1]*f[12]+jacobgeo_inv[0]*f[9]+jacobgeo_inv[3]*f[8]+jacobgeo_inv[2]*f[4]); 
  J_inv_f[10] = 0.5*(jacobgeo_inv[3]*f[15]+jacobgeo_inv[2]*f[14]+jacobgeo_inv[1]*f[13]+jacobgeo_inv[0]*f[10]); 
  J_inv_f[11] = 0.5*(jacobgeo_inv[0]*f[11]+jacobgeo_inv[1]*f[7]+jacobgeo_inv[2]*f[6]+f[3]*jacobgeo_inv[3]); 
  J_inv_f[12] = 0.5*(jacobgeo_inv[0]*f[12]+jacobgeo_inv[1]*f[9]+jacobgeo_inv[2]*f[8]+jacobgeo_inv[3]*f[4]); 
  J_inv_f[13] = 0.5*(jacobgeo_inv[2]*f[15]+jacobgeo_inv[3]*f[14]+jacobgeo_inv[0]*f[13]+jacobgeo_inv[1]*f[10]); 
  J_inv_f[14] = 0.5*(jacobgeo_inv[1]*f[15]+jacobgeo_inv[0]*f[14]+jacobgeo_inv[3]*f[13]+jacobgeo_inv[2]*f[10]); 
  J_inv_f[15] = 0.5*(jacobgeo_inv[0]*f[15]+jacobgeo_inv[1]*f[14]+jacobgeo_inv[2]*f[13]+jacobgeo_inv[3]*f[10]); 
  J_inv_f[16] = 0.5*jacobgeo_inv[3]*f[20]+0.5000000000000001*(jacobgeo_inv[2]*f[18]+jacobgeo_inv[1]*f[17])+0.5*jacobgeo_inv[0]*f[16]; 
  J_inv_f[17] = 0.5000000000000001*jacobgeo_inv[2]*f[20]+0.5*(jacobgeo_inv[3]*f[18]+jacobgeo_inv[0]*f[17])+0.5000000000000001*jacobgeo_inv[1]*f[16]; 
  J_inv_f[18] = 0.5000000000000001*jacobgeo_inv[1]*f[20]+0.5*(jacobgeo_inv[0]*f[18]+jacobgeo_inv[3]*f[17])+0.5000000000000001*jacobgeo_inv[2]*f[16]; 
  J_inv_f[19] = 0.5*jacobgeo_inv[3]*f[23]+0.5000000000000001*(jacobgeo_inv[2]*f[22]+jacobgeo_inv[1]*f[21])+0.5*jacobgeo_inv[0]*f[19]; 
  J_inv_f[20] = 0.5*jacobgeo_inv[0]*f[20]+0.5000000000000001*(jacobgeo_inv[1]*f[18]+jacobgeo_inv[2]*f[17])+0.5*jacobgeo_inv[3]*f[16]; 
  J_inv_f[21] = 0.5000000000000001*jacobgeo_inv[2]*f[23]+0.5*(jacobgeo_inv[3]*f[22]+jacobgeo_inv[0]*f[21])+0.5000000000000001*jacobgeo_inv[1]*f[19]; 
  J_inv_f[22] = 0.5000000000000001*jacobgeo_inv[1]*f[23]+0.5*(jacobgeo_inv[0]*f[22]+jacobgeo_inv[3]*f[21])+0.5000000000000001*jacobgeo_inv[2]*f[19]; 
  J_inv_f[23] = 0.5*jacobgeo_inv[0]*f[23]+0.5000000000000001*(jacobgeo_inv[1]*f[22]+jacobgeo_inv[2]*f[21])+0.5*jacobgeo_inv[3]*f[19]; 
  double tmp[4]; 
  tmp[0] = (4.0*J_inv_f[0]*wx2)/m_+(1.154700538379252*J_inv_f[4]*dv2)/m_; 
  tmp[1] = (4.0*J_inv_f[1]*wx2)/m_+(1.154700538379252*J_inv_f[8]*dv2)/m_; 
  tmp[2] = (4.0*J_inv_f[2]*wx2)/m_+(1.154700538379252*J_inv_f[9]*dv2)/m_; 
  tmp[3] = (4.0*J_inv_f[5]*wx2)/m_+(1.154700538379252*J_inv_f[12]*dv2)/m_; 
 
  out[0] += 4.0*J_inv_f[0]*volFact; 
  out[1] += volFact*(4.0*J_inv_f[0]*wx1+1.154700538379252*J_inv_f[3]*dv1); 
  out[2] += volFact*(4.0*J_inv_f[0]*wx1_sq+2.309401076758503*J_inv_f[3]*dv1*wx1+0.2981423969999719*J_inv_f[16]*dv1_sq+0.3333333333333333*J_inv_f[0]*dv1_sq); 
  out[3] += (bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
