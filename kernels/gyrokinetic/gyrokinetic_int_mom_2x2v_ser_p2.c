#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_2x2v_ser_p2(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  double tmp[8]; 
  tmp[0] = (4.0*f[0]*wx2)/m_+(1.154700538379252*f[4]*dv2)/m_; 
  tmp[1] = (4.0*f[1]*wx2)/m_+(1.154700538379252*f[8]*dv2)/m_; 
  tmp[2] = (4.0*f[2]*wx2)/m_+(1.154700538379252*f[9]*dv2)/m_; 
  tmp[3] = (4.0*f[5]*wx2)/m_+(1.154700538379252*f[16]*dv2)/m_; 
  tmp[4] = (4.0*f[11]*wx2)/m_+(1.154700538379251*f[25]*dv2)/m_; 
  tmp[5] = (4.0*f[12]*wx2)/m_+(1.154700538379251*f[26]*dv2)/m_; 
  tmp[6] = (4.0*f[19]*wx2)/m_+(1.154700538379251*f[35]*dv2)/m_; 
  tmp[7] = (4.0*f[20]*wx2)/m_+(1.154700538379251*f[36]*dv2)/m_; 
 
  out[0] += 4.0*f[0]*volFact; 
  out[1] += volFact*(4.0*f[0]*wx1+1.154700538379252*f[3]*dv1); 
  out[2] += volFact*(4.0*f[0]*wx1_sq+2.309401076758503*f[3]*dv1*wx1+0.2981423969999719*f[13]*dv1_sq+0.3333333333333333*f[0]*dv1_sq); 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
