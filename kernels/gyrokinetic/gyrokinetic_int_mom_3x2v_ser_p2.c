#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_3x2v_ser_p2(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]*dxv[4]/4; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
 
  double tmp[20]; 
  tmp[0] = (4.0*f[0]*wx2)/m_+(1.154700538379252*f[5]*dv2)/m_; 
  tmp[1] = (4.0*f[1]*wx2)/m_+(1.154700538379252*f[12]*dv2)/m_; 
  tmp[2] = (4.0*f[2]*wx2)/m_+(1.154700538379252*f[13]*dv2)/m_; 
  tmp[3] = (4.0*f[3]*wx2)/m_+(1.154700538379252*f[14]*dv2)/m_; 
  tmp[4] = (4.0*f[6]*wx2)/m_+(1.154700538379252*f[25]*dv2)/m_; 
  tmp[5] = (4.0*f[7]*wx2)/m_+(1.154700538379252*f[26]*dv2)/m_; 
  tmp[6] = (4.0*f[8]*wx2)/m_+(1.154700538379252*f[27]*dv2)/m_; 
  tmp[7] = (4.0*f[16]*wx2)/m_+(1.154700538379251*f[43]*dv2)/m_; 
  tmp[8] = (4.0*f[17]*wx2)/m_+(1.154700538379251*f[44]*dv2)/m_; 
  tmp[9] = (4.0*f[18]*wx2)/m_+(1.154700538379251*f[45]*dv2)/m_; 
  tmp[10] = (4.0*f[21]*wx2)/m_+(1.154700538379252*f[52]*dv2)/m_; 
  tmp[11] = (4.0*f[31]*wx2)/m_+(1.154700538379251*f[68]*dv2)/m_; 
  tmp[12] = (4.0*f[32]*wx2)/m_+(1.154700538379251*f[69]*dv2)/m_; 
  tmp[13] = (4.0*f[33]*wx2)/m_+(1.154700538379251*f[70]*dv2)/m_; 
  tmp[14] = (4.0*f[34]*wx2)/m_+(1.154700538379251*f[71]*dv2)/m_; 
  tmp[15] = (4.0*f[35]*wx2)/m_+(1.154700538379251*f[72]*dv2)/m_; 
  tmp[16] = (4.0*f[36]*wx2)/m_+(1.154700538379251*f[73]*dv2)/m_; 
  tmp[17] = (4.0*f[56]*wx2)/m_+(1.154700538379251*f[91]*dv2)/m_; 
  tmp[18] = (4.0*f[57]*wx2)/m_+(1.154700538379251*f[92]*dv2)/m_; 
  tmp[19] = (4.0*f[58]*wx2)/m_+(1.154700538379251*f[93]*dv2)/m_; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += volFact*(5.656854249492382*f[0]*wx1+1.632993161855453*f[4]*dv1); 
  out[2] += volFact*(5.656854249492382*f[0]*wx1_sq+3.265986323710906*f[4]*dv1*wx1+0.421637021355784*f[19]*dv1_sq+0.4714045207910317*f[0]*dv1_sq); 
  out[3] += (bmag[19]*tmp[19]+bmag[18]*tmp[18]+bmag[17]*tmp[17]+bmag[16]*tmp[16]+bmag[15]*tmp[15]+bmag[14]*tmp[14]+bmag[13]*tmp[13]+bmag[12]*tmp[12]+bmag[11]*tmp[11]+bmag[10]*tmp[10]+bmag[9]*tmp[9]+bmag[8]*tmp[8]+bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
