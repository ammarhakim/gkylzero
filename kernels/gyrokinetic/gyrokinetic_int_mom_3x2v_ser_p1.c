#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_3x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  const double wx1 = w[3], dv1 = dxv[3]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[4], dv2 = dxv[4]; 
 
  double tmp[8]; 
  tmp[0] = (4.0*f[0]*wx2)/m_+(1.154700538379252*f[5]*dv2)/m_; 
  tmp[1] = (4.0*f[1]*wx2)/m_+(1.154700538379252*f[12]*dv2)/m_; 
  tmp[2] = (4.0*f[2]*wx2)/m_+(1.154700538379252*f[13]*dv2)/m_; 
  tmp[3] = (4.0*f[3]*wx2)/m_+(1.154700538379252*f[14]*dv2)/m_; 
  tmp[4] = (4.0*f[6]*wx2)/m_+(1.154700538379252*f[20]*dv2)/m_; 
  tmp[5] = (4.0*f[7]*wx2)/m_+(1.154700538379252*f[21]*dv2)/m_; 
  tmp[6] = (4.0*f[8]*wx2)/m_+(1.154700538379252*f[22]*dv2)/m_; 
  tmp[7] = (4.0*f[16]*wx2)/m_+(1.154700538379252*f[27]*dv2)/m_; 
 
  out[0] += 5.656854249492383*f[0]*volFact; 
  out[1] += volFact*(5.656854249492383*f[0]*wx1+1.632993161855453*f[4]*dv1); 
  out[2] += volFact*(5.656854249492383*f[0]*wx1_sq+3.265986323710906*f[4]*dv1*wx1+0.421637021355784*f[32]*dv1_sq+0.4714045207910317*f[0]*dv1_sq); 
  out[3] += (bmag[7]*tmp[7]+bmag[6]*tmp[6]+bmag[5]*tmp[5]+bmag[4]*tmp[4]+bmag[3]*tmp[3]+bmag[2]*tmp[2]+bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
