#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_int_mom_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[0]*dxv[1]*dxv[2]*0.125; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
 
  double tmp[2]; 
  tmp[0] = (4.0*f[0]*wx2)/m_+(1.154700538379252*f[3]*dv2)/m_; 
  tmp[1] = (4.0*f[1]*wx2)/m_+(1.154700538379252*f[5]*dv2)/m_; 
 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += volFact*(2.828427124746191*f[0]*wx1+0.8164965809277261*f[2]*dv1); 
  out[2] += volFact*(2.828427124746191*f[0]*wx1_sq+1.632993161855453*f[2]*dv1*wx1+0.210818510677892*f[8]*dv1_sq+0.2357022603955158*f[0]*dv1_sq); 
  out[3] += (bmag[1]*tmp[1]+bmag[0]*tmp[0])*volFact; 
} 
