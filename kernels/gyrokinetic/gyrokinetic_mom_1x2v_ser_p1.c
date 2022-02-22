#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M1_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
} 
GKYL_CU_DH void gyrokinetic_M2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  out[0] += (2.0*(0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M2_par_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M2_perp_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  out[0] += ((0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M3_par_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[2]*dv1*wx1_sq+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[4]*dv1*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M3_perp_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += (volFact*(1.414213562373095*Bmag[1]*f[1]*wx1*wx2+1.414213562373095*Bmag[0]*f[0]*wx1*wx2+0.408248290463863*Bmag[1]*f[4]*dv1*wx2+0.408248290463863*Bmag[0]*f[2]*dv1*wx2+0.408248290463863*Bmag[1]*f[5]*dv2*wx1+0.408248290463863*Bmag[0]*f[3]*dv2*wx1+0.1178511301977579*Bmag[1]*f[7]*dv1*dv2+0.1178511301977579*Bmag[0]*f[6]*dv1*dv2))/m_; 
  out[1] += (volFact*(1.414213562373095*Bmag[0]*f[1]*wx1*wx2+1.414213562373095*f[0]*Bmag[1]*wx1*wx2+0.408248290463863*Bmag[0]*f[4]*dv1*wx2+0.408248290463863*Bmag[1]*f[2]*dv1*wx2+0.408248290463863*Bmag[0]*f[5]*dv2*wx1+0.408248290463863*Bmag[1]*f[3]*dv2*wx1+0.1178511301977579*Bmag[0]*f[7]*dv1*dv2+0.1178511301977579*Bmag[1]*f[6]*dv1*dv2))/m_; 
} 
GKYL_CU_DH void gyrokinetic_ThreeMoments_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT outM0, double* GKYL_RESTRICT outM1, double* GKYL_RESTRICT outM2) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  outM0[0] += 2.0*f[0]*volFact; 
  outM0[1] += 2.0*f[1]*volFact; 
  outM1[0] += 0.3333333333333333*volFact*(6.0*f[0]*wx1+1.732050807568877*f[2]*dv1); 
  outM1[1] += 0.3333333333333333*volFact*(6.0*f[1]*wx1+1.732050807568877*f[4]*dv1); 
  outM2[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  outM2[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  double tmp[2]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2; 
  outM2[0] += (2.0*(0.7071067811865475*Bmag[1]*tmp[1]+0.7071067811865475*Bmag[0]*tmp[0])*volFact)/m_; 
  outM2[1] += (2.0*(0.7071067811865475*Bmag[0]*tmp[1]+0.7071067811865475*tmp[0]*Bmag[1])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M0_step1_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[1]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[3]*volFact; 
  out[3] += 1.414213562373095*f[5]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M0_step2_1x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *Bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
} 
