#include <gkyl_mom_gyrokinetic_kernels.h> 
GKYL_CU_DH void gyrokinetic_M0_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M1_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[2] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[3] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
} 
GKYL_CU_DH void gyrokinetic_M2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  out[0] += (2.0*(0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += (2.0*(0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += (2.0*(0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += (2.0*(0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M2_par_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M2_perp_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  out[0] += ((0.5*bmag[3]*tmp[3]+0.5*bmag[2]*tmp[2]+0.5*bmag[1]*tmp[1]+0.5*bmag[0]*tmp[0])*volFact)/m_; 
  out[1] += ((0.5*bmag[2]*tmp[3]+0.5*tmp[2]*bmag[3]+0.5*bmag[0]*tmp[1]+0.5*tmp[0]*bmag[1])*volFact)/m_; 
  out[2] += ((0.5*bmag[1]*tmp[3]+0.5*tmp[1]*bmag[3]+0.5*bmag[0]*tmp[2]+0.5*tmp[0]*bmag[2])*volFact)/m_; 
  out[3] += ((0.5*bmag[0]*tmp[3]+0.5*tmp[0]*bmag[3]+0.5*bmag[1]*tmp[2]+0.5*tmp[1]*bmag[2])*volFact)/m_; 
} 
GKYL_CU_DH void gyrokinetic_M3_par_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[3]*dv1*wx1_sq+0.4472135954999579*f[16]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[3]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[6]*dv1*wx1_sq+0.447213595499958*f[17]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[6]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[2]*wx1*wx1_sq+1.732050807568877*f[7]*dv1*wx1_sq+0.447213595499958*f[18]*dv1_sq*wx1+0.5*f[2]*dv1_sq*wx1+0.08660254037844387*f[7]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[5]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.4472135954999579*f[20]*dv1_sq*wx1+0.5*f[5]*dv1_sq*wx1+0.08660254037844387*f[11]*dv1*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M3_perp_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  out[0] += (volFact*(bmag[3]*f[5]*wx1*wx2+bmag[2]*f[2]*wx1*wx2+bmag[1]*f[1]*wx1*wx2+bmag[0]*f[0]*wx1*wx2+0.2886751345948129*bmag[3]*f[11]*dv1*wx2+0.2886751345948129*bmag[2]*f[7]*dv1*wx2+0.2886751345948129*bmag[1]*f[6]*dv1*wx2+0.2886751345948129*bmag[0]*f[3]*dv1*wx2+0.2886751345948129*bmag[3]*f[12]*dv2*wx1+0.2886751345948129*bmag[2]*f[9]*dv2*wx1+0.2886751345948129*bmag[1]*f[8]*dv2*wx1+0.2886751345948129*bmag[0]*f[4]*dv2*wx1+0.08333333333333333*bmag[3]*f[15]*dv1*dv2+0.08333333333333333*bmag[2]*f[14]*dv1*dv2+0.08333333333333333*bmag[1]*f[13]*dv1*dv2+0.08333333333333333*bmag[0]*f[10]*dv1*dv2))/m_; 
  out[1] += (volFact*(bmag[2]*f[5]*wx1*wx2+f[2]*bmag[3]*wx1*wx2+bmag[0]*f[1]*wx1*wx2+f[0]*bmag[1]*wx1*wx2+0.2886751345948129*bmag[2]*f[11]*dv1*wx2+0.2886751345948129*bmag[3]*f[7]*dv1*wx2+0.2886751345948129*bmag[0]*f[6]*dv1*wx2+0.2886751345948129*bmag[1]*f[3]*dv1*wx2+0.2886751345948129*bmag[2]*f[12]*dv2*wx1+0.2886751345948129*bmag[3]*f[9]*dv2*wx1+0.2886751345948129*bmag[0]*f[8]*dv2*wx1+0.2886751345948129*bmag[1]*f[4]*dv2*wx1+0.08333333333333333*bmag[2]*f[15]*dv1*dv2+0.08333333333333333*bmag[3]*f[14]*dv1*dv2+0.08333333333333333*bmag[0]*f[13]*dv1*dv2+0.08333333333333333*bmag[1]*f[10]*dv1*dv2))/m_; 
  out[2] += (volFact*(bmag[1]*f[5]*wx1*wx2+f[1]*bmag[3]*wx1*wx2+bmag[0]*f[2]*wx1*wx2+f[0]*bmag[2]*wx1*wx2+0.2886751345948129*bmag[1]*f[11]*dv1*wx2+0.2886751345948129*bmag[0]*f[7]*dv1*wx2+0.2886751345948129*bmag[3]*f[6]*dv1*wx2+0.2886751345948129*bmag[2]*f[3]*dv1*wx2+0.2886751345948129*bmag[1]*f[12]*dv2*wx1+0.2886751345948129*bmag[0]*f[9]*dv2*wx1+0.2886751345948129*bmag[3]*f[8]*dv2*wx1+0.2886751345948129*bmag[2]*f[4]*dv2*wx1+0.08333333333333333*bmag[1]*f[15]*dv1*dv2+0.08333333333333333*bmag[0]*f[14]*dv1*dv2+0.08333333333333333*bmag[3]*f[13]*dv1*dv2+0.08333333333333333*bmag[2]*f[10]*dv1*dv2))/m_; 
  out[3] += (volFact*(bmag[0]*f[5]*wx1*wx2+f[0]*bmag[3]*wx1*wx2+bmag[1]*f[2]*wx1*wx2+f[1]*bmag[2]*wx1*wx2+0.2886751345948129*bmag[0]*f[11]*dv1*wx2+0.2886751345948129*bmag[1]*f[7]*dv1*wx2+0.2886751345948129*bmag[2]*f[6]*dv1*wx2+0.2886751345948129*bmag[3]*f[3]*dv1*wx2+0.2886751345948129*bmag[0]*f[12]*dv2*wx1+0.2886751345948129*bmag[1]*f[9]*dv2*wx1+0.2886751345948129*bmag[2]*f[8]*dv2*wx1+0.2886751345948129*bmag[3]*f[4]*dv2*wx1+0.08333333333333333*bmag[0]*f[15]*dv1*dv2+0.08333333333333333*bmag[1]*f[14]*dv1*dv2+0.08333333333333333*bmag[2]*f[13]*dv1*dv2+0.08333333333333333*bmag[3]*f[10]*dv1*dv2))/m_; 
} 
GKYL_CU_DH void gyrokinetic_three_moments_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[2]*dxv[3]/4; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  double tmp[4]; 
  tmp[0] = 2.0*f[0]*wx2+0.5773502691896258*f[4]*dv2; 
  tmp[1] = 2.0*f[1]*wx2+0.5773502691896258*f[8]*dv2; 
  tmp[2] = 2.0*f[2]*wx2+0.5773502691896258*f[9]*dv2; 
  tmp[3] = 2.0*f[5]*wx2+0.5773502691896258*f[12]*dv2; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[2]*volFact; 
  out[3] += 2.0*f[5]*volFact; 
  out[4] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[3]*dv1); 
  out[5] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[6]*dv1); 
  out[6] += volFact*(2.0*f[2]*wx1+0.5773502691896258*f[7]*dv1); 
  out[7] += volFact*(2.0*f[5]*wx1+0.5773502691896258*f[11]*dv1); 
  out[8] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[3]*dv1*wx1+(bmag[3]*tmp[3])/m_+(bmag[2]*tmp[2])/m_+(bmag[1]*tmp[1])/m_+(bmag[0]*tmp[0])/m_+0.149071198499986*f[16]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[9] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[6]*dv1*wx1+(bmag[2]*tmp[3])/m_+(tmp[2]*bmag[3])/m_+(bmag[0]*tmp[1])/m_+(tmp[0]*bmag[1])/m_+0.149071198499986*f[17]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[10] += volFact*(2.0*f[2]*wx1_sq+1.154700538379252*f[7]*dv1*wx1+(bmag[1]*tmp[3])/m_+(tmp[1]*bmag[3])/m_+(bmag[0]*tmp[2])/m_+(tmp[0]*bmag[2])/m_+0.149071198499986*f[18]*dv1_sq+0.1666666666666667*f[2]*dv1_sq); 
  out[11] += volFact*(2.0*f[5]*wx1_sq+1.154700538379252*f[11]*dv1*wx1+(bmag[0]*tmp[3])/m_+(tmp[0]*bmag[3])/m_+(bmag[1]*tmp[2])/m_+(tmp[1]*bmag[2])/m_+0.149071198499986*f[20]*dv1_sq+0.1666666666666667*f[5]*dv1_sq); 
} 
GKYL_CU_DH void gyrokinetic_M0_step1_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[2]/2; 
  out[0] += 1.414213562373095*f[0]*volFact; 
  out[1] += 1.414213562373095*f[1]*volFact; 
  out[2] += 1.414213562373095*f[2]*volFact; 
  out[3] += 1.414213562373095*f[4]*volFact; 
  out[4] += 1.414213562373095*f[5]*volFact; 
  out[5] += 1.414213562373095*f[8]*volFact; 
  out[6] += 1.414213562373095*f[9]*volFact; 
  out[7] += 1.414213562373095*f[12]*volFact; 
} 
GKYL_CU_DH void gyrokinetic_M0_step2_2x2v_ser_p1(const double *w, const double *dxv, const int *idx, double m_, const double *bmag, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = 2.0*M_PI/m_*dxv[3]/2; 
  out[0] += 2.828427124746191*f[0]*volFact; 
  out[1] += 2.828427124746191*f[1]*volFact; 
  out[2] += 2.828427124746191*f[2]*volFact; 
  out[3] += 2.828427124746191*f[4]*volFact; 
} 
