#include <gkyl_vlasov_mom_kernels.h>

void vlasov_mom_1x2v_m0_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
} 
void vlasov_mom_1x2v_m0_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
void vlasov_mom_1x2v_m1i_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
  out[2] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2); 
  out[3] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2); 
} 
void vlasov_mom_1x2v_m1i_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
  out[2] += volFact*(2.0*f[7]*wx1+0.5773502691896257*f[11]*dv1); 
  out[3] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2); 
  out[4] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2); 
  out[5] += volFact*(2.0*f[7]*wx2+0.5773502691896257*f[13]*dv2); 
} 
void vlasov_mom_1x2v_m2ij_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[2]*dv1*wx2+0.5773502691896258*f[3]*dv2*wx1+0.1666666666666667*f[6]*dv1*dv2); 
  out[3] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[4]*dv1*wx2+0.5773502691896258*f[5]*dv2*wx1+0.1666666666666667*f[7]*dv1*dv2); 
  out[4] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+0.1666666666666667*f[0]*dv2_sq); 
  out[5] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+0.1666666666666667*f[1]*dv2_sq); 
} 
void vlasov_mom_1x2v_m2ij_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[12]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx1_sq+1.154700538379251*f[11]*dv1*wx1+0.1666666666666667*f[7]*dv1_sq); 
  out[3] += volFact*(2.0*f[0]*wx1*wx2+0.5773502691896258*f[2]*dv1*wx2+0.5773502691896258*f[3]*dv2*wx1+0.1666666666666667*f[6]*dv1*dv2); 
  out[4] += volFact*(2.0*f[1]*wx1*wx2+0.5773502691896258*f[4]*dv1*wx2+0.5773502691896258*f[5]*dv2*wx1+0.1666666666666667*f[10]*dv1*dv2); 
  out[5] += volFact*(2.0*f[7]*wx1*wx2+0.5773502691896257*f[11]*dv1*wx2+0.5773502691896257*f[13]*dv2*wx1+0.1666666666666667*f[17]*dv1*dv2); 
  out[6] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+0.149071198499986*f[9]*dv2_sq+0.1666666666666667*f[0]*dv2_sq); 
  out[7] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+0.149071198499986*f[15]*dv2_sq+0.1666666666666667*f[1]*dv2_sq); 
  out[8] += volFact*(2.0*f[7]*wx2_sq+1.154700538379251*f[13]*dv2*wx2+0.1666666666666667*f[7]*dv2_sq); 
} 
void vlasov_mom_1x2v_m2_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.1666666666666667*f[0]*dv2_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.1666666666666667*f[1]*dv2_sq+0.1666666666666667*f[1]*dv1_sq); 
} 
void vlasov_mom_1x2v_m2_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[9]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[15]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[12]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx2_sq+1.154700538379251*f[13]*dv2*wx2+2.0*f[7]*wx1_sq+1.154700538379251*f[11]*dv1*wx1+0.1666666666666667*f[7]*dv2_sq+0.1666666666666667*f[7]*dv1_sq); 
} 
void vlasov_mom_1x2v_m3i_p1(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.154700538379252*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[2]*dv1*wx1_sq+0.1666666666666667*f[0]*dv2_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04811252243246882*f[2]*dv1*dv2_sq+0.08660254037844387*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.154700538379252*f[5]*dv2*wx1*wx2+0.3333333333333333*f[7]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.1666666666666667*f[1]*dv2_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04811252243246882*f[4]*dv1*dv2_sq+0.08660254037844387*f[4]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[3]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[2]*dv1*wx1*wx2+0.5*f[0]*dv2_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.08660254037844387*f[3]*dv2*dv2_sq+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[3] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[5]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[4]*dv1*wx1*wx2+0.5*f[1]*dv2_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[7]*dv1*dv2*wx1+0.08660254037844387*f[5]*dv2*dv2_sq+0.04811252243246882*f[5]*dv1_sq*dv2); 
} 
void vlasov_mom_1x2v_m3i_p2(const double *w, const double *dxv, const int *idx, const double *restrict f, double* restrict out) 
{ 
  const double volFact = dxv[1]*dxv[2]/4; 
  const double wx1 = w[1], dv1 = dxv[1]; 
  const double wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const double wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const double wx2 = w[2], dv2 = dxv[2]; 
  const double wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const double wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.154700538379252*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[2]*dv1*wx1_sq+0.149071198499986*f[9]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.4472135954999579*f[8]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04303314829119351*f[16]*dv1*dv2_sq+0.04811252243246882*f[2]*dv1*dv2_sq+0.08660254037844387*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.154700538379252*f[5]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.149071198499986*f[15]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.447213595499958*f[12]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04303314829119353*f[19]*dv1*dv2_sq+0.04811252243246882*f[4]*dv1*dv2_sq+0.08660254037844387*f[4]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx1*wx2_sq+0.5773502691896257*f[11]*dv1*wx2_sq+1.154700538379251*f[13]*dv2*wx1*wx2+0.3333333333333333*f[17]*dv1*dv2*wx2+2.0*f[7]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.1666666666666667*f[7]*dv2_sq*wx1+0.5*f[7]*dv1_sq*wx1+0.04811252243246881*f[11]*dv1*dv2_sq+0.08660254037844385*f[11]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[3]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[2]*dv1*wx1*wx2+0.4472135954999579*f[9]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.149071198499986*f[8]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.08660254037844387*f[3]*dv2*dv2_sq+0.04303314829119351*f[14]*dv1_sq*dv2+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[4] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[5]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[4]*dv1*wx1*wx2+0.447213595499958*f[15]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.149071198499986*f[12]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.08660254037844387*f[5]*dv2*dv2_sq+0.04303314829119353*f[18]*dv1_sq*dv2+0.04811252243246882*f[5]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[7]*wx2*wx2_sq+1.732050807568877*f[13]*dv2*wx2_sq+2.0*f[7]*wx1_sq*wx2+1.154700538379251*f[11]*dv1*wx1*wx2+0.5*f[7]*dv2_sq*wx2+0.1666666666666667*f[7]*dv1_sq*wx2+0.5773502691896257*f[13]*dv2*wx1_sq+0.3333333333333333*f[17]*dv1*dv2*wx1+0.08660254037844385*f[13]*dv2*dv2_sq+0.04811252243246881*f[13]*dv1_sq*dv2); 
} 
