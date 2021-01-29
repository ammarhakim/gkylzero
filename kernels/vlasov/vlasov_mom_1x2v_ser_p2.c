#include <gkyl_vlasov_mom_kernels.h> 
void vlasov_M0_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  out[0] += 2.0*f[0]*volFact; 
  out[1] += 2.0*f[1]*volFact; 
  out[2] += 2.0*f[7]*volFact; 
} 
void vlasov_M1i_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  out[0] += volFact*(2.0*f[0]*wx1+0.5773502691896258*f[2]*dv1); 
  out[1] += volFact*(2.0*f[1]*wx1+0.5773502691896258*f[4]*dv1); 
  out[2] += volFact*(2.0*f[7]*wx1+0.5773502691896257*f[11]*dv1); 
  out[3] += volFact*(2.0*f[0]*wx2+0.5773502691896258*f[3]*dv2); 
  out[4] += volFact*(2.0*f[1]*wx2+0.5773502691896258*f[5]*dv2); 
  out[5] += volFact*(2.0*f[7]*wx2+0.5773502691896257*f[13]*dv2); 
} 
void vlasov_M2_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx2_sq+1.154700538379252*f[3]*dv2*wx2+2.0*f[0]*wx1_sq+1.154700538379252*f[2]*dv1*wx1+0.149071198499986*f[9]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx2_sq+1.154700538379252*f[5]*dv2*wx2+2.0*f[1]*wx1_sq+1.154700538379252*f[4]*dv1*wx1+0.149071198499986*f[15]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[12]*dv1_sq+0.1666666666666667*f[1]*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx2_sq+1.154700538379251*f[13]*dv2*wx2+2.0*f[7]*wx1_sq+1.154700538379251*f[11]*dv1*wx1+0.1666666666666667*f[7]*dv2_sq+0.1666666666666667*f[7]*dv1_sq); 
} 
void vlasov_FiveMoments_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict outM0, gkyl_real* restrict outM1i, gkyl_real* restrict outM2) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  gkyl_real tempM0[3], tempM1i[6]; 

  tempM0[0] = 2.0*f[0]*volFact; 
  tempM0[1] = 2.0*f[1]*volFact; 
  tempM0[2] = 2.0*f[7]*volFact; 

  tempM1i[0] = tempM0[0]*wx1+0.5773502691896258*f[2]*dv1*volFact; 
  tempM1i[1] = tempM0[1]*wx1+0.5773502691896258*f[4]*dv1*volFact; 
  tempM1i[2] = tempM0[2]*wx1+0.5773502691896257*f[11]*dv1*volFact; 
  tempM1i[3] = tempM0[0]*wx2+0.5773502691896258*f[3]*dv2*volFact; 
  tempM1i[4] = tempM0[1]*wx2+0.5773502691896258*f[5]*dv2*volFact; 
  tempM1i[5] = tempM0[2]*wx2+0.5773502691896257*f[13]*dv2*volFact; 

  outM0[0] += tempM0[0]; 
  outM0[1] += tempM0[1]; 
  outM0[2] += tempM0[2]; 
  outM1i[0] += tempM1i[0]; 
  outM1i[1] += tempM1i[1]; 
  outM1i[2] += tempM1i[2]; 
  outM1i[3] += tempM1i[3]; 
  outM1i[4] += tempM1i[4]; 
  outM1i[5] += tempM1i[5]; 
  outM2[0] += tempM0[0]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[3]*wx2+2.0*tempM1i[0]*wx1+(0.149071198499986*f[9]*dv2_sq+0.1666666666666667*f[0]*dv2_sq+0.149071198499986*f[8]*dv1_sq+0.1666666666666667*f[0]*dv1_sq)*volFact; 
  outM2[1] += tempM0[1]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[4]*wx2+2.0*tempM1i[1]*wx1+(0.149071198499986*f[15]*dv2_sq+0.1666666666666667*f[1]*dv2_sq+0.149071198499986*f[12]*dv1_sq+0.1666666666666667*f[1]*dv1_sq)*volFact; 
  outM2[2] += tempM0[2]*((-1.0*wx2_sq)-1.0*wx1_sq)+2.0*tempM1i[5]*wx2+2.0*tempM1i[2]*wx1+(0.1666666666666667*f[7]*dv2_sq+0.1666666666666667*f[7]*dv1_sq)*volFact; 
} 
void vlasov_M2ij_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
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
void vlasov_M3i_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const gkyl_real wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.154700538379252*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[2]*dv1*wx1_sq+0.149071198499986*f[9]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.4472135954999579*f[8]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.04303314829119351*f[16]*dv1*dv2_sq+0.04811252243246882*f[2]*dv1*dv2_sq+0.08660254037844387*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.154700538379252*f[5]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.149071198499986*f[15]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.447213595499958*f[12]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.04303314829119353*f[19]*dv1*dv2_sq+0.04811252243246882*f[4]*dv1*dv2_sq+0.08660254037844387*f[4]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx1*wx2_sq+0.5773502691896257*f[11]*dv1*wx2_sq+1.154700538379251*f[13]*dv2*wx1*wx2+0.3333333333333333*f[17]*dv1*dv2*wx2+2.0*f[7]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.1666666666666667*f[7]*dv2_sq*wx1+0.5*f[7]*dv1_sq*wx1+0.04811252243246881*f[11]*dv1*dv2_sq+0.08660254037844385*f[11]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[3]*dv2*wx2_sq+2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[2]*dv1*wx1*wx2+0.4472135954999579*f[9]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.149071198499986*f[8]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.08660254037844387*f[3]*dv2*dv2_sq+0.04303314829119351*f[14]*dv1_sq*dv2+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[4] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[5]*dv2*wx2_sq+2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[4]*dv1*wx1*wx2+0.447213595499958*f[15]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.149071198499986*f[12]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.08660254037844387*f[5]*dv2*dv2_sq+0.04303314829119353*f[18]*dv1_sq*dv2+0.04811252243246882*f[5]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[7]*wx2*wx2_sq+1.732050807568877*f[13]*dv2*wx2_sq+2.0*f[7]*wx1_sq*wx2+1.154700538379251*f[11]*dv1*wx1*wx2+0.5*f[7]*dv2_sq*wx2+0.1666666666666667*f[7]*dv1_sq*wx2+0.5773502691896257*f[13]*dv2*wx1_sq+0.3333333333333333*f[17]*dv1*dv2*wx1+0.08660254037844385*f[13]*dv2*dv2_sq+0.04811252243246881*f[13]*dv1_sq*dv2); 
} 
void vlasov_M3ijk_1x2v_ser_p2(const gkyl_real *w, const gkyl_real *dxv, const int *idx, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  const gkyl_real volFact = dxv[1]*dxv[2]/4; 
  const gkyl_real wx1 = w[1], dv1 = dxv[1]; 
  const gkyl_real wx1_sq = wx1*wx1, dv1_sq = dv1*dv1; 
  const gkyl_real wx1_cu = wx1*wx1*wx1, dv1_cu = dv1*dv1*dv1; 
  const gkyl_real wx2 = w[2], dv2 = dxv[2]; 
  const gkyl_real wx2_sq = wx2*wx2, dv2_sq = dv2*dv2; 
  const gkyl_real wx2_cu = wx2*wx2*wx2, dv2_cu = dv2*dv2*dv2; 
  out[0] += volFact*(2.0*f[0]*wx1*wx1_sq+1.732050807568877*f[2]*dv1*wx1_sq+0.4472135954999579*f[8]*dv1_sq*wx1+0.5*f[0]*dv1_sq*wx1+0.08660254037844387*f[2]*dv1*dv1_sq); 
  out[1] += volFact*(2.0*f[1]*wx1*wx1_sq+1.732050807568877*f[4]*dv1*wx1_sq+0.447213595499958*f[12]*dv1_sq*wx1+0.5*f[1]*dv1_sq*wx1+0.08660254037844387*f[4]*dv1*dv1_sq); 
  out[2] += volFact*(2.0*f[7]*wx1*wx1_sq+1.732050807568877*f[11]*dv1*wx1_sq+0.5*f[7]*dv1_sq*wx1+0.08660254037844385*f[11]*dv1*dv1_sq); 
  out[3] += volFact*(2.0*f[0]*wx1_sq*wx2+1.154700538379252*f[2]*dv1*wx1*wx2+0.149071198499986*f[8]*dv1_sq*wx2+0.1666666666666667*f[0]*dv1_sq*wx2+0.5773502691896258*f[3]*dv2*wx1_sq+0.3333333333333333*f[6]*dv1*dv2*wx1+0.04303314829119351*f[14]*dv1_sq*dv2+0.04811252243246882*f[3]*dv1_sq*dv2); 
  out[4] += volFact*(2.0*f[1]*wx1_sq*wx2+1.154700538379252*f[4]*dv1*wx1*wx2+0.149071198499986*f[12]*dv1_sq*wx2+0.1666666666666667*f[1]*dv1_sq*wx2+0.5773502691896258*f[5]*dv2*wx1_sq+0.3333333333333333*f[10]*dv1*dv2*wx1+0.04303314829119353*f[18]*dv1_sq*dv2+0.04811252243246882*f[5]*dv1_sq*dv2); 
  out[5] += volFact*(2.0*f[7]*wx1_sq*wx2+1.154700538379251*f[11]*dv1*wx1*wx2+0.1666666666666667*f[7]*dv1_sq*wx2+0.5773502691896257*f[13]*dv2*wx1_sq+0.3333333333333333*f[17]*dv1*dv2*wx1+0.04811252243246881*f[13]*dv1_sq*dv2); 
  out[6] += volFact*(2.0*f[0]*wx1*wx2_sq+0.5773502691896258*f[2]*dv1*wx2_sq+1.154700538379252*f[3]*dv2*wx1*wx2+0.3333333333333333*f[6]*dv1*dv2*wx2+0.149071198499986*f[9]*dv2_sq*wx1+0.1666666666666667*f[0]*dv2_sq*wx1+0.04303314829119351*f[16]*dv1*dv2_sq+0.04811252243246882*f[2]*dv1*dv2_sq); 
  out[7] += volFact*(2.0*f[1]*wx1*wx2_sq+0.5773502691896258*f[4]*dv1*wx2_sq+1.154700538379252*f[5]*dv2*wx1*wx2+0.3333333333333333*f[10]*dv1*dv2*wx2+0.149071198499986*f[15]*dv2_sq*wx1+0.1666666666666667*f[1]*dv2_sq*wx1+0.04303314829119353*f[19]*dv1*dv2_sq+0.04811252243246882*f[4]*dv1*dv2_sq); 
  out[8] += volFact*(2.0*f[7]*wx1*wx2_sq+0.5773502691896257*f[11]*dv1*wx2_sq+1.154700538379251*f[13]*dv2*wx1*wx2+0.3333333333333333*f[17]*dv1*dv2*wx2+0.1666666666666667*f[7]*dv2_sq*wx1+0.04811252243246881*f[11]*dv1*dv2_sq); 
  out[9] += volFact*(2.0*f[0]*wx2*wx2_sq+1.732050807568877*f[3]*dv2*wx2_sq+0.4472135954999579*f[9]*dv2_sq*wx2+0.5*f[0]*dv2_sq*wx2+0.08660254037844387*f[3]*dv2*dv2_sq); 
  out[10] += volFact*(2.0*f[1]*wx2*wx2_sq+1.732050807568877*f[5]*dv2*wx2_sq+0.447213595499958*f[15]*dv2_sq*wx2+0.5*f[1]*dv2_sq*wx2+0.08660254037844387*f[5]*dv2*dv2_sq); 
  out[11] += volFact*(2.0*f[7]*wx2*wx2_sq+1.732050807568877*f[13]*dv2*wx2_sq+0.5*f[7]*dv2_sq*wx2+0.08660254037844385*f[13]*dv2*dv2_sq); 
} 
