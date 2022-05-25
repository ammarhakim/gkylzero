#include <gkyl_mom_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_int_mom_2x3v_ser_p1(const double *w, const double *dxv, const int *idx, const double *p_over_gamma, const double *f, double* GKYL_RESTRICT out) 
{ 
  const double volFact = dxv[0]*dxv[1]*dxv[2]*dxv[3]*dxv[4]*0.03125; 
  const double wx1 = w[2], dv1 = dxv[2]; 
  const double *p0_over_gamma = &p_over_gamma[0]; 
  const double wx2 = w[3], dv2 = dxv[3]; 
  const double *p1_over_gamma = &p_over_gamma[8]; 
  const double wx3 = w[4], dv3 = dxv[4]; 
  const double *p2_over_gamma = &p_over_gamma[16]; 
 
  out[0] += 5.656854249492382*f[0]*volFact; 
  out[1] += (2.0*p0_over_gamma[7]*f[25]+2.0*p0_over_gamma[6]*f[15]+2.0*p0_over_gamma[5]*f[14]+2.0*p0_over_gamma[4]*f[11]+2.0*p0_over_gamma[3]*f[5]+2.0*p0_over_gamma[2]*f[4]+2.0*p0_over_gamma[1]*f[3]+2.0*f[0]*p0_over_gamma[0])*volFact; 
  out[2] += (2.0*p1_over_gamma[7]*f[25]+2.0*p1_over_gamma[6]*f[15]+2.0*p1_over_gamma[5]*f[14]+2.0*p1_over_gamma[4]*f[11]+2.0*p1_over_gamma[3]*f[5]+2.0*p1_over_gamma[2]*f[4]+2.0*p1_over_gamma[1]*f[3]+2.0*f[0]*p1_over_gamma[0])*volFact; 
  out[3] += (2.0*p2_over_gamma[7]*f[25]+2.0*p2_over_gamma[6]*f[15]+2.0*p2_over_gamma[5]*f[14]+2.0*p2_over_gamma[4]*f[11]+2.0*p2_over_gamma[3]*f[5]+2.0*p2_over_gamma[2]*f[4]+2.0*p2_over_gamma[1]*f[3]+2.0*f[0]*p2_over_gamma[0])*volFact; 
  out[4] += volFact*(2.0*p2_over_gamma[7]*f[25]*wx3+2.0*p2_over_gamma[6]*f[15]*wx3+2.0*p2_over_gamma[5]*f[14]*wx3+2.0*p2_over_gamma[4]*f[11]*wx3+2.0*p2_over_gamma[3]*f[5]*wx3+2.0*p2_over_gamma[2]*f[4]*wx3+2.0*p2_over_gamma[1]*f[3]*wx3+2.0*f[0]*p2_over_gamma[0]*wx3+2.0*p1_over_gamma[7]*f[25]*wx2+2.0*p1_over_gamma[6]*f[15]*wx2+2.0*p1_over_gamma[5]*f[14]*wx2+2.0*p1_over_gamma[4]*f[11]*wx2+2.0*p1_over_gamma[3]*f[5]*wx2+2.0*p1_over_gamma[2]*f[4]*wx2+2.0*p1_over_gamma[1]*f[3]*wx2+2.0*f[0]*p1_over_gamma[0]*wx2+2.0*p0_over_gamma[7]*f[25]*wx1+2.0*p0_over_gamma[6]*f[15]*wx1+2.0*p0_over_gamma[5]*f[14]*wx1+2.0*p0_over_gamma[4]*f[11]*wx1+2.0*p0_over_gamma[3]*f[5]*wx1+2.0*p0_over_gamma[2]*f[4]*wx1+2.0*p0_over_gamma[1]*f[3]*wx1+2.0*f[0]*p0_over_gamma[0]*wx1+0.5773502691896258*p2_over_gamma[4]*f[25]*dv3+0.5773502691896258*p2_over_gamma[2]*f[15]*dv3+0.5773502691896258*p2_over_gamma[1]*f[14]*dv3+0.5773502691896258*p2_over_gamma[7]*f[11]*dv3+0.5773502691896258*f[4]*p2_over_gamma[6]*dv3+0.5773502691896258*f[3]*p2_over_gamma[5]*dv3+0.5773502691896258*p2_over_gamma[0]*f[5]*dv3+0.5773502691896258*f[0]*p2_over_gamma[3]*dv3+0.5773502691896258*p1_over_gamma[5]*f[25]*dv2+0.5773502691896258*p1_over_gamma[3]*f[15]*dv2+0.5773502691896258*p1_over_gamma[7]*f[14]*dv2+0.5773502691896258*p1_over_gamma[1]*f[11]*dv2+0.5773502691896258*f[5]*p1_over_gamma[6]*dv2+0.5773502691896258*f[3]*p1_over_gamma[4]*dv2+0.5773502691896258*p1_over_gamma[0]*f[4]*dv2+0.5773502691896258*f[0]*p1_over_gamma[2]*dv2+0.5773502691896258*p0_over_gamma[6]*f[25]*dv1+0.5773502691896258*p0_over_gamma[7]*f[15]*dv1+0.5773502691896258*p0_over_gamma[3]*f[14]*dv1+0.5773502691896258*p0_over_gamma[2]*f[11]*dv1+0.5773502691896258*f[5]*p0_over_gamma[5]*dv1+0.5773502691896258*f[4]*p0_over_gamma[4]*dv1+0.5773502691896258*p0_over_gamma[0]*f[3]*dv1+0.5773502691896258*f[0]*p0_over_gamma[1]*dv1); 
} 
