#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_vol_1x1v_ser_p2(const double *w, const double *dxv, const double *field, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // field:      q/m*EM fields.
  // cot_vec:   Only used in gen geo.
  // f:         Input distribution function.
  // out:       Incremented output.
  double dv0dx0 = dxv[1]/dxv[0]; 
  double w0dx0 = w[1]/dxv[0]; 
  const double dv10 = 2/dxv[1]; 
  const double *E0 = &field[0]; 
  const double dv1 = dxv[1], wv1 = w[1]; 

  double cflFreq_mid = 0.0; 
  double alpha_cdim[8] = {0.0}; 
  double alpha_vdim[8] = {0.0}; 

  alpha_cdim[0] = 4.0*w0dx0; 
  alpha_cdim[2] = 1.154700538379252*dv0dx0; 
  cflFreq_mid += 5.0*(fabs(w0dx0)+0.5*dv0dx0); 

  alpha_vdim[0] = 1.414213562373095*E0[0]*dv10; 
  alpha_vdim[1] = 1.414213562373095*E0[1]*dv10; 
  alpha_vdim[4] = 1.414213562373095*E0[2]*dv10; 
  cflFreq_mid += 5.0*fabs(0.25*alpha_vdim[0]-0.2795084971874737*alpha_vdim[4]); 

  out[1] += 0.8660254037844386*(alpha_cdim[2]*f[2]+alpha_cdim[0]*f[0]); 
  out[2] += 0.8660254037844386*(alpha_vdim[4]*f[4]+alpha_vdim[1]*f[1]+alpha_vdim[0]*f[0]); 
  out[3] += 0.7745966692414833*(alpha_cdim[2]*f[5]+alpha_vdim[1]*f[4]+f[1]*alpha_vdim[4])+0.8660254037844386*(alpha_cdim[0]*f[2]+f[0]*alpha_cdim[2]+alpha_vdim[0]*f[1]+f[0]*alpha_vdim[1]); 
  out[4] += 1.936491673103709*(alpha_cdim[2]*f[3]+alpha_cdim[0]*f[1]); 
  out[5] += 1.936491673103709*(alpha_vdim[4]*f[6]+alpha_vdim[1]*f[3]+alpha_vdim[0]*f[2]); 
  out[6] += 1.732050807568877*alpha_cdim[2]*f[7]+0.5532833351724881*alpha_vdim[4]*f[4]+0.8660254037844386*(alpha_vdim[0]*f[4]+f[0]*alpha_vdim[4])+1.936491673103709*alpha_cdim[0]*f[3]+f[1]*(1.936491673103709*alpha_cdim[2]+0.7745966692414833*alpha_vdim[1]); 
  out[7] += 1.732050807568877*alpha_vdim[1]*f[6]+0.8660254037844386*alpha_cdim[0]*f[5]+f[3]*(1.732050807568877*alpha_vdim[4]+1.936491673103709*alpha_vdim[0])+(0.7745966692414833*alpha_cdim[2]+1.936491673103709*alpha_vdim[1])*f[2]; 

  return cflFreq_mid; 
} 
