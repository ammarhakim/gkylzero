#include <gkyl_vlasov_kernels.h> 
gkyl_real vlasov_stream_vol_1x2v_ser_p1(const gkyl_real *w, const gkyl_real *dxv, const gkyl_real *f, gkyl_real* restrict out) 
{ 
  // w[NDIM]: Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // f: Input distribution function.
  // out: Incremented output.
  gkyl_real w1Ddx0  = w[1]/dxv[0]; 
  gkyl_real dv1Ddx0 = dxv[1]/dxv[0]; 

  out[1] += 3.464101615137754*f[0]*w1Ddx0+f[2]*dv1Ddx0; 
  out[4] += 3.464101615137754*f[2]*w1Ddx0+f[0]*dv1Ddx0; 
  out[5] += 3.464101615137754*f[3]*w1Ddx0+f[6]*dv1Ddx0; 
  out[7] += 3.464101615137754*f[6]*w1Ddx0+f[3]*dv1Ddx0; 

  return fabs(w1Ddx0)+0.5*dv1Ddx0;
} 
