#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_2x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[3]:      Cell-center coordinates. 
  // dxv[3]:    Cell spacing. 
  // nu:     collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvpar = 2.0/dxv[2]; 
  const double dvpar = dxv[2], wvpar = w[2]; 
  double alphaDrag[12]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -1.414213562373095*nu[2]*rdvpar*wvpar; 
  alphaDrag[3] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[4] = -1.414213562373095*nu[3]*rdvpar*wvpar; 
  alphaDrag[5] = -0.408248290463863*nu[1]*dvpar*rdvpar; 
  alphaDrag[6] = -0.408248290463863*nu[2]*dvpar*rdvpar; 
  alphaDrag[7] = -0.408248290463863*nu[3]*dvpar*rdvpar; 

  out[3] += 0.6123724356957944*(alphaDrag[7]*f[7]+alphaDrag[6]*f[6]+alphaDrag[5]*f[5]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[5] += 0.6123724356957944*(alphaDrag[6]*f[7]+f[6]*alphaDrag[7]+alphaDrag[3]*f[5]+f[3]*alphaDrag[5]+alphaDrag[2]*f[4]+f[2]*alphaDrag[4]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[6] += 0.6123724356957944*(alphaDrag[5]*f[7]+f[5]*alphaDrag[7]+alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[7] += 0.6123724356957944*(alphaDrag[3]*f[7]+f[3]*alphaDrag[7]+alphaDrag[5]*f[6]+f[5]*alphaDrag[6]+alphaDrag[0]*f[4]+f[0]*alphaDrag[4]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[8] += 1.224744871391589*(alphaDrag[7]*f[11]+alphaDrag[6]*f[10]+alphaDrag[5]*f[9]+alphaDrag[3]*f[8])+1.369306393762915*(alphaDrag[4]*f[7]+f[4]*alphaDrag[7]+alphaDrag[2]*f[6]+f[2]*alphaDrag[6]+alphaDrag[1]*f[5]+f[1]*alphaDrag[5]+alphaDrag[0]*f[3]+f[0]*alphaDrag[3]); 
  out[9] += 1.224744871391589*(alphaDrag[6]*f[11]+alphaDrag[7]*f[10]+alphaDrag[3]*f[9]+alphaDrag[5]*f[8])+1.369306393762915*(alphaDrag[2]*f[7]+f[2]*alphaDrag[7]+alphaDrag[4]*f[6]+f[4]*alphaDrag[6]+alphaDrag[0]*f[5]+f[0]*alphaDrag[5]+alphaDrag[1]*f[3]+f[1]*alphaDrag[3]); 
  out[10] += 1.224744871391589*(alphaDrag[5]*f[11]+alphaDrag[3]*f[10]+alphaDrag[7]*f[9]+alphaDrag[6]*f[8])+1.369306393762915*(alphaDrag[1]*f[7]+f[1]*alphaDrag[7]+alphaDrag[0]*f[6]+f[0]*alphaDrag[6]+alphaDrag[4]*f[5]+f[4]*alphaDrag[5]+alphaDrag[2]*f[3]+f[2]*alphaDrag[3]); 
  out[11] += 1.224744871391589*(alphaDrag[3]*f[11]+alphaDrag[5]*f[10]+alphaDrag[6]*f[9]+alphaDrag[7]*f[8])+1.369306393762915*(alphaDrag[0]*f[7]+f[0]*alphaDrag[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[2]*f[5]+f[2]*alphaDrag[5]+alphaDrag[3]*f[4]+f[3]*alphaDrag[4]); 

  return fabs(0.8838834764831842*alphaDrag[0]); 

} 
