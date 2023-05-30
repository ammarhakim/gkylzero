#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_1x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:   Cell-center coordinates. 
  // dxv[2]: Cell spacing. 
  // nu:      collisionality. 
  // f:       Input distribution function.
  // out:     Incremented output 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0 = &f[0]; 
  const double *G_1 = &f[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 
  double alphaDrag[6]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[3] = -0.408248290463863*nu[1]*dvpar*rdvpar; 

  out_F_0[2] += 0.8660254037844386*F_0[3]*alphaDrag[3]+0.8660254037844386*F_0[2]*alphaDrag[2]+0.8660254037844386*F_0[1]*alphaDrag[1]+0.8660254037844386*F_0[0]*alphaDrag[0]; 
  out_F_0[3] += 0.8660254037844386*F_0[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*F_0[3]+0.8660254037844386*F_0[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*F_0[1]; 
  out_F_0[4] += 1.732050807568877*alphaDrag[3]*F_0[5]+1.732050807568877*alphaDrag[2]*F_0[4]+1.936491673103709*F_0[1]*alphaDrag[3]+1.936491673103709*alphaDrag[1]*F_0[3]+1.936491673103709*F_0[0]*alphaDrag[2]+1.936491673103709*alphaDrag[0]*F_0[2]; 
  out_F_0[5] += 1.732050807568877*alphaDrag[2]*F_0[5]+1.732050807568877*alphaDrag[3]*F_0[4]+1.936491673103709*F_0[0]*alphaDrag[3]+1.936491673103709*alphaDrag[0]*F_0[3]+1.936491673103709*F_0[1]*alphaDrag[2]+1.936491673103709*alphaDrag[1]*F_0[2]; 
  out_G_1[2] += 0.8660254037844386*G_1[3]*alphaDrag[3]+0.8660254037844386*G_1[2]*alphaDrag[2]+0.8660254037844386*G_1[1]*alphaDrag[1]+0.8660254037844386*G_1[0]*alphaDrag[0]; 
  out_G_1[3] += 0.8660254037844386*G_1[2]*alphaDrag[3]+0.8660254037844386*alphaDrag[2]*G_1[3]+0.8660254037844386*G_1[0]*alphaDrag[1]+0.8660254037844386*alphaDrag[0]*G_1[1]; 
  out_G_1[4] += 1.732050807568877*alphaDrag[3]*G_1[5]+1.732050807568877*alphaDrag[2]*G_1[4]+1.936491673103709*G_1[1]*alphaDrag[3]+1.936491673103709*alphaDrag[1]*G_1[3]+1.936491673103709*G_1[0]*alphaDrag[2]+1.936491673103709*alphaDrag[0]*G_1[2]; 
  out_G_1[5] += 1.732050807568877*alphaDrag[2]*G_1[5]+1.732050807568877*alphaDrag[3]*G_1[4]+1.936491673103709*G_1[0]*alphaDrag[3]+1.936491673103709*alphaDrag[0]*G_1[3]+1.936491673103709*G_1[1]*alphaDrag[2]+1.936491673103709*alphaDrag[1]*G_1[2]; 

  return fabs(1.25*alphaDrag[0]); 

} 
