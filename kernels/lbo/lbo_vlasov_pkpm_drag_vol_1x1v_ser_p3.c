#include <gkyl_lbo_vlasov_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_vol_1x1v_ser_p3(const double *w, const double *dxv, const double *nu, const double *f, double* GKYL_RESTRICT out) 
{ 
  // w[2]:      Cell-center coordinates. 
  // dxv[2]:    Cell spacing. 
  // nu:     collisionality. 
  // f:         Input distribution function.
  // out:       Incremented output 
  const double rdvpar = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  double alphaDrag[12]; 
  // Expand rdvpar*(nu*vx) in phase basis.
  alphaDrag[0] = -1.414213562373095*nu[0]*rdvpar*wvpar; 
  alphaDrag[1] = -1.414213562373095*nu[1]*rdvpar*wvpar; 
  alphaDrag[2] = -0.408248290463863*nu[0]*dvpar*rdvpar; 
  alphaDrag[3] = -0.408248290463863*nu[1]*dvpar*rdvpar; 
  alphaDrag[4] = -1.414213562373095*nu[2]*rdvpar*wvpar; 
  alphaDrag[6] = -0.408248290463863*nu[2]*dvpar*rdvpar; 
  alphaDrag[8] = -1.414213562373095*nu[3]*rdvpar*wvpar; 
  alphaDrag[10] = -0.408248290463863*nu[3]*dvpar*rdvpar; 

  out[2] += 0.8660254037844386*(alphaDrag[10]*f[10]+alphaDrag[8]*f[8]+alphaDrag[6]*f[6]+alphaDrag[4]*f[4]+alphaDrag[3]*f[3]+alphaDrag[2]*f[2]+alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[3] += 0.7606388292556648*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10]+alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+0.7745966692414833*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+0.8660254037844386*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3]+alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 
  out[5] += 1.936491673103709*(alphaDrag[8]*f[10]+f[8]*alphaDrag[10])+1.732050807568877*alphaDrag[3]*f[7]+1.936491673103709*(alphaDrag[4]*f[6]+f[4]*alphaDrag[6])+1.732050807568877*alphaDrag[2]*f[5]+1.936491673103709*(alphaDrag[1]*f[3]+f[1]*alphaDrag[3]+alphaDrag[0]*f[2]+f[0]*alphaDrag[2]); 
  out[6] += 0.5163977794943223*alphaDrag[10]*f[10]+0.7606388292556648*(alphaDrag[3]*f[10]+f[3]*alphaDrag[10])+0.5163977794943223*alphaDrag[8]*f[8]+0.7606388292556648*(alphaDrag[1]*f[8]+f[1]*alphaDrag[8])+0.5532833351724881*alphaDrag[6]*f[6]+0.8660254037844386*(alphaDrag[2]*f[6]+f[2]*alphaDrag[6])+0.5532833351724881*alphaDrag[4]*f[4]+0.8660254037844386*(alphaDrag[0]*f[4]+f[0]*alphaDrag[4])+0.7745966692414833*(alphaDrag[3]*f[3]+alphaDrag[1]*f[1]); 
  out[7] += 1.700840128541522*(alphaDrag[4]*f[10]+f[4]*alphaDrag[10]+alphaDrag[6]*f[8]+f[6]*alphaDrag[8])+1.549193338482967*alphaDrag[6]*f[7]+1.732050807568877*(alphaDrag[2]*f[7]+alphaDrag[1]*f[6]+f[1]*alphaDrag[6]+alphaDrag[3]*(f[5]+f[4])+f[3]*alphaDrag[4])+1.936491673103709*(alphaDrag[0]*f[3]+f[0]*alphaDrag[3]+alphaDrag[1]*f[2]+f[1]*alphaDrag[2]); 
  out[9] += 2.598076211353316*alphaDrag[3]*f[11]+3.968626966596886*alphaDrag[10]*f[10]+2.598076211353316*alphaDrag[2]*f[9]+1.322875655532295*alphaDrag[8]*f[8]+2.958039891549809*alphaDrag[1]*f[7]+3.968626966596886*alphaDrag[6]*f[6]+2.958039891549809*alphaDrag[0]*f[5]+1.322875655532295*alphaDrag[4]*f[4]+3.968626966596886*(alphaDrag[3]*f[3]+alphaDrag[2]*f[2])+1.322875655532295*(alphaDrag[1]*f[1]+alphaDrag[0]*f[0]); 
  out[10] += (0.5163977794943223*alphaDrag[6]+0.8660254037844386*alphaDrag[2])*f[10]+(0.5163977794943223*f[6]+0.8660254037844386*f[2])*alphaDrag[10]+(0.5163977794943223*alphaDrag[4]+0.8660254037844386*alphaDrag[0])*f[8]+(0.5163977794943223*f[4]+0.8660254037844386*f[0])*alphaDrag[8]+0.7606388292556648*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6]+alphaDrag[1]*f[4]+f[1]*alphaDrag[4]); 
  out[11] += (2.32379000772445*alphaDrag[6]+2.598076211353316*alphaDrag[2])*f[11]+3.485685011586674*(alphaDrag[6]*f[10]+f[6]*alphaDrag[10])+2.598076211353316*alphaDrag[3]*f[9]+1.161895003862225*(alphaDrag[4]*f[8]+f[4]*alphaDrag[8])+(2.645751311064591*alphaDrag[4]+2.958039891549809*alphaDrag[0])*f[7]+3.54964786985977*(alphaDrag[3]*f[6]+f[3]*alphaDrag[6])+2.958039891549809*alphaDrag[1]*f[5]+1.183215956619923*(alphaDrag[1]*f[4]+f[1]*alphaDrag[4])+3.968626966596886*(alphaDrag[2]*f[3]+f[2]*alphaDrag[3])+1.322875655532295*(alphaDrag[0]*f[1]+f[0]*alphaDrag[1]); 

  return fabs(1.75*alphaDrag[0]-1.956559480312316*alphaDrag[4]); 

} 
