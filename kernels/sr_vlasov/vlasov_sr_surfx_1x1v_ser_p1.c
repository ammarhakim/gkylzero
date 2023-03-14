#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_sr_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *p_over_gamma, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // p_over_gamma: p/gamma (velocity).
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double *p0_over_gamma = &p_over_gamma[0]; 
  double alpha[2] = {0.0}; 
  alpha[0] = p0_over_gamma[0]; 
  alpha[1] = p0_over_gamma[1]; 

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[2]; 
  double Ghat_l[2]; 
  if (wv>0) { 

  Ghat_r[0] = alpha[1]*(0.8660254037844386*fc[3]+0.5*fc[2])+alpha[0]*(0.8660254037844386*fc[1]+0.5*fc[0]); 
  Ghat_r[1] = alpha[0]*(0.8660254037844386*fc[3]+0.5*fc[2])+alpha[1]*(0.8660254037844386*fc[1]+0.5*fc[0]); 

  Ghat_l[0] = alpha[1]*(0.8660254037844386*fl[3]+0.5*fl[2])+alpha[0]*(0.8660254037844386*fl[1]+0.5*fl[0]); 
  Ghat_l[1] = alpha[0]*(0.8660254037844386*fl[3]+0.5*fl[2])+alpha[1]*(0.8660254037844386*fl[1]+0.5*fl[0]); 

  } else { 

  Ghat_r[0] = -0.5*(alpha[1]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[0]*(1.732050807568877*fr[1]-1.0*fr[0])); 
  Ghat_r[1] = -0.5*(alpha[0]*(1.732050807568877*fr[3]-1.0*fr[2])+alpha[1]*(1.732050807568877*fr[1]-1.0*fr[0])); 

  Ghat_l[0] = -0.5*(alpha[1]*(1.732050807568877*fc[3]-1.0*fc[2])+alpha[0]*(1.732050807568877*fc[1]-1.0*fc[0])); 
  Ghat_l[1] = -0.5*(alpha[0]*(1.732050807568877*fc[3]-1.0*fc[2])+alpha[1]*(1.732050807568877*fc[1]-1.0*fc[0])); 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
} 
