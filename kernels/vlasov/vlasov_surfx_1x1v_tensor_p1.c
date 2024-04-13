#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfx_1x1v_tensor_p1(const double *w, const double *dxv, const double *alpha_geo, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // alpha_geo: Fields used only for general geometry.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[3]; 
  double Ghat_l[3]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[3]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[3]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[2] = (0.3162277660168379*fc[3]+0.1825741858350554*fc[2])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[3]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[3]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[2] = (0.3162277660168379*fl[3]+0.1825741858350554*fl[2])*dv; 

  } else { 

  Ghat_r[0] = -0.1178511301977579*((10.39230484541326*fr[1]-6.0*fr[0])*wv+(3.0*fr[3]-1.732050807568877*fr[2])*dv); 
  Ghat_r[1] = -0.1178511301977579*((10.39230484541326*fr[3]-6.0*fr[2])*wv+(3.0*fr[1]-1.732050807568877*fr[0])*dv); 
  Ghat_r[2] = -0.04714045207910316*(6.708203932499369*fr[3]-3.872983346207417*fr[2])*dv; 

  Ghat_l[0] = -0.1178511301977579*((10.39230484541326*fc[1]-6.0*fc[0])*wv+(3.0*fc[3]-1.732050807568877*fc[2])*dv); 
  Ghat_l[1] = -0.1178511301977579*((10.39230484541326*fc[3]-6.0*fc[2])*wv+(3.0*fc[1]-1.732050807568877*fc[0])*dv); 
  Ghat_l[2] = -0.04714045207910316*(6.708203932499369*fc[3]-3.872983346207417*fc[2])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 

  return 0.;

} 
