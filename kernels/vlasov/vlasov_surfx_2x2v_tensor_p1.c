#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfx_2x2v_tensor_p1(const double *w, const double *dxv, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // alpha_surf_l: Surface expansion of phase space flux on the left (used by general geometry version).
  // alpha_surf_r: Surface expansion of phase space flux on the right (used by general geometry version).
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points (used by general geometry version).
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points (used by general geometry version).
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1 (used by general geometry version).
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1 (used by general geometry version).
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[2], wv = w[2]; 
  double Ghat_r[16]; 
  double Ghat_l[16]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[6]+0.2041241452319315*fc[3])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[5]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[11]+0.2041241452319315*fc[7])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[6]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[8]+0.7071067811865475*fc[4])*wv+(0.3535533905932737*fc[13]+0.2041241452319315*fc[10])*dv; 
  Ghat_r[4] = (1.224744871391589*fc[11]+0.7071067811865475*fc[7])*wv+(0.3535533905932737*fc[5]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[5] = (1.224744871391589*fc[12]+0.7071067811865475*fc[9])*wv+(0.3535533905932737*fc[15]+0.2041241452319315*fc[14])*dv; 
  Ghat_r[6] = (1.224744871391589*fc[13]+0.7071067811865475*fc[10])*wv+(0.3535533905932737*fc[8]+0.2041241452319315*fc[4])*dv; 
  Ghat_r[7] = (1.224744871391589*fc[15]+0.7071067811865475*fc[14])*wv+(0.3535533905932737*fc[12]+0.2041241452319315*fc[9])*dv; 
  Ghat_r[8] = (0.3162277660168379*fc[6]+0.1825741858350554*fc[3])*dv; 
  Ghat_r[9] = (0.3162277660168379*fc[11]+0.1825741858350554*fc[7])*dv; 
  Ghat_r[10] = (0.3162277660168379*fc[13]+0.1825741858350554*fc[10])*dv; 
  Ghat_r[11] = (0.3162277660168379*fc[15]+0.1825741858350554*fc[14])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[6]+0.2041241452319315*fl[3])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[5]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[11]+0.2041241452319315*fl[7])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[6]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[8]+0.7071067811865475*fl[4])*wv+(0.3535533905932737*fl[13]+0.2041241452319315*fl[10])*dv; 
  Ghat_l[4] = (1.224744871391589*fl[11]+0.7071067811865475*fl[7])*wv+(0.3535533905932737*fl[5]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[5] = (1.224744871391589*fl[12]+0.7071067811865475*fl[9])*wv+(0.3535533905932737*fl[15]+0.2041241452319315*fl[14])*dv; 
  Ghat_l[6] = (1.224744871391589*fl[13]+0.7071067811865475*fl[10])*wv+(0.3535533905932737*fl[8]+0.2041241452319315*fl[4])*dv; 
  Ghat_l[7] = (1.224744871391589*fl[15]+0.7071067811865475*fl[14])*wv+(0.3535533905932737*fl[12]+0.2041241452319315*fl[9])*dv; 
  Ghat_l[8] = (0.3162277660168379*fl[6]+0.1825741858350554*fl[3])*dv; 
  Ghat_l[9] = (0.3162277660168379*fl[11]+0.1825741858350554*fl[7])*dv; 
  Ghat_l[10] = (0.3162277660168379*fl[13]+0.1825741858350554*fl[10])*dv; 
  Ghat_l[11] = (0.3162277660168379*fl[15]+0.1825741858350554*fl[14])*dv; 

  } else { 

  Ghat_r[0] = -0.1178511301977579*((10.39230484541326*fr[1]-6.0*fr[0])*wv+(3.0*fr[6]-1.732050807568877*fr[3])*dv); 
  Ghat_r[1] = -0.1178511301977579*((10.39230484541326*fr[5]-6.0*fr[2])*wv+(3.0*fr[11]-1.732050807568877*fr[7])*dv); 
  Ghat_r[2] = -0.1178511301977579*((10.39230484541326*fr[6]-6.0*fr[3])*wv+(3.0*fr[1]-1.732050807568877*fr[0])*dv); 
  Ghat_r[3] = -0.1178511301977579*((10.39230484541326*fr[8]-6.0*fr[4])*wv+(3.0*fr[13]-1.732050807568877*fr[10])*dv); 
  Ghat_r[4] = -0.1178511301977579*((10.39230484541326*fr[11]-6.0*fr[7])*wv+(3.0*fr[5]-1.732050807568877*fr[2])*dv); 
  Ghat_r[5] = -0.1178511301977579*((10.39230484541326*fr[12]-6.0*fr[9])*wv+(3.0*fr[15]-1.732050807568877*fr[14])*dv); 
  Ghat_r[6] = -0.1178511301977579*((10.39230484541326*fr[13]-6.0*fr[10])*wv+(3.0*fr[8]-1.732050807568877*fr[4])*dv); 
  Ghat_r[7] = -0.1178511301977579*((10.39230484541326*fr[15]-6.0*fr[14])*wv+(3.0*fr[12]-1.732050807568877*fr[9])*dv); 
  Ghat_r[8] = -0.04714045207910316*(6.708203932499369*fr[6]-3.872983346207417*fr[3])*dv; 
  Ghat_r[9] = -0.04714045207910316*(6.708203932499369*fr[11]-3.872983346207417*fr[7])*dv; 
  Ghat_r[10] = -0.04714045207910316*(6.708203932499369*fr[13]-3.872983346207417*fr[10])*dv; 
  Ghat_r[11] = -0.04714045207910316*(6.708203932499369*fr[15]-3.872983346207417*fr[14])*dv; 

  Ghat_l[0] = -0.1178511301977579*((10.39230484541326*fc[1]-6.0*fc[0])*wv+(3.0*fc[6]-1.732050807568877*fc[3])*dv); 
  Ghat_l[1] = -0.1178511301977579*((10.39230484541326*fc[5]-6.0*fc[2])*wv+(3.0*fc[11]-1.732050807568877*fc[7])*dv); 
  Ghat_l[2] = -0.1178511301977579*((10.39230484541326*fc[6]-6.0*fc[3])*wv+(3.0*fc[1]-1.732050807568877*fc[0])*dv); 
  Ghat_l[3] = -0.1178511301977579*((10.39230484541326*fc[8]-6.0*fc[4])*wv+(3.0*fc[13]-1.732050807568877*fc[10])*dv); 
  Ghat_l[4] = -0.1178511301977579*((10.39230484541326*fc[11]-6.0*fc[7])*wv+(3.0*fc[5]-1.732050807568877*fc[2])*dv); 
  Ghat_l[5] = -0.1178511301977579*((10.39230484541326*fc[12]-6.0*fc[9])*wv+(3.0*fc[15]-1.732050807568877*fc[14])*dv); 
  Ghat_l[6] = -0.1178511301977579*((10.39230484541326*fc[13]-6.0*fc[10])*wv+(3.0*fc[8]-1.732050807568877*fc[4])*dv); 
  Ghat_l[7] = -0.1178511301977579*((10.39230484541326*fc[15]-6.0*fc[14])*wv+(3.0*fc[12]-1.732050807568877*fc[9])*dv); 
  Ghat_l[8] = -0.04714045207910316*(6.708203932499369*fc[6]-3.872983346207417*fc[3])*dv; 
  Ghat_l[9] = -0.04714045207910316*(6.708203932499369*fc[11]-3.872983346207417*fc[7])*dv; 
  Ghat_l[10] = -0.04714045207910316*(6.708203932499369*fc[13]-3.872983346207417*fc[10])*dv; 
  Ghat_l[11] = -0.04714045207910316*(6.708203932499369*fc[15]-3.872983346207417*fc[14])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[7] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dx10; 
  out[8] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
  out[9] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dx10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dx10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dx10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dx10; 
  out[13] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dx10; 
  out[14] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dx10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dx10; 

  return 0.;

} 
