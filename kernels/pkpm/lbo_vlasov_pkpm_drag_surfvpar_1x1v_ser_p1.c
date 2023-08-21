#include <gkyl_lbo_vlasov_pkpm_kernels.h> 
GKYL_CU_DH double lbo_vlasov_pkpm_drag_surfvpar_1x1v_ser_p1(const double *w, const double *dxv, const double *nu, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:    cell-center coordinates. 
  // dxv[2]:  cell spacing. 
  // nu:       collisionalities added (self and cross species collisionalities). 
  // fl/fc/fr: Input Distribution function [F_0, T_perp G = T_perp (F_1 - F_0)] in left/center/right cells 
  // out:      incremented distribution function in cell 

  const double dv1par = 2.0/dxv[1]; 
  const double dvpar = dxv[1], wvpar = w[1]; 
  const double *F_0l = &fl[0]; 
  const double *G_1l = &fl[6]; 
  const double *F_0c = &fc[0]; 
  const double *G_1c = &fc[6]; 
  const double *F_0r = &fr[0]; 
  const double *G_1r = &fr[6]; 
  double *out_F_0 = &out[0]; 
  double *out_G_1 = &out[6]; 

  double alphaDrSurf_l[2] = {0.0}; 
  alphaDrSurf_l[0] = nu[0]*wvpar-0.5*nu[0]*dvpar; 
  alphaDrSurf_l[1] = nu[1]*wvpar-0.5*nu[1]*dvpar; 

  double alphaDrSurf_r[2] = {0.0}; 
  alphaDrSurf_r[0] = nu[0]*wvpar+0.5*nu[0]*dvpar; 
  alphaDrSurf_r[1] = nu[1]*wvpar+0.5*nu[1]*dvpar; 

  double Ghat_F_0_r[2]; 
  double Ghat_F_0_l[2]; 
  double Ghat_G_1_r[2]; 
  double Ghat_G_1_l[2]; 
  if (wvpar>0) { 

  Ghat_F_0_r[0] = 1.118033988749895*(alphaDrSurf_r[1]*F_0r[5]+alphaDrSurf_r[0]*F_0r[4])-0.8660254037844386*(alphaDrSurf_r[1]*F_0r[3]+alphaDrSurf_r[0]*F_0r[2])+0.5*(F_0r[1]*alphaDrSurf_r[1]+F_0r[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 1.118033988749895*(alphaDrSurf_r[0]*F_0r[5]+alphaDrSurf_r[1]*F_0r[4])-0.8660254037844386*(alphaDrSurf_r[0]*F_0r[3]+alphaDrSurf_r[1]*F_0r[2])+0.5*(F_0r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0r[1]); 
  Ghat_G_1_r[0] = 1.118033988749895*(alphaDrSurf_r[1]*G_1r[5]+alphaDrSurf_r[0]*G_1r[4])-0.8660254037844386*(alphaDrSurf_r[1]*G_1r[3]+alphaDrSurf_r[0]*G_1r[2])+0.5*(G_1r[1]*alphaDrSurf_r[1]+G_1r[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 1.118033988749895*(alphaDrSurf_r[0]*G_1r[5]+alphaDrSurf_r[1]*G_1r[4])-0.8660254037844386*(alphaDrSurf_r[0]*G_1r[3]+alphaDrSurf_r[1]*G_1r[2])+0.5*(G_1r[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1r[1]); 

  Ghat_F_0_l[0] = 1.118033988749895*(alphaDrSurf_l[1]*F_0c[5]+alphaDrSurf_l[0]*F_0c[4])-0.8660254037844386*(alphaDrSurf_l[1]*F_0c[3]+alphaDrSurf_l[0]*F_0c[2])+0.5*(F_0c[1]*alphaDrSurf_l[1]+F_0c[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 1.118033988749895*(alphaDrSurf_l[0]*F_0c[5]+alphaDrSurf_l[1]*F_0c[4])-0.8660254037844386*(alphaDrSurf_l[0]*F_0c[3]+alphaDrSurf_l[1]*F_0c[2])+0.5*(F_0c[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0c[1]); 
  Ghat_G_1_l[0] = 1.118033988749895*(alphaDrSurf_l[1]*G_1c[5]+alphaDrSurf_l[0]*G_1c[4])-0.8660254037844386*(alphaDrSurf_l[1]*G_1c[3]+alphaDrSurf_l[0]*G_1c[2])+0.5*(G_1c[1]*alphaDrSurf_l[1]+G_1c[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 1.118033988749895*(alphaDrSurf_l[0]*G_1c[5]+alphaDrSurf_l[1]*G_1c[4])-0.8660254037844386*(alphaDrSurf_l[0]*G_1c[3]+alphaDrSurf_l[1]*G_1c[2])+0.5*(G_1c[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1c[1]); 

  } else { 

  Ghat_F_0_r[0] = 1.118033988749895*(alphaDrSurf_r[1]*F_0c[5]+alphaDrSurf_r[0]*F_0c[4])+0.8660254037844386*(alphaDrSurf_r[1]*F_0c[3]+alphaDrSurf_r[0]*F_0c[2])+0.5*(F_0c[1]*alphaDrSurf_r[1]+F_0c[0]*alphaDrSurf_r[0]); 
  Ghat_F_0_r[1] = 1.118033988749895*(alphaDrSurf_r[0]*F_0c[5]+alphaDrSurf_r[1]*F_0c[4])+0.8660254037844386*(alphaDrSurf_r[0]*F_0c[3]+alphaDrSurf_r[1]*F_0c[2])+0.5*(F_0c[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*F_0c[1]); 
  Ghat_G_1_r[0] = 1.118033988749895*(alphaDrSurf_r[1]*G_1c[5]+alphaDrSurf_r[0]*G_1c[4])+0.8660254037844386*(alphaDrSurf_r[1]*G_1c[3]+alphaDrSurf_r[0]*G_1c[2])+0.5*(G_1c[1]*alphaDrSurf_r[1]+G_1c[0]*alphaDrSurf_r[0]); 
  Ghat_G_1_r[1] = 1.118033988749895*(alphaDrSurf_r[0]*G_1c[5]+alphaDrSurf_r[1]*G_1c[4])+0.8660254037844386*(alphaDrSurf_r[0]*G_1c[3]+alphaDrSurf_r[1]*G_1c[2])+0.5*(G_1c[0]*alphaDrSurf_r[1]+alphaDrSurf_r[0]*G_1c[1]); 

  Ghat_F_0_l[0] = 1.118033988749895*(alphaDrSurf_l[1]*F_0l[5]+alphaDrSurf_l[0]*F_0l[4])+0.8660254037844386*(alphaDrSurf_l[1]*F_0l[3]+alphaDrSurf_l[0]*F_0l[2])+0.5*(F_0l[1]*alphaDrSurf_l[1]+F_0l[0]*alphaDrSurf_l[0]); 
  Ghat_F_0_l[1] = 1.118033988749895*(alphaDrSurf_l[0]*F_0l[5]+alphaDrSurf_l[1]*F_0l[4])+0.8660254037844386*(alphaDrSurf_l[0]*F_0l[3]+alphaDrSurf_l[1]*F_0l[2])+0.5*(F_0l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*F_0l[1]); 
  Ghat_G_1_l[0] = 1.118033988749895*(alphaDrSurf_l[1]*G_1l[5]+alphaDrSurf_l[0]*G_1l[4])+0.8660254037844386*(alphaDrSurf_l[1]*G_1l[3]+alphaDrSurf_l[0]*G_1l[2])+0.5*(G_1l[1]*alphaDrSurf_l[1]+G_1l[0]*alphaDrSurf_l[0]); 
  Ghat_G_1_l[1] = 1.118033988749895*(alphaDrSurf_l[0]*G_1l[5]+alphaDrSurf_l[1]*G_1l[4])+0.8660254037844386*(alphaDrSurf_l[0]*G_1l[3]+alphaDrSurf_l[1]*G_1l[2])+0.5*(G_1l[0]*alphaDrSurf_l[1]+alphaDrSurf_l[0]*G_1l[1]); 

  } 
  out_F_0[0] += (0.7071067811865475*Ghat_F_0_r[0]-0.7071067811865475*Ghat_F_0_l[0])*dv1par; 
  out_F_0[1] += (0.7071067811865475*Ghat_F_0_r[1]-0.7071067811865475*Ghat_F_0_l[1])*dv1par; 
  out_F_0[2] += 1.224744871391589*(Ghat_F_0_r[0]+Ghat_F_0_l[0])*dv1par; 
  out_F_0[3] += 1.224744871391589*(Ghat_F_0_r[1]+Ghat_F_0_l[1])*dv1par; 
  out_F_0[4] += (1.58113883008419*Ghat_F_0_r[0]-1.58113883008419*Ghat_F_0_l[0])*dv1par; 
  out_F_0[5] += (1.58113883008419*Ghat_F_0_r[1]-1.58113883008419*Ghat_F_0_l[1])*dv1par; 
  out_G_1[0] += (0.7071067811865475*Ghat_G_1_r[0]-0.7071067811865475*Ghat_G_1_l[0])*dv1par; 
  out_G_1[1] += (0.7071067811865475*Ghat_G_1_r[1]-0.7071067811865475*Ghat_G_1_l[1])*dv1par; 
  out_G_1[2] += 1.224744871391589*(Ghat_G_1_r[0]+Ghat_G_1_l[0])*dv1par; 
  out_G_1[3] += 1.224744871391589*(Ghat_G_1_r[1]+Ghat_G_1_l[1])*dv1par; 
  out_G_1[4] += (1.58113883008419*Ghat_G_1_r[0]-1.58113883008419*Ghat_G_1_l[0])*dv1par; 
  out_G_1[5] += (1.58113883008419*Ghat_G_1_r[1]-1.58113883008419*Ghat_G_1_l[1])*dv1par; 

  return 0.;

} 
