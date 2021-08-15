#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *E0 = &qmem[0]; 

  double alpha[2]; 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  double fUpwindQuad_l[2];
  double fUpwindQuad_r[2];
  fUpwindQuad_l[0] = copysignf(1.0,alpha[0]-alpha[1])*((-0.8660254037844386*(fl[3]+fc[3]))+0.8660254037844386*(fl[2]+fc[2])-0.5*fl[1]+0.5*(fc[1]+fl[0])-0.5*fc[0])-0.8660254037844386*fl[3]+0.8660254037844386*(fc[3]+fl[2])-0.8660254037844386*fc[2]-0.5*(fl[1]+fc[1])+0.5*(fl[0]+fc[0]); 
  fUpwindQuad_r[0] = copysignf(1.0,alpha[0]-alpha[1])*((-0.8660254037844386*(fr[3]+fc[3]))+0.8660254037844386*(fr[2]+fc[2])+0.5*fr[1]-0.5*(fc[1]+fr[0])+0.5*fc[0])+0.8660254037844386*fr[3]-0.8660254037844386*(fc[3]+fr[2])+0.8660254037844386*fc[2]-0.5*(fr[1]+fc[1])+0.5*(fr[0]+fc[0]); 
  fUpwindQuad_l[1] = copysignf(1.0,alpha[1]+alpha[0])*(0.8660254037844386*(fl[3]+fc[3]+fl[2]+fc[2])+0.5*(fl[1]+fl[0])-0.5*(fc[1]+fc[0]))+0.8660254037844386*(fl[3]+fl[2])-0.8660254037844386*(fc[3]+fc[2])+0.5*(fl[1]+fc[1]+fl[0]+fc[0]); 
  fUpwindQuad_r[1] = copysignf(1.0,alpha[1]+alpha[0])*(0.8660254037844386*(fr[3]+fc[3]+fr[2]+fc[2])-0.5*(fr[1]+fr[0])+0.5*(fc[1]+fc[0]))-0.8660254037844386*(fr[3]+fr[2])+0.8660254037844386*(fc[3]+fc[2])+0.5*(fr[1]+fc[1]+fr[0]+fc[0]); 
  double fUpwind_l[2];
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  double fUpwind_r[2];
  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  double Ghat_l[2]; 
  double Ghat_r[2]; 
  Ghat_l[0] = 0.7071067811865475*(alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] = 0.7071067811865475*(alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 

  Ghat_r[0] = 0.7071067811865475*(alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] = 0.7071067811865475*(alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 

} 
