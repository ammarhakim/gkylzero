#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* restrict out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv10 = 2/dxv[1]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double *E0 = &qmem[0]; 

  double Ghat_r[3]; 
  double favg_r[3]; 
  double Ghat_l[3]; 
  double favg_l[3]; 
  double alpha[3]; 

  favg_r[0] = 1.58113883008419*(fr[5]+fc[5])-1.224744871391589*fr[2]+1.224744871391589*fc[2]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = 1.58113883008419*(fr[7]+fc[7])-1.224744871391589*fr[3]+1.224744871391589*fc[3]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[6])+1.224744871391589*fc[6]+0.7071067811865475*(fr[4]+fc[4]); 

  favg_l[0] = 1.58113883008419*(fl[5]+fc[5])+1.224744871391589*fl[2]-1.224744871391589*fc[2]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.58113883008419*(fl[7]+fc[7])+1.224744871391589*fl[3]-1.224744871391589*fc[3]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[6]-1.224744871391589*fc[6]+0.7071067811865475*(fl[4]+fc[4]); 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 
  alpha[2] = E0[2]; 

  double amid = 0.7071067811865475*alpha[0]-0.7905694150420947*alpha[2]; 

  Ghat_r[0] = ((-0.7905694150420947*fr[5])+0.7905694150420947*fc[5]+0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*fr[0])*amax+0.3535533905932737*(fc[0]*amax+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = ((-0.7905694150420947*fr[7])+0.7905694150420947*fc[7]+0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.3162277660168379*(alpha[1]*favg_r[2]+favg_r[1]*alpha[2])+0.3535533905932737*(alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[6]+fc[6])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4])*amax+0.2258769757263128*alpha[2]*favg_r[2]+0.3535533905932737*(alpha[0]*favg_r[2]+favg_r[0]*alpha[2])+0.3162277660168379*alpha[1]*favg_r[1]; 

  Ghat_l[0] = (0.7905694150420947*fl[5]-0.7905694150420947*fc[5]+0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.3535533905932737*(alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.7905694150420947*fl[7]-0.7905694150420947*fc[7]+0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.3162277660168379*(alpha[1]*favg_l[2]+favg_l[1]*alpha[2])+0.3535533905932737*(alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[6]+fc[6])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4])*amax+0.2258769757263128*alpha[2]*favg_l[2]+0.3535533905932737*(alpha[0]*favg_l[2]+favg_l[0]*alpha[2])+0.3162277660168379*alpha[1]*favg_l[1]; 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[5] += (1.58113883008419*Ghat_l[0]-1.58113883008419*Ghat_r[0])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += (1.58113883008419*Ghat_l[1]-1.58113883008419*Ghat_r[1])*dv10; 

  return fabs(amid); 
} 
