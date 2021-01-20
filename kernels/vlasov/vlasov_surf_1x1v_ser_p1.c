#include <gkyl_vlasov_kernels.h> 
void vlasov_surfx_1x1v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* restrict out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[2]; 
  double Ghat_l[2]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[3]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[3]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[3]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[3]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 

  } else { 

  Ghat_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[1])*wv+(0.2041241452319315*fr[2]-0.3535533905932737*fr[3])*dv; 
  Ghat_r[1] = (0.7071067811865475*fr[2]-1.224744871391589*fr[3])*wv+(0.2041241452319315*fr[0]-0.3535533905932737*fr[1])*dv; 

  Ghat_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[1])*wv+(0.2041241452319315*fc[2]-0.3535533905932737*fc[3])*dv; 
  Ghat_l[1] = (0.7071067811865475*fc[2]-1.224744871391589*fc[3])*wv+(0.2041241452319315*fc[0]-0.3535533905932737*fc[1])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
} 
double vlasov_surfvx_1x1v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* restrict out) 
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

  double Ghat_r[2]; 
  double favg_r[2]; 
  double Ghat_l[2]; 
  double favg_l[2]; 
  double alpha[2]; 

  favg_r[0] = (-1.224744871391589*fr[2])+1.224744871391589*fc[2]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[3])+1.224744871391589*fc[3]+0.7071067811865475*(fr[1]+fc[1]); 

  favg_l[0] = 1.224744871391589*fl[2]-1.224744871391589*fc[2]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[3]-1.224744871391589*fc[3]+0.7071067811865475*(fl[1]+fc[1]); 

  alpha[0] = E0[0]; 
  alpha[1] = E0[1]; 

  double amid = 0.7071067811865475*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*fr[0])*amax+0.3535533905932737*(fc[0]*amax+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*fr[1])*amax+0.3535533905932737*(fc[1]*amax+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 

  Ghat_l[0] = (0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.3535533905932737*(alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.3535533905932737*(alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 

  return fabs(amid); 
} 
