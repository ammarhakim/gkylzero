#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvx_1x3v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* restrict out) 
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
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *E0 = &qmem[0]; 
  const double *B0 = &qmem[6]; 
  const double *B1 = &qmem[8]; 
  const double *B2 = &qmem[10]; 

  double Ghat_r[8]; 
  double favg_r[8]; 
  double Ghat_l[8]; 
  double favg_l[8]; 
  double alpha[8]; 

  favg_r[0] = (-1.224744871391589*fr[2])+1.224744871391589*fc[2]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[5])+1.224744871391589*fc[5]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[7])+1.224744871391589*fc[7]+0.7071067811865475*(fr[3]+fc[3]); 
  favg_r[3] = (-1.224744871391589*fr[9])+1.224744871391589*fc[9]+0.7071067811865475*(fr[4]+fc[4]); 
  favg_r[4] = (-1.224744871391589*fr[11])+1.224744871391589*fc[11]+0.7071067811865475*(fr[6]+fc[6]); 
  favg_r[5] = (-1.224744871391589*fr[12])+1.224744871391589*fc[12]+0.7071067811865475*(fr[8]+fc[8]); 
  favg_r[6] = (-1.224744871391589*fr[14])+1.224744871391589*fc[14]+0.7071067811865475*(fr[10]+fc[10]); 
  favg_r[7] = (-1.224744871391589*fr[15])+1.224744871391589*fc[15]+0.7071067811865475*(fr[13]+fc[13]); 

  favg_l[0] = 1.224744871391589*fl[2]-1.224744871391589*fc[2]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[5]-1.224744871391589*fc[5]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[7]-1.224744871391589*fc[7]+0.7071067811865475*(fl[3]+fc[3]); 
  favg_l[3] = 1.224744871391589*fl[9]-1.224744871391589*fc[9]+0.7071067811865475*(fl[4]+fc[4]); 
  favg_l[4] = 1.224744871391589*fl[11]-1.224744871391589*fc[11]+0.7071067811865475*(fl[6]+fc[6]); 
  favg_l[5] = 1.224744871391589*fl[12]-1.224744871391589*fc[12]+0.7071067811865475*(fl[8]+fc[8]); 
  favg_l[6] = 1.224744871391589*fl[14]-1.224744871391589*fc[14]+0.7071067811865475*(fl[10]+fc[10]); 
  favg_l[7] = 1.224744871391589*fl[15]-1.224744871391589*fc[15]+0.7071067811865475*(fl[13]+fc[13]); 

  alpha[0] = 2.0*(B2[0]*wv2+E0[0])-2.0*B1[0]*wv3; 
  alpha[1] = 2.0*(B2[1]*wv2+E0[1])-2.0*B1[1]*wv3; 
  alpha[2] = 0.5773502691896258*B2[0]*dv2; 
  alpha[3] = -0.5773502691896258*B1[0]*dv3; 
  alpha[4] = 0.5773502691896258*B2[1]*dv2; 
  alpha[5] = -0.5773502691896258*B1[1]*dv3; 

  double amid = 0.3535533905932737*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*fr[0]+0.3535533905932737*fc[0])*amax+0.1767766952966368*(alpha[5]*favg_r[5]+alpha[4]*favg_r[4]+alpha[3]*favg_r[3]+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[5]+fc[5])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.1767766952966368*(alpha[3]*favg_r[5]+favg_r[3]*alpha[5]+alpha[2]*favg_r[4]+favg_r[2]*alpha[4]+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[7]+fc[7])-0.3535533905932737*fr[3]+0.3535533905932737*fc[3])*amax+0.1767766952966368*(alpha[5]*favg_r[7]+alpha[3]*favg_r[6]+alpha[1]*favg_r[4]+favg_r[1]*alpha[4]+alpha[0]*favg_r[2]+favg_r[0]*alpha[2]); 
  Ghat_r[3] = (0.6123724356957944*(fr[9]+fc[9])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4])*amax+0.1767766952966368*(alpha[4]*favg_r[7]+alpha[2]*favg_r[6]+alpha[1]*favg_r[5]+favg_r[1]*alpha[5]+alpha[0]*favg_r[3]+favg_r[0]*alpha[3]); 
  Ghat_r[4] = (0.6123724356957944*(fr[11]+fc[11])-0.3535533905932737*fr[6]+0.3535533905932737*fc[6])*amax+0.1767766952966368*(alpha[3]*favg_r[7]+alpha[5]*favg_r[6]+alpha[0]*favg_r[4]+favg_r[0]*alpha[4]+alpha[1]*favg_r[2]+favg_r[1]*alpha[2]); 
  Ghat_r[5] = (0.6123724356957944*(fr[12]+fc[12])-0.3535533905932737*fr[8]+0.3535533905932737*fc[8])*amax+0.1767766952966368*(alpha[2]*favg_r[7]+alpha[4]*favg_r[6]+alpha[0]*favg_r[5]+favg_r[0]*alpha[5]+alpha[1]*favg_r[3]+favg_r[1]*alpha[3]); 
  Ghat_r[6] = (0.6123724356957944*(fr[14]+fc[14])-0.3535533905932737*fr[10]+0.3535533905932737*fc[10])*amax+0.1767766952966368*(alpha[1]*favg_r[7]+alpha[0]*favg_r[6]+alpha[4]*favg_r[5]+favg_r[4]*alpha[5]+alpha[2]*favg_r[3]+favg_r[2]*alpha[3]); 
  Ghat_r[7] = (0.6123724356957944*(fr[15]+fc[15])-0.3535533905932737*fr[13]+0.3535533905932737*fc[13])*amax+0.1767766952966368*(alpha[0]*favg_r[7]+alpha[1]*favg_r[6]+alpha[2]*favg_r[5]+favg_r[2]*alpha[5]+alpha[3]*favg_r[4]+favg_r[3]*alpha[4]); 

  Ghat_l[0] = (0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.1767766952966368*(alpha[5]*favg_l[5]+alpha[4]*favg_l[4]+alpha[3]*favg_l[3]+alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[5]+fc[5])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.1767766952966368*(alpha[3]*favg_l[5]+favg_l[3]*alpha[5]+alpha[2]*favg_l[4]+favg_l[2]*alpha[4]+alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[7]+fc[7])+0.3535533905932737*fl[3]-0.3535533905932737*fc[3])*amax+0.1767766952966368*(alpha[5]*favg_l[7]+alpha[3]*favg_l[6]+alpha[1]*favg_l[4]+favg_l[1]*alpha[4]+alpha[0]*favg_l[2]+favg_l[0]*alpha[2]); 
  Ghat_l[3] = (0.6123724356957944*(fl[9]+fc[9])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4])*amax+0.1767766952966368*(alpha[4]*favg_l[7]+alpha[2]*favg_l[6]+alpha[1]*favg_l[5]+favg_l[1]*alpha[5]+alpha[0]*favg_l[3]+favg_l[0]*alpha[3]); 
  Ghat_l[4] = (0.6123724356957944*(fl[11]+fc[11])+0.3535533905932737*fl[6]-0.3535533905932737*fc[6])*amax+0.1767766952966368*(alpha[3]*favg_l[7]+alpha[5]*favg_l[6]+alpha[0]*favg_l[4]+favg_l[0]*alpha[4]+alpha[1]*favg_l[2]+favg_l[1]*alpha[2]); 
  Ghat_l[5] = (0.6123724356957944*(fl[12]+fc[12])+0.3535533905932737*fl[8]-0.3535533905932737*fc[8])*amax+0.1767766952966368*(alpha[2]*favg_l[7]+alpha[4]*favg_l[6]+alpha[0]*favg_l[5]+favg_l[0]*alpha[5]+alpha[1]*favg_l[3]+favg_l[1]*alpha[3]); 
  Ghat_l[6] = (0.6123724356957944*(fl[14]+fc[14])+0.3535533905932737*fl[10]-0.3535533905932737*fc[10])*amax+0.1767766952966368*(alpha[1]*favg_l[7]+alpha[0]*favg_l[6]+alpha[4]*favg_l[5]+favg_l[4]*alpha[5]+alpha[2]*favg_l[3]+favg_l[2]*alpha[3]); 
  Ghat_l[7] = (0.6123724356957944*(fl[15]+fc[15])+0.3535533905932737*fl[13]-0.3535533905932737*fc[13])*amax+0.1767766952966368*(alpha[0]*favg_l[7]+alpha[1]*favg_l[6]+alpha[2]*favg_l[5]+favg_l[2]*alpha[5]+alpha[3]*favg_l[4]+favg_l[3]*alpha[4]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[9] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[12] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 

  return fabs(amid); 
} 
