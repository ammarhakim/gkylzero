#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH double vlasov_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* restrict out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv11 = 2/dxv[3]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double *E1 = &qmem[4]; 
  const double *B2 = &qmem[20]; 

  double Ghat_r[8]; 
  double favg_r[8]; 
  double Ghat_l[8]; 
  double favg_l[8]; 
  double alpha[8]; 

  favg_r[0] = (-1.224744871391589*fr[4])+1.224744871391589*fc[4]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[8])+1.224744871391589*fc[8]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[9])+1.224744871391589*fc[9]+0.7071067811865475*(fr[2]+fc[2]); 
  favg_r[3] = (-1.224744871391589*fr[10])+1.224744871391589*fc[10]+0.7071067811865475*(fr[3]+fc[3]); 
  favg_r[4] = (-1.224744871391589*fr[12])+1.224744871391589*fc[12]+0.7071067811865475*(fr[5]+fc[5]); 
  favg_r[5] = (-1.224744871391589*fr[13])+1.224744871391589*fc[13]+0.7071067811865475*(fr[6]+fc[6]); 
  favg_r[6] = (-1.224744871391589*fr[14])+1.224744871391589*fc[14]+0.7071067811865475*(fr[7]+fc[7]); 
  favg_r[7] = (-1.224744871391589*fr[15])+1.224744871391589*fc[15]+0.7071067811865475*(fr[11]+fc[11]); 

  favg_l[0] = 1.224744871391589*fl[4]-1.224744871391589*fc[4]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[8]-1.224744871391589*fc[8]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[9]-1.224744871391589*fc[9]+0.7071067811865475*(fl[2]+fc[2]); 
  favg_l[3] = 1.224744871391589*fl[10]-1.224744871391589*fc[10]+0.7071067811865475*(fl[3]+fc[3]); 
  favg_l[4] = 1.224744871391589*fl[12]-1.224744871391589*fc[12]+0.7071067811865475*(fl[5]+fc[5]); 
  favg_l[5] = 1.224744871391589*fl[13]-1.224744871391589*fc[13]+0.7071067811865475*(fl[6]+fc[6]); 
  favg_l[6] = 1.224744871391589*fl[14]-1.224744871391589*fc[14]+0.7071067811865475*(fl[7]+fc[7]); 
  favg_l[7] = 1.224744871391589*fl[15]-1.224744871391589*fc[15]+0.7071067811865475*(fl[11]+fc[11]); 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = 1.414213562373095*E1[2]-1.414213562373095*B2[2]*wv1; 
  alpha[3] = -0.408248290463863*B2[0]*dv1; 
  alpha[4] = 1.414213562373095*E1[3]-1.414213562373095*B2[3]*wv1; 
  alpha[5] = -0.408248290463863*B2[1]*dv1; 
  alpha[6] = -0.408248290463863*B2[2]*dv1; 
  alpha[7] = -0.408248290463863*B2[3]*dv1; 

  double amid = 0.3535533905932737*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[4]+fc[4])-0.3535533905932737*fr[0]+0.3535533905932737*fc[0])*amax+0.1767766952966368*(alpha[7]*favg_r[7]+alpha[6]*favg_r[6]+alpha[5]*favg_r[5]+alpha[4]*favg_r[4]+alpha[3]*favg_r[3]+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[8]+fc[8])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.1767766952966368*(alpha[6]*favg_r[7]+favg_r[6]*alpha[7]+alpha[3]*favg_r[5]+favg_r[3]*alpha[5]+alpha[2]*favg_r[4]+favg_r[2]*alpha[4]+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[9]+fc[9])-0.3535533905932737*fr[2]+0.3535533905932737*fc[2])*amax+0.1767766952966368*(alpha[5]*favg_r[7]+favg_r[5]*alpha[7]+alpha[3]*favg_r[6]+favg_r[3]*alpha[6]+alpha[1]*favg_r[4]+favg_r[1]*alpha[4]+alpha[0]*favg_r[2]+favg_r[0]*alpha[2]); 
  Ghat_r[3] = (0.6123724356957944*(fr[10]+fc[10])-0.3535533905932737*fr[3]+0.3535533905932737*fc[3])*amax+0.1767766952966368*(alpha[4]*favg_r[7]+favg_r[4]*alpha[7]+alpha[2]*favg_r[6]+favg_r[2]*alpha[6]+alpha[1]*favg_r[5]+favg_r[1]*alpha[5]+alpha[0]*favg_r[3]+favg_r[0]*alpha[3]); 
  Ghat_r[4] = (0.6123724356957944*(fr[12]+fc[12])-0.3535533905932737*fr[5]+0.3535533905932737*fc[5])*amax+0.1767766952966368*(alpha[3]*favg_r[7]+favg_r[3]*alpha[7]+alpha[5]*favg_r[6]+favg_r[5]*alpha[6]+alpha[0]*favg_r[4]+favg_r[0]*alpha[4]+alpha[1]*favg_r[2]+favg_r[1]*alpha[2]); 
  Ghat_r[5] = (0.6123724356957944*(fr[13]+fc[13])-0.3535533905932737*fr[6]+0.3535533905932737*fc[6])*amax+0.1767766952966368*(alpha[2]*favg_r[7]+favg_r[2]*alpha[7]+alpha[4]*favg_r[6]+favg_r[4]*alpha[6]+alpha[0]*favg_r[5]+favg_r[0]*alpha[5]+alpha[1]*favg_r[3]+favg_r[1]*alpha[3]); 
  Ghat_r[6] = (0.6123724356957944*(fr[14]+fc[14])-0.3535533905932737*fr[7]+0.3535533905932737*fc[7])*amax+0.1767766952966368*(alpha[1]*favg_r[7]+favg_r[1]*alpha[7]+alpha[0]*favg_r[6]+favg_r[0]*alpha[6]+alpha[4]*favg_r[5]+favg_r[4]*alpha[5]+alpha[2]*favg_r[3]+favg_r[2]*alpha[3]); 
  Ghat_r[7] = (0.6123724356957944*(fr[15]+fc[15])-0.3535533905932737*fr[11]+0.3535533905932737*fc[11])*amax+0.1767766952966368*(alpha[0]*favg_r[7]+favg_r[0]*alpha[7]+alpha[1]*favg_r[6]+favg_r[1]*alpha[6]+alpha[2]*favg_r[5]+favg_r[2]*alpha[5]+alpha[3]*favg_r[4]+favg_r[3]*alpha[4]); 

  Ghat_l[0] = (0.6123724356957944*(fl[4]+fc[4])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.1767766952966368*(alpha[7]*favg_l[7]+alpha[6]*favg_l[6]+alpha[5]*favg_l[5]+alpha[4]*favg_l[4]+alpha[3]*favg_l[3]+alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[8]+fc[8])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.1767766952966368*(alpha[6]*favg_l[7]+favg_l[6]*alpha[7]+alpha[3]*favg_l[5]+favg_l[3]*alpha[5]+alpha[2]*favg_l[4]+favg_l[2]*alpha[4]+alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[9]+fc[9])+0.3535533905932737*fl[2]-0.3535533905932737*fc[2])*amax+0.1767766952966368*(alpha[5]*favg_l[7]+favg_l[5]*alpha[7]+alpha[3]*favg_l[6]+favg_l[3]*alpha[6]+alpha[1]*favg_l[4]+favg_l[1]*alpha[4]+alpha[0]*favg_l[2]+favg_l[0]*alpha[2]); 
  Ghat_l[3] = (0.6123724356957944*(fl[10]+fc[10])+0.3535533905932737*fl[3]-0.3535533905932737*fc[3])*amax+0.1767766952966368*(alpha[4]*favg_l[7]+favg_l[4]*alpha[7]+alpha[2]*favg_l[6]+favg_l[2]*alpha[6]+alpha[1]*favg_l[5]+favg_l[1]*alpha[5]+alpha[0]*favg_l[3]+favg_l[0]*alpha[3]); 
  Ghat_l[4] = (0.6123724356957944*(fl[12]+fc[12])+0.3535533905932737*fl[5]-0.3535533905932737*fc[5])*amax+0.1767766952966368*(alpha[3]*favg_l[7]+favg_l[3]*alpha[7]+alpha[5]*favg_l[6]+favg_l[5]*alpha[6]+alpha[0]*favg_l[4]+favg_l[0]*alpha[4]+alpha[1]*favg_l[2]+favg_l[1]*alpha[2]); 
  Ghat_l[5] = (0.6123724356957944*(fl[13]+fc[13])+0.3535533905932737*fl[6]-0.3535533905932737*fc[6])*amax+0.1767766952966368*(alpha[2]*favg_l[7]+favg_l[2]*alpha[7]+alpha[4]*favg_l[6]+favg_l[4]*alpha[6]+alpha[0]*favg_l[5]+favg_l[0]*alpha[5]+alpha[1]*favg_l[3]+favg_l[1]*alpha[3]); 
  Ghat_l[6] = (0.6123724356957944*(fl[14]+fc[14])+0.3535533905932737*fl[7]-0.3535533905932737*fc[7])*amax+0.1767766952966368*(alpha[1]*favg_l[7]+favg_l[1]*alpha[7]+alpha[0]*favg_l[6]+favg_l[0]*alpha[6]+alpha[4]*favg_l[5]+favg_l[4]*alpha[5]+alpha[2]*favg_l[3]+favg_l[2]*alpha[3]); 
  Ghat_l[7] = (0.6123724356957944*(fl[15]+fc[15])+0.3535533905932737*fl[11]-0.3535533905932737*fc[11])*amax+0.1767766952966368*(alpha[0]*favg_l[7]+favg_l[0]*alpha[7]+alpha[1]*favg_l[6]+favg_l[1]*alpha[6]+alpha[2]*favg_l[5]+favg_l[2]*alpha[5]+alpha[3]*favg_l[4]+favg_l[3]*alpha[4]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[4] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[7] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[8] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[9] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[11] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[12] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 

  return fabs(amid); 
} 
