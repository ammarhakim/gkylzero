#include <gkyl_vlasov_kernels.h> 
GKYL_CU_DH void vlasov_surfx_1x2v_ser_p1(const double *w, const double *dxv, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]:   Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fl/fc/fr:  Input Distribution function in left/center/right cells.
  // out:       Incremented distribution function in center cell.
  const double dx10 = 2/dxv[0]; 
  const double dv = dxv[1], wv = w[1]; 
  double Ghat_r[4]; 
  double Ghat_l[4]; 
  if (wv>0) { 

  Ghat_r[0] = (1.224744871391589*fc[1]+0.7071067811865475*fc[0])*wv+(0.3535533905932737*fc[4]+0.2041241452319315*fc[2])*dv; 
  Ghat_r[1] = (1.224744871391589*fc[4]+0.7071067811865475*fc[2])*wv+(0.3535533905932737*fc[1]+0.2041241452319315*fc[0])*dv; 
  Ghat_r[2] = (1.224744871391589*fc[5]+0.7071067811865475*fc[3])*wv+(0.3535533905932737*fc[7]+0.2041241452319315*fc[6])*dv; 
  Ghat_r[3] = (1.224744871391589*fc[7]+0.7071067811865475*fc[6])*wv+(0.3535533905932737*fc[5]+0.2041241452319315*fc[3])*dv; 

  Ghat_l[0] = (1.224744871391589*fl[1]+0.7071067811865475*fl[0])*wv+(0.3535533905932737*fl[4]+0.2041241452319315*fl[2])*dv; 
  Ghat_l[1] = (1.224744871391589*fl[4]+0.7071067811865475*fl[2])*wv+(0.3535533905932737*fl[1]+0.2041241452319315*fl[0])*dv; 
  Ghat_l[2] = (1.224744871391589*fl[5]+0.7071067811865475*fl[3])*wv+(0.3535533905932737*fl[7]+0.2041241452319315*fl[6])*dv; 
  Ghat_l[3] = (1.224744871391589*fl[7]+0.7071067811865475*fl[6])*wv+(0.3535533905932737*fl[5]+0.2041241452319315*fl[3])*dv; 

  } else { 

  Ghat_r[0] = (0.7071067811865475*fr[0]-1.224744871391589*fr[1])*wv+(0.2041241452319315*fr[2]-0.3535533905932737*fr[4])*dv; 
  Ghat_r[1] = (0.7071067811865475*fr[2]-1.224744871391589*fr[4])*wv+(0.2041241452319315*fr[0]-0.3535533905932737*fr[1])*dv; 
  Ghat_r[2] = (0.7071067811865475*fr[3]-1.224744871391589*fr[5])*wv+(0.2041241452319315*fr[6]-0.3535533905932737*fr[7])*dv; 
  Ghat_r[3] = (0.7071067811865475*fr[6]-1.224744871391589*fr[7])*wv+(0.2041241452319315*fr[3]-0.3535533905932737*fr[5])*dv; 

  Ghat_l[0] = (0.7071067811865475*fc[0]-1.224744871391589*fc[1])*wv+(0.2041241452319315*fc[2]-0.3535533905932737*fc[4])*dv; 
  Ghat_l[1] = (0.7071067811865475*fc[2]-1.224744871391589*fc[4])*wv+(0.2041241452319315*fc[0]-0.3535533905932737*fc[1])*dv; 
  Ghat_l[2] = (0.7071067811865475*fc[3]-1.224744871391589*fc[5])*wv+(0.2041241452319315*fc[6]-0.3535533905932737*fc[7])*dv; 
  Ghat_l[3] = (0.7071067811865475*fc[6]-1.224744871391589*fc[7])*wv+(0.2041241452319315*fc[3]-0.3535533905932737*fc[5])*dv; 

  } 
  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dx10; 
  out[1] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dx10; 
  out[2] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dx10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dx10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dx10; 
  out[5] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dx10; 
  out[6] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dx10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dx10; 
} 
GKYL_CU_DH double vlasov_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
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
  const double *E0 = &qmem[0]; 
  const double *B2 = &qmem[10]; 

  double Ghat_r[4]; 
  double favg_r[4]; 
  double Ghat_l[4]; 
  double favg_l[4]; 
  double alpha[4]; 

  favg_r[0] = (-1.224744871391589*fr[2])+1.224744871391589*fc[2]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[4])+1.224744871391589*fc[4]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[6])+1.224744871391589*fc[6]+0.7071067811865475*(fr[3]+fc[3]); 
  favg_r[3] = (-1.224744871391589*fr[7])+1.224744871391589*fc[7]+0.7071067811865475*(fr[5]+fc[5]); 

  favg_l[0] = 1.224744871391589*fl[2]-1.224744871391589*fc[2]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[4]-1.224744871391589*fc[4]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[6]-1.224744871391589*fc[6]+0.7071067811865475*(fl[3]+fc[3]); 
  favg_l[3] = 1.224744871391589*fl[7]-1.224744871391589*fc[7]+0.7071067811865475*(fl[5]+fc[5]); 

  alpha[0] = 1.414213562373095*(B2[0]*wv2+E0[0]); 
  alpha[1] = 1.414213562373095*(B2[1]*wv2+E0[1]); 
  alpha[2] = 0.408248290463863*B2[0]*dv2; 
  alpha[3] = 0.408248290463863*B2[1]*dv2; 

  double amid = 0.5*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*fr[0]+0.3535533905932737*fc[0])*amax+0.25*(alpha[3]*favg_r[3]+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[4]+fc[4])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.25*(alpha[2]*favg_r[3]+favg_r[2]*alpha[3]+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[6]+fc[6])-0.3535533905932737*fr[3]+0.3535533905932737*fc[3])*amax+0.25*(alpha[1]*favg_r[3]+favg_r[1]*alpha[3]+alpha[0]*favg_r[2]+favg_r[0]*alpha[2]); 
  Ghat_r[3] = (0.6123724356957944*(fr[7]+fc[7])-0.3535533905932737*fr[5]+0.3535533905932737*fc[5])*amax+0.25*(alpha[0]*favg_r[3]+favg_r[0]*alpha[3]+alpha[1]*favg_r[2]+favg_r[1]*alpha[2]); 

  Ghat_l[0] = (0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.25*(alpha[3]*favg_l[3]+alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[4]+fc[4])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.25*(alpha[2]*favg_l[3]+favg_l[2]*alpha[3]+alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[6]+fc[6])+0.3535533905932737*fl[3]-0.3535533905932737*fc[3])*amax+0.25*(alpha[1]*favg_l[3]+favg_l[1]*alpha[3]+alpha[0]*favg_l[2]+favg_l[0]*alpha[2]); 
  Ghat_l[3] = (0.6123724356957944*(fl[7]+fc[7])+0.3535533905932737*fl[5]-0.3535533905932737*fc[5])*amax+0.25*(alpha[0]*favg_l[3]+favg_l[0]*alpha[3]+alpha[1]*favg_l[2]+favg_l[1]*alpha[2]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[3] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[4] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 

  return fabs(amid); 
} 
GKYL_CU_DH double vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double amax, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // amax:      amax in global lax flux.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  // returns abs(amid) for use in determining amax in cfl and global lax flux 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 

  double Ghat_r[4]; 
  double favg_r[4]; 
  double Ghat_l[4]; 
  double favg_l[4]; 
  double alpha[4]; 

  favg_r[0] = (-1.224744871391589*fr[3])+1.224744871391589*fc[3]+0.7071067811865475*(fr[0]+fc[0]); 
  favg_r[1] = (-1.224744871391589*fr[5])+1.224744871391589*fc[5]+0.7071067811865475*(fr[1]+fc[1]); 
  favg_r[2] = (-1.224744871391589*fr[6])+1.224744871391589*fc[6]+0.7071067811865475*(fr[2]+fc[2]); 
  favg_r[3] = (-1.224744871391589*fr[7])+1.224744871391589*fc[7]+0.7071067811865475*(fr[4]+fc[4]); 

  favg_l[0] = 1.224744871391589*fl[3]-1.224744871391589*fc[3]+0.7071067811865475*(fl[0]+fc[0]); 
  favg_l[1] = 1.224744871391589*fl[5]-1.224744871391589*fc[5]+0.7071067811865475*(fl[1]+fc[1]); 
  favg_l[2] = 1.224744871391589*fl[6]-1.224744871391589*fc[6]+0.7071067811865475*(fl[2]+fc[2]); 
  favg_l[3] = 1.224744871391589*fl[7]-1.224744871391589*fc[7]+0.7071067811865475*(fl[4]+fc[4]); 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 

  double amid = 0.5*alpha[0]; 

  Ghat_r[0] = (0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*fr[0]+0.3535533905932737*fc[0])*amax+0.25*(alpha[3]*favg_r[3]+alpha[2]*favg_r[2]+alpha[1]*favg_r[1]+alpha[0]*favg_r[0]); 
  Ghat_r[1] = (0.6123724356957944*(fr[5]+fc[5])-0.3535533905932737*fr[1]+0.3535533905932737*fc[1])*amax+0.25*(alpha[2]*favg_r[3]+favg_r[2]*alpha[3]+alpha[0]*favg_r[1]+favg_r[0]*alpha[1]); 
  Ghat_r[2] = (0.6123724356957944*(fr[6]+fc[6])-0.3535533905932737*fr[2]+0.3535533905932737*fc[2])*amax+0.25*(alpha[1]*favg_r[3]+favg_r[1]*alpha[3]+alpha[0]*favg_r[2]+favg_r[0]*alpha[2]); 
  Ghat_r[3] = (0.6123724356957944*(fr[7]+fc[7])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4])*amax+0.25*(alpha[0]*favg_r[3]+favg_r[0]*alpha[3]+alpha[1]*favg_r[2]+favg_r[1]*alpha[2]); 

  Ghat_l[0] = (0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*fl[0]-0.3535533905932737*fc[0])*amax+0.25*(alpha[3]*favg_l[3]+alpha[2]*favg_l[2]+alpha[1]*favg_l[1]+alpha[0]*favg_l[0]); 
  Ghat_l[1] = (0.6123724356957944*(fl[5]+fc[5])+0.3535533905932737*fl[1]-0.3535533905932737*fc[1])*amax+0.25*(alpha[2]*favg_l[3]+favg_l[2]*alpha[3]+alpha[0]*favg_l[1]+favg_l[0]*alpha[1]); 
  Ghat_l[2] = (0.6123724356957944*(fl[6]+fc[6])+0.3535533905932737*fl[2]-0.3535533905932737*fc[2])*amax+0.25*(alpha[1]*favg_l[3]+favg_l[1]*alpha[3]+alpha[0]*favg_l[2]+favg_l[0]*alpha[2]); 
  Ghat_l[3] = (0.6123724356957944*(fl[7]+fc[7])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4])*amax+0.25*(alpha[0]*favg_l[3]+favg_l[0]*alpha[3]+alpha[1]*favg_l[2]+favg_l[1]*alpha[2]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 

  return fabs(amid); 
} 
