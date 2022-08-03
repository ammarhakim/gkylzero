#include <gkyl_lbo_gyrokinetic_kernels.h> 
GKYL_CU_DH void lbo_gyrokinetic_drag_surfmu_1x2v_ser_p2(const double *w, const double *dxv, const double m_, const double *bmag_inv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:     cell-center coordinates. 
  // dxv[3]:   cell spacing. 
  // m_:        species mass.
  // bmag_inv:  1/(magnetic field magnitude). 
  // nuSum:     collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:  distribution function in cells 
  // out:       incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 

  double Ghat_r[8] = {0.0}; 
  double Ghat_l[8] = {0.0}; 
  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = 2.828427124746191*nuSum[0]*w[2]-1.414213562373095*nuSum[0]*dxv[2]; 
  alphaDrSurf_l[1] = 2.828427124746191*nuSum[1]*w[2]-1.414213562373095*nuSum[1]*dxv[2]; 
  alphaDrSurf_l[4] = 2.828427124746191*nuSum[2]*w[2]-1.414213562373095*dxv[2]*nuSum[2]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 2.828427124746191*nuSum[0]*w[2]+1.414213562373095*nuSum[0]*dxv[2]; 
  alphaDrSurf_r[1] = 2.828427124746191*nuSum[1]*w[2]+1.414213562373095*nuSum[1]*dxv[2]; 
  alphaDrSurf_r[4] = 2.828427124746191*nuSum[2]*w[2]+1.414213562373095*dxv[2]*nuSum[2]; 

  Ghat_l[0] = 0.7905694150420948*alphaDrSurf_l[1]*fc[15]-0.6123724356957944*alphaDrSurf_l[4]*fc[13]+0.7905694150420947*alphaDrSurf_l[0]*fc[9]+0.3535533905932737*alphaDrSurf_l[4]*fc[7]-0.6123724356957944*alphaDrSurf_l[1]*fc[5]-0.6123724356957944*alphaDrSurf_l[0]*fc[3]+0.3535533905932737*alphaDrSurf_l[1]*fc[1]+0.3535533905932737*alphaDrSurf_l[0]*fc[0]; 
  Ghat_l[1] = 0.7071067811865475*alphaDrSurf_l[4]*fc[15]+0.7905694150420948*alphaDrSurf_l[0]*fc[15]-0.5477225575051661*alphaDrSurf_l[1]*fc[13]+0.7905694150420947*alphaDrSurf_l[1]*fc[9]+0.3162277660168379*alphaDrSurf_l[1]*fc[7]-0.5477225575051661*alphaDrSurf_l[4]*fc[5]-0.6123724356957944*alphaDrSurf_l[0]*fc[5]+0.3162277660168379*fc[1]*alphaDrSurf_l[4]-0.6123724356957944*alphaDrSurf_l[1]*fc[3]+0.3535533905932737*alphaDrSurf_l[0]*fc[1]+0.3535533905932737*fc[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = 0.7905694150420947*alphaDrSurf_l[1]*fc[19]-0.6123724356957944*alphaDrSurf_l[4]*fc[17]+0.7905694150420948*alphaDrSurf_l[0]*fc[16]+0.3535533905932737*alphaDrSurf_l[4]*fc[11]-0.6123724356957944*alphaDrSurf_l[1]*fc[10]-0.6123724356957944*alphaDrSurf_l[0]*fc[6]+0.3535533905932737*alphaDrSurf_l[1]*fc[4]+0.3535533905932737*alphaDrSurf_l[0]*fc[2]; 
  Ghat_l[3] = 0.7071067811865475*alphaDrSurf_l[4]*fc[19]+0.7905694150420947*alphaDrSurf_l[0]*fc[19]-0.5477225575051661*alphaDrSurf_l[1]*fc[17]+0.7905694150420948*alphaDrSurf_l[1]*fc[16]+0.3162277660168379*alphaDrSurf_l[1]*fc[11]-0.5477225575051661*alphaDrSurf_l[4]*fc[10]-0.6123724356957944*alphaDrSurf_l[0]*fc[10]-0.6123724356957944*alphaDrSurf_l[1]*fc[6]+0.3162277660168379*alphaDrSurf_l[4]*fc[4]+0.3535533905932737*alphaDrSurf_l[0]*fc[4]+0.3535533905932737*alphaDrSurf_l[1]*fc[2]; 
  Ghat_l[4] = 0.7071067811865475*alphaDrSurf_l[1]*fc[15]-0.3912303982179757*alphaDrSurf_l[4]*fc[13]-0.6123724356957944*alphaDrSurf_l[0]*fc[13]+0.7905694150420947*alphaDrSurf_l[4]*fc[9]+0.2258769757263128*alphaDrSurf_l[4]*fc[7]+0.3535533905932737*alphaDrSurf_l[0]*fc[7]-0.5477225575051661*alphaDrSurf_l[1]*fc[5]-0.6123724356957944*fc[3]*alphaDrSurf_l[4]+0.3535533905932737*fc[0]*alphaDrSurf_l[4]+0.3162277660168379*alphaDrSurf_l[1]*fc[1]; 
  Ghat_l[5] = (-0.6123724356957944*alphaDrSurf_l[1]*fc[18])-0.6123724356957944*alphaDrSurf_l[0]*fc[14]+0.3535533905932737*alphaDrSurf_l[1]*fc[12]+0.3535533905932737*alphaDrSurf_l[0]*fc[8]; 
  Ghat_l[6] = 0.7071067811865475*alphaDrSurf_l[1]*fc[19]-0.3912303982179757*alphaDrSurf_l[4]*fc[17]-0.6123724356957944*alphaDrSurf_l[0]*fc[17]+0.7905694150420947*alphaDrSurf_l[4]*fc[16]+0.2258769757263128*alphaDrSurf_l[4]*fc[11]+0.3535533905932737*alphaDrSurf_l[0]*fc[11]-0.5477225575051661*alphaDrSurf_l[1]*fc[10]-0.6123724356957944*alphaDrSurf_l[4]*fc[6]+0.3162277660168379*alphaDrSurf_l[1]*fc[4]+0.3535533905932737*fc[2]*alphaDrSurf_l[4]; 
  Ghat_l[7] = (-0.5477225575051661*alphaDrSurf_l[4]*fc[18])-0.6123724356957944*alphaDrSurf_l[0]*fc[18]-0.6123724356957944*alphaDrSurf_l[1]*fc[14]+0.3162277660168379*alphaDrSurf_l[4]*fc[12]+0.3535533905932737*alphaDrSurf_l[0]*fc[12]+0.3535533905932737*alphaDrSurf_l[1]*fc[8]; 

  Ghat_r[0] = 0.7905694150420948*alphaDrSurf_r[1]*fr[15]-0.6123724356957944*alphaDrSurf_r[4]*fr[13]+0.7905694150420947*alphaDrSurf_r[0]*fr[9]+0.3535533905932737*alphaDrSurf_r[4]*fr[7]-0.6123724356957944*alphaDrSurf_r[1]*fr[5]-0.6123724356957944*alphaDrSurf_r[0]*fr[3]+0.3535533905932737*alphaDrSurf_r[1]*fr[1]+0.3535533905932737*alphaDrSurf_r[0]*fr[0]; 
  Ghat_r[1] = 0.7071067811865475*alphaDrSurf_r[4]*fr[15]+0.7905694150420948*alphaDrSurf_r[0]*fr[15]-0.5477225575051661*alphaDrSurf_r[1]*fr[13]+0.7905694150420947*alphaDrSurf_r[1]*fr[9]+0.3162277660168379*alphaDrSurf_r[1]*fr[7]-0.5477225575051661*alphaDrSurf_r[4]*fr[5]-0.6123724356957944*alphaDrSurf_r[0]*fr[5]+0.3162277660168379*fr[1]*alphaDrSurf_r[4]-0.6123724356957944*alphaDrSurf_r[1]*fr[3]+0.3535533905932737*alphaDrSurf_r[0]*fr[1]+0.3535533905932737*fr[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = 0.7905694150420947*alphaDrSurf_r[1]*fr[19]-0.6123724356957944*alphaDrSurf_r[4]*fr[17]+0.7905694150420948*alphaDrSurf_r[0]*fr[16]+0.3535533905932737*alphaDrSurf_r[4]*fr[11]-0.6123724356957944*alphaDrSurf_r[1]*fr[10]-0.6123724356957944*alphaDrSurf_r[0]*fr[6]+0.3535533905932737*alphaDrSurf_r[1]*fr[4]+0.3535533905932737*alphaDrSurf_r[0]*fr[2]; 
  Ghat_r[3] = 0.7071067811865475*alphaDrSurf_r[4]*fr[19]+0.7905694150420947*alphaDrSurf_r[0]*fr[19]-0.5477225575051661*alphaDrSurf_r[1]*fr[17]+0.7905694150420948*alphaDrSurf_r[1]*fr[16]+0.3162277660168379*alphaDrSurf_r[1]*fr[11]-0.5477225575051661*alphaDrSurf_r[4]*fr[10]-0.6123724356957944*alphaDrSurf_r[0]*fr[10]-0.6123724356957944*alphaDrSurf_r[1]*fr[6]+0.3162277660168379*alphaDrSurf_r[4]*fr[4]+0.3535533905932737*alphaDrSurf_r[0]*fr[4]+0.3535533905932737*alphaDrSurf_r[1]*fr[2]; 
  Ghat_r[4] = 0.7071067811865475*alphaDrSurf_r[1]*fr[15]-0.3912303982179757*alphaDrSurf_r[4]*fr[13]-0.6123724356957944*alphaDrSurf_r[0]*fr[13]+0.7905694150420947*alphaDrSurf_r[4]*fr[9]+0.2258769757263128*alphaDrSurf_r[4]*fr[7]+0.3535533905932737*alphaDrSurf_r[0]*fr[7]-0.5477225575051661*alphaDrSurf_r[1]*fr[5]-0.6123724356957944*fr[3]*alphaDrSurf_r[4]+0.3535533905932737*fr[0]*alphaDrSurf_r[4]+0.3162277660168379*alphaDrSurf_r[1]*fr[1]; 
  Ghat_r[5] = (-0.6123724356957944*alphaDrSurf_r[1]*fr[18])-0.6123724356957944*alphaDrSurf_r[0]*fr[14]+0.3535533905932737*alphaDrSurf_r[1]*fr[12]+0.3535533905932737*alphaDrSurf_r[0]*fr[8]; 
  Ghat_r[6] = 0.7071067811865475*alphaDrSurf_r[1]*fr[19]-0.3912303982179757*alphaDrSurf_r[4]*fr[17]-0.6123724356957944*alphaDrSurf_r[0]*fr[17]+0.7905694150420947*alphaDrSurf_r[4]*fr[16]+0.2258769757263128*alphaDrSurf_r[4]*fr[11]+0.3535533905932737*alphaDrSurf_r[0]*fr[11]-0.5477225575051661*alphaDrSurf_r[1]*fr[10]-0.6123724356957944*alphaDrSurf_r[4]*fr[6]+0.3162277660168379*alphaDrSurf_r[1]*fr[4]+0.3535533905932737*fr[2]*alphaDrSurf_r[4]; 
  Ghat_r[7] = (-0.5477225575051661*alphaDrSurf_r[4]*fr[18])-0.6123724356957944*alphaDrSurf_r[0]*fr[18]-0.6123724356957944*alphaDrSurf_r[1]*fr[14]+0.3162277660168379*alphaDrSurf_r[4]*fr[12]+0.3535533905932737*alphaDrSurf_r[0]*fr[12]+0.3535533905932737*alphaDrSurf_r[1]*fr[8]; 

  out[0] += (0.7071067811865475*Ghat_r[0]-0.7071067811865475*Ghat_l[0])*rdv2; 
  out[1] += (0.7071067811865475*Ghat_r[1]-0.7071067811865475*Ghat_l[1])*rdv2; 
  out[2] += (0.7071067811865475*Ghat_r[2]-0.7071067811865475*Ghat_l[2])*rdv2; 
  out[3] += 1.224744871391589*(Ghat_r[0]+Ghat_l[0])*rdv2; 
  out[4] += (0.7071067811865475*Ghat_r[3]-0.7071067811865475*Ghat_l[3])*rdv2; 
  out[5] += 1.224744871391589*(Ghat_r[1]+Ghat_l[1])*rdv2; 
  out[6] += 1.224744871391589*(Ghat_r[2]+Ghat_l[2])*rdv2; 
  out[7] += (0.7071067811865475*Ghat_r[4]-0.7071067811865475*Ghat_l[4])*rdv2; 
  out[8] += (0.7071067811865475*Ghat_r[5]-0.7071067811865475*Ghat_l[5])*rdv2; 
  out[9] += (1.58113883008419*Ghat_r[0]-1.58113883008419*Ghat_l[0])*rdv2; 
  out[10] += 1.224744871391589*(Ghat_r[3]+Ghat_l[3])*rdv2; 
  out[11] += (0.7071067811865475*Ghat_r[6]-0.7071067811865475*Ghat_l[6])*rdv2; 
  out[12] += (0.7071067811865475*Ghat_r[7]-0.7071067811865475*Ghat_l[7])*rdv2; 
  out[13] += 1.224744871391589*(Ghat_r[4]+Ghat_l[4])*rdv2; 
  out[14] += 1.224744871391589*(Ghat_r[5]+Ghat_l[5])*rdv2; 
  out[15] += (1.58113883008419*Ghat_r[1]-1.58113883008419*Ghat_l[1])*rdv2; 
  out[16] += (1.58113883008419*Ghat_r[2]-1.58113883008419*Ghat_l[2])*rdv2; 
  out[17] += 1.224744871391589*(Ghat_r[6]+Ghat_l[6])*rdv2; 
  out[18] += 1.224744871391589*(Ghat_r[7]+Ghat_l[7])*rdv2; 
  out[19] += (1.58113883008419*Ghat_r[3]-1.58113883008419*Ghat_l[3])*rdv2; 
} 
