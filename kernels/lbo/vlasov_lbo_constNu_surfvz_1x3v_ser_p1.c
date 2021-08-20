#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_surfvz_1x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *sumNuUz_l = &nuUSum_l[4]; 
  const double *sumNuUz_r = &nuUSum_r[4]; 

  double alphaDrSurf_l[8]; 
  alphaDrSurf_l[0] = (2.82842712474619*w[3]+1.414213562373095*dxv[3])*nuSum-2.0*sumNuUz_l[0]; 
  alphaDrSurf_l[1] = -2.0*sumNuUz_l[1]; 

  double alphaDrSurf_r[8]; 
  alphaDrSurf_r[0] = ((-2.82842712474619*w[3])-1.414213562373095*dxv[3])*nuSum+2.0*sumNuUz_r[0]; 
  alphaDrSurf_r[1] = 2.0*sumNuUz_r[1]; 

  double fUpwindQuad_l[8];
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = (-0.4330127018922193*fl[15])+0.4330127018922193*(fl[14]+fl[13]+fl[12])-0.25*fl[11]-0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*fl[4]-0.25*(fl[3]+fl[2]+fl[1])+0.25*fl[0]; 
  } else { 
    fUpwindQuad_l[0] = 0.4330127018922193*fc[15]-0.4330127018922193*(fc[14]+fc[13]+fc[12])-0.25*fc[11]+0.4330127018922193*(fc[10]+fc[9]+fc[8])+0.25*(fc[7]+fc[6]+fc[5])-0.4330127018922193*fc[4]-0.25*(fc[3]+fc[2]+fc[1])+0.25*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = 0.4330127018922193*(fl[15]+fl[14])-0.4330127018922193*(fl[13]+fl[12])+0.25*fl[11]-0.4330127018922193*(fl[10]+fl[9])+0.4330127018922193*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.4330127018922193*fl[4]-0.25*(fl[3]+fl[2])+0.25*(fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[1] = (-0.4330127018922193*(fc[15]+fc[14]))+0.4330127018922193*(fc[13]+fc[12])+0.25*fc[11]+0.4330127018922193*(fc[10]+fc[9])-0.4330127018922193*fc[8]+0.25*fc[7]-0.25*(fc[6]+fc[5])-0.4330127018922193*fc[4]-0.25*(fc[3]+fc[2])+0.25*(fc[1]+fc[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[2] = 0.4330127018922193*fl[15]-0.4330127018922193*fl[14]+0.4330127018922193*fl[13]-0.4330127018922193*fl[12]+0.25*fl[11]-0.4330127018922193*fl[10]+0.4330127018922193*fl[9]-0.4330127018922193*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.4330127018922193*fl[4]-0.25*fl[3]+0.25*fl[2]-0.25*fl[1]+0.25*fl[0]; 
  } else { 
    fUpwindQuad_l[2] = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.4330127018922193*fc[13]+0.4330127018922193*fc[12]+0.25*fc[11]+0.4330127018922193*fc[10]-0.4330127018922193*fc[9]+0.4330127018922193*fc[8]-0.25*fc[7]+0.25*fc[6]-0.25*fc[5]-0.4330127018922193*fc[4]-0.25*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = (-0.4330127018922193*(fl[15]+fl[14]+fl[13]))+0.4330127018922193*fl[12]-0.25*fl[11]-0.4330127018922193*fl[10]+0.4330127018922193*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]+0.4330127018922193*fl[4]-0.25*fl[3]+0.25*(fl[2]+fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[3] = 0.4330127018922193*(fc[15]+fc[14]+fc[13])-0.4330127018922193*fc[12]-0.25*fc[11]+0.4330127018922193*fc[10]-0.4330127018922193*(fc[9]+fc[8])-0.25*(fc[7]+fc[6])+0.25*fc[5]-0.4330127018922193*fc[4]-0.25*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[4] = 0.4330127018922193*fl[15]-0.4330127018922193*(fl[14]+fl[13])+0.4330127018922193*fl[12]+0.25*fl[11]+0.4330127018922193*fl[10]-0.4330127018922193*(fl[9]+fl[8])-0.25*(fl[7]+fl[6])+0.25*fl[5]+0.4330127018922193*fl[4]+0.25*fl[3]-0.25*(fl[2]+fl[1])+0.25*fl[0]; 
  } else { 
    fUpwindQuad_l[4] = (-0.4330127018922193*fc[15])+0.4330127018922193*(fc[14]+fc[13])-0.4330127018922193*fc[12]+0.25*fc[11]-0.4330127018922193*fc[10]+0.4330127018922193*(fc[9]+fc[8])-0.25*(fc[7]+fc[6])+0.25*fc[5]-0.4330127018922193*fc[4]+0.25*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = (-0.4330127018922193*(fl[15]+fl[14]))+0.4330127018922193*fl[13]-0.4330127018922193*fl[12]-0.25*fl[11]+0.4330127018922193*fl[10]-0.4330127018922193*fl[9]+0.4330127018922193*fl[8]-0.25*fl[7]+0.25*fl[6]-0.25*fl[5]+0.4330127018922193*fl[4]+0.25*fl[3]-0.25*fl[2]+0.25*(fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[5] = 0.4330127018922193*(fc[15]+fc[14])-0.4330127018922193*fc[13]+0.4330127018922193*fc[12]-0.25*fc[11]-0.4330127018922193*fc[10]+0.4330127018922193*fc[9]-0.4330127018922193*fc[8]-0.25*fc[7]+0.25*fc[6]-0.25*fc[5]-0.4330127018922193*fc[4]+0.25*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[6] = (-0.4330127018922193*fl[15])+0.4330127018922193*fl[14]-0.4330127018922193*(fl[13]+fl[12])-0.25*fl[11]+0.4330127018922193*(fl[10]+fl[9])-0.4330127018922193*fl[8]+0.25*fl[7]-0.25*(fl[6]+fl[5])+0.4330127018922193*fl[4]+0.25*(fl[3]+fl[2])-0.25*fl[1]+0.25*fl[0]; 
  } else { 
    fUpwindQuad_l[6] = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.4330127018922193*(fc[13]+fc[12])-0.25*fc[11]-0.4330127018922193*(fc[10]+fc[9])+0.4330127018922193*fc[8]+0.25*fc[7]-0.25*(fc[6]+fc[5])-0.4330127018922193*fc[4]+0.25*(fc[3]+fc[2])-0.25*fc[1]+0.25*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = 0.4330127018922193*(fl[15]+fl[14]+fl[13]+fl[12])+0.25*fl[11]+0.4330127018922193*(fl[10]+fl[9]+fl[8])+0.25*(fl[7]+fl[6]+fl[5])+0.4330127018922193*fl[4]+0.25*(fl[3]+fl[2]+fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[7] = (-0.4330127018922193*(fc[15]+fc[14]+fc[13]+fc[12]))+0.25*fc[11]-0.4330127018922193*(fc[10]+fc[9]+fc[8])+0.25*(fc[7]+fc[6]+fc[5])-0.4330127018922193*fc[4]+0.25*(fc[3]+fc[2]+fc[1]+fc[0]); 
  } 

  double fUpwind_l[8];
  fUpwind_l[0] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[5] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.3535533905932737*(fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.3535533905932737*(fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 

  double fUpwindQuad_r[8];
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = (-0.4330127018922193*fc[15])+0.4330127018922193*(fc[14]+fc[13]+fc[12])-0.25*fc[11]-0.4330127018922193*(fc[10]+fc[9]+fc[8])+0.25*(fc[7]+fc[6]+fc[5])+0.4330127018922193*fc[4]-0.25*(fc[3]+fc[2]+fc[1])+0.25*fc[0]; 
  } else { 
    fUpwindQuad_r[0] = 0.4330127018922193*fr[15]-0.4330127018922193*(fr[14]+fr[13]+fr[12])-0.25*fr[11]+0.4330127018922193*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.4330127018922193*fr[4]-0.25*(fr[3]+fr[2]+fr[1])+0.25*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = 0.4330127018922193*(fc[15]+fc[14])-0.4330127018922193*(fc[13]+fc[12])+0.25*fc[11]-0.4330127018922193*(fc[10]+fc[9])+0.4330127018922193*fc[8]+0.25*fc[7]-0.25*(fc[6]+fc[5])+0.4330127018922193*fc[4]-0.25*(fc[3]+fc[2])+0.25*(fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[1] = (-0.4330127018922193*(fr[15]+fr[14]))+0.4330127018922193*(fr[13]+fr[12])+0.25*fr[11]+0.4330127018922193*(fr[10]+fr[9])-0.4330127018922193*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])-0.4330127018922193*fr[4]-0.25*(fr[3]+fr[2])+0.25*(fr[1]+fr[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[2] = 0.4330127018922193*fc[15]-0.4330127018922193*fc[14]+0.4330127018922193*fc[13]-0.4330127018922193*fc[12]+0.25*fc[11]-0.4330127018922193*fc[10]+0.4330127018922193*fc[9]-0.4330127018922193*fc[8]-0.25*fc[7]+0.25*fc[6]-0.25*fc[5]+0.4330127018922193*fc[4]-0.25*fc[3]+0.25*fc[2]-0.25*fc[1]+0.25*fc[0]; 
  } else { 
    fUpwindQuad_r[2] = (-0.4330127018922193*fr[15])+0.4330127018922193*fr[14]-0.4330127018922193*fr[13]+0.4330127018922193*fr[12]+0.25*fr[11]+0.4330127018922193*fr[10]-0.4330127018922193*fr[9]+0.4330127018922193*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]-0.4330127018922193*fr[4]-0.25*fr[3]+0.25*fr[2]-0.25*fr[1]+0.25*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = (-0.4330127018922193*(fc[15]+fc[14]+fc[13]))+0.4330127018922193*fc[12]-0.25*fc[11]-0.4330127018922193*fc[10]+0.4330127018922193*(fc[9]+fc[8])-0.25*(fc[7]+fc[6])+0.25*fc[5]+0.4330127018922193*fc[4]-0.25*fc[3]+0.25*(fc[2]+fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[3] = 0.4330127018922193*(fr[15]+fr[14]+fr[13])-0.4330127018922193*fr[12]-0.25*fr[11]+0.4330127018922193*fr[10]-0.4330127018922193*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.4330127018922193*fr[4]-0.25*fr[3]+0.25*(fr[2]+fr[1]+fr[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[4] = 0.4330127018922193*fc[15]-0.4330127018922193*(fc[14]+fc[13])+0.4330127018922193*fc[12]+0.25*fc[11]+0.4330127018922193*fc[10]-0.4330127018922193*(fc[9]+fc[8])-0.25*(fc[7]+fc[6])+0.25*fc[5]+0.4330127018922193*fc[4]+0.25*fc[3]-0.25*(fc[2]+fc[1])+0.25*fc[0]; 
  } else { 
    fUpwindQuad_r[4] = (-0.4330127018922193*fr[15])+0.4330127018922193*(fr[14]+fr[13])-0.4330127018922193*fr[12]+0.25*fr[11]-0.4330127018922193*fr[10]+0.4330127018922193*(fr[9]+fr[8])-0.25*(fr[7]+fr[6])+0.25*fr[5]-0.4330127018922193*fr[4]+0.25*fr[3]-0.25*(fr[2]+fr[1])+0.25*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = (-0.4330127018922193*(fc[15]+fc[14]))+0.4330127018922193*fc[13]-0.4330127018922193*fc[12]-0.25*fc[11]+0.4330127018922193*fc[10]-0.4330127018922193*fc[9]+0.4330127018922193*fc[8]-0.25*fc[7]+0.25*fc[6]-0.25*fc[5]+0.4330127018922193*fc[4]+0.25*fc[3]-0.25*fc[2]+0.25*(fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[5] = 0.4330127018922193*(fr[15]+fr[14])-0.4330127018922193*fr[13]+0.4330127018922193*fr[12]-0.25*fr[11]-0.4330127018922193*fr[10]+0.4330127018922193*fr[9]-0.4330127018922193*fr[8]-0.25*fr[7]+0.25*fr[6]-0.25*fr[5]-0.4330127018922193*fr[4]+0.25*fr[3]-0.25*fr[2]+0.25*(fr[1]+fr[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[6] = (-0.4330127018922193*fc[15])+0.4330127018922193*fc[14]-0.4330127018922193*(fc[13]+fc[12])-0.25*fc[11]+0.4330127018922193*(fc[10]+fc[9])-0.4330127018922193*fc[8]+0.25*fc[7]-0.25*(fc[6]+fc[5])+0.4330127018922193*fc[4]+0.25*(fc[3]+fc[2])-0.25*fc[1]+0.25*fc[0]; 
  } else { 
    fUpwindQuad_r[6] = 0.4330127018922193*fr[15]-0.4330127018922193*fr[14]+0.4330127018922193*(fr[13]+fr[12])-0.25*fr[11]-0.4330127018922193*(fr[10]+fr[9])+0.4330127018922193*fr[8]+0.25*fr[7]-0.25*(fr[6]+fr[5])-0.4330127018922193*fr[4]+0.25*(fr[3]+fr[2])-0.25*fr[1]+0.25*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = 0.4330127018922193*(fc[15]+fc[14]+fc[13]+fc[12])+0.25*fc[11]+0.4330127018922193*(fc[10]+fc[9]+fc[8])+0.25*(fc[7]+fc[6]+fc[5])+0.4330127018922193*fc[4]+0.25*(fc[3]+fc[2]+fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[7] = (-0.4330127018922193*(fr[15]+fr[14]+fr[13]+fr[12]))+0.25*fr[11]-0.4330127018922193*(fr[10]+fr[9]+fr[8])+0.25*(fr[7]+fr[6]+fr[5])-0.4330127018922193*fr[4]+0.25*(fr[3]+fr[2]+fr[1]+fr[0]); 
  } 

  double fUpwind_r[8];
  fUpwind_r[0] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[5] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.3535533905932737*(fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.3535533905932737*(fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 

  double Gdiff_l[8]; 
  double Gdiff_r[8]; 
  double Ghat_l[8]; 
  double Ghat_r[8]; 
  double Gdiff2_l[8]; 
  double Gdiff2_r[8]; 


  Gdiff2_l[0] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[8]-3.464101615137754*nuVtSqSum_l[1]*fc[8]+3.464101615137754*nuVtSqSum_l[0]*fl[4]-3.464101615137754*nuVtSqSum_l[0]*fc[4]+(3.0*fl[1]+3.0*fc[1])*nuVtSqSum_l[1]+(3.0*fl[0]+3.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[8]-3.464101615137754*nuVtSqSum_l[0]*fc[8]+3.464101615137754*nuVtSqSum_l[1]*fl[4]-3.464101615137754*nuVtSqSum_l[1]*fc[4]+(3.0*fl[0]+3.0*fc[0])*nuVtSqSum_l[1]+3.0*nuVtSqSum_l[0]*fl[1]+3.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff2_l[2] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[12]-3.464101615137754*nuVtSqSum_l[1]*fc[12]+3.464101615137754*nuVtSqSum_l[0]*fl[9]-3.464101615137754*nuVtSqSum_l[0]*fc[9]+3.0*nuVtSqSum_l[1]*fl[5]+3.0*nuVtSqSum_l[1]*fc[5]+3.0*nuVtSqSum_l[0]*fl[2]+3.0*nuVtSqSum_l[0]*fc[2]); 
  Gdiff2_l[3] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[13]-3.464101615137754*nuVtSqSum_l[1]*fc[13]+3.464101615137754*nuVtSqSum_l[0]*fl[10]-3.464101615137754*nuVtSqSum_l[0]*fc[10]+3.0*nuVtSqSum_l[1]*fl[6]+3.0*nuVtSqSum_l[1]*fc[6]+3.0*nuVtSqSum_l[0]*fl[3]+3.0*nuVtSqSum_l[0]*fc[3]); 
  Gdiff2_l[4] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[12]-3.464101615137754*nuVtSqSum_l[0]*fc[12]+3.464101615137754*nuVtSqSum_l[1]*fl[9]-3.464101615137754*nuVtSqSum_l[1]*fc[9]+3.0*nuVtSqSum_l[0]*fl[5]+3.0*nuVtSqSum_l[0]*fc[5]+3.0*nuVtSqSum_l[1]*fl[2]+3.0*nuVtSqSum_l[1]*fc[2]); 
  Gdiff2_l[5] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[13]-3.464101615137754*nuVtSqSum_l[0]*fc[13]+3.464101615137754*nuVtSqSum_l[1]*fl[10]-3.464101615137754*nuVtSqSum_l[1]*fc[10]+3.0*nuVtSqSum_l[0]*fl[6]+3.0*nuVtSqSum_l[0]*fc[6]+3.0*nuVtSqSum_l[1]*fl[3]+3.0*nuVtSqSum_l[1]*fc[3]); 
  Gdiff2_l[6] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[15]-3.464101615137754*nuVtSqSum_l[1]*fc[15]+3.464101615137754*nuVtSqSum_l[0]*fl[14]-3.464101615137754*nuVtSqSum_l[0]*fc[14]+3.0*nuVtSqSum_l[1]*fl[11]+3.0*nuVtSqSum_l[1]*fc[11]+3.0*nuVtSqSum_l[0]*fl[7]+3.0*nuVtSqSum_l[0]*fc[7]); 
  Gdiff2_l[7] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[15]-3.464101615137754*nuVtSqSum_l[0]*fc[15]+3.464101615137754*nuVtSqSum_l[1]*fl[14]-3.464101615137754*nuVtSqSum_l[1]*fc[14]+3.0*nuVtSqSum_l[0]*fl[11]+3.0*nuVtSqSum_l[0]*fc[11]+3.0*nuVtSqSum_l[1]*fl[7]+3.0*nuVtSqSum_l[1]*fc[7]); 


  Gdiff2_r[0] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[8]-3.464101615137754*nuVtSqSum_r[1]*fc[8]+3.464101615137754*nuVtSqSum_r[0]*fr[4]-3.464101615137754*nuVtSqSum_r[0]*fc[4]+((-3.0*fr[1])-3.0*fc[1])*nuVtSqSum_r[1]+((-3.0*fr[0])-3.0*fc[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[8]-3.464101615137754*nuVtSqSum_r[0]*fc[8]+3.464101615137754*nuVtSqSum_r[1]*fr[4]-3.464101615137754*nuVtSqSum_r[1]*fc[4]+((-3.0*fr[0])-3.0*fc[0])*nuVtSqSum_r[1]-3.0*nuVtSqSum_r[0]*fr[1]-3.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff2_r[2] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[12]-3.464101615137754*nuVtSqSum_r[1]*fc[12]+3.464101615137754*nuVtSqSum_r[0]*fr[9]-3.464101615137754*nuVtSqSum_r[0]*fc[9]-3.0*nuVtSqSum_r[1]*fr[5]-3.0*nuVtSqSum_r[1]*fc[5]-3.0*nuVtSqSum_r[0]*fr[2]-3.0*nuVtSqSum_r[0]*fc[2]); 
  Gdiff2_r[3] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[13]-3.464101615137754*nuVtSqSum_r[1]*fc[13]+3.464101615137754*nuVtSqSum_r[0]*fr[10]-3.464101615137754*nuVtSqSum_r[0]*fc[10]-3.0*nuVtSqSum_r[1]*fr[6]-3.0*nuVtSqSum_r[1]*fc[6]-3.0*nuVtSqSum_r[0]*fr[3]-3.0*nuVtSqSum_r[0]*fc[3]); 
  Gdiff2_r[4] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[12]-3.464101615137754*nuVtSqSum_r[0]*fc[12]+3.464101615137754*nuVtSqSum_r[1]*fr[9]-3.464101615137754*nuVtSqSum_r[1]*fc[9]-3.0*nuVtSqSum_r[0]*fr[5]-3.0*nuVtSqSum_r[0]*fc[5]-3.0*nuVtSqSum_r[1]*fr[2]-3.0*nuVtSqSum_r[1]*fc[2]); 
  Gdiff2_r[5] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[13]-3.464101615137754*nuVtSqSum_r[0]*fc[13]+3.464101615137754*nuVtSqSum_r[1]*fr[10]-3.464101615137754*nuVtSqSum_r[1]*fc[10]-3.0*nuVtSqSum_r[0]*fr[6]-3.0*nuVtSqSum_r[0]*fc[6]-3.0*nuVtSqSum_r[1]*fr[3]-3.0*nuVtSqSum_r[1]*fc[3]); 
  Gdiff2_r[6] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[15]-3.464101615137754*nuVtSqSum_r[1]*fc[15]+3.464101615137754*nuVtSqSum_r[0]*fr[14]-3.464101615137754*nuVtSqSum_r[0]*fc[14]-3.0*nuVtSqSum_r[1]*fr[11]-3.0*nuVtSqSum_r[1]*fc[11]-3.0*nuVtSqSum_r[0]*fr[7]-3.0*nuVtSqSum_r[0]*fc[7]); 
  Gdiff2_r[7] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[15]-3.464101615137754*nuVtSqSum_r[0]*fc[15]+3.464101615137754*nuVtSqSum_r[1]*fr[14]-3.464101615137754*nuVtSqSum_r[1]*fc[14]-3.0*nuVtSqSum_r[0]*fr[11]-3.0*nuVtSqSum_r[0]*fc[11]-3.0*nuVtSqSum_r[1]*fr[7]-3.0*nuVtSqSum_r[1]*fc[7]); 


  Gdiff_l[0] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[8]+8.660254037844386*nuVtSqSum_l[1]*fc[8]+8.660254037844386*nuVtSqSum_l[0]*fl[4]+8.660254037844386*nuVtSqSum_l[0]*fc[4]+(9.0*fl[1]-9.0*fc[1])*nuVtSqSum_l[1]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[8]+8.660254037844386*nuVtSqSum_l[0]*fc[8]+8.660254037844386*nuVtSqSum_l[1]*fl[4]+8.660254037844386*nuVtSqSum_l[1]*fc[4]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSum_l[1]+9.0*nuVtSqSum_l[0]*fl[1]-9.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff_l[2] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[12]+8.660254037844386*nuVtSqSum_l[1]*fc[12]+8.660254037844386*nuVtSqSum_l[0]*fl[9]+8.660254037844386*nuVtSqSum_l[0]*fc[9]+9.0*nuVtSqSum_l[1]*fl[5]-9.0*nuVtSqSum_l[1]*fc[5]+9.0*nuVtSqSum_l[0]*fl[2]-9.0*nuVtSqSum_l[0]*fc[2]); 
  Gdiff_l[3] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[13]+8.660254037844386*nuVtSqSum_l[1]*fc[13]+8.660254037844386*nuVtSqSum_l[0]*fl[10]+8.660254037844386*nuVtSqSum_l[0]*fc[10]+9.0*nuVtSqSum_l[1]*fl[6]-9.0*nuVtSqSum_l[1]*fc[6]+9.0*nuVtSqSum_l[0]*fl[3]-9.0*nuVtSqSum_l[0]*fc[3]); 
  Gdiff_l[4] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[12]+8.660254037844386*nuVtSqSum_l[0]*fc[12]+8.660254037844386*nuVtSqSum_l[1]*fl[9]+8.660254037844386*nuVtSqSum_l[1]*fc[9]+9.0*nuVtSqSum_l[0]*fl[5]-9.0*nuVtSqSum_l[0]*fc[5]+9.0*nuVtSqSum_l[1]*fl[2]-9.0*nuVtSqSum_l[1]*fc[2]); 
  Gdiff_l[5] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[13]+8.660254037844386*nuVtSqSum_l[0]*fc[13]+8.660254037844386*nuVtSqSum_l[1]*fl[10]+8.660254037844386*nuVtSqSum_l[1]*fc[10]+9.0*nuVtSqSum_l[0]*fl[6]-9.0*nuVtSqSum_l[0]*fc[6]+9.0*nuVtSqSum_l[1]*fl[3]-9.0*nuVtSqSum_l[1]*fc[3]); 
  Gdiff_l[6] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[15]+8.660254037844386*nuVtSqSum_l[1]*fc[15]+8.660254037844386*nuVtSqSum_l[0]*fl[14]+8.660254037844386*nuVtSqSum_l[0]*fc[14]+9.0*nuVtSqSum_l[1]*fl[11]-9.0*nuVtSqSum_l[1]*fc[11]+9.0*nuVtSqSum_l[0]*fl[7]-9.0*nuVtSqSum_l[0]*fc[7]); 
  Gdiff_l[7] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[15]+8.660254037844386*nuVtSqSum_l[0]*fc[15]+8.660254037844386*nuVtSqSum_l[1]*fl[14]+8.660254037844386*nuVtSqSum_l[1]*fc[14]+9.0*nuVtSqSum_l[0]*fl[11]-9.0*nuVtSqSum_l[0]*fc[11]+9.0*nuVtSqSum_l[1]*fl[7]-9.0*nuVtSqSum_l[1]*fc[7]); 


  Gdiff_r[0] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[8]+8.660254037844386*nuVtSqSum_r[1]*fc[8]+8.660254037844386*nuVtSqSum_r[0]*fr[4]+8.660254037844386*nuVtSqSum_r[0]*fc[4]+(9.0*fc[1]-9.0*fr[1])*nuVtSqSum_r[1]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[8]+8.660254037844386*nuVtSqSum_r[0]*fc[8]+8.660254037844386*nuVtSqSum_r[1]*fr[4]+8.660254037844386*nuVtSqSum_r[1]*fc[4]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSum_r[1]-9.0*nuVtSqSum_r[0]*fr[1]+9.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff_r[2] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[12]+8.660254037844386*nuVtSqSum_r[1]*fc[12]+8.660254037844386*nuVtSqSum_r[0]*fr[9]+8.660254037844386*nuVtSqSum_r[0]*fc[9]-9.0*nuVtSqSum_r[1]*fr[5]+9.0*nuVtSqSum_r[1]*fc[5]-9.0*nuVtSqSum_r[0]*fr[2]+9.0*nuVtSqSum_r[0]*fc[2]); 
  Gdiff_r[3] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[13]+8.660254037844386*nuVtSqSum_r[1]*fc[13]+8.660254037844386*nuVtSqSum_r[0]*fr[10]+8.660254037844386*nuVtSqSum_r[0]*fc[10]-9.0*nuVtSqSum_r[1]*fr[6]+9.0*nuVtSqSum_r[1]*fc[6]-9.0*nuVtSqSum_r[0]*fr[3]+9.0*nuVtSqSum_r[0]*fc[3]); 
  Gdiff_r[4] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[12]+8.660254037844386*nuVtSqSum_r[0]*fc[12]+8.660254037844386*nuVtSqSum_r[1]*fr[9]+8.660254037844386*nuVtSqSum_r[1]*fc[9]-9.0*nuVtSqSum_r[0]*fr[5]+9.0*nuVtSqSum_r[0]*fc[5]-9.0*nuVtSqSum_r[1]*fr[2]+9.0*nuVtSqSum_r[1]*fc[2]); 
  Gdiff_r[5] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[13]+8.660254037844386*nuVtSqSum_r[0]*fc[13]+8.660254037844386*nuVtSqSum_r[1]*fr[10]+8.660254037844386*nuVtSqSum_r[1]*fc[10]-9.0*nuVtSqSum_r[0]*fr[6]+9.0*nuVtSqSum_r[0]*fc[6]-9.0*nuVtSqSum_r[1]*fr[3]+9.0*nuVtSqSum_r[1]*fc[3]); 
  Gdiff_r[6] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[15]+8.660254037844386*nuVtSqSum_r[1]*fc[15]+8.660254037844386*nuVtSqSum_r[0]*fr[14]+8.660254037844386*nuVtSqSum_r[0]*fc[14]-9.0*nuVtSqSum_r[1]*fr[11]+9.0*nuVtSqSum_r[1]*fc[11]-9.0*nuVtSqSum_r[0]*fr[7]+9.0*nuVtSqSum_r[0]*fc[7]); 
  Gdiff_r[7] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[15]+8.660254037844386*nuVtSqSum_r[0]*fc[15]+8.660254037844386*nuVtSqSum_r[1]*fr[14]+8.660254037844386*nuVtSqSum_r[1]*fc[14]-9.0*nuVtSqSum_r[0]*fr[11]+9.0*nuVtSqSum_r[0]*fc[11]-9.0*nuVtSqSum_r[1]*fr[7]+9.0*nuVtSqSum_r[1]*fc[7]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]; 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = (-1.0*Gdiff_r[4]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[5] = (-1.0*Gdiff_r[5]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[6] = (-1.0*Gdiff_r[6]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]; 
  Ghat_r[7] = (-1.0*Gdiff_r[7]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[4] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_l[4]*rdv2-0.7071067811865475*Ghat_r[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_l[5]*rdv2-0.7071067811865475*Ghat_r[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_l[6]*rdv2-0.7071067811865475*Ghat_r[6]*rdv2; 
  out[8] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[9] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[10] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_l[7]*rdv2-0.7071067811865475*Ghat_r[7]*rdv2; 
  out[12] += 1.224744871391589*Gdiff2_r[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2-1.224744871391589*Ghat_l[4]*rdv2; 
  out[13] += 1.224744871391589*Gdiff2_r[5]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2-1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2-1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2-1.224744871391589*Ghat_l[7]*rdv2; 
} 
