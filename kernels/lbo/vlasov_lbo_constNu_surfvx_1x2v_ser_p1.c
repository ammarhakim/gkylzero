#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum_l, const double *nuUSum_r, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx_l = &nuUSum_l[0]; 
  const double *sumNuUx_r = &nuUSum_r[0]; 

  double alphaDrSurf_l[4]; 
  alphaDrSurf_l[0] = (2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_l[0]; 
  alphaDrSurf_l[1] = -1.414213562373095*sumNuUx_l[1]; 

  double alphaDrSurf_r[4]; 
  alphaDrSurf_r[0] = ((-2.0*w[1])-1.0*dxv[1])*nuSum+1.414213562373095*sumNuUx_r[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*sumNuUx_r[1]; 

  double fUpwindQuad_l[4];
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = 0.6123724356957944*fl[7]-0.6123724356957944*fl[6]+0.3535533905932737*fl[5]-0.6123724356957944*fl[4]-0.3535533905932737*fl[3]+0.6123724356957944*fl[2]-0.3535533905932737*fl[1]+0.3535533905932737*fl[0]; 
  } else { 
    fUpwindQuad_l[0] = (-0.6123724356957944*fc[7])+0.6123724356957944*fc[6]+0.3535533905932737*fc[5]+0.6123724356957944*fc[4]-0.3535533905932737*fc[3]-0.6123724356957944*fc[2]-0.3535533905932737*fc[1]+0.3535533905932737*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = (-0.6123724356957944*(fl[7]+fl[6]))-0.3535533905932737*fl[5]+0.6123724356957944*fl[4]-0.3535533905932737*fl[3]+0.6123724356957944*fl[2]+0.3535533905932737*(fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[1] = 0.6123724356957944*(fc[7]+fc[6])-0.3535533905932737*fc[5]-0.6123724356957944*fc[4]-0.3535533905932737*fc[3]-0.6123724356957944*fc[2]+0.3535533905932737*(fc[1]+fc[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[2] = (-0.6123724356957944*fl[7])+0.6123724356957944*fl[6]-0.3535533905932737*fl[5]-0.6123724356957944*fl[4]+0.3535533905932737*fl[3]+0.6123724356957944*fl[2]-0.3535533905932737*fl[1]+0.3535533905932737*fl[0]; 
  } else { 
    fUpwindQuad_l[2] = 0.6123724356957944*fc[7]-0.6123724356957944*fc[6]-0.3535533905932737*fc[5]+0.6123724356957944*fc[4]+0.3535533905932737*fc[3]-0.6123724356957944*fc[2]-0.3535533905932737*fc[1]+0.3535533905932737*fc[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = 0.6123724356957944*(fl[7]+fl[6])+0.3535533905932737*fl[5]+0.6123724356957944*fl[4]+0.3535533905932737*fl[3]+0.6123724356957944*fl[2]+0.3535533905932737*(fl[1]+fl[0]); 
  } else { 
    fUpwindQuad_l[3] = (-0.6123724356957944*(fc[7]+fc[6]))+0.3535533905932737*fc[5]-0.6123724356957944*fc[4]+0.3535533905932737*fc[3]-0.6123724356957944*fc[2]+0.3535533905932737*(fc[1]+fc[0]); 
  } 

  double fUpwind_l[4];
  fUpwind_l[0] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5*(fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.5*(fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  double fUpwindQuad_r[4];
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = 0.6123724356957944*fc[7]-0.6123724356957944*fc[6]+0.3535533905932737*fc[5]-0.6123724356957944*fc[4]-0.3535533905932737*fc[3]+0.6123724356957944*fc[2]-0.3535533905932737*fc[1]+0.3535533905932737*fc[0]; 
  } else { 
    fUpwindQuad_r[0] = (-0.6123724356957944*fr[7])+0.6123724356957944*fr[6]+0.3535533905932737*fr[5]+0.6123724356957944*fr[4]-0.3535533905932737*fr[3]-0.6123724356957944*fr[2]-0.3535533905932737*fr[1]+0.3535533905932737*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = (-0.6123724356957944*(fc[7]+fc[6]))-0.3535533905932737*fc[5]+0.6123724356957944*fc[4]-0.3535533905932737*fc[3]+0.6123724356957944*fc[2]+0.3535533905932737*(fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[1] = 0.6123724356957944*(fr[7]+fr[6])-0.3535533905932737*fr[5]-0.6123724356957944*fr[4]-0.3535533905932737*fr[3]-0.6123724356957944*fr[2]+0.3535533905932737*(fr[1]+fr[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[2] = (-0.6123724356957944*fc[7])+0.6123724356957944*fc[6]-0.3535533905932737*fc[5]-0.6123724356957944*fc[4]+0.3535533905932737*fc[3]+0.6123724356957944*fc[2]-0.3535533905932737*fc[1]+0.3535533905932737*fc[0]; 
  } else { 
    fUpwindQuad_r[2] = 0.6123724356957944*fr[7]-0.6123724356957944*fr[6]-0.3535533905932737*fr[5]+0.6123724356957944*fr[4]+0.3535533905932737*fr[3]-0.6123724356957944*fr[2]-0.3535533905932737*fr[1]+0.3535533905932737*fr[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = 0.6123724356957944*(fc[7]+fc[6])+0.3535533905932737*fc[5]+0.6123724356957944*fc[4]+0.3535533905932737*fc[3]+0.6123724356957944*fc[2]+0.3535533905932737*(fc[1]+fc[0]); 
  } else { 
    fUpwindQuad_r[3] = (-0.6123724356957944*(fr[7]+fr[6]))+0.3535533905932737*fr[5]-0.6123724356957944*fr[4]+0.3535533905932737*fr[3]-0.6123724356957944*fr[2]+0.3535533905932737*(fr[1]+fr[0]); 
  } 

  double fUpwind_r[4];
  fUpwind_r[0] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5*(fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.5*(fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

  double Gdiff_l[4]; 
  double Gdiff_r[4]; 
  double Ghat_l[4]; 
  double Ghat_r[4]; 
  double Gdiff2_l[4]; 
  double Gdiff2_r[4]; 


  Gdiff2_l[0] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[4]-3.464101615137754*nuVtSqSum_l[1]*fc[4]+3.464101615137754*nuVtSqSum_l[0]*fl[2]-3.464101615137754*nuVtSqSum_l[0]*fc[2]+(3.0*fl[1]+3.0*fc[1])*nuVtSqSum_l[1]+(3.0*fl[0]+3.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[4]-3.464101615137754*nuVtSqSum_l[0]*fc[4]+3.464101615137754*nuVtSqSum_l[1]*fl[2]-3.464101615137754*nuVtSqSum_l[1]*fc[2]+(3.0*fl[0]+3.0*fc[0])*nuVtSqSum_l[1]+3.0*nuVtSqSum_l[0]*fl[1]+3.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff2_l[2] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fl[7]-3.464101615137754*nuVtSqSum_l[1]*fc[7]+3.464101615137754*nuVtSqSum_l[0]*fl[6]-3.464101615137754*nuVtSqSum_l[0]*fc[6]+3.0*nuVtSqSum_l[1]*fl[5]+3.0*nuVtSqSum_l[1]*fc[5]+3.0*nuVtSqSum_l[0]*fl[3]+3.0*nuVtSqSum_l[0]*fc[3]); 
  Gdiff2_l[3] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fl[7]-3.464101615137754*nuVtSqSum_l[0]*fc[7]+3.464101615137754*nuVtSqSum_l[1]*fl[6]-3.464101615137754*nuVtSqSum_l[1]*fc[6]+3.0*nuVtSqSum_l[0]*fl[5]+3.0*nuVtSqSum_l[0]*fc[5]+3.0*nuVtSqSum_l[1]*fl[3]+3.0*nuVtSqSum_l[1]*fc[3]); 


  Gdiff2_r[0] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[4]-3.464101615137754*nuVtSqSum_r[1]*fc[4]+3.464101615137754*nuVtSqSum_r[0]*fr[2]-3.464101615137754*nuVtSqSum_r[0]*fc[2]+((-3.0*fr[1])-3.0*fc[1])*nuVtSqSum_r[1]+((-3.0*fr[0])-3.0*fc[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[4]-3.464101615137754*nuVtSqSum_r[0]*fc[4]+3.464101615137754*nuVtSqSum_r[1]*fr[2]-3.464101615137754*nuVtSqSum_r[1]*fc[2]+((-3.0*fr[0])-3.0*fc[0])*nuVtSqSum_r[1]-3.0*nuVtSqSum_r[0]*fr[1]-3.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff2_r[2] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fr[7]-3.464101615137754*nuVtSqSum_r[1]*fc[7]+3.464101615137754*nuVtSqSum_r[0]*fr[6]-3.464101615137754*nuVtSqSum_r[0]*fc[6]-3.0*nuVtSqSum_r[1]*fr[5]-3.0*nuVtSqSum_r[1]*fc[5]-3.0*nuVtSqSum_r[0]*fr[3]-3.0*nuVtSqSum_r[0]*fc[3]); 
  Gdiff2_r[3] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fr[7]-3.464101615137754*nuVtSqSum_r[0]*fc[7]+3.464101615137754*nuVtSqSum_r[1]*fr[6]-3.464101615137754*nuVtSqSum_r[1]*fc[6]-3.0*nuVtSqSum_r[0]*fr[5]-3.0*nuVtSqSum_r[0]*fc[5]-3.0*nuVtSqSum_r[1]*fr[3]-3.0*nuVtSqSum_r[1]*fc[3]); 


  Gdiff_l[0] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[4]+8.660254037844386*nuVtSqSum_l[1]*fc[4]+8.660254037844386*nuVtSqSum_l[0]*fl[2]+8.660254037844386*nuVtSqSum_l[0]*fc[2]+(9.0*fl[1]-9.0*fc[1])*nuVtSqSum_l[1]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[4]+8.660254037844386*nuVtSqSum_l[0]*fc[4]+8.660254037844386*nuVtSqSum_l[1]*fl[2]+8.660254037844386*nuVtSqSum_l[1]*fc[2]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSum_l[1]+9.0*nuVtSqSum_l[0]*fl[1]-9.0*nuVtSqSum_l[0]*fc[1]); 
  Gdiff_l[2] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fl[7]+8.660254037844386*nuVtSqSum_l[1]*fc[7]+8.660254037844386*nuVtSqSum_l[0]*fl[6]+8.660254037844386*nuVtSqSum_l[0]*fc[6]+9.0*nuVtSqSum_l[1]*fl[5]-9.0*nuVtSqSum_l[1]*fc[5]+9.0*nuVtSqSum_l[0]*fl[3]-9.0*nuVtSqSum_l[0]*fc[3]); 
  Gdiff_l[3] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fl[7]+8.660254037844386*nuVtSqSum_l[0]*fc[7]+8.660254037844386*nuVtSqSum_l[1]*fl[6]+8.660254037844386*nuVtSqSum_l[1]*fc[6]+9.0*nuVtSqSum_l[0]*fl[5]-9.0*nuVtSqSum_l[0]*fc[5]+9.0*nuVtSqSum_l[1]*fl[3]-9.0*nuVtSqSum_l[1]*fc[3]); 


  Gdiff_r[0] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[4]+8.660254037844386*nuVtSqSum_r[1]*fc[4]+8.660254037844386*nuVtSqSum_r[0]*fr[2]+8.660254037844386*nuVtSqSum_r[0]*fc[2]+(9.0*fc[1]-9.0*fr[1])*nuVtSqSum_r[1]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[4]+8.660254037844386*nuVtSqSum_r[0]*fc[4]+8.660254037844386*nuVtSqSum_r[1]*fr[2]+8.660254037844386*nuVtSqSum_r[1]*fc[2]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSum_r[1]-9.0*nuVtSqSum_r[0]*fr[1]+9.0*nuVtSqSum_r[0]*fc[1]); 
  Gdiff_r[2] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fr[7]+8.660254037844386*nuVtSqSum_r[1]*fc[7]+8.660254037844386*nuVtSqSum_r[0]*fr[6]+8.660254037844386*nuVtSqSum_r[0]*fc[6]-9.0*nuVtSqSum_r[1]*fr[5]+9.0*nuVtSqSum_r[1]*fc[5]-9.0*nuVtSqSum_r[0]*fr[3]+9.0*nuVtSqSum_r[0]*fc[3]); 
  Gdiff_r[3] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fr[7]+8.660254037844386*nuVtSqSum_r[0]*fc[7]+8.660254037844386*nuVtSqSum_r[1]*fr[6]+8.660254037844386*nuVtSqSum_r[1]*fc[6]-9.0*nuVtSqSum_r[0]*fr[5]+9.0*nuVtSqSum_r[0]*fc[5]-9.0*nuVtSqSum_r[1]*fr[3]+9.0*nuVtSqSum_r[1]*fc[3]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[4] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
} 
