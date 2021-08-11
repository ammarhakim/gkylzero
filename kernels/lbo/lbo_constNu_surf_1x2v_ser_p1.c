#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vx_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out) 
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

  const double *sumNuUxL = &nuUSumL[0]; 
  const double *sumNuUxR = &nuUSumR[0]; 

  double alphaDrSurfL[4]; 
  alphaDrSurfL[0] = (2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUxL[0]; 
  alphaDrSurfL[1] = -1.414213562373095*sumNuUxL[1]; 

  double alphaDrSurfR[4]; 
  alphaDrSurfR[0] = ((-2.0*w[1])-1.0*dxv[1])*nuSum+1.414213562373095*sumNuUxR[0]; 
  alphaDrSurfR[1] = 1.414213562373095*sumNuUxR[1]; 

  double fUpwindQuadL[4];
  fUpwindQuadL[0] = 0.5*(0.6123724356957944*fl[7]-0.6123724356957944*(fc[7]+fl[6])+0.6123724356957944*fc[6]+0.3535533905932737*(fl[5]+fc[5])-0.6123724356957944*fl[4]+0.6123724356957944*fc[4]-0.3535533905932737*(fl[3]+fc[3])+0.6123724356957944*fl[2]-0.6123724356957944*fc[2]-0.3535533905932737*(fl[1]+fc[1])+0.3535533905932737*(fl[0]+fc[0]))-(0.5*(alphaDrSurfL[0]-alphaDrSurfL[1])*(0.6123724356957944*(fl[7]+fc[7])-0.6123724356957944*(fl[6]+fc[6])+0.3535533905932737*fl[5]-0.3535533905932737*fc[5]-0.6123724356957944*(fl[4]+fc[4])-0.3535533905932737*fl[3]+0.3535533905932737*fc[3]+0.6123724356957944*(fl[2]+fc[2])-0.3535533905932737*fl[1]+0.3535533905932737*(fc[1]+fl[0])-0.3535533905932737*fc[0]))/fabs(alphaDrSurfL[0]-alphaDrSurfL[1]); 
  fUpwindQuadL[1] = 0.5*((-0.6123724356957944*(fl[7]+fl[6]))+0.6123724356957944*(fc[7]+fc[6])-0.3535533905932737*(fl[5]+fc[5])+0.6123724356957944*fl[4]-0.6123724356957944*fc[4]-0.3535533905932737*(fl[3]+fc[3])+0.6123724356957944*fl[2]-0.6123724356957944*fc[2]+0.3535533905932737*(fl[1]+fc[1]+fl[0]+fc[0]))-(0.5*(alphaDrSurfL[1]+alphaDrSurfL[0])*((-0.6123724356957944*(fl[7]+fc[7]+fl[6]+fc[6]))-0.3535533905932737*fl[5]+0.3535533905932737*fc[5]+0.6123724356957944*(fl[4]+fc[4])-0.3535533905932737*fl[3]+0.3535533905932737*fc[3]+0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*(fl[1]+fl[0])-0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfL[1]+alphaDrSurfL[0]); 
  fUpwindQuadL[2] = 0.5*((-0.6123724356957944*fl[7])+0.6123724356957944*(fc[7]+fl[6])-0.6123724356957944*fc[6]-0.3535533905932737*(fl[5]+fc[5])-0.6123724356957944*fl[4]+0.6123724356957944*fc[4]+0.3535533905932737*(fl[3]+fc[3])+0.6123724356957944*fl[2]-0.6123724356957944*fc[2]-0.3535533905932737*(fl[1]+fc[1])+0.3535533905932737*(fl[0]+fc[0]))-(0.5*(alphaDrSurfL[0]-alphaDrSurfL[1])*((-0.6123724356957944*(fl[7]+fc[7]))+0.6123724356957944*(fl[6]+fc[6])-0.3535533905932737*fl[5]+0.3535533905932737*fc[5]-0.6123724356957944*(fl[4]+fc[4])+0.3535533905932737*fl[3]-0.3535533905932737*fc[3]+0.6123724356957944*(fl[2]+fc[2])-0.3535533905932737*fl[1]+0.3535533905932737*(fc[1]+fl[0])-0.3535533905932737*fc[0]))/fabs(alphaDrSurfL[0]-alphaDrSurfL[1]); 
  fUpwindQuadL[3] = 0.5*(0.6123724356957944*(fl[7]+fl[6])-0.6123724356957944*(fc[7]+fc[6])+0.3535533905932737*(fl[5]+fc[5])+0.6123724356957944*fl[4]-0.6123724356957944*fc[4]+0.3535533905932737*(fl[3]+fc[3])+0.6123724356957944*fl[2]-0.6123724356957944*fc[2]+0.3535533905932737*(fl[1]+fc[1]+fl[0]+fc[0]))-(0.5*(alphaDrSurfL[1]+alphaDrSurfL[0])*(0.6123724356957944*(fl[7]+fc[7]+fl[6]+fc[6])+0.3535533905932737*fl[5]-0.3535533905932737*fc[5]+0.6123724356957944*(fl[4]+fc[4])+0.3535533905932737*fl[3]-0.3535533905932737*fc[3]+0.6123724356957944*(fl[2]+fc[2])+0.3535533905932737*(fl[1]+fl[0])-0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfL[1]+alphaDrSurfL[0]); 

  double fUpwindL[4];
  fUpwindL[0] = 0.5*(fUpwindQuadL[3]+fUpwindQuadL[2]+fUpwindQuadL[1]+fUpwindQuadL[0]); 
  fUpwindL[1] = 0.5*(fUpwindQuadL[3]-1.0*fUpwindQuadL[2]+fUpwindQuadL[1]-1.0*fUpwindQuadL[0]); 
  fUpwindL[2] = 0.5*(fUpwindQuadL[3]+fUpwindQuadL[2]-1.0*(fUpwindQuadL[1]+fUpwindQuadL[0])); 
  fUpwindL[3] = 0.5*(fUpwindQuadL[3]-1.0*(fUpwindQuadL[2]+fUpwindQuadL[1])+fUpwindQuadL[0]); 

  double fUpwindQuadR[4];
  fUpwindQuadR[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fc[7]+fr[6])-0.6123724356957944*fc[6]+0.3535533905932737*(fr[5]+fc[5])+0.6123724356957944*fr[4]-0.6123724356957944*fc[4]-0.3535533905932737*(fr[3]+fc[3])-0.6123724356957944*fr[2]+0.6123724356957944*fc[2]-0.3535533905932737*(fr[1]+fc[1])+0.3535533905932737*(fr[0]+fc[0]))-(0.5*(alphaDrSurfR[0]-alphaDrSurfR[1])*(0.6123724356957944*(fr[7]+fc[7])-0.6123724356957944*(fr[6]+fc[6])-0.3535533905932737*fr[5]+0.3535533905932737*fc[5]-0.6123724356957944*(fr[4]+fc[4])+0.3535533905932737*fr[3]-0.3535533905932737*fc[3]+0.6123724356957944*(fr[2]+fc[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fc[1]+fr[0])+0.3535533905932737*fc[0]))/fabs(alphaDrSurfR[0]-alphaDrSurfR[1]); 
  fUpwindQuadR[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fc[7]+fc[6])-0.3535533905932737*(fr[5]+fc[5])-0.6123724356957944*fr[4]+0.6123724356957944*fc[4]-0.3535533905932737*(fr[3]+fc[3])-0.6123724356957944*fr[2]+0.6123724356957944*fc[2]+0.3535533905932737*(fr[1]+fc[1]+fr[0]+fc[0]))-(0.5*(alphaDrSurfR[1]+alphaDrSurfR[0])*((-0.6123724356957944*(fr[7]+fc[7]+fr[6]+fc[6]))+0.3535533905932737*fr[5]-0.3535533905932737*fc[5]+0.6123724356957944*(fr[4]+fc[4])+0.3535533905932737*fr[3]-0.3535533905932737*fc[3]+0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfR[1]+alphaDrSurfR[0]); 
  fUpwindQuadR[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fc[7]+fr[6])+0.6123724356957944*fc[6]-0.3535533905932737*(fr[5]+fc[5])+0.6123724356957944*fr[4]-0.6123724356957944*fc[4]+0.3535533905932737*(fr[3]+fc[3])-0.6123724356957944*fr[2]+0.6123724356957944*fc[2]-0.3535533905932737*(fr[1]+fc[1])+0.3535533905932737*(fr[0]+fc[0]))-(0.5*(alphaDrSurfR[0]-alphaDrSurfR[1])*((-0.6123724356957944*(fr[7]+fc[7]))+0.6123724356957944*(fr[6]+fc[6])+0.3535533905932737*fr[5]-0.3535533905932737*fc[5]-0.6123724356957944*(fr[4]+fc[4])-0.3535533905932737*fr[3]+0.3535533905932737*fc[3]+0.6123724356957944*(fr[2]+fc[2])+0.3535533905932737*fr[1]-0.3535533905932737*(fc[1]+fr[0])+0.3535533905932737*fc[0]))/fabs(alphaDrSurfR[0]-alphaDrSurfR[1]); 
  fUpwindQuadR[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]))+0.6123724356957944*(fc[7]+fc[6])+0.3535533905932737*(fr[5]+fc[5])-0.6123724356957944*fr[4]+0.6123724356957944*fc[4]+0.3535533905932737*(fr[3]+fc[3])-0.6123724356957944*fr[2]+0.6123724356957944*fc[2]+0.3535533905932737*(fr[1]+fc[1]+fr[0]+fc[0]))-(0.5*(alphaDrSurfR[1]+alphaDrSurfR[0])*(0.6123724356957944*(fr[7]+fc[7]+fr[6]+fc[6])-0.3535533905932737*fr[5]+0.3535533905932737*fc[5]+0.6123724356957944*(fr[4]+fc[4])-0.3535533905932737*fr[3]+0.3535533905932737*fc[3]+0.6123724356957944*(fr[2]+fc[2])-0.3535533905932737*(fr[1]+fr[0])+0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfR[1]+alphaDrSurfR[0]); 

  double fUpwindR[4];
  fUpwindR[0] = 0.5*(fUpwindQuadR[3]+fUpwindQuadR[2]+fUpwindQuadR[1]+fUpwindQuadR[0]); 
  fUpwindR[1] = 0.5*(fUpwindQuadR[3]-1.0*fUpwindQuadR[2]+fUpwindQuadR[1]-1.0*fUpwindQuadR[0]); 
  fUpwindR[2] = 0.5*(fUpwindQuadR[3]+fUpwindQuadR[2]-1.0*(fUpwindQuadR[1]+fUpwindQuadR[0])); 
  fUpwindR[3] = 0.5*(fUpwindQuadR[3]-1.0*(fUpwindQuadR[2]+fUpwindQuadR[1])+fUpwindQuadR[0]); 

  double GdiffL[8]; 
  double GdiffR[8]; 
  double GhatL[8]; 
  double GhatR[8]; 
  double incr2_l[8]; 
  double incr2_r[8]; 


  incr2_l[2] = -0.1767766952966367*(2.0*nuVtSqSumL[1]*fl[4]-2.0*nuVtSqSumL[1]*fc[4]+2.0*nuVtSqSumL[0]*fl[2]-2.0*nuVtSqSumL[0]*fc[2]+(1.732050807568877*fl[1]+1.732050807568877*fc[1])*nuVtSqSumL[1]+(1.732050807568877*fl[0]+1.732050807568877*fc[0])*nuVtSqSumL[0]); 
  incr2_l[4] = -0.1767766952966367*(2.0*nuVtSqSumL[0]*fl[4]-2.0*nuVtSqSumL[0]*fc[4]+2.0*nuVtSqSumL[1]*fl[2]-2.0*nuVtSqSumL[1]*fc[2]+(1.732050807568877*fl[0]+1.732050807568877*fc[0])*nuVtSqSumL[1]+1.732050807568877*nuVtSqSumL[0]*fl[1]+1.732050807568877*nuVtSqSumL[0]*fc[1]); 
  incr2_l[6] = -0.1767766952966367*(2.0*nuVtSqSumL[1]*fl[7]-2.0*nuVtSqSumL[1]*fc[7]+2.0*nuVtSqSumL[0]*fl[6]-2.0*nuVtSqSumL[0]*fc[6]+1.732050807568877*nuVtSqSumL[1]*fl[5]+1.732050807568877*nuVtSqSumL[1]*fc[5]+1.732050807568877*nuVtSqSumL[0]*fl[3]+1.732050807568877*nuVtSqSumL[0]*fc[3]); 
  incr2_l[7] = -0.1767766952966367*(2.0*nuVtSqSumL[0]*fl[7]-2.0*nuVtSqSumL[0]*fc[7]+2.0*nuVtSqSumL[1]*fl[6]-2.0*nuVtSqSumL[1]*fc[6]+1.732050807568877*nuVtSqSumL[0]*fl[5]+1.732050807568877*nuVtSqSumL[0]*fc[5]+1.732050807568877*nuVtSqSumL[1]*fl[3]+1.732050807568877*nuVtSqSumL[1]*fc[3]); 


  incr2_r[2] = -0.1767766952966367*(2.0*nuVtSqSumR[1]*fr[4]-2.0*nuVtSqSumR[1]*fc[4]+2.0*nuVtSqSumR[0]*fr[2]-2.0*nuVtSqSumR[0]*fc[2]+((-1.732050807568877*fr[1])-1.732050807568877*fc[1])*nuVtSqSumR[1]+((-1.732050807568877*fr[0])-1.732050807568877*fc[0])*nuVtSqSumR[0]); 
  incr2_r[4] = -0.1767766952966367*(2.0*nuVtSqSumR[0]*fr[4]-2.0*nuVtSqSumR[0]*fc[4]+2.0*nuVtSqSumR[1]*fr[2]-2.0*nuVtSqSumR[1]*fc[2]+((-1.732050807568877*fr[0])-1.732050807568877*fc[0])*nuVtSqSumR[1]-1.732050807568877*nuVtSqSumR[0]*fr[1]-1.732050807568877*nuVtSqSumR[0]*fc[1]); 
  incr2_r[6] = -0.1767766952966367*(2.0*nuVtSqSumR[1]*fr[7]-2.0*nuVtSqSumR[1]*fc[7]+2.0*nuVtSqSumR[0]*fr[6]-2.0*nuVtSqSumR[0]*fc[6]-1.732050807568877*nuVtSqSumR[1]*fr[5]-1.732050807568877*nuVtSqSumR[1]*fc[5]-1.732050807568877*nuVtSqSumR[0]*fr[3]-1.732050807568877*nuVtSqSumR[0]*fc[3]); 
  incr2_r[7] = -0.1767766952966367*(2.0*nuVtSqSumR[0]*fr[7]-2.0*nuVtSqSumR[0]*fc[7]+2.0*nuVtSqSumR[1]*fr[6]-2.0*nuVtSqSumR[1]*fc[6]-1.732050807568877*nuVtSqSumR[0]*fr[5]-1.732050807568877*nuVtSqSumR[0]*fc[5]-1.732050807568877*nuVtSqSumR[1]*fr[3]-1.732050807568877*nuVtSqSumR[1]*fc[3]); 


  GdiffL[0] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[1]*fl[4]+8.660254037844386*nuVtSqSumL[1]*fc[4]+8.660254037844386*nuVtSqSumL[0]*fl[2]+8.660254037844386*nuVtSqSumL[0]*fc[2]+(9.0*fl[1]-9.0*fc[1])*nuVtSqSumL[1]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSumL[0]); 
  GdiffL[1] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[0]*fl[4]+8.660254037844386*nuVtSqSumL[0]*fc[4]+8.660254037844386*nuVtSqSumL[1]*fl[2]+8.660254037844386*nuVtSqSumL[1]*fc[2]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSumL[1]+9.0*nuVtSqSumL[0]*fl[1]-9.0*nuVtSqSumL[0]*fc[1]); 
  GdiffL[3] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[1]*fl[7]+8.660254037844386*nuVtSqSumL[1]*fc[7]+8.660254037844386*nuVtSqSumL[0]*fl[6]+8.660254037844386*nuVtSqSumL[0]*fc[6]+9.0*nuVtSqSumL[1]*fl[5]-9.0*nuVtSqSumL[1]*fc[5]+9.0*nuVtSqSumL[0]*fl[3]-9.0*nuVtSqSumL[0]*fc[3]); 
  GdiffL[5] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[0]*fl[7]+8.660254037844386*nuVtSqSumL[0]*fc[7]+8.660254037844386*nuVtSqSumL[1]*fl[6]+8.660254037844386*nuVtSqSumL[1]*fc[6]+9.0*nuVtSqSumL[0]*fl[5]-9.0*nuVtSqSumL[0]*fc[5]+9.0*nuVtSqSumL[1]*fl[3]-9.0*nuVtSqSumL[1]*fc[3]); 


  GdiffR[0] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[1]*fr[4]+8.660254037844386*nuVtSqSumR[1]*fc[4]+8.660254037844386*nuVtSqSumR[0]*fr[2]+8.660254037844386*nuVtSqSumR[0]*fc[2]+(9.0*fc[1]-9.0*fr[1])*nuVtSqSumR[1]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSumR[0]); 
  GdiffR[1] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[0]*fr[4]+8.660254037844386*nuVtSqSumR[0]*fc[4]+8.660254037844386*nuVtSqSumR[1]*fr[2]+8.660254037844386*nuVtSqSumR[1]*fc[2]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSumR[1]-9.0*nuVtSqSumR[0]*fr[1]+9.0*nuVtSqSumR[0]*fc[1]); 
  GdiffR[3] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[1]*fr[7]+8.660254037844386*nuVtSqSumR[1]*fc[7]+8.660254037844386*nuVtSqSumR[0]*fr[6]+8.660254037844386*nuVtSqSumR[0]*fc[6]-9.0*nuVtSqSumR[1]*fr[5]+9.0*nuVtSqSumR[1]*fc[5]-9.0*nuVtSqSumR[0]*fr[3]+9.0*nuVtSqSumR[0]*fc[3]); 
  GdiffR[5] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[0]*fr[7]+8.660254037844386*nuVtSqSumR[0]*fc[7]+8.660254037844386*nuVtSqSumR[1]*fr[6]+8.660254037844386*nuVtSqSumR[1]*fc[6]-9.0*nuVtSqSumR[0]*fr[5]+9.0*nuVtSqSumR[0]*fc[5]-9.0*nuVtSqSumR[1]*fr[3]+9.0*nuVtSqSumR[1]*fc[3]); 

  GhatL[0] = GdiffL[0]*rdv2+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[1]+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[0]; 
  GhatL[1] = GdiffL[1]*rdv2+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[1]+0.7071067811865475*fUpwindL[0]*alphaDrSurfL[1]; 
  GhatL[3] = GdiffL[3]*rdv2+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[3]+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[2]; 
  GhatL[5] = GdiffL[5]*rdv2+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[3]+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[2]; 

  GhatR[0] = (-1.0*GdiffR[0]*rdv2)+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[1]+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[0]; 
  GhatR[1] = (-1.0*GdiffR[1]*rdv2)+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[1]+0.7071067811865475*fUpwindR[0]*alphaDrSurfR[1]; 
  GhatR[3] = (-1.0*GdiffR[3]*rdv2)+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[3]+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[2]; 
  GhatR[5] = (-1.0*GdiffR[5]*rdv2)+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[3]+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[2]; 

  out[0] += 0.5*GhatL[0]*rdv2-0.5*GhatR[0]*rdv2; 
  out[1] += 0.5*GhatL[1]*rdv2-0.5*GhatR[1]*rdv2; 
  out[2] += incr2_r[2]*rdvSq4+incr2_l[2]*rdvSq4-0.8660254037844386*GhatR[0]*rdv2-0.8660254037844386*GhatL[0]*rdv2; 
  out[3] += 0.5*GhatL[3]*rdv2-0.5*GhatR[3]*rdv2; 
  out[4] += incr2_r[4]*rdvSq4+incr2_l[4]*rdvSq4-0.8660254037844386*GhatR[1]*rdv2-0.8660254037844386*GhatL[1]*rdv2; 
  out[5] += 0.5*GhatL[5]*rdv2-0.5*GhatR[5]*rdv2; 
  out[6] += incr2_r[6]*rdvSq4+incr2_l[6]*rdvSq4-0.8660254037844386*GhatR[3]*rdv2-0.8660254037844386*GhatL[3]*rdv2; 
  out[7] += incr2_r[7]*rdvSq4+incr2_l[7]*rdvSq4-0.8660254037844386*GhatR[5]*rdv2-0.8660254037844386*GhatL[5]*rdv2; 
} 
GKYL_CU_DH void lbo_constNu_surf_1x2v_ser_vy_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double *out) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUyL = &nuUSumL[2]; 
  const double *sumNuUyR = &nuUSumR[2]; 

  double alphaDrSurfL[4]; 
  alphaDrSurfL[0] = (2.0*w[2]+dxv[2])*nuSum-1.414213562373095*sumNuUyL[0]; 
  alphaDrSurfL[1] = -1.414213562373095*sumNuUyL[1]; 

  double alphaDrSurfR[4]; 
  alphaDrSurfR[0] = ((-2.0*w[2])-1.0*dxv[2])*nuSum+1.414213562373095*sumNuUyR[0]; 
  alphaDrSurfR[1] = 1.414213562373095*sumNuUyR[1]; 

  double fUpwindQuadL[4];
  fUpwindQuadL[0] = 0.5*(0.6123724356957944*fl[7]-0.6123724356957944*(fc[7]+fl[6]+fl[5])+0.6123724356957944*(fc[6]+fc[5])+0.3535533905932737*(fl[4]+fc[4])+0.6123724356957944*fl[3]-0.6123724356957944*fc[3]-0.3535533905932737*(fl[2]+fc[2]+fl[1]+fc[1])+0.3535533905932737*(fl[0]+fc[0]))-(0.5*(alphaDrSurfL[0]-alphaDrSurfL[1])*(0.6123724356957944*(fl[7]+fc[7])-0.6123724356957944*(fl[6]+fc[6]+fl[5]+fc[5])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4]+0.6123724356957944*(fl[3]+fc[3])-0.3535533905932737*(fl[2]+fl[1])+0.3535533905932737*(fc[2]+fc[1]+fl[0])-0.3535533905932737*fc[0]))/fabs(alphaDrSurfL[0]-alphaDrSurfL[1]); 
  fUpwindQuadL[1] = 0.5*((-0.6123724356957944*(fl[7]+fl[6]))+0.6123724356957944*(fc[7]+fc[6]+fl[5])-0.6123724356957944*fc[5]-0.3535533905932737*(fl[4]+fc[4])+0.6123724356957944*fl[3]-0.6123724356957944*fc[3]-0.3535533905932737*(fl[2]+fc[2])+0.3535533905932737*(fl[1]+fc[1]+fl[0]+fc[0]))-(0.5*(alphaDrSurfL[1]+alphaDrSurfL[0])*((-0.6123724356957944*(fl[7]+fc[7]+fl[6]+fc[6]))+0.6123724356957944*(fl[5]+fc[5])-0.3535533905932737*fl[4]+0.3535533905932737*fc[4]+0.6123724356957944*(fl[3]+fc[3])-0.3535533905932737*fl[2]+0.3535533905932737*(fc[2]+fl[1]+fl[0])-0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfL[1]+alphaDrSurfL[0]); 
  fUpwindQuadL[2] = 0.5*((-0.6123724356957944*fl[7])+0.6123724356957944*(fc[7]+fl[6])-0.6123724356957944*(fc[6]+fl[5])+0.6123724356957944*fc[5]-0.3535533905932737*(fl[4]+fc[4])+0.6123724356957944*fl[3]-0.6123724356957944*fc[3]+0.3535533905932737*(fl[2]+fc[2])-0.3535533905932737*(fl[1]+fc[1])+0.3535533905932737*(fl[0]+fc[0]))-(0.5*(alphaDrSurfL[0]-alphaDrSurfL[1])*((-0.6123724356957944*(fl[7]+fc[7]))+0.6123724356957944*(fl[6]+fc[6])-0.6123724356957944*(fl[5]+fc[5])-0.3535533905932737*fl[4]+0.3535533905932737*fc[4]+0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*fl[2]-0.3535533905932737*(fc[2]+fl[1])+0.3535533905932737*(fc[1]+fl[0])-0.3535533905932737*fc[0]))/fabs(alphaDrSurfL[0]-alphaDrSurfL[1]); 
  fUpwindQuadL[3] = 0.5*(0.6123724356957944*(fl[7]+fl[6]+fl[5])-0.6123724356957944*(fc[7]+fc[6]+fc[5])+0.3535533905932737*(fl[4]+fc[4])+0.6123724356957944*fl[3]-0.6123724356957944*fc[3]+0.3535533905932737*(fl[2]+fc[2]+fl[1]+fc[1]+fl[0]+fc[0]))-(0.5*(alphaDrSurfL[1]+alphaDrSurfL[0])*(0.6123724356957944*(fl[7]+fc[7]+fl[6]+fc[6]+fl[5]+fc[5])+0.3535533905932737*fl[4]-0.3535533905932737*fc[4]+0.6123724356957944*(fl[3]+fc[3])+0.3535533905932737*(fl[2]+fl[1]+fl[0])-0.3535533905932737*(fc[2]+fc[1]+fc[0])))/fabs(alphaDrSurfL[1]+alphaDrSurfL[0]); 

  double fUpwindL[4];
  fUpwindL[0] = 0.5*(fUpwindQuadL[3]+fUpwindQuadL[2]+fUpwindQuadL[1]+fUpwindQuadL[0]); 
  fUpwindL[1] = 0.5*(fUpwindQuadL[3]-1.0*fUpwindQuadL[2]+fUpwindQuadL[1]-1.0*fUpwindQuadL[0]); 
  fUpwindL[2] = 0.5*(fUpwindQuadL[3]+fUpwindQuadL[2]-1.0*(fUpwindQuadL[1]+fUpwindQuadL[0])); 
  fUpwindL[3] = 0.5*(fUpwindQuadL[3]-1.0*(fUpwindQuadL[2]+fUpwindQuadL[1])+fUpwindQuadL[0]); 

  double fUpwindQuadR[4];
  fUpwindQuadR[0] = 0.5*((-0.6123724356957944*fr[7])+0.6123724356957944*(fc[7]+fr[6]+fr[5])-0.6123724356957944*(fc[6]+fc[5])+0.3535533905932737*(fr[4]+fc[4])-0.6123724356957944*fr[3]+0.6123724356957944*fc[3]-0.3535533905932737*(fr[2]+fc[2]+fr[1]+fc[1])+0.3535533905932737*(fr[0]+fc[0]))-(0.5*(alphaDrSurfR[0]-alphaDrSurfR[1])*(0.6123724356957944*(fr[7]+fc[7])-0.6123724356957944*(fr[6]+fc[6]+fr[5]+fc[5])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4]+0.6123724356957944*(fr[3]+fc[3])+0.3535533905932737*(fr[2]+fr[1])-0.3535533905932737*(fc[2]+fc[1]+fr[0])+0.3535533905932737*fc[0]))/fabs(alphaDrSurfR[0]-alphaDrSurfR[1]); 
  fUpwindQuadR[1] = 0.5*(0.6123724356957944*(fr[7]+fr[6])-0.6123724356957944*(fc[7]+fc[6]+fr[5])+0.6123724356957944*fc[5]-0.3535533905932737*(fr[4]+fc[4])-0.6123724356957944*fr[3]+0.6123724356957944*fc[3]-0.3535533905932737*(fr[2]+fc[2])+0.3535533905932737*(fr[1]+fc[1]+fr[0]+fc[0]))-(0.5*(alphaDrSurfR[1]+alphaDrSurfR[0])*((-0.6123724356957944*(fr[7]+fc[7]+fr[6]+fc[6]))+0.6123724356957944*(fr[5]+fc[5])+0.3535533905932737*fr[4]-0.3535533905932737*fc[4]+0.6123724356957944*(fr[3]+fc[3])+0.3535533905932737*fr[2]-0.3535533905932737*(fc[2]+fr[1]+fr[0])+0.3535533905932737*(fc[1]+fc[0])))/fabs(alphaDrSurfR[1]+alphaDrSurfR[0]); 
  fUpwindQuadR[2] = 0.5*(0.6123724356957944*fr[7]-0.6123724356957944*(fc[7]+fr[6])+0.6123724356957944*(fc[6]+fr[5])-0.6123724356957944*fc[5]-0.3535533905932737*(fr[4]+fc[4])-0.6123724356957944*fr[3]+0.6123724356957944*fc[3]+0.3535533905932737*(fr[2]+fc[2])-0.3535533905932737*(fr[1]+fc[1])+0.3535533905932737*(fr[0]+fc[0]))-(0.5*(alphaDrSurfR[0]-alphaDrSurfR[1])*((-0.6123724356957944*(fr[7]+fc[7]))+0.6123724356957944*(fr[6]+fc[6])-0.6123724356957944*(fr[5]+fc[5])+0.3535533905932737*fr[4]-0.3535533905932737*fc[4]+0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*fr[2]+0.3535533905932737*(fc[2]+fr[1])-0.3535533905932737*(fc[1]+fr[0])+0.3535533905932737*fc[0]))/fabs(alphaDrSurfR[0]-alphaDrSurfR[1]); 
  fUpwindQuadR[3] = 0.5*((-0.6123724356957944*(fr[7]+fr[6]+fr[5]))+0.6123724356957944*(fc[7]+fc[6]+fc[5])+0.3535533905932737*(fr[4]+fc[4])-0.6123724356957944*fr[3]+0.6123724356957944*fc[3]+0.3535533905932737*(fr[2]+fc[2]+fr[1]+fc[1]+fr[0]+fc[0]))-(0.5*(alphaDrSurfR[1]+alphaDrSurfR[0])*(0.6123724356957944*(fr[7]+fc[7]+fr[6]+fc[6]+fr[5]+fc[5])-0.3535533905932737*fr[4]+0.3535533905932737*fc[4]+0.6123724356957944*(fr[3]+fc[3])-0.3535533905932737*(fr[2]+fr[1]+fr[0])+0.3535533905932737*(fc[2]+fc[1]+fc[0])))/fabs(alphaDrSurfR[1]+alphaDrSurfR[0]); 

  double fUpwindR[4];
  fUpwindR[0] = 0.5*(fUpwindQuadR[3]+fUpwindQuadR[2]+fUpwindQuadR[1]+fUpwindQuadR[0]); 
  fUpwindR[1] = 0.5*(fUpwindQuadR[3]-1.0*fUpwindQuadR[2]+fUpwindQuadR[1]-1.0*fUpwindQuadR[0]); 
  fUpwindR[2] = 0.5*(fUpwindQuadR[3]+fUpwindQuadR[2]-1.0*(fUpwindQuadR[1]+fUpwindQuadR[0])); 
  fUpwindR[3] = 0.5*(fUpwindQuadR[3]-1.0*(fUpwindQuadR[2]+fUpwindQuadR[1])+fUpwindQuadR[0]); 

  double GdiffL[8]; 
  double GdiffR[8]; 
  double GhatL[8]; 
  double GhatR[8]; 
  double incr2_l[8]; 
  double incr2_r[8]; 


  incr2_l[3] = -0.1767766952966367*(2.0*nuVtSqSumL[1]*fl[5]-2.0*nuVtSqSumL[1]*fc[5]+2.0*nuVtSqSumL[0]*fl[3]-2.0*nuVtSqSumL[0]*fc[3]+(1.732050807568877*fl[1]+1.732050807568877*fc[1])*nuVtSqSumL[1]+(1.732050807568877*fl[0]+1.732050807568877*fc[0])*nuVtSqSumL[0]); 
  incr2_l[5] = -0.1767766952966367*(2.0*nuVtSqSumL[0]*fl[5]-2.0*nuVtSqSumL[0]*fc[5]+2.0*nuVtSqSumL[1]*fl[3]-2.0*nuVtSqSumL[1]*fc[3]+(1.732050807568877*fl[0]+1.732050807568877*fc[0])*nuVtSqSumL[1]+1.732050807568877*nuVtSqSumL[0]*fl[1]+1.732050807568877*nuVtSqSumL[0]*fc[1]); 
  incr2_l[6] = -0.1767766952966367*(2.0*nuVtSqSumL[1]*fl[7]-2.0*nuVtSqSumL[1]*fc[7]+2.0*nuVtSqSumL[0]*fl[6]-2.0*nuVtSqSumL[0]*fc[6]+1.732050807568877*nuVtSqSumL[1]*fl[4]+1.732050807568877*nuVtSqSumL[1]*fc[4]+1.732050807568877*nuVtSqSumL[0]*fl[2]+1.732050807568877*nuVtSqSumL[0]*fc[2]); 
  incr2_l[7] = -0.1767766952966367*(2.0*nuVtSqSumL[0]*fl[7]-2.0*nuVtSqSumL[0]*fc[7]+2.0*nuVtSqSumL[1]*fl[6]-2.0*nuVtSqSumL[1]*fc[6]+1.732050807568877*nuVtSqSumL[0]*fl[4]+1.732050807568877*nuVtSqSumL[0]*fc[4]+1.732050807568877*nuVtSqSumL[1]*fl[2]+1.732050807568877*nuVtSqSumL[1]*fc[2]); 


  incr2_r[3] = -0.1767766952966367*(2.0*nuVtSqSumR[1]*fr[5]-2.0*nuVtSqSumR[1]*fc[5]+2.0*nuVtSqSumR[0]*fr[3]-2.0*nuVtSqSumR[0]*fc[3]+((-1.732050807568877*fr[1])-1.732050807568877*fc[1])*nuVtSqSumR[1]+((-1.732050807568877*fr[0])-1.732050807568877*fc[0])*nuVtSqSumR[0]); 
  incr2_r[5] = -0.1767766952966367*(2.0*nuVtSqSumR[0]*fr[5]-2.0*nuVtSqSumR[0]*fc[5]+2.0*nuVtSqSumR[1]*fr[3]-2.0*nuVtSqSumR[1]*fc[3]+((-1.732050807568877*fr[0])-1.732050807568877*fc[0])*nuVtSqSumR[1]-1.732050807568877*nuVtSqSumR[0]*fr[1]-1.732050807568877*nuVtSqSumR[0]*fc[1]); 
  incr2_r[6] = -0.1767766952966367*(2.0*nuVtSqSumR[1]*fr[7]-2.0*nuVtSqSumR[1]*fc[7]+2.0*nuVtSqSumR[0]*fr[6]-2.0*nuVtSqSumR[0]*fc[6]-1.732050807568877*nuVtSqSumR[1]*fr[4]-1.732050807568877*nuVtSqSumR[1]*fc[4]-1.732050807568877*nuVtSqSumR[0]*fr[2]-1.732050807568877*nuVtSqSumR[0]*fc[2]); 
  incr2_r[7] = -0.1767766952966367*(2.0*nuVtSqSumR[0]*fr[7]-2.0*nuVtSqSumR[0]*fc[7]+2.0*nuVtSqSumR[1]*fr[6]-2.0*nuVtSqSumR[1]*fc[6]-1.732050807568877*nuVtSqSumR[0]*fr[4]-1.732050807568877*nuVtSqSumR[0]*fc[4]-1.732050807568877*nuVtSqSumR[1]*fr[2]-1.732050807568877*nuVtSqSumR[1]*fc[2]); 


  GdiffL[0] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[1]*fl[5]+8.660254037844386*nuVtSqSumL[1]*fc[5]+8.660254037844386*nuVtSqSumL[0]*fl[3]+8.660254037844386*nuVtSqSumL[0]*fc[3]+(9.0*fl[1]-9.0*fc[1])*nuVtSqSumL[1]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSumL[0]); 
  GdiffL[1] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[0]*fl[5]+8.660254037844386*nuVtSqSumL[0]*fc[5]+8.660254037844386*nuVtSqSumL[1]*fl[3]+8.660254037844386*nuVtSqSumL[1]*fc[3]+(9.0*fl[0]-9.0*fc[0])*nuVtSqSumL[1]+9.0*nuVtSqSumL[0]*fl[1]-9.0*nuVtSqSumL[0]*fc[1]); 
  GdiffL[2] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[1]*fl[7]+8.660254037844386*nuVtSqSumL[1]*fc[7]+8.660254037844386*nuVtSqSumL[0]*fl[6]+8.660254037844386*nuVtSqSumL[0]*fc[6]+9.0*nuVtSqSumL[1]*fl[4]-9.0*nuVtSqSumL[1]*fc[4]+9.0*nuVtSqSumL[0]*fl[2]-9.0*nuVtSqSumL[0]*fc[2]); 
  GdiffL[4] = -0.08838834764831838*(8.660254037844386*nuVtSqSumL[0]*fl[7]+8.660254037844386*nuVtSqSumL[0]*fc[7]+8.660254037844386*nuVtSqSumL[1]*fl[6]+8.660254037844386*nuVtSqSumL[1]*fc[6]+9.0*nuVtSqSumL[0]*fl[4]-9.0*nuVtSqSumL[0]*fc[4]+9.0*nuVtSqSumL[1]*fl[2]-9.0*nuVtSqSumL[1]*fc[2]); 


  GdiffR[0] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[1]*fr[5]+8.660254037844386*nuVtSqSumR[1]*fc[5]+8.660254037844386*nuVtSqSumR[0]*fr[3]+8.660254037844386*nuVtSqSumR[0]*fc[3]+(9.0*fc[1]-9.0*fr[1])*nuVtSqSumR[1]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSumR[0]); 
  GdiffR[1] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[0]*fr[5]+8.660254037844386*nuVtSqSumR[0]*fc[5]+8.660254037844386*nuVtSqSumR[1]*fr[3]+8.660254037844386*nuVtSqSumR[1]*fc[3]+(9.0*fc[0]-9.0*fr[0])*nuVtSqSumR[1]-9.0*nuVtSqSumR[0]*fr[1]+9.0*nuVtSqSumR[0]*fc[1]); 
  GdiffR[2] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[1]*fr[7]+8.660254037844386*nuVtSqSumR[1]*fc[7]+8.660254037844386*nuVtSqSumR[0]*fr[6]+8.660254037844386*nuVtSqSumR[0]*fc[6]-9.0*nuVtSqSumR[1]*fr[4]+9.0*nuVtSqSumR[1]*fc[4]-9.0*nuVtSqSumR[0]*fr[2]+9.0*nuVtSqSumR[0]*fc[2]); 
  GdiffR[4] = -0.08838834764831838*(8.660254037844386*nuVtSqSumR[0]*fr[7]+8.660254037844386*nuVtSqSumR[0]*fc[7]+8.660254037844386*nuVtSqSumR[1]*fr[6]+8.660254037844386*nuVtSqSumR[1]*fc[6]-9.0*nuVtSqSumR[0]*fr[4]+9.0*nuVtSqSumR[0]*fc[4]-9.0*nuVtSqSumR[1]*fr[2]+9.0*nuVtSqSumR[1]*fc[2]); 

  GhatL[0] = GdiffL[0]*rdv2+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[1]+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[0]; 
  GhatL[1] = GdiffL[1]*rdv2+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[1]+0.7071067811865475*fUpwindL[0]*alphaDrSurfL[1]; 
  GhatL[2] = GdiffL[2]*rdv2+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[3]+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[2]; 
  GhatL[4] = GdiffL[4]*rdv2+0.7071067811865475*alphaDrSurfL[0]*fUpwindL[3]+0.7071067811865475*alphaDrSurfL[1]*fUpwindL[2]; 

  GhatR[0] = (-1.0*GdiffR[0]*rdv2)+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[1]+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[0]; 
  GhatR[1] = (-1.0*GdiffR[1]*rdv2)+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[1]+0.7071067811865475*fUpwindR[0]*alphaDrSurfR[1]; 
  GhatR[2] = (-1.0*GdiffR[2]*rdv2)+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[3]+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[2]; 
  GhatR[4] = (-1.0*GdiffR[4]*rdv2)+0.7071067811865475*alphaDrSurfR[0]*fUpwindR[3]+0.7071067811865475*alphaDrSurfR[1]*fUpwindR[2]; 

  out[0] += 0.5*GhatL[0]*rdv2-0.5*GhatR[0]*rdv2; 
  out[1] += 0.5*GhatL[1]*rdv2-0.5*GhatR[1]*rdv2; 
  out[2] += 0.5*GhatL[2]*rdv2-0.5*GhatR[2]*rdv2; 
  out[3] += incr2_r[3]*rdvSq4+incr2_l[3]*rdvSq4-0.8660254037844386*GhatR[0]*rdv2-0.8660254037844386*GhatL[0]*rdv2; 
  out[4] += 0.5*GhatL[4]*rdv2-0.5*GhatR[4]*rdv2; 
  out[5] += incr2_r[5]*rdvSq4+incr2_l[5]*rdvSq4-0.8660254037844386*GhatR[1]*rdv2-0.8660254037844386*GhatL[1]*rdv2; 
  out[6] += incr2_r[6]*rdvSq4+incr2_l[6]*rdvSq4-0.8660254037844386*GhatR[2]*rdv2-0.8660254037844386*GhatL[2]*rdv2; 
  out[7] += incr2_r[7]*rdvSq4+incr2_l[7]*rdvSq4-0.8660254037844386*GhatR[4]*rdv2-0.8660254037844386*GhatL[4]*rdv2; 
} 
