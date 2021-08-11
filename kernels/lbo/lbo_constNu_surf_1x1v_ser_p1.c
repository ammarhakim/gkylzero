#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void lbo_constNu_surf_1x1v_ser_vx_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSumL, const double *nuUSumR, const double *nuVtSqSumL, const double *nuVtSqSumR, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:          Cell-center coordinates. 
  // dxv[2]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[2]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUxL = &nuUSumL[0]; 
  const double *sumNuUxR = &nuUSumR[0]; 

  double alphaDrSurfL[2]; 
  alphaDrSurfL[0] = 0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUxL[0]); 
  alphaDrSurfL[1] = -1.0*sumNuUxL[1]; 

  double alphaDrSurfR[2]; 
  alphaDrSurfR[0] = -0.7071067811865475*((2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUxR[0]); 
  alphaDrSurfR[1] = sumNuUxR[1]; 

  double fUpwindQuadL[2];
  fUpwindQuadL[0] = 0.5*((-0.8660254037844386*fl[3])+0.8660254037844386*(fc[3]+fl[2])-0.8660254037844386*fc[2]-0.5*(fl[1]+fc[1])+0.5*(fl[0]+fc[0]))-(0.5*(alphaDrSurfL[0]-alphaDrSurfL[1])*((-0.8660254037844386*(fl[3]+fc[3]))+0.8660254037844386*(fl[2]+fc[2])-0.5*fl[1]+0.5*(fc[1]+fl[0])-0.5*fc[0]))/fabs(alphaDrSurfL[0]-alphaDrSurfL[1]); 
  fUpwindQuadL[1] = 0.5*(0.8660254037844386*(fl[3]+fl[2])-0.8660254037844386*(fc[3]+fc[2])+0.5*(fl[1]+fc[1]+fl[0]+fc[0]))-(0.5*(alphaDrSurfL[1]+alphaDrSurfL[0])*(0.8660254037844386*(fl[3]+fc[3]+fl[2]+fc[2])+0.5*(fl[1]+fl[0])-0.5*(fc[1]+fc[0])))/fabs(alphaDrSurfL[1]+alphaDrSurfL[0]); 

  double fUpwindL[2];
  fUpwindL[0] = 0.7071067811865475*(fUpwindQuadL[1]+fUpwindQuadL[0]); 
  fUpwindL[1] = 0.7071067811865475*(fUpwindQuadL[1]-1.0*fUpwindQuadL[0]); 

  double fUpwindQuadR[2];
  fUpwindQuadR[0] = 0.5*(0.8660254037844386*fr[3]-0.8660254037844386*(fc[3]+fr[2])+0.8660254037844386*fc[2]-0.5*(fr[1]+fc[1])+0.5*(fr[0]+fc[0]))-(0.5*(alphaDrSurfR[0]-alphaDrSurfR[1])*((-0.8660254037844386*(fr[3]+fc[3]))+0.8660254037844386*(fr[2]+fc[2])+0.5*fr[1]-0.5*(fc[1]+fr[0])+0.5*fc[0]))/fabs(alphaDrSurfR[0]-alphaDrSurfR[1]); 
  fUpwindQuadR[1] = 0.5*((-0.8660254037844386*(fr[3]+fr[2]))+0.8660254037844386*(fc[3]+fc[2])+0.5*(fr[1]+fc[1]+fr[0]+fc[0]))-(0.5*(alphaDrSurfR[1]+alphaDrSurfR[0])*(0.8660254037844386*(fr[3]+fc[3]+fr[2]+fc[2])-0.5*(fr[1]+fr[0])+0.5*(fc[1]+fc[0])))/fabs(alphaDrSurfR[1]+alphaDrSurfR[0]); 

  double fUpwindR[2];
  fUpwindR[0] = 0.7071067811865475*(fUpwindQuadR[1]+fUpwindQuadR[0]); 
  fUpwindR[1] = 0.7071067811865475*(fUpwindQuadR[1]-1.0*fUpwindQuadR[0]); 

  double GdiffL[4]; 
  double GdiffR[4]; 
  double GhatL[4]; 
  double GhatR[4]; 
  double incr2_l[4]; 
  double incr2_r[4]; 


  incr2_l[2] = -0.125*(2.82842712474619*nuVtSqSumL[1]*fl[3]-2.82842712474619*nuVtSqSumL[1]*fc[3]+2.82842712474619*nuVtSqSumL[0]*fl[2]-2.82842712474619*nuVtSqSumL[0]*fc[2]+(2.449489742783178*fl[1]+2.449489742783178*fc[1])*nuVtSqSumL[1]+(2.449489742783178*fl[0]+2.449489742783178*fc[0])*nuVtSqSumL[0]); 
  incr2_l[3] = -0.125*(2.82842712474619*nuVtSqSumL[0]*fl[3]-2.82842712474619*nuVtSqSumL[0]*fc[3]+2.82842712474619*nuVtSqSumL[1]*fl[2]-2.82842712474619*nuVtSqSumL[1]*fc[2]+(2.449489742783178*fl[0]+2.449489742783178*fc[0])*nuVtSqSumL[1]+2.449489742783178*nuVtSqSumL[0]*fl[1]+2.449489742783178*nuVtSqSumL[0]*fc[1]); 


  incr2_r[2] = -0.125*(2.82842712474619*nuVtSqSumR[1]*fr[3]-2.82842712474619*nuVtSqSumR[1]*fc[3]+2.82842712474619*nuVtSqSumR[0]*fr[2]-2.82842712474619*nuVtSqSumR[0]*fc[2]+((-2.449489742783178*fr[1])-2.449489742783178*fc[1])*nuVtSqSumR[1]+((-2.449489742783178*fr[0])-2.449489742783178*fc[0])*nuVtSqSumR[0]); 
  incr2_r[3] = -0.125*(2.82842712474619*nuVtSqSumR[0]*fr[3]-2.82842712474619*nuVtSqSumR[0]*fc[3]+2.82842712474619*nuVtSqSumR[1]*fr[2]-2.82842712474619*nuVtSqSumR[1]*fc[2]+((-2.449489742783178*fr[0])-2.449489742783178*fc[0])*nuVtSqSumR[1]-2.449489742783178*nuVtSqSumR[0]*fr[1]-2.449489742783178*nuVtSqSumR[0]*fc[1]); 


  GdiffL[0] = -0.0625*(12.24744871391589*nuVtSqSumL[1]*fl[3]+12.24744871391589*nuVtSqSumL[1]*fc[3]+12.24744871391589*nuVtSqSumL[0]*fl[2]+12.24744871391589*nuVtSqSumL[0]*fc[2]+(12.72792206135785*fl[1]-12.72792206135785*fc[1])*nuVtSqSumL[1]+(12.72792206135785*fl[0]-12.72792206135785*fc[0])*nuVtSqSumL[0]); 
  GdiffL[1] = -0.0625*(12.24744871391589*nuVtSqSumL[0]*fl[3]+12.24744871391589*nuVtSqSumL[0]*fc[3]+12.24744871391589*nuVtSqSumL[1]*fl[2]+12.24744871391589*nuVtSqSumL[1]*fc[2]+(12.72792206135785*fl[0]-12.72792206135785*fc[0])*nuVtSqSumL[1]+12.72792206135785*nuVtSqSumL[0]*fl[1]-12.72792206135785*nuVtSqSumL[0]*fc[1]); 


  GdiffR[0] = -0.0625*(12.24744871391589*nuVtSqSumR[1]*fr[3]+12.24744871391589*nuVtSqSumR[1]*fc[3]+12.24744871391589*nuVtSqSumR[0]*fr[2]+12.24744871391589*nuVtSqSumR[0]*fc[2]+(12.72792206135785*fc[1]-12.72792206135785*fr[1])*nuVtSqSumR[1]+(12.72792206135785*fc[0]-12.72792206135785*fr[0])*nuVtSqSumR[0]); 
  GdiffR[1] = -0.0625*(12.24744871391589*nuVtSqSumR[0]*fr[3]+12.24744871391589*nuVtSqSumR[0]*fc[3]+12.24744871391589*nuVtSqSumR[1]*fr[2]+12.24744871391589*nuVtSqSumR[1]*fc[2]+(12.72792206135785*fc[0]-12.72792206135785*fr[0])*nuVtSqSumR[1]-12.72792206135785*nuVtSqSumR[0]*fr[1]+12.72792206135785*nuVtSqSumR[0]*fc[1]); 

  GhatL[0] = GdiffL[0]*rdv2+alphaDrSurfL[1]*fUpwindL[1]+alphaDrSurfL[0]*fUpwindL[0]; 
  GhatL[1] = GdiffL[1]*rdv2+alphaDrSurfL[0]*fUpwindL[1]+fUpwindL[0]*alphaDrSurfL[1]; 

  GhatR[0] = (-1.0*GdiffR[0]*rdv2)+alphaDrSurfR[1]*fUpwindR[1]+alphaDrSurfR[0]*fUpwindR[0]; 
  GhatR[1] = (-1.0*GdiffR[1]*rdv2)+alphaDrSurfR[0]*fUpwindR[1]+fUpwindR[0]*alphaDrSurfR[1]; 

  out[0] += 0.5*GhatL[0]*rdv2-0.5*GhatR[0]*rdv2; 
  out[1] += 0.5*GhatL[1]*rdv2-0.5*GhatR[1]*rdv2; 
  out[2] += incr2_r[2]*rdvSq4+incr2_l[2]*rdvSq4-0.8660254037844386*GhatR[0]*rdv2-0.8660254037844386*GhatL[0]*rdv2; 
  out[3] += incr2_r[3]*rdvSq4+incr2_l[3]*rdvSq4-0.8660254037844386*GhatR[1]*rdv2-0.8660254037844386*GhatL[1]*rdv2; 
} 
