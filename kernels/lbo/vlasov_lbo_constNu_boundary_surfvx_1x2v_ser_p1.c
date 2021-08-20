#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvx_1x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[3]:          Cell-center coordinates. 
  // dxv[3]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[4]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[2]:  sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  double Gdiff2_l[4]; 
  double Gdiff2_r[4]; 


  Gdiff2_l[0] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[4]-3.464101615137754*nuVtSqSum_l[1]*fEdge[4]+3.464101615137754*nuVtSqSum_l[0]*fSkin[2]-3.464101615137754*nuVtSqSum_l[0]*fEdge[2]+((-3.0*fSkin[1])-3.0*fEdge[1])*nuVtSqSum_l[1]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[4]-3.464101615137754*nuVtSqSum_l[0]*fEdge[4]+3.464101615137754*nuVtSqSum_l[1]*fSkin[2]-3.464101615137754*nuVtSqSum_l[1]*fEdge[2]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[1]-3.0*nuVtSqSum_l[0]*fSkin[1]-3.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff2_l[2] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[1]*fSkin[7]-3.464101615137754*nuVtSqSum_l[1]*fEdge[7]+3.464101615137754*nuVtSqSum_l[0]*fSkin[6]-3.464101615137754*nuVtSqSum_l[0]*fEdge[6]-3.0*nuVtSqSum_l[1]*fSkin[5]-3.0*nuVtSqSum_l[1]*fEdge[5]-3.0*nuVtSqSum_l[0]*fSkin[3]-3.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff2_l[3] = -0.08333333333333333*(3.464101615137754*nuVtSqSum_l[0]*fSkin[7]-3.464101615137754*nuVtSqSum_l[0]*fEdge[7]+3.464101615137754*nuVtSqSum_l[1]*fSkin[6]-3.464101615137754*nuVtSqSum_l[1]*fEdge[6]-3.0*nuVtSqSum_l[0]*fSkin[5]-3.0*nuVtSqSum_l[0]*fEdge[5]-3.0*nuVtSqSum_l[1]*fSkin[3]-3.0*nuVtSqSum_l[1]*fEdge[3]); 


  Gdiff2_r[0] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[4]-3.464101615137754*nuVtSqSum_r[1]*fEdge[4]+3.464101615137754*nuVtSqSum_r[0]*fSkin[2]-3.464101615137754*nuVtSqSum_r[0]*fEdge[2]+(3.0*fSkin[1]+3.0*fEdge[1])*nuVtSqSum_r[1]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[4]-3.464101615137754*nuVtSqSum_r[0]*fEdge[4]+3.464101615137754*nuVtSqSum_r[1]*fSkin[2]-3.464101615137754*nuVtSqSum_r[1]*fEdge[2]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[1]+3.0*nuVtSqSum_r[0]*fSkin[1]+3.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff2_r[2] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[1]*fSkin[7]-3.464101615137754*nuVtSqSum_r[1]*fEdge[7]+3.464101615137754*nuVtSqSum_r[0]*fSkin[6]-3.464101615137754*nuVtSqSum_r[0]*fEdge[6]+3.0*nuVtSqSum_r[1]*fSkin[5]+3.0*nuVtSqSum_r[1]*fEdge[5]+3.0*nuVtSqSum_r[0]*fSkin[3]+3.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff2_r[3] = 0.08333333333333333*(3.464101615137754*nuVtSqSum_r[0]*fSkin[7]-3.464101615137754*nuVtSqSum_r[0]*fEdge[7]+3.464101615137754*nuVtSqSum_r[1]*fSkin[6]-3.464101615137754*nuVtSqSum_r[1]*fEdge[6]+3.0*nuVtSqSum_r[0]*fSkin[5]+3.0*nuVtSqSum_r[0]*fEdge[5]+3.0*nuVtSqSum_r[1]*fSkin[3]+3.0*nuVtSqSum_r[1]*fEdge[3]); 

  if (edge == -1) { 

  const double *sumNuUx_r = &nuUSum[0]; 

  double alphaDrSurf_r[4]; 
  alphaDrSurf_r[0] = ((-2.0*w[1])-1.0*dxv[1])*nuSum+1.414213562373095*sumNuUx_r[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*sumNuUx_r[1]; 

  double fUpwindQuad_r[4];
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[0] = 0.6123724356957944*fSkin[7]-0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 
    fUpwindQuad_r[0] = (-0.6123724356957944*fEdge[7])+0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = (-0.6123724356957944*(fSkin[7]+fSkin[6]))-0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[1] = 0.6123724356957944*(fEdge[7]+fEdge[6])-0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } 
  if (alphaDrSurf_r[0]-alphaDrSurf_r[1] > 0) { 
    fUpwindQuad_r[2] = (-0.6123724356957944*fSkin[7])+0.6123724356957944*fSkin[6]-0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } else { 
    fUpwindQuad_r[2] = 0.6123724356957944*fEdge[7]-0.6123724356957944*fEdge[6]-0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } 
  if (alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = 0.6123724356957944*(fSkin[7]+fSkin[6])+0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]+0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[3] = (-0.6123724356957944*(fEdge[7]+fEdge[6]))+0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]-0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } 

  double fUpwind_r[4];
  fUpwind_r[0] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.5*(fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.5*(fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.5*(fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

  double Gdiff_r[4]; 
  double Ghat_r[4]; 


  Gdiff_r[0] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[4]+8.660254037844386*nuVtSqSum_r[1]*fEdge[4]+8.660254037844386*nuVtSqSum_r[0]*fSkin[2]+8.660254037844386*nuVtSqSum_r[0]*fEdge[2]+(9.0*fSkin[1]-9.0*fEdge[1])*nuVtSqSum_r[1]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[4]+8.660254037844386*nuVtSqSum_r[0]*fEdge[4]+8.660254037844386*nuVtSqSum_r[1]*fSkin[2]+8.660254037844386*nuVtSqSum_r[1]*fEdge[2]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[1]+9.0*nuVtSqSum_r[0]*fSkin[1]-9.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff_r[2] = -0.0625*(8.660254037844386*nuVtSqSum_r[1]*fSkin[7]+8.660254037844386*nuVtSqSum_r[1]*fEdge[7]+8.660254037844386*nuVtSqSum_r[0]*fSkin[6]+8.660254037844386*nuVtSqSum_r[0]*fEdge[6]+9.0*nuVtSqSum_r[1]*fSkin[5]-9.0*nuVtSqSum_r[1]*fEdge[5]+9.0*nuVtSqSum_r[0]*fSkin[3]-9.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff_r[3] = -0.0625*(8.660254037844386*nuVtSqSum_r[0]*fSkin[7]+8.660254037844386*nuVtSqSum_r[0]*fEdge[7]+8.660254037844386*nuVtSqSum_r[1]*fSkin[6]+8.660254037844386*nuVtSqSum_r[1]*fEdge[6]+9.0*nuVtSqSum_r[0]*fSkin[5]-9.0*nuVtSqSum_r[0]*fEdge[5]+9.0*nuVtSqSum_r[1]*fSkin[3]-9.0*nuVtSqSum_r[1]*fEdge[3]); 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 

  out[0] += -0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[0]*rdvSq4+1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2; 
  out[3] += -0.7071067811865475*Ghat_r[2]*rdv2; 
  out[4] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum_l[1]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2; 
  out[5] += -0.7071067811865475*Ghat_r[3]*rdv2; 
  out[6] += 1.060660171779821*nuVtSqSum_l[1]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum_l[0]*fSkin[6]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[5]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[3]*rdvSq4+1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2; 
  out[7] += 1.060660171779821*nuVtSqSum_l[0]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum_l[1]*fSkin[6]*rdvSq4-0.6123724356957944*nuVtSqSum_l[0]*fSkin[5]*rdvSq4-0.6123724356957944*nuVtSqSum_l[1]*fSkin[3]*rdvSq4+1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2; 

  } else {

  const double *sumNuUx_l = &nuUSum[0]; 

  double alphaDrSurf_l[4]; 
  alphaDrSurf_l[0] = (2.0*w[1]+dxv[1])*nuSum-1.414213562373095*sumNuUx_l[0]; 
  alphaDrSurf_l[1] = -1.414213562373095*sumNuUx_l[1]; 

  double fUpwindQuad_l[4];
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[0] = 0.6123724356957944*fEdge[7]-0.6123724356957944*fEdge[6]+0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 
    fUpwindQuad_l[0] = (-0.6123724356957944*fSkin[7])+0.6123724356957944*fSkin[6]+0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = (-0.6123724356957944*(fEdge[7]+fEdge[6]))-0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]-0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[1] = 0.6123724356957944*(fSkin[7]+fSkin[6])-0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]-0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } 
  if (alphaDrSurf_l[0]-alphaDrSurf_l[1] > 0) { 
    fUpwindQuad_l[2] = (-0.6123724356957944*fEdge[7])+0.6123724356957944*fEdge[6]-0.3535533905932737*fEdge[5]-0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]-0.3535533905932737*fEdge[1]+0.3535533905932737*fEdge[0]; 
  } else { 
    fUpwindQuad_l[2] = 0.6123724356957944*fSkin[7]-0.6123724356957944*fSkin[6]-0.3535533905932737*fSkin[5]+0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]-0.3535533905932737*fSkin[1]+0.3535533905932737*fSkin[0]; 
  } 
  if (alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = 0.6123724356957944*(fEdge[7]+fEdge[6])+0.3535533905932737*fEdge[5]+0.6123724356957944*fEdge[4]+0.3535533905932737*fEdge[3]+0.6123724356957944*fEdge[2]+0.3535533905932737*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[3] = (-0.6123724356957944*(fSkin[7]+fSkin[6]))+0.3535533905932737*fSkin[5]-0.6123724356957944*fSkin[4]+0.3535533905932737*fSkin[3]-0.6123724356957944*fSkin[2]+0.3535533905932737*(fSkin[1]+fSkin[0]); 
  } 

  double fUpwind_l[4];
  fUpwind_l[0] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.5*(fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.5*(fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.5*(fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  double Gdiff_l[4]; 
  double Ghat_l[4]; 


  Gdiff_l[0] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[4]+8.660254037844386*nuVtSqSum_l[1]*fEdge[4]+8.660254037844386*nuVtSqSum_l[0]*fSkin[2]+8.660254037844386*nuVtSqSum_l[0]*fEdge[2]+(9.0*fEdge[1]-9.0*fSkin[1])*nuVtSqSum_l[1]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[4]+8.660254037844386*nuVtSqSum_l[0]*fEdge[4]+8.660254037844386*nuVtSqSum_l[1]*fSkin[2]+8.660254037844386*nuVtSqSum_l[1]*fEdge[2]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[1]-9.0*nuVtSqSum_l[0]*fSkin[1]+9.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff_l[2] = -0.0625*(8.660254037844386*nuVtSqSum_l[1]*fSkin[7]+8.660254037844386*nuVtSqSum_l[1]*fEdge[7]+8.660254037844386*nuVtSqSum_l[0]*fSkin[6]+8.660254037844386*nuVtSqSum_l[0]*fEdge[6]-9.0*nuVtSqSum_l[1]*fSkin[5]+9.0*nuVtSqSum_l[1]*fEdge[5]-9.0*nuVtSqSum_l[0]*fSkin[3]+9.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff_l[3] = -0.0625*(8.660254037844386*nuVtSqSum_l[0]*fSkin[7]+8.660254037844386*nuVtSqSum_l[0]*fEdge[7]+8.660254037844386*nuVtSqSum_l[1]*fSkin[6]+8.660254037844386*nuVtSqSum_l[1]*fEdge[6]-9.0*nuVtSqSum_l[0]*fSkin[5]+9.0*nuVtSqSum_l[0]*fEdge[5]-9.0*nuVtSqSum_l[1]*fSkin[3]+9.0*nuVtSqSum_l[1]*fEdge[3]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[2]*rdv2; 
  out[4] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum_r[1]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_l[3]*rdv2; 
  out[6] += 1.060660171779821*nuVtSqSum_r[1]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum_r[0]*fSkin[6]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[5]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[3]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 1.060660171779821*nuVtSqSum_r[0]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum_r[1]*fSkin[6]*rdvSq4+0.6123724356957944*nuVtSqSum_r[0]*fSkin[5]*rdvSq4+0.6123724356957944*nuVtSqSum_r[1]*fSkin[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_l[3]*rdv2; 
  } 
} 
