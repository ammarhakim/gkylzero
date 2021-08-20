#include <gkyl_lbo_kernels.h> 
GKYL_CU_DH void vlasov_lbo_constNu_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum_l, const double *nuVtSqSum_r, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:          Cell-center coordinates. 
  // dxv[4]:        Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:     sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]:  sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:      Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  double Gdiff2_l[8]; 
  double Gdiff2_r[8]; 


  Gdiff2_l[0] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[3]*fSkin[12]-3.464101615137754*nuVtSqSum_l[3]*fEdge[12]+3.464101615137754*nuVtSqSum_l[2]*fSkin[9]-3.464101615137754*nuVtSqSum_l[2]*fEdge[9]+3.464101615137754*nuVtSqSum_l[1]*fSkin[8]-3.464101615137754*nuVtSqSum_l[1]*fEdge[8]-3.0*nuVtSqSum_l[3]*fSkin[5]-3.0*nuVtSqSum_l[3]*fEdge[5]+3.464101615137754*nuVtSqSum_l[0]*fSkin[4]-3.464101615137754*nuVtSqSum_l[0]*fEdge[4]+((-3.0*fSkin[2])-3.0*fEdge[2])*nuVtSqSum_l[2]+((-3.0*fSkin[1])-3.0*fEdge[1])*nuVtSqSum_l[1]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[0]); 
  Gdiff2_l[1] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[2]*fSkin[12]-3.464101615137754*nuVtSqSum_l[2]*fEdge[12]+3.464101615137754*nuVtSqSum_l[3]*fSkin[9]-3.464101615137754*nuVtSqSum_l[3]*fEdge[9]+3.464101615137754*nuVtSqSum_l[0]*fSkin[8]-3.464101615137754*nuVtSqSum_l[0]*fEdge[8]-3.0*nuVtSqSum_l[2]*fSkin[5]-3.0*nuVtSqSum_l[2]*fEdge[5]+3.464101615137754*nuVtSqSum_l[1]*fSkin[4]-3.464101615137754*nuVtSqSum_l[1]*fEdge[4]+((-3.0*fSkin[2])-3.0*fEdge[2])*nuVtSqSum_l[3]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[1]-3.0*nuVtSqSum_l[0]*fSkin[1]-3.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff2_l[2] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[1]*fSkin[12]-3.464101615137754*nuVtSqSum_l[1]*fEdge[12]+3.464101615137754*nuVtSqSum_l[0]*fSkin[9]-3.464101615137754*nuVtSqSum_l[0]*fEdge[9]+3.464101615137754*nuVtSqSum_l[3]*fSkin[8]-3.464101615137754*nuVtSqSum_l[3]*fEdge[8]-3.0*nuVtSqSum_l[1]*fSkin[5]-3.0*nuVtSqSum_l[1]*fEdge[5]+3.464101615137754*nuVtSqSum_l[2]*fSkin[4]-3.464101615137754*nuVtSqSum_l[2]*fEdge[4]+((-3.0*fSkin[1])-3.0*fEdge[1])*nuVtSqSum_l[3]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[2]-3.0*nuVtSqSum_l[0]*fSkin[2]-3.0*nuVtSqSum_l[0]*fEdge[2]); 
  Gdiff2_l[3] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[3]*fSkin[15]-3.464101615137754*nuVtSqSum_l[3]*fEdge[15]+3.464101615137754*nuVtSqSum_l[2]*fSkin[14]-3.464101615137754*nuVtSqSum_l[2]*fEdge[14]+3.464101615137754*nuVtSqSum_l[1]*fSkin[13]-3.464101615137754*nuVtSqSum_l[1]*fEdge[13]-3.0*nuVtSqSum_l[3]*fSkin[11]-3.0*nuVtSqSum_l[3]*fEdge[11]+3.464101615137754*nuVtSqSum_l[0]*fSkin[10]-3.464101615137754*nuVtSqSum_l[0]*fEdge[10]-3.0*nuVtSqSum_l[2]*fSkin[7]-3.0*nuVtSqSum_l[2]*fEdge[7]-3.0*nuVtSqSum_l[1]*fSkin[6]-3.0*nuVtSqSum_l[1]*fEdge[6]-3.0*nuVtSqSum_l[0]*fSkin[3]-3.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff2_l[4] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[0]*fSkin[12]-3.464101615137754*nuVtSqSum_l[0]*fEdge[12]+3.464101615137754*nuVtSqSum_l[1]*fSkin[9]-3.464101615137754*nuVtSqSum_l[1]*fEdge[9]+3.464101615137754*nuVtSqSum_l[2]*fSkin[8]-3.464101615137754*nuVtSqSum_l[2]*fEdge[8]-3.0*nuVtSqSum_l[0]*fSkin[5]-3.0*nuVtSqSum_l[0]*fEdge[5]+3.464101615137754*nuVtSqSum_l[3]*fSkin[4]-3.464101615137754*nuVtSqSum_l[3]*fEdge[4]+((-3.0*fSkin[0])-3.0*fEdge[0])*nuVtSqSum_l[3]+((-3.0*fSkin[1])-3.0*fEdge[1])*nuVtSqSum_l[2]-3.0*nuVtSqSum_l[1]*fSkin[2]-3.0*nuVtSqSum_l[1]*fEdge[2]); 
  Gdiff2_l[5] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[2]*fSkin[15]-3.464101615137754*nuVtSqSum_l[2]*fEdge[15]+3.464101615137754*nuVtSqSum_l[3]*fSkin[14]-3.464101615137754*nuVtSqSum_l[3]*fEdge[14]+3.464101615137754*nuVtSqSum_l[0]*fSkin[13]-3.464101615137754*nuVtSqSum_l[0]*fEdge[13]-3.0*nuVtSqSum_l[2]*fSkin[11]-3.0*nuVtSqSum_l[2]*fEdge[11]+3.464101615137754*nuVtSqSum_l[1]*fSkin[10]-3.464101615137754*nuVtSqSum_l[1]*fEdge[10]-3.0*nuVtSqSum_l[3]*fSkin[7]-3.0*nuVtSqSum_l[3]*fEdge[7]-3.0*nuVtSqSum_l[0]*fSkin[6]-3.0*nuVtSqSum_l[0]*fEdge[6]-3.0*nuVtSqSum_l[1]*fSkin[3]-3.0*nuVtSqSum_l[1]*fEdge[3]); 
  Gdiff2_l[6] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[1]*fSkin[15]-3.464101615137754*nuVtSqSum_l[1]*fEdge[15]+3.464101615137754*nuVtSqSum_l[0]*fSkin[14]-3.464101615137754*nuVtSqSum_l[0]*fEdge[14]+3.464101615137754*nuVtSqSum_l[3]*fSkin[13]-3.464101615137754*nuVtSqSum_l[3]*fEdge[13]-3.0*nuVtSqSum_l[1]*fSkin[11]-3.0*nuVtSqSum_l[1]*fEdge[11]+3.464101615137754*nuVtSqSum_l[2]*fSkin[10]-3.464101615137754*nuVtSqSum_l[2]*fEdge[10]-3.0*nuVtSqSum_l[0]*fSkin[7]-3.0*nuVtSqSum_l[0]*fEdge[7]-3.0*nuVtSqSum_l[3]*fSkin[6]-3.0*nuVtSqSum_l[3]*fEdge[6]-3.0*nuVtSqSum_l[2]*fSkin[3]-3.0*nuVtSqSum_l[2]*fEdge[3]); 
  Gdiff2_l[7] = -0.05892556509887892*(3.464101615137754*nuVtSqSum_l[0]*fSkin[15]-3.464101615137754*nuVtSqSum_l[0]*fEdge[15]+3.464101615137754*nuVtSqSum_l[1]*fSkin[14]-3.464101615137754*nuVtSqSum_l[1]*fEdge[14]+3.464101615137754*nuVtSqSum_l[2]*fSkin[13]-3.464101615137754*nuVtSqSum_l[2]*fEdge[13]-3.0*nuVtSqSum_l[0]*fSkin[11]-3.0*nuVtSqSum_l[0]*fEdge[11]+3.464101615137754*nuVtSqSum_l[3]*fSkin[10]-3.464101615137754*nuVtSqSum_l[3]*fEdge[10]-3.0*nuVtSqSum_l[1]*fSkin[7]-3.0*nuVtSqSum_l[1]*fEdge[7]-3.0*nuVtSqSum_l[2]*fSkin[6]-3.0*nuVtSqSum_l[2]*fEdge[6]+((-3.0*fSkin[3])-3.0*fEdge[3])*nuVtSqSum_l[3]); 


  Gdiff2_r[0] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[3]*fSkin[12]-3.464101615137754*nuVtSqSum_r[3]*fEdge[12]+3.464101615137754*nuVtSqSum_r[2]*fSkin[9]-3.464101615137754*nuVtSqSum_r[2]*fEdge[9]+3.464101615137754*nuVtSqSum_r[1]*fSkin[8]-3.464101615137754*nuVtSqSum_r[1]*fEdge[8]+3.0*nuVtSqSum_r[3]*fSkin[5]+3.0*nuVtSqSum_r[3]*fEdge[5]+3.464101615137754*nuVtSqSum_r[0]*fSkin[4]-3.464101615137754*nuVtSqSum_r[0]*fEdge[4]+(3.0*fSkin[2]+3.0*fEdge[2])*nuVtSqSum_r[2]+(3.0*fSkin[1]+3.0*fEdge[1])*nuVtSqSum_r[1]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff2_r[1] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[2]*fSkin[12]-3.464101615137754*nuVtSqSum_r[2]*fEdge[12]+3.464101615137754*nuVtSqSum_r[3]*fSkin[9]-3.464101615137754*nuVtSqSum_r[3]*fEdge[9]+3.464101615137754*nuVtSqSum_r[0]*fSkin[8]-3.464101615137754*nuVtSqSum_r[0]*fEdge[8]+3.0*nuVtSqSum_r[2]*fSkin[5]+3.0*nuVtSqSum_r[2]*fEdge[5]+3.464101615137754*nuVtSqSum_r[1]*fSkin[4]-3.464101615137754*nuVtSqSum_r[1]*fEdge[4]+(3.0*fSkin[2]+3.0*fEdge[2])*nuVtSqSum_r[3]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[1]+3.0*nuVtSqSum_r[0]*fSkin[1]+3.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff2_r[2] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[1]*fSkin[12]-3.464101615137754*nuVtSqSum_r[1]*fEdge[12]+3.464101615137754*nuVtSqSum_r[0]*fSkin[9]-3.464101615137754*nuVtSqSum_r[0]*fEdge[9]+3.464101615137754*nuVtSqSum_r[3]*fSkin[8]-3.464101615137754*nuVtSqSum_r[3]*fEdge[8]+3.0*nuVtSqSum_r[1]*fSkin[5]+3.0*nuVtSqSum_r[1]*fEdge[5]+3.464101615137754*nuVtSqSum_r[2]*fSkin[4]-3.464101615137754*nuVtSqSum_r[2]*fEdge[4]+(3.0*fSkin[1]+3.0*fEdge[1])*nuVtSqSum_r[3]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[2]+3.0*nuVtSqSum_r[0]*fSkin[2]+3.0*nuVtSqSum_r[0]*fEdge[2]); 
  Gdiff2_r[3] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[3]*fSkin[15]-3.464101615137754*nuVtSqSum_r[3]*fEdge[15]+3.464101615137754*nuVtSqSum_r[2]*fSkin[14]-3.464101615137754*nuVtSqSum_r[2]*fEdge[14]+3.464101615137754*nuVtSqSum_r[1]*fSkin[13]-3.464101615137754*nuVtSqSum_r[1]*fEdge[13]+3.0*nuVtSqSum_r[3]*fSkin[11]+3.0*nuVtSqSum_r[3]*fEdge[11]+3.464101615137754*nuVtSqSum_r[0]*fSkin[10]-3.464101615137754*nuVtSqSum_r[0]*fEdge[10]+3.0*nuVtSqSum_r[2]*fSkin[7]+3.0*nuVtSqSum_r[2]*fEdge[7]+3.0*nuVtSqSum_r[1]*fSkin[6]+3.0*nuVtSqSum_r[1]*fEdge[6]+3.0*nuVtSqSum_r[0]*fSkin[3]+3.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff2_r[4] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[0]*fSkin[12]-3.464101615137754*nuVtSqSum_r[0]*fEdge[12]+3.464101615137754*nuVtSqSum_r[1]*fSkin[9]-3.464101615137754*nuVtSqSum_r[1]*fEdge[9]+3.464101615137754*nuVtSqSum_r[2]*fSkin[8]-3.464101615137754*nuVtSqSum_r[2]*fEdge[8]+3.0*nuVtSqSum_r[0]*fSkin[5]+3.0*nuVtSqSum_r[0]*fEdge[5]+3.464101615137754*nuVtSqSum_r[3]*fSkin[4]-3.464101615137754*nuVtSqSum_r[3]*fEdge[4]+(3.0*fSkin[0]+3.0*fEdge[0])*nuVtSqSum_r[3]+(3.0*fSkin[1]+3.0*fEdge[1])*nuVtSqSum_r[2]+3.0*nuVtSqSum_r[1]*fSkin[2]+3.0*nuVtSqSum_r[1]*fEdge[2]); 
  Gdiff2_r[5] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[2]*fSkin[15]-3.464101615137754*nuVtSqSum_r[2]*fEdge[15]+3.464101615137754*nuVtSqSum_r[3]*fSkin[14]-3.464101615137754*nuVtSqSum_r[3]*fEdge[14]+3.464101615137754*nuVtSqSum_r[0]*fSkin[13]-3.464101615137754*nuVtSqSum_r[0]*fEdge[13]+3.0*nuVtSqSum_r[2]*fSkin[11]+3.0*nuVtSqSum_r[2]*fEdge[11]+3.464101615137754*nuVtSqSum_r[1]*fSkin[10]-3.464101615137754*nuVtSqSum_r[1]*fEdge[10]+3.0*nuVtSqSum_r[3]*fSkin[7]+3.0*nuVtSqSum_r[3]*fEdge[7]+3.0*nuVtSqSum_r[0]*fSkin[6]+3.0*nuVtSqSum_r[0]*fEdge[6]+3.0*nuVtSqSum_r[1]*fSkin[3]+3.0*nuVtSqSum_r[1]*fEdge[3]); 
  Gdiff2_r[6] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[1]*fSkin[15]-3.464101615137754*nuVtSqSum_r[1]*fEdge[15]+3.464101615137754*nuVtSqSum_r[0]*fSkin[14]-3.464101615137754*nuVtSqSum_r[0]*fEdge[14]+3.464101615137754*nuVtSqSum_r[3]*fSkin[13]-3.464101615137754*nuVtSqSum_r[3]*fEdge[13]+3.0*nuVtSqSum_r[1]*fSkin[11]+3.0*nuVtSqSum_r[1]*fEdge[11]+3.464101615137754*nuVtSqSum_r[2]*fSkin[10]-3.464101615137754*nuVtSqSum_r[2]*fEdge[10]+3.0*nuVtSqSum_r[0]*fSkin[7]+3.0*nuVtSqSum_r[0]*fEdge[7]+3.0*nuVtSqSum_r[3]*fSkin[6]+3.0*nuVtSqSum_r[3]*fEdge[6]+3.0*nuVtSqSum_r[2]*fSkin[3]+3.0*nuVtSqSum_r[2]*fEdge[3]); 
  Gdiff2_r[7] = 0.05892556509887892*(3.464101615137754*nuVtSqSum_r[0]*fSkin[15]-3.464101615137754*nuVtSqSum_r[0]*fEdge[15]+3.464101615137754*nuVtSqSum_r[1]*fSkin[14]-3.464101615137754*nuVtSqSum_r[1]*fEdge[14]+3.464101615137754*nuVtSqSum_r[2]*fSkin[13]-3.464101615137754*nuVtSqSum_r[2]*fEdge[13]+3.0*nuVtSqSum_r[0]*fSkin[11]+3.0*nuVtSqSum_r[0]*fEdge[11]+3.464101615137754*nuVtSqSum_r[3]*fSkin[10]-3.464101615137754*nuVtSqSum_r[3]*fEdge[10]+3.0*nuVtSqSum_r[1]*fSkin[7]+3.0*nuVtSqSum_r[1]*fEdge[7]+3.0*nuVtSqSum_r[2]*fSkin[6]+3.0*nuVtSqSum_r[2]*fEdge[6]+(3.0*fSkin[3]+3.0*fEdge[3])*nuVtSqSum_r[3]); 

  if (edge == -1) { 

  const double *sumNuUy_r = &nuUSum[4]; 

  double alphaDrSurf_r[8]; 
  alphaDrSurf_r[0] = -0.7071067811865475*((4.0*w[3]+2.0*dxv[3])*nuSum-2.0*sumNuUy_r[0]); 
  alphaDrSurf_r[1] = 1.414213562373095*sumNuUy_r[1]; 
  alphaDrSurf_r[2] = 1.414213562373095*sumNuUy_r[2]; 
  alphaDrSurf_r[4] = 1.414213562373095*sumNuUy_r[3]; 

  double fUpwindQuad_r[8];
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*(fSkin[14]+fSkin[13]+fSkin[12])-0.25*fSkin[11]-0.4330127018922193*(fSkin[10]+fSkin[9]+fSkin[8])+0.25*(fSkin[7]+fSkin[6]+fSkin[5])+0.4330127018922193*fSkin[4]-0.25*(fSkin[3]+fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[0] = 0.4330127018922193*fEdge[15]-0.4330127018922193*(fEdge[14]+fEdge[13]+fEdge[12])-0.25*fEdge[11]+0.4330127018922193*(fEdge[10]+fEdge[9]+fEdge[8])+0.25*(fEdge[7]+fEdge[6]+fEdge[5])-0.4330127018922193*fEdge[4]-0.25*(fEdge[3]+fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } 
  if ((-alphaDrSurf_r[4])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.4330127018922193*(fSkin[13]+fSkin[12])+0.25*fSkin[11]-0.4330127018922193*(fSkin[10]+fSkin[9])+0.4330127018922193*fSkin[8]+0.25*fSkin[7]-0.25*(fSkin[6]+fSkin[5])+0.4330127018922193*fSkin[4]-0.25*(fSkin[3]+fSkin[2])+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[1] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.4330127018922193*(fEdge[13]+fEdge[12])+0.25*fEdge[11]+0.4330127018922193*(fEdge[10]+fEdge[9])-0.4330127018922193*fEdge[8]+0.25*fEdge[7]-0.25*(fEdge[6]+fEdge[5])-0.4330127018922193*fEdge[4]-0.25*(fEdge[3]+fEdge[2])+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alphaDrSurf_r[4])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.4330127018922193*fSkin[13]-0.4330127018922193*fSkin[12]+0.25*fSkin[11]-0.4330127018922193*fSkin[10]+0.4330127018922193*fSkin[9]-0.4330127018922193*fSkin[8]-0.25*fSkin[7]+0.25*fSkin[6]-0.25*fSkin[5]+0.4330127018922193*fSkin[4]-0.25*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[2] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.4330127018922193*fEdge[13]+0.4330127018922193*fEdge[12]+0.25*fEdge[11]+0.4330127018922193*fEdge[10]-0.4330127018922193*fEdge[9]+0.4330127018922193*fEdge[8]-0.25*fEdge[7]+0.25*fEdge[6]-0.25*fEdge[5]-0.4330127018922193*fEdge[4]-0.25*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = (-0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13]))+0.4330127018922193*fSkin[12]-0.25*fSkin[11]-0.4330127018922193*fSkin[10]+0.4330127018922193*(fSkin[9]+fSkin[8])-0.25*(fSkin[7]+fSkin[6])+0.25*fSkin[5]+0.4330127018922193*fSkin[4]-0.25*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[3] = 0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13])-0.4330127018922193*fEdge[12]-0.25*fEdge[11]+0.4330127018922193*fEdge[10]-0.4330127018922193*(fEdge[9]+fEdge[8])-0.25*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.4330127018922193*fEdge[4]-0.25*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } 
  if (alphaDrSurf_r[4]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = 0.4330127018922193*fSkin[15]-0.4330127018922193*(fSkin[14]+fSkin[13])+0.4330127018922193*fSkin[12]+0.25*fSkin[11]+0.4330127018922193*fSkin[10]-0.4330127018922193*(fSkin[9]+fSkin[8])-0.25*(fSkin[7]+fSkin[6])+0.25*fSkin[5]+0.4330127018922193*fSkin[4]+0.25*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[4] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*(fEdge[14]+fEdge[13])-0.4330127018922193*fEdge[12]+0.25*fEdge[11]-0.4330127018922193*fEdge[10]+0.4330127018922193*(fEdge[9]+fEdge[8])-0.25*(fEdge[7]+fEdge[6])+0.25*fEdge[5]-0.4330127018922193*fEdge[4]+0.25*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } 
  if ((-alphaDrSurf_r[4])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.4330127018922193*fSkin[13]-0.4330127018922193*fSkin[12]-0.25*fSkin[11]+0.4330127018922193*fSkin[10]-0.4330127018922193*fSkin[9]+0.4330127018922193*fSkin[8]-0.25*fSkin[7]+0.25*fSkin[6]-0.25*fSkin[5]+0.4330127018922193*fSkin[4]+0.25*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[5] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.4330127018922193*fEdge[13]+0.4330127018922193*fEdge[12]-0.25*fEdge[11]-0.4330127018922193*fEdge[10]+0.4330127018922193*fEdge[9]-0.4330127018922193*fEdge[8]-0.25*fEdge[7]+0.25*fEdge[6]-0.25*fEdge[5]-0.4330127018922193*fEdge[4]+0.25*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } 
  if ((-alphaDrSurf_r[4])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.4330127018922193*(fSkin[13]+fSkin[12])-0.25*fSkin[11]+0.4330127018922193*(fSkin[10]+fSkin[9])-0.4330127018922193*fSkin[8]+0.25*fSkin[7]-0.25*(fSkin[6]+fSkin[5])+0.4330127018922193*fSkin[4]+0.25*(fSkin[3]+fSkin[2])-0.25*fSkin[1]+0.25*fSkin[0]; 
  } else { 
    fUpwindQuad_r[6] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.4330127018922193*(fEdge[13]+fEdge[12])-0.25*fEdge[11]-0.4330127018922193*(fEdge[10]+fEdge[9])+0.4330127018922193*fEdge[8]+0.25*fEdge[7]-0.25*(fEdge[6]+fEdge[5])-0.4330127018922193*fEdge[4]+0.25*(fEdge[3]+fEdge[2])-0.25*fEdge[1]+0.25*fEdge[0]; 
  } 
  if (alphaDrSurf_r[4]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = 0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13]+fSkin[12])+0.25*fSkin[11]+0.4330127018922193*(fSkin[10]+fSkin[9]+fSkin[8])+0.25*(fSkin[7]+fSkin[6]+fSkin[5])+0.4330127018922193*fSkin[4]+0.25*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
  } else { 
    fUpwindQuad_r[7] = (-0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13]+fEdge[12]))+0.25*fEdge[11]-0.4330127018922193*(fEdge[10]+fEdge[9]+fEdge[8])+0.25*(fEdge[7]+fEdge[6]+fEdge[5])-0.4330127018922193*fEdge[4]+0.25*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
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

  double Gdiff_r[8]; 
  double Ghat_r[8]; 


  Gdiff_r[0] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[3]*fSkin[12]+5.0*nuVtSqSum_r[3]*fEdge[12]+5.0*nuVtSqSum_r[2]*fSkin[9]+5.0*nuVtSqSum_r[2]*fEdge[9]+5.0*nuVtSqSum_r[1]*fSkin[8]+5.0*nuVtSqSum_r[1]*fEdge[8]+5.0*nuVtSqSum_r[0]*fSkin[4]+5.0*nuVtSqSum_r[0]*fEdge[4])+9.0*nuVtSqSum_r[3]*fSkin[5]-9.0*nuVtSqSum_r[3]*fEdge[5]+(9.0*fSkin[2]-9.0*fEdge[2])*nuVtSqSum_r[2]+(9.0*fSkin[1]-9.0*fEdge[1])*nuVtSqSum_r[1]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[0]); 
  Gdiff_r[1] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[2]*fSkin[12]+5.0*nuVtSqSum_r[2]*fEdge[12]+5.0*nuVtSqSum_r[3]*fSkin[9]+5.0*nuVtSqSum_r[3]*fEdge[9]+5.0*nuVtSqSum_r[0]*fSkin[8]+5.0*nuVtSqSum_r[0]*fEdge[8]+5.0*nuVtSqSum_r[1]*fSkin[4]+5.0*nuVtSqSum_r[1]*fEdge[4])+9.0*nuVtSqSum_r[2]*fSkin[5]-9.0*nuVtSqSum_r[2]*fEdge[5]+(9.0*fSkin[2]-9.0*fEdge[2])*nuVtSqSum_r[3]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[1]+9.0*nuVtSqSum_r[0]*fSkin[1]-9.0*nuVtSqSum_r[0]*fEdge[1]); 
  Gdiff_r[2] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[1]*fSkin[12]+5.0*nuVtSqSum_r[1]*fEdge[12]+5.0*nuVtSqSum_r[0]*fSkin[9]+5.0*nuVtSqSum_r[0]*fEdge[9]+5.0*nuVtSqSum_r[3]*fSkin[8]+5.0*nuVtSqSum_r[3]*fEdge[8]+5.0*nuVtSqSum_r[2]*fSkin[4]+5.0*nuVtSqSum_r[2]*fEdge[4])+9.0*nuVtSqSum_r[1]*fSkin[5]-9.0*nuVtSqSum_r[1]*fEdge[5]+(9.0*fSkin[1]-9.0*fEdge[1])*nuVtSqSum_r[3]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[2]+9.0*nuVtSqSum_r[0]*fSkin[2]-9.0*nuVtSqSum_r[0]*fEdge[2]); 
  Gdiff_r[3] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[3]*fSkin[15]+5.0*nuVtSqSum_r[3]*fEdge[15]+5.0*nuVtSqSum_r[2]*fSkin[14]+5.0*nuVtSqSum_r[2]*fEdge[14]+5.0*nuVtSqSum_r[1]*fSkin[13]+5.0*nuVtSqSum_r[1]*fEdge[13]+5.0*nuVtSqSum_r[0]*fSkin[10]+5.0*nuVtSqSum_r[0]*fEdge[10])+9.0*nuVtSqSum_r[3]*fSkin[11]-9.0*nuVtSqSum_r[3]*fEdge[11]+9.0*nuVtSqSum_r[2]*fSkin[7]-9.0*nuVtSqSum_r[2]*fEdge[7]+9.0*nuVtSqSum_r[1]*fSkin[6]-9.0*nuVtSqSum_r[1]*fEdge[6]+9.0*nuVtSqSum_r[0]*fSkin[3]-9.0*nuVtSqSum_r[0]*fEdge[3]); 
  Gdiff_r[4] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[0]*fSkin[12]+5.0*nuVtSqSum_r[0]*fEdge[12]+5.0*nuVtSqSum_r[1]*fSkin[9]+5.0*nuVtSqSum_r[1]*fEdge[9]+5.0*nuVtSqSum_r[2]*fSkin[8]+5.0*nuVtSqSum_r[2]*fEdge[8]+5.0*nuVtSqSum_r[3]*fSkin[4]+5.0*nuVtSqSum_r[3]*fEdge[4])+9.0*nuVtSqSum_r[0]*fSkin[5]-9.0*nuVtSqSum_r[0]*fEdge[5]+(9.0*fSkin[0]-9.0*fEdge[0])*nuVtSqSum_r[3]+(9.0*fSkin[1]-9.0*fEdge[1])*nuVtSqSum_r[2]+9.0*nuVtSqSum_r[1]*fSkin[2]-9.0*nuVtSqSum_r[1]*fEdge[2]); 
  Gdiff_r[5] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[2]*fSkin[15]+5.0*nuVtSqSum_r[2]*fEdge[15]+5.0*nuVtSqSum_r[3]*fSkin[14]+5.0*nuVtSqSum_r[3]*fEdge[14]+5.0*nuVtSqSum_r[0]*fSkin[13]+5.0*nuVtSqSum_r[0]*fEdge[13]+5.0*nuVtSqSum_r[1]*fSkin[10]+5.0*nuVtSqSum_r[1]*fEdge[10])+9.0*nuVtSqSum_r[2]*fSkin[11]-9.0*nuVtSqSum_r[2]*fEdge[11]+9.0*nuVtSqSum_r[3]*fSkin[7]-9.0*nuVtSqSum_r[3]*fEdge[7]+9.0*nuVtSqSum_r[0]*fSkin[6]-9.0*nuVtSqSum_r[0]*fEdge[6]+9.0*nuVtSqSum_r[1]*fSkin[3]-9.0*nuVtSqSum_r[1]*fEdge[3]); 
  Gdiff_r[6] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[1]*fSkin[15]+5.0*nuVtSqSum_r[1]*fEdge[15]+5.0*nuVtSqSum_r[0]*fSkin[14]+5.0*nuVtSqSum_r[0]*fEdge[14]+5.0*nuVtSqSum_r[3]*fSkin[13]+5.0*nuVtSqSum_r[3]*fEdge[13]+5.0*nuVtSqSum_r[2]*fSkin[10]+5.0*nuVtSqSum_r[2]*fEdge[10])+9.0*nuVtSqSum_r[1]*fSkin[11]-9.0*nuVtSqSum_r[1]*fEdge[11]+9.0*nuVtSqSum_r[0]*fSkin[7]-9.0*nuVtSqSum_r[0]*fEdge[7]+9.0*nuVtSqSum_r[3]*fSkin[6]-9.0*nuVtSqSum_r[3]*fEdge[6]+9.0*nuVtSqSum_r[2]*fSkin[3]-9.0*nuVtSqSum_r[2]*fEdge[3]); 
  Gdiff_r[7] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_r[0]*fSkin[15]+5.0*nuVtSqSum_r[0]*fEdge[15]+5.0*nuVtSqSum_r[1]*fSkin[14]+5.0*nuVtSqSum_r[1]*fEdge[14]+5.0*nuVtSqSum_r[2]*fSkin[13]+5.0*nuVtSqSum_r[2]*fEdge[13]+5.0*nuVtSqSum_r[3]*fSkin[10]+5.0*nuVtSqSum_r[3]*fEdge[10])+9.0*nuVtSqSum_r[0]*fSkin[11]-9.0*nuVtSqSum_r[0]*fEdge[11]+9.0*nuVtSqSum_r[1]*fSkin[7]-9.0*nuVtSqSum_r[1]*fEdge[7]+9.0*nuVtSqSum_r[2]*fSkin[6]-9.0*nuVtSqSum_r[2]*fEdge[6]+(9.0*fSkin[3]-9.0*fEdge[3])*nuVtSqSum_r[3]); 

  Ghat_r[0] = (-1.0*Gdiff_r[0]*rdv2)+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[4]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[2]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[1]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = (-1.0*Gdiff_r[1]*rdv2)+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[2]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[1]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = (-1.0*Gdiff_r[2]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = (-1.0*Gdiff_r[3]*rdv2)+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = (-1.0*Gdiff_r[4]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[4]+0.3535533905932737*fUpwind_r[0]*alphaDrSurf_r[4]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[2]+0.3535533905932737*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Ghat_r[5] = (-1.0*Gdiff_r[5]*rdv2)+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[6] = (-1.0*Gdiff_r[6]*rdv2)+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[4]*fUpwind_r[5]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Ghat_r[7] = (-1.0*Gdiff_r[7]*rdv2)+0.3535533905932737*alphaDrSurf_r[0]*fUpwind_r[7]+0.3535533905932737*alphaDrSurf_r[1]*fUpwind_r[6]+0.3535533905932737*alphaDrSurf_r[2]*fUpwind_r[5]+0.3535533905932737*fUpwind_r[3]*alphaDrSurf_r[4]; 

  out[0] += -0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat_r[3]*rdv2; 
  out[4] += 0.75*nuVtSqSum_l[3]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum_l[3]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum_l[2]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum_l[1]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum_l[0]*rdvSq4+1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat_r[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat_r[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat_r[6]*rdv2; 
  out[8] += 0.75*nuVtSqSum_l[2]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum_l[2]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum_l[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum_l[1]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2; 
  out[9] += 0.75*nuVtSqSum_l[1]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum_l[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum_l[2]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[2]*rdvSq4+1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2; 
  out[10] += 0.75*nuVtSqSum_l[3]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum_l[3]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum_l[2]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[6]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[3]*rdvSq4+1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat_r[7]*rdv2; 
  out[12] += 0.75*nuVtSqSum_l[0]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[4]*rdvSq4+1.224744871391589*Gdiff2_r[4]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum_l[3]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum_l[2]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2; 
  out[13] += 0.75*nuVtSqSum_l[2]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum_l[2]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum_l[3]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2_r[5]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2; 
  out[14] += 0.75*nuVtSqSum_l[1]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_l[0]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum_l[3]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2_r[6]*rdvSq4-0.4330127018922193*nuVtSqSum_l[2]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2; 
  out[15] += 0.75*nuVtSqSum_l[0]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_l[1]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_l[2]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum_l[0]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_l[3]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum_l[1]*fSkin[7]*rdvSq4+1.224744871391589*Gdiff2_r[7]*rdvSq4-0.4330127018922193*nuVtSqSum_l[2]*fSkin[6]*rdvSq4-0.4330127018922193*fSkin[3]*nuVtSqSum_l[3]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2; 

  } else {

  const double *sumNuUy_l = &nuUSum[4]; 

  double alphaDrSurf_l[8]; 
  alphaDrSurf_l[0] = 0.7071067811865475*((4.0*w[3]+2.0*dxv[3])*nuSum-2.0*sumNuUy_l[0]); 
  alphaDrSurf_l[1] = -1.414213562373095*sumNuUy_l[1]; 
  alphaDrSurf_l[2] = -1.414213562373095*sumNuUy_l[2]; 
  alphaDrSurf_l[4] = -1.414213562373095*sumNuUy_l[3]; 

  double fUpwindQuad_l[8];
  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*(fEdge[14]+fEdge[13]+fEdge[12])-0.25*fEdge[11]-0.4330127018922193*(fEdge[10]+fEdge[9]+fEdge[8])+0.25*(fEdge[7]+fEdge[6]+fEdge[5])+0.4330127018922193*fEdge[4]-0.25*(fEdge[3]+fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[0] = 0.4330127018922193*fSkin[15]-0.4330127018922193*(fSkin[14]+fSkin[13]+fSkin[12])-0.25*fSkin[11]+0.4330127018922193*(fSkin[10]+fSkin[9]+fSkin[8])+0.25*(fSkin[7]+fSkin[6]+fSkin[5])-0.4330127018922193*fSkin[4]-0.25*(fSkin[3]+fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } 
  if ((-alphaDrSurf_l[4])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = 0.4330127018922193*(fEdge[15]+fEdge[14])-0.4330127018922193*(fEdge[13]+fEdge[12])+0.25*fEdge[11]-0.4330127018922193*(fEdge[10]+fEdge[9])+0.4330127018922193*fEdge[8]+0.25*fEdge[7]-0.25*(fEdge[6]+fEdge[5])+0.4330127018922193*fEdge[4]-0.25*(fEdge[3]+fEdge[2])+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[1] = (-0.4330127018922193*(fSkin[15]+fSkin[14]))+0.4330127018922193*(fSkin[13]+fSkin[12])+0.25*fSkin[11]+0.4330127018922193*(fSkin[10]+fSkin[9])-0.4330127018922193*fSkin[8]+0.25*fSkin[7]-0.25*(fSkin[6]+fSkin[5])-0.4330127018922193*fSkin[4]-0.25*(fSkin[3]+fSkin[2])+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alphaDrSurf_l[4])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = 0.4330127018922193*fEdge[15]-0.4330127018922193*fEdge[14]+0.4330127018922193*fEdge[13]-0.4330127018922193*fEdge[12]+0.25*fEdge[11]-0.4330127018922193*fEdge[10]+0.4330127018922193*fEdge[9]-0.4330127018922193*fEdge[8]-0.25*fEdge[7]+0.25*fEdge[6]-0.25*fEdge[5]+0.4330127018922193*fEdge[4]-0.25*fEdge[3]+0.25*fEdge[2]-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[2] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*fSkin[14]-0.4330127018922193*fSkin[13]+0.4330127018922193*fSkin[12]+0.25*fSkin[11]+0.4330127018922193*fSkin[10]-0.4330127018922193*fSkin[9]+0.4330127018922193*fSkin[8]-0.25*fSkin[7]+0.25*fSkin[6]-0.25*fSkin[5]-0.4330127018922193*fSkin[4]-0.25*fSkin[3]+0.25*fSkin[2]-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = (-0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13]))+0.4330127018922193*fEdge[12]-0.25*fEdge[11]-0.4330127018922193*fEdge[10]+0.4330127018922193*(fEdge[9]+fEdge[8])-0.25*(fEdge[7]+fEdge[6])+0.25*fEdge[5]+0.4330127018922193*fEdge[4]-0.25*fEdge[3]+0.25*(fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[3] = 0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13])-0.4330127018922193*fSkin[12]-0.25*fSkin[11]+0.4330127018922193*fSkin[10]-0.4330127018922193*(fSkin[9]+fSkin[8])-0.25*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.4330127018922193*fSkin[4]-0.25*fSkin[3]+0.25*(fSkin[2]+fSkin[1]+fSkin[0]); 
  } 
  if (alphaDrSurf_l[4]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = 0.4330127018922193*fEdge[15]-0.4330127018922193*(fEdge[14]+fEdge[13])+0.4330127018922193*fEdge[12]+0.25*fEdge[11]+0.4330127018922193*fEdge[10]-0.4330127018922193*(fEdge[9]+fEdge[8])-0.25*(fEdge[7]+fEdge[6])+0.25*fEdge[5]+0.4330127018922193*fEdge[4]+0.25*fEdge[3]-0.25*(fEdge[2]+fEdge[1])+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[4] = (-0.4330127018922193*fSkin[15])+0.4330127018922193*(fSkin[14]+fSkin[13])-0.4330127018922193*fSkin[12]+0.25*fSkin[11]-0.4330127018922193*fSkin[10]+0.4330127018922193*(fSkin[9]+fSkin[8])-0.25*(fSkin[7]+fSkin[6])+0.25*fSkin[5]-0.4330127018922193*fSkin[4]+0.25*fSkin[3]-0.25*(fSkin[2]+fSkin[1])+0.25*fSkin[0]; 
  } 
  if ((-alphaDrSurf_l[4])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = (-0.4330127018922193*(fEdge[15]+fEdge[14]))+0.4330127018922193*fEdge[13]-0.4330127018922193*fEdge[12]-0.25*fEdge[11]+0.4330127018922193*fEdge[10]-0.4330127018922193*fEdge[9]+0.4330127018922193*fEdge[8]-0.25*fEdge[7]+0.25*fEdge[6]-0.25*fEdge[5]+0.4330127018922193*fEdge[4]+0.25*fEdge[3]-0.25*fEdge[2]+0.25*(fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[5] = 0.4330127018922193*(fSkin[15]+fSkin[14])-0.4330127018922193*fSkin[13]+0.4330127018922193*fSkin[12]-0.25*fSkin[11]-0.4330127018922193*fSkin[10]+0.4330127018922193*fSkin[9]-0.4330127018922193*fSkin[8]-0.25*fSkin[7]+0.25*fSkin[6]-0.25*fSkin[5]-0.4330127018922193*fSkin[4]+0.25*fSkin[3]-0.25*fSkin[2]+0.25*(fSkin[1]+fSkin[0]); 
  } 
  if ((-alphaDrSurf_l[4])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = (-0.4330127018922193*fEdge[15])+0.4330127018922193*fEdge[14]-0.4330127018922193*(fEdge[13]+fEdge[12])-0.25*fEdge[11]+0.4330127018922193*(fEdge[10]+fEdge[9])-0.4330127018922193*fEdge[8]+0.25*fEdge[7]-0.25*(fEdge[6]+fEdge[5])+0.4330127018922193*fEdge[4]+0.25*(fEdge[3]+fEdge[2])-0.25*fEdge[1]+0.25*fEdge[0]; 
  } else { 
    fUpwindQuad_l[6] = 0.4330127018922193*fSkin[15]-0.4330127018922193*fSkin[14]+0.4330127018922193*(fSkin[13]+fSkin[12])-0.25*fSkin[11]-0.4330127018922193*(fSkin[10]+fSkin[9])+0.4330127018922193*fSkin[8]+0.25*fSkin[7]-0.25*(fSkin[6]+fSkin[5])-0.4330127018922193*fSkin[4]+0.25*(fSkin[3]+fSkin[2])-0.25*fSkin[1]+0.25*fSkin[0]; 
  } 
  if (alphaDrSurf_l[4]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = 0.4330127018922193*(fEdge[15]+fEdge[14]+fEdge[13]+fEdge[12])+0.25*fEdge[11]+0.4330127018922193*(fEdge[10]+fEdge[9]+fEdge[8])+0.25*(fEdge[7]+fEdge[6]+fEdge[5])+0.4330127018922193*fEdge[4]+0.25*(fEdge[3]+fEdge[2]+fEdge[1]+fEdge[0]); 
  } else { 
    fUpwindQuad_l[7] = (-0.4330127018922193*(fSkin[15]+fSkin[14]+fSkin[13]+fSkin[12]))+0.25*fSkin[11]-0.4330127018922193*(fSkin[10]+fSkin[9]+fSkin[8])+0.25*(fSkin[7]+fSkin[6]+fSkin[5])-0.4330127018922193*fSkin[4]+0.25*(fSkin[3]+fSkin[2]+fSkin[1]+fSkin[0]); 
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

  double Gdiff_l[8]; 
  double Ghat_l[8]; 


  Gdiff_l[0] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[3]*fSkin[12]+5.0*nuVtSqSum_l[3]*fEdge[12]+5.0*nuVtSqSum_l[2]*fSkin[9]+5.0*nuVtSqSum_l[2]*fEdge[9]+5.0*nuVtSqSum_l[1]*fSkin[8]+5.0*nuVtSqSum_l[1]*fEdge[8]+5.0*nuVtSqSum_l[0]*fSkin[4]+5.0*nuVtSqSum_l[0]*fEdge[4])-9.0*nuVtSqSum_l[3]*fSkin[5]+9.0*nuVtSqSum_l[3]*fEdge[5]+(9.0*fEdge[2]-9.0*fSkin[2])*nuVtSqSum_l[2]+(9.0*fEdge[1]-9.0*fSkin[1])*nuVtSqSum_l[1]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[0]); 
  Gdiff_l[1] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[2]*fSkin[12]+5.0*nuVtSqSum_l[2]*fEdge[12]+5.0*nuVtSqSum_l[3]*fSkin[9]+5.0*nuVtSqSum_l[3]*fEdge[9]+5.0*nuVtSqSum_l[0]*fSkin[8]+5.0*nuVtSqSum_l[0]*fEdge[8]+5.0*nuVtSqSum_l[1]*fSkin[4]+5.0*nuVtSqSum_l[1]*fEdge[4])-9.0*nuVtSqSum_l[2]*fSkin[5]+9.0*nuVtSqSum_l[2]*fEdge[5]+(9.0*fEdge[2]-9.0*fSkin[2])*nuVtSqSum_l[3]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[1]-9.0*nuVtSqSum_l[0]*fSkin[1]+9.0*nuVtSqSum_l[0]*fEdge[1]); 
  Gdiff_l[2] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[1]*fSkin[12]+5.0*nuVtSqSum_l[1]*fEdge[12]+5.0*nuVtSqSum_l[0]*fSkin[9]+5.0*nuVtSqSum_l[0]*fEdge[9]+5.0*nuVtSqSum_l[3]*fSkin[8]+5.0*nuVtSqSum_l[3]*fEdge[8]+5.0*nuVtSqSum_l[2]*fSkin[4]+5.0*nuVtSqSum_l[2]*fEdge[4])-9.0*nuVtSqSum_l[1]*fSkin[5]+9.0*nuVtSqSum_l[1]*fEdge[5]+(9.0*fEdge[1]-9.0*fSkin[1])*nuVtSqSum_l[3]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[2]-9.0*nuVtSqSum_l[0]*fSkin[2]+9.0*nuVtSqSum_l[0]*fEdge[2]); 
  Gdiff_l[3] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[3]*fSkin[15]+5.0*nuVtSqSum_l[3]*fEdge[15]+5.0*nuVtSqSum_l[2]*fSkin[14]+5.0*nuVtSqSum_l[2]*fEdge[14]+5.0*nuVtSqSum_l[1]*fSkin[13]+5.0*nuVtSqSum_l[1]*fEdge[13]+5.0*nuVtSqSum_l[0]*fSkin[10]+5.0*nuVtSqSum_l[0]*fEdge[10])-9.0*nuVtSqSum_l[3]*fSkin[11]+9.0*nuVtSqSum_l[3]*fEdge[11]-9.0*nuVtSqSum_l[2]*fSkin[7]+9.0*nuVtSqSum_l[2]*fEdge[7]-9.0*nuVtSqSum_l[1]*fSkin[6]+9.0*nuVtSqSum_l[1]*fEdge[6]-9.0*nuVtSqSum_l[0]*fSkin[3]+9.0*nuVtSqSum_l[0]*fEdge[3]); 
  Gdiff_l[4] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[0]*fSkin[12]+5.0*nuVtSqSum_l[0]*fEdge[12]+5.0*nuVtSqSum_l[1]*fSkin[9]+5.0*nuVtSqSum_l[1]*fEdge[9]+5.0*nuVtSqSum_l[2]*fSkin[8]+5.0*nuVtSqSum_l[2]*fEdge[8]+5.0*nuVtSqSum_l[3]*fSkin[4]+5.0*nuVtSqSum_l[3]*fEdge[4])-9.0*nuVtSqSum_l[0]*fSkin[5]+9.0*nuVtSqSum_l[0]*fEdge[5]+(9.0*fEdge[0]-9.0*fSkin[0])*nuVtSqSum_l[3]+(9.0*fEdge[1]-9.0*fSkin[1])*nuVtSqSum_l[2]-9.0*nuVtSqSum_l[1]*fSkin[2]+9.0*nuVtSqSum_l[1]*fEdge[2]); 
  Gdiff_l[5] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[2]*fSkin[15]+5.0*nuVtSqSum_l[2]*fEdge[15]+5.0*nuVtSqSum_l[3]*fSkin[14]+5.0*nuVtSqSum_l[3]*fEdge[14]+5.0*nuVtSqSum_l[0]*fSkin[13]+5.0*nuVtSqSum_l[0]*fEdge[13]+5.0*nuVtSqSum_l[1]*fSkin[10]+5.0*nuVtSqSum_l[1]*fEdge[10])-9.0*nuVtSqSum_l[2]*fSkin[11]+9.0*nuVtSqSum_l[2]*fEdge[11]-9.0*nuVtSqSum_l[3]*fSkin[7]+9.0*nuVtSqSum_l[3]*fEdge[7]-9.0*nuVtSqSum_l[0]*fSkin[6]+9.0*nuVtSqSum_l[0]*fEdge[6]-9.0*nuVtSqSum_l[1]*fSkin[3]+9.0*nuVtSqSum_l[1]*fEdge[3]); 
  Gdiff_l[6] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[1]*fSkin[15]+5.0*nuVtSqSum_l[1]*fEdge[15]+5.0*nuVtSqSum_l[0]*fSkin[14]+5.0*nuVtSqSum_l[0]*fEdge[14]+5.0*nuVtSqSum_l[3]*fSkin[13]+5.0*nuVtSqSum_l[3]*fEdge[13]+5.0*nuVtSqSum_l[2]*fSkin[10]+5.0*nuVtSqSum_l[2]*fEdge[10])-9.0*nuVtSqSum_l[1]*fSkin[11]+9.0*nuVtSqSum_l[1]*fEdge[11]-9.0*nuVtSqSum_l[0]*fSkin[7]+9.0*nuVtSqSum_l[0]*fEdge[7]-9.0*nuVtSqSum_l[3]*fSkin[6]+9.0*nuVtSqSum_l[3]*fEdge[6]-9.0*nuVtSqSum_l[2]*fSkin[3]+9.0*nuVtSqSum_l[2]*fEdge[3]); 
  Gdiff_l[7] = -0.04419417382415917*(1.732050807568877*(5.0*nuVtSqSum_l[0]*fSkin[15]+5.0*nuVtSqSum_l[0]*fEdge[15]+5.0*nuVtSqSum_l[1]*fSkin[14]+5.0*nuVtSqSum_l[1]*fEdge[14]+5.0*nuVtSqSum_l[2]*fSkin[13]+5.0*nuVtSqSum_l[2]*fEdge[13]+5.0*nuVtSqSum_l[3]*fSkin[10]+5.0*nuVtSqSum_l[3]*fEdge[10])-9.0*nuVtSqSum_l[0]*fSkin[11]+9.0*nuVtSqSum_l[0]*fEdge[11]-9.0*nuVtSqSum_l[1]*fSkin[7]+9.0*nuVtSqSum_l[1]*fEdge[7]-9.0*nuVtSqSum_l[2]*fSkin[6]+9.0*nuVtSqSum_l[2]*fEdge[6]+(9.0*fEdge[3]-9.0*fSkin[3])*nuVtSqSum_l[3]); 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[4]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[2]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[1]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[2]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[1]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[4]+0.3535533905932737*fUpwind_l[0]*alphaDrSurf_l[4]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[2]+0.3535533905932737*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[4]*fUpwind_l[5]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.3535533905932737*alphaDrSurf_l[0]*fUpwind_l[7]+0.3535533905932737*alphaDrSurf_l[1]*fUpwind_l[6]+0.3535533905932737*alphaDrSurf_l[2]*fUpwind_l[5]+0.3535533905932737*fUpwind_l[3]*alphaDrSurf_l[4]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[3]*rdv2; 
  out[4] += 0.75*nuVtSqSum_r[3]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum_r[3]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum_r[2]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum_r[1]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_l[6]*rdv2; 
  out[8] += 0.75*nuVtSqSum_r[2]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum_r[2]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum_r[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum_r[1]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_l[1]*rdv2; 
  out[9] += 0.75*nuVtSqSum_r[1]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum_r[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum_r[2]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_l[2]*rdv2; 
  out[10] += 0.75*nuVtSqSum_r[3]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum_r[3]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum_r[2]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[6]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_l[7]*rdv2; 
  out[12] += 0.75*nuVtSqSum_r[0]*fSkin[12]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[9]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[5]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum_r[3]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum_r[2]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat_l[4]*rdv2; 
  out[13] += 0.75*nuVtSqSum_r[2]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum_r[2]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum_r[3]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_l[5]*rdv2; 
  out[14] += 0.75*nuVtSqSum_r[1]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_r[0]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum_r[3]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4+0.4330127018922193*nuVtSqSum_r[2]*fSkin[3]*rdvSq4-1.224744871391589*Ghat_l[6]*rdv2; 
  out[15] += 0.75*nuVtSqSum_r[0]*fSkin[15]*rdvSq4+0.75*nuVtSqSum_r[1]*fSkin[14]*rdvSq4+0.75*nuVtSqSum_r[2]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum_r[0]*fSkin[11]*rdvSq4+0.75*nuVtSqSum_r[3]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum_r[1]*fSkin[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4+0.4330127018922193*nuVtSqSum_r[2]*fSkin[6]*rdvSq4+0.4330127018922193*fSkin[3]*nuVtSqSum_r[3]*rdvSq4-1.224744871391589*Ghat_l[7]*rdv2; 
  } 
} 
