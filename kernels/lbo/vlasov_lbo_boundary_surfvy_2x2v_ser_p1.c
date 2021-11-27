#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_2x2v_p1_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvy_2x2v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[4]:         Cell-center coordinates. 
  // dxv[4]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[8]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *sumNuUy = &nuUSum[4]; 

  double alphaDrSurf[8] = {0.0}; 
  double fUpwindQuad[8] = {0.0};
  double fUpwind[8] = {0.0};;
  double Ghat[8] = {0.0}; 
  double Gdiff[8] = {0.0}; 
  double Gdiff2[8] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x2v_p1_surfvy_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_2x2v_p1_surfvy_quad_0(-1, fEdge); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x2v_p1_surfvy_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_2x2v_p1_surfvy_quad_1(-1, fEdge); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x2v_p1_surfvy_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_2x2v_p1_surfvy_quad_2(-1, fEdge); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x2v_p1_surfvy_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_2x2v_p1_surfvy_quad_3(-1, fEdge); 
  } 
  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_2x2v_p1_surfvy_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_2x2v_p1_surfvy_quad_4(-1, fEdge); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x2v_p1_surfvy_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_2x2v_p1_surfvy_quad_5(-1, fEdge); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x2v_p1_surfvy_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_2x2v_p1_surfvy_quad_6(-1, fEdge); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x2v_p1_surfvy_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_2x2v_p1_surfvy_quad_7(-1, fEdge); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Gdiff2[0] = 0.2041241452319315*nuVtSqSum[3]*fSkin[12]-0.2041241452319315*nuVtSqSum[3]*fEdge[12]+0.2041241452319315*nuVtSqSum[2]*fSkin[9]-0.2041241452319315*nuVtSqSum[2]*fEdge[9]+0.2041241452319315*nuVtSqSum[1]*fSkin[8]-0.2041241452319315*nuVtSqSum[1]*fEdge[8]+0.1767766952966368*nuVtSqSum[3]*fSkin[5]+0.1767766952966368*nuVtSqSum[3]*fEdge[5]+0.2041241452319315*nuVtSqSum[0]*fSkin[4]-0.2041241452319315*nuVtSqSum[0]*fEdge[4]+0.1767766952966368*fSkin[2]*nuVtSqSum[2]+0.1767766952966368*fEdge[2]*nuVtSqSum[2]+0.1767766952966368*fSkin[1]*nuVtSqSum[1]+0.1767766952966368*fEdge[1]*nuVtSqSum[1]+0.1767766952966368*fSkin[0]*nuVtSqSum[0]+0.1767766952966368*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.2041241452319315*nuVtSqSum[2]*fSkin[12]-0.2041241452319315*nuVtSqSum[2]*fEdge[12]+0.2041241452319315*nuVtSqSum[3]*fSkin[9]-0.2041241452319315*nuVtSqSum[3]*fEdge[9]+0.2041241452319315*nuVtSqSum[0]*fSkin[8]-0.2041241452319315*nuVtSqSum[0]*fEdge[8]+0.1767766952966368*nuVtSqSum[2]*fSkin[5]+0.1767766952966368*nuVtSqSum[2]*fEdge[5]+0.2041241452319315*nuVtSqSum[1]*fSkin[4]-0.2041241452319315*nuVtSqSum[1]*fEdge[4]+0.1767766952966368*fSkin[2]*nuVtSqSum[3]+0.1767766952966368*fEdge[2]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[1]+0.1767766952966368*fEdge[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fSkin[1]+0.1767766952966368*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.2041241452319315*nuVtSqSum[1]*fSkin[12]-0.2041241452319315*nuVtSqSum[1]*fEdge[12]+0.2041241452319315*nuVtSqSum[0]*fSkin[9]-0.2041241452319315*nuVtSqSum[0]*fEdge[9]+0.2041241452319315*nuVtSqSum[3]*fSkin[8]-0.2041241452319315*nuVtSqSum[3]*fEdge[8]+0.1767766952966368*nuVtSqSum[1]*fSkin[5]+0.1767766952966368*nuVtSqSum[1]*fEdge[5]+0.2041241452319315*nuVtSqSum[2]*fSkin[4]-0.2041241452319315*nuVtSqSum[2]*fEdge[4]+0.1767766952966368*fSkin[1]*nuVtSqSum[3]+0.1767766952966368*fEdge[1]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[2]+0.1767766952966368*fEdge[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fSkin[2]+0.1767766952966368*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = 0.2041241452319315*nuVtSqSum[3]*fSkin[15]-0.2041241452319315*nuVtSqSum[3]*fEdge[15]+0.2041241452319315*nuVtSqSum[2]*fSkin[14]-0.2041241452319315*nuVtSqSum[2]*fEdge[14]+0.2041241452319315*nuVtSqSum[1]*fSkin[13]-0.2041241452319315*nuVtSqSum[1]*fEdge[13]+0.1767766952966368*nuVtSqSum[3]*fSkin[11]+0.1767766952966368*nuVtSqSum[3]*fEdge[11]+0.2041241452319315*nuVtSqSum[0]*fSkin[10]-0.2041241452319315*nuVtSqSum[0]*fEdge[10]+0.1767766952966368*nuVtSqSum[2]*fSkin[7]+0.1767766952966368*nuVtSqSum[2]*fEdge[7]+0.1767766952966368*nuVtSqSum[1]*fSkin[6]+0.1767766952966368*nuVtSqSum[1]*fEdge[6]+0.1767766952966368*nuVtSqSum[0]*fSkin[3]+0.1767766952966368*nuVtSqSum[0]*fEdge[3]; 
  Gdiff2[4] = 0.2041241452319315*nuVtSqSum[0]*fSkin[12]-0.2041241452319315*nuVtSqSum[0]*fEdge[12]+0.2041241452319315*nuVtSqSum[1]*fSkin[9]-0.2041241452319315*nuVtSqSum[1]*fEdge[9]+0.2041241452319315*nuVtSqSum[2]*fSkin[8]-0.2041241452319315*nuVtSqSum[2]*fEdge[8]+0.1767766952966368*nuVtSqSum[0]*fSkin[5]+0.1767766952966368*nuVtSqSum[0]*fEdge[5]+0.2041241452319315*nuVtSqSum[3]*fSkin[4]-0.2041241452319315*nuVtSqSum[3]*fEdge[4]+0.1767766952966368*fSkin[0]*nuVtSqSum[3]+0.1767766952966368*fEdge[0]*nuVtSqSum[3]+0.1767766952966368*fSkin[1]*nuVtSqSum[2]+0.1767766952966368*fEdge[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fSkin[2]+0.1767766952966368*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[5] = 0.2041241452319315*nuVtSqSum[2]*fSkin[15]-0.2041241452319315*nuVtSqSum[2]*fEdge[15]+0.2041241452319315*nuVtSqSum[3]*fSkin[14]-0.2041241452319315*nuVtSqSum[3]*fEdge[14]+0.2041241452319315*nuVtSqSum[0]*fSkin[13]-0.2041241452319315*nuVtSqSum[0]*fEdge[13]+0.1767766952966368*nuVtSqSum[2]*fSkin[11]+0.1767766952966368*nuVtSqSum[2]*fEdge[11]+0.2041241452319315*nuVtSqSum[1]*fSkin[10]-0.2041241452319315*nuVtSqSum[1]*fEdge[10]+0.1767766952966368*nuVtSqSum[3]*fSkin[7]+0.1767766952966368*nuVtSqSum[3]*fEdge[7]+0.1767766952966368*nuVtSqSum[0]*fSkin[6]+0.1767766952966368*nuVtSqSum[0]*fEdge[6]+0.1767766952966368*nuVtSqSum[1]*fSkin[3]+0.1767766952966368*nuVtSqSum[1]*fEdge[3]; 
  Gdiff2[6] = 0.2041241452319315*nuVtSqSum[1]*fSkin[15]-0.2041241452319315*nuVtSqSum[1]*fEdge[15]+0.2041241452319315*nuVtSqSum[0]*fSkin[14]-0.2041241452319315*nuVtSqSum[0]*fEdge[14]+0.2041241452319315*nuVtSqSum[3]*fSkin[13]-0.2041241452319315*nuVtSqSum[3]*fEdge[13]+0.1767766952966368*nuVtSqSum[1]*fSkin[11]+0.1767766952966368*nuVtSqSum[1]*fEdge[11]+0.2041241452319315*nuVtSqSum[2]*fSkin[10]-0.2041241452319315*nuVtSqSum[2]*fEdge[10]+0.1767766952966368*nuVtSqSum[0]*fSkin[7]+0.1767766952966368*nuVtSqSum[0]*fEdge[7]+0.1767766952966368*nuVtSqSum[3]*fSkin[6]+0.1767766952966368*nuVtSqSum[3]*fEdge[6]+0.1767766952966368*nuVtSqSum[2]*fSkin[3]+0.1767766952966368*nuVtSqSum[2]*fEdge[3]; 
  Gdiff2[7] = 0.2041241452319315*nuVtSqSum[0]*fSkin[15]-0.2041241452319315*nuVtSqSum[0]*fEdge[15]+0.2041241452319315*nuVtSqSum[1]*fSkin[14]-0.2041241452319315*nuVtSqSum[1]*fEdge[14]+0.2041241452319315*nuVtSqSum[2]*fSkin[13]-0.2041241452319315*nuVtSqSum[2]*fEdge[13]+0.1767766952966368*nuVtSqSum[0]*fSkin[11]+0.1767766952966368*nuVtSqSum[0]*fEdge[11]+0.2041241452319315*nuVtSqSum[3]*fSkin[10]-0.2041241452319315*nuVtSqSum[3]*fEdge[10]+0.1767766952966368*nuVtSqSum[1]*fSkin[7]+0.1767766952966368*nuVtSqSum[1]*fEdge[7]+0.1767766952966368*nuVtSqSum[2]*fSkin[6]+0.1767766952966368*nuVtSqSum[2]*fEdge[6]+0.1767766952966368*fSkin[3]*nuVtSqSum[3]+0.1767766952966368*fEdge[3]*nuVtSqSum[3]; 

  Gdiff[0] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[12])-0.3827327723098713*nuVtSqSum[3]*fEdge[12]-0.3827327723098713*nuVtSqSum[2]*fSkin[9]-0.3827327723098713*nuVtSqSum[2]*fEdge[9]-0.3827327723098713*nuVtSqSum[1]*fSkin[8]-0.3827327723098713*nuVtSqSum[1]*fEdge[8]-0.3977475644174328*nuVtSqSum[3]*fSkin[5]+0.3977475644174328*nuVtSqSum[3]*fEdge[5]-0.3827327723098713*nuVtSqSum[0]*fSkin[4]-0.3827327723098713*nuVtSqSum[0]*fEdge[4]-0.3977475644174328*fSkin[2]*nuVtSqSum[2]+0.3977475644174328*fEdge[2]*nuVtSqSum[2]-0.3977475644174328*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.3977475644174328*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[12])-0.3827327723098713*nuVtSqSum[2]*fEdge[12]-0.3827327723098713*nuVtSqSum[3]*fSkin[9]-0.3827327723098713*nuVtSqSum[3]*fEdge[9]-0.3827327723098713*nuVtSqSum[0]*fSkin[8]-0.3827327723098713*nuVtSqSum[0]*fEdge[8]-0.3977475644174328*nuVtSqSum[2]*fSkin[5]+0.3977475644174328*nuVtSqSum[2]*fEdge[5]-0.3827327723098713*nuVtSqSum[1]*fSkin[4]-0.3827327723098713*nuVtSqSum[1]*fEdge[4]-0.3977475644174328*fSkin[2]*nuVtSqSum[3]+0.3977475644174328*fEdge[2]*nuVtSqSum[3]-0.3977475644174328*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.3977475644174328*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[12])-0.3827327723098713*nuVtSqSum[1]*fEdge[12]-0.3827327723098713*nuVtSqSum[0]*fSkin[9]-0.3827327723098713*nuVtSqSum[0]*fEdge[9]-0.3827327723098713*nuVtSqSum[3]*fSkin[8]-0.3827327723098713*nuVtSqSum[3]*fEdge[8]-0.3977475644174328*nuVtSqSum[1]*fSkin[5]+0.3977475644174328*nuVtSqSum[1]*fEdge[5]-0.3827327723098713*nuVtSqSum[2]*fSkin[4]-0.3827327723098713*nuVtSqSum[2]*fEdge[4]-0.3977475644174328*fSkin[1]*nuVtSqSum[3]+0.3977475644174328*fEdge[1]*nuVtSqSum[3]-0.3977475644174328*fSkin[0]*nuVtSqSum[2]+0.3977475644174328*fEdge[0]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[0]*fSkin[2]+0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[15])-0.3827327723098713*nuVtSqSum[3]*fEdge[15]-0.3827327723098713*nuVtSqSum[2]*fSkin[14]-0.3827327723098713*nuVtSqSum[2]*fEdge[14]-0.3827327723098713*nuVtSqSum[1]*fSkin[13]-0.3827327723098713*nuVtSqSum[1]*fEdge[13]-0.3977475644174328*nuVtSqSum[3]*fSkin[11]+0.3977475644174328*nuVtSqSum[3]*fEdge[11]-0.3827327723098713*nuVtSqSum[0]*fSkin[10]-0.3827327723098713*nuVtSqSum[0]*fEdge[10]-0.3977475644174328*nuVtSqSum[2]*fSkin[7]+0.3977475644174328*nuVtSqSum[2]*fEdge[7]-0.3977475644174328*nuVtSqSum[1]*fSkin[6]+0.3977475644174328*nuVtSqSum[1]*fEdge[6]-0.3977475644174328*nuVtSqSum[0]*fSkin[3]+0.3977475644174328*nuVtSqSum[0]*fEdge[3]; 
  Gdiff[4] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[12])-0.3827327723098713*nuVtSqSum[0]*fEdge[12]-0.3827327723098713*nuVtSqSum[1]*fSkin[9]-0.3827327723098713*nuVtSqSum[1]*fEdge[9]-0.3827327723098713*nuVtSqSum[2]*fSkin[8]-0.3827327723098713*nuVtSqSum[2]*fEdge[8]-0.3977475644174328*nuVtSqSum[0]*fSkin[5]+0.3977475644174328*nuVtSqSum[0]*fEdge[5]-0.3827327723098713*nuVtSqSum[3]*fSkin[4]-0.3827327723098713*nuVtSqSum[3]*fEdge[4]-0.3977475644174328*fSkin[0]*nuVtSqSum[3]+0.3977475644174328*fEdge[0]*nuVtSqSum[3]-0.3977475644174328*fSkin[1]*nuVtSqSum[2]+0.3977475644174328*fEdge[1]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[1]*fSkin[2]+0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[5] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[15])-0.3827327723098713*nuVtSqSum[2]*fEdge[15]-0.3827327723098713*nuVtSqSum[3]*fSkin[14]-0.3827327723098713*nuVtSqSum[3]*fEdge[14]-0.3827327723098713*nuVtSqSum[0]*fSkin[13]-0.3827327723098713*nuVtSqSum[0]*fEdge[13]-0.3977475644174328*nuVtSqSum[2]*fSkin[11]+0.3977475644174328*nuVtSqSum[2]*fEdge[11]-0.3827327723098713*nuVtSqSum[1]*fSkin[10]-0.3827327723098713*nuVtSqSum[1]*fEdge[10]-0.3977475644174328*nuVtSqSum[3]*fSkin[7]+0.3977475644174328*nuVtSqSum[3]*fEdge[7]-0.3977475644174328*nuVtSqSum[0]*fSkin[6]+0.3977475644174328*nuVtSqSum[0]*fEdge[6]-0.3977475644174328*nuVtSqSum[1]*fSkin[3]+0.3977475644174328*nuVtSqSum[1]*fEdge[3]; 
  Gdiff[6] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[15])-0.3827327723098713*nuVtSqSum[1]*fEdge[15]-0.3827327723098713*nuVtSqSum[0]*fSkin[14]-0.3827327723098713*nuVtSqSum[0]*fEdge[14]-0.3827327723098713*nuVtSqSum[3]*fSkin[13]-0.3827327723098713*nuVtSqSum[3]*fEdge[13]-0.3977475644174328*nuVtSqSum[1]*fSkin[11]+0.3977475644174328*nuVtSqSum[1]*fEdge[11]-0.3827327723098713*nuVtSqSum[2]*fSkin[10]-0.3827327723098713*nuVtSqSum[2]*fEdge[10]-0.3977475644174328*nuVtSqSum[0]*fSkin[7]+0.3977475644174328*nuVtSqSum[0]*fEdge[7]-0.3977475644174328*nuVtSqSum[3]*fSkin[6]+0.3977475644174328*nuVtSqSum[3]*fEdge[6]-0.3977475644174328*nuVtSqSum[2]*fSkin[3]+0.3977475644174328*nuVtSqSum[2]*fEdge[3]; 
  Gdiff[7] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[15])-0.3827327723098713*nuVtSqSum[0]*fEdge[15]-0.3827327723098713*nuVtSqSum[1]*fSkin[14]-0.3827327723098713*nuVtSqSum[1]*fEdge[14]-0.3827327723098713*nuVtSqSum[2]*fSkin[13]-0.3827327723098713*nuVtSqSum[2]*fEdge[13]-0.3977475644174328*nuVtSqSum[0]*fSkin[11]+0.3977475644174328*nuVtSqSum[0]*fEdge[11]-0.3827327723098713*nuVtSqSum[3]*fSkin[10]-0.3827327723098713*nuVtSqSum[3]*fEdge[10]-0.3977475644174328*nuVtSqSum[1]*fSkin[7]+0.3977475644174328*nuVtSqSum[1]*fEdge[7]-0.3977475644174328*nuVtSqSum[2]*fSkin[6]+0.3977475644174328*nuVtSqSum[2]*fEdge[6]-0.3977475644174328*fSkin[3]*nuVtSqSum[3]+0.3977475644174328*fEdge[3]*nuVtSqSum[3]; 

  Ghat[0] = (-1.0*Gdiff[0]*rdv2)+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = (-1.0*Gdiff[1]*rdv2)+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = (-1.0*Gdiff[2]*rdv2)+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = (-1.0*Gdiff[3]*rdv2)+0.3535533905932737*alphaDrSurf[4]*fUpwind[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = (-1.0*Gdiff[4]*rdv2)+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = (-1.0*Gdiff[5]*rdv2)+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = (-1.0*Gdiff[6]*rdv2)+0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = (-1.0*Gdiff[7]*rdv2)+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 0.75*nuVtSqSum[3]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[1]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[0]*rdvSq4+1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.75*nuVtSqSum[2]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[1]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 0.75*nuVtSqSum[1]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[4]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[2]*rdvSq4+1.224744871391589*Gdiff2[2]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 0.75*nuVtSqSum[3]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[6]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[3]*rdvSq4+1.224744871391589*Gdiff2[3]*rdvSq4-1.224744871391589*Ghat[3]*rdv2; 
  out[11] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 0.75*nuVtSqSum[0]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[8]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[4]*rdvSq4+1.224744871391589*Gdiff2[4]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[3]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 0.75*nuVtSqSum[2]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2[5]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 0.75*nuVtSqSum[1]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2[6]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[3]*rdvSq4-1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 0.75*nuVtSqSum[0]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[7]*rdvSq4+1.224744871391589*Gdiff2[7]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[6]*rdvSq4-0.4330127018922193*fSkin[3]*nuVtSqSum[3]*rdvSq4-1.224744871391589*Ghat[7]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.7071067811865475*(nuSum[0]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[0]); 
  alphaDrSurf[1] = 0.7071067811865475*(nuSum[1]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[1]); 
  alphaDrSurf[2] = 0.7071067811865475*(nuSum[2]*(2.0*w[3]+dxv[3])-2.0*sumNuUy[2]); 
  alphaDrSurf[4] = 0.7071067811865475*(2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]); 

  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_2x2v_p1_surfvy_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x2v_p1_surfvy_quad_0(-1, fSkin); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[1] = ser_2x2v_p1_surfvy_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x2v_p1_surfvy_quad_1(-1, fSkin); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_2x2v_p1_surfvy_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x2v_p1_surfvy_quad_2(-1, fSkin); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[3] = ser_2x2v_p1_surfvy_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x2v_p1_surfvy_quad_3(-1, fSkin); 
  } 
  if (alphaDrSurf[4]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[4] = ser_2x2v_p1_surfvy_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_2x2v_p1_surfvy_quad_4(-1, fSkin); 
  } 
  if ((-alphaDrSurf[4])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[5] = ser_2x2v_p1_surfvy_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_2x2v_p1_surfvy_quad_5(-1, fSkin); 
  } 
  if ((-alphaDrSurf[4])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[6] = ser_2x2v_p1_surfvy_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_2x2v_p1_surfvy_quad_6(-1, fSkin); 
  } 
  if (alphaDrSurf[4]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] < 0) { 
    fUpwindQuad[7] = ser_2x2v_p1_surfvy_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_2x2v_p1_surfvy_quad_7(-1, fSkin); 
  } 

  fUpwind[0] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[5] = 0.3535533905932737*(fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[6] = 0.3535533905932737*(fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.3535533905932737*(fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 

  Gdiff2[0] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[12])+0.2041241452319315*nuVtSqSum[3]*fEdge[12]-0.2041241452319315*nuVtSqSum[2]*fSkin[9]+0.2041241452319315*nuVtSqSum[2]*fEdge[9]-0.2041241452319315*nuVtSqSum[1]*fSkin[8]+0.2041241452319315*nuVtSqSum[1]*fEdge[8]+0.1767766952966368*nuVtSqSum[3]*fSkin[5]+0.1767766952966368*nuVtSqSum[3]*fEdge[5]-0.2041241452319315*nuVtSqSum[0]*fSkin[4]+0.2041241452319315*nuVtSqSum[0]*fEdge[4]+0.1767766952966368*fSkin[2]*nuVtSqSum[2]+0.1767766952966368*fEdge[2]*nuVtSqSum[2]+0.1767766952966368*fSkin[1]*nuVtSqSum[1]+0.1767766952966368*fEdge[1]*nuVtSqSum[1]+0.1767766952966368*fSkin[0]*nuVtSqSum[0]+0.1767766952966368*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[12])+0.2041241452319315*nuVtSqSum[2]*fEdge[12]-0.2041241452319315*nuVtSqSum[3]*fSkin[9]+0.2041241452319315*nuVtSqSum[3]*fEdge[9]-0.2041241452319315*nuVtSqSum[0]*fSkin[8]+0.2041241452319315*nuVtSqSum[0]*fEdge[8]+0.1767766952966368*nuVtSqSum[2]*fSkin[5]+0.1767766952966368*nuVtSqSum[2]*fEdge[5]-0.2041241452319315*nuVtSqSum[1]*fSkin[4]+0.2041241452319315*nuVtSqSum[1]*fEdge[4]+0.1767766952966368*fSkin[2]*nuVtSqSum[3]+0.1767766952966368*fEdge[2]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[1]+0.1767766952966368*fEdge[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fSkin[1]+0.1767766952966368*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[12])+0.2041241452319315*nuVtSqSum[1]*fEdge[12]-0.2041241452319315*nuVtSqSum[0]*fSkin[9]+0.2041241452319315*nuVtSqSum[0]*fEdge[9]-0.2041241452319315*nuVtSqSum[3]*fSkin[8]+0.2041241452319315*nuVtSqSum[3]*fEdge[8]+0.1767766952966368*nuVtSqSum[1]*fSkin[5]+0.1767766952966368*nuVtSqSum[1]*fEdge[5]-0.2041241452319315*nuVtSqSum[2]*fSkin[4]+0.2041241452319315*nuVtSqSum[2]*fEdge[4]+0.1767766952966368*fSkin[1]*nuVtSqSum[3]+0.1767766952966368*fEdge[1]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[2]+0.1767766952966368*fEdge[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fSkin[2]+0.1767766952966368*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[15])+0.2041241452319315*nuVtSqSum[3]*fEdge[15]-0.2041241452319315*nuVtSqSum[2]*fSkin[14]+0.2041241452319315*nuVtSqSum[2]*fEdge[14]-0.2041241452319315*nuVtSqSum[1]*fSkin[13]+0.2041241452319315*nuVtSqSum[1]*fEdge[13]+0.1767766952966368*nuVtSqSum[3]*fSkin[11]+0.1767766952966368*nuVtSqSum[3]*fEdge[11]-0.2041241452319315*nuVtSqSum[0]*fSkin[10]+0.2041241452319315*nuVtSqSum[0]*fEdge[10]+0.1767766952966368*nuVtSqSum[2]*fSkin[7]+0.1767766952966368*nuVtSqSum[2]*fEdge[7]+0.1767766952966368*nuVtSqSum[1]*fSkin[6]+0.1767766952966368*nuVtSqSum[1]*fEdge[6]+0.1767766952966368*nuVtSqSum[0]*fSkin[3]+0.1767766952966368*nuVtSqSum[0]*fEdge[3]; 
  Gdiff2[4] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[12])+0.2041241452319315*nuVtSqSum[0]*fEdge[12]-0.2041241452319315*nuVtSqSum[1]*fSkin[9]+0.2041241452319315*nuVtSqSum[1]*fEdge[9]-0.2041241452319315*nuVtSqSum[2]*fSkin[8]+0.2041241452319315*nuVtSqSum[2]*fEdge[8]+0.1767766952966368*nuVtSqSum[0]*fSkin[5]+0.1767766952966368*nuVtSqSum[0]*fEdge[5]-0.2041241452319315*nuVtSqSum[3]*fSkin[4]+0.2041241452319315*nuVtSqSum[3]*fEdge[4]+0.1767766952966368*fSkin[0]*nuVtSqSum[3]+0.1767766952966368*fEdge[0]*nuVtSqSum[3]+0.1767766952966368*fSkin[1]*nuVtSqSum[2]+0.1767766952966368*fEdge[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fSkin[2]+0.1767766952966368*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[5] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[15])+0.2041241452319315*nuVtSqSum[2]*fEdge[15]-0.2041241452319315*nuVtSqSum[3]*fSkin[14]+0.2041241452319315*nuVtSqSum[3]*fEdge[14]-0.2041241452319315*nuVtSqSum[0]*fSkin[13]+0.2041241452319315*nuVtSqSum[0]*fEdge[13]+0.1767766952966368*nuVtSqSum[2]*fSkin[11]+0.1767766952966368*nuVtSqSum[2]*fEdge[11]-0.2041241452319315*nuVtSqSum[1]*fSkin[10]+0.2041241452319315*nuVtSqSum[1]*fEdge[10]+0.1767766952966368*nuVtSqSum[3]*fSkin[7]+0.1767766952966368*nuVtSqSum[3]*fEdge[7]+0.1767766952966368*nuVtSqSum[0]*fSkin[6]+0.1767766952966368*nuVtSqSum[0]*fEdge[6]+0.1767766952966368*nuVtSqSum[1]*fSkin[3]+0.1767766952966368*nuVtSqSum[1]*fEdge[3]; 
  Gdiff2[6] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[15])+0.2041241452319315*nuVtSqSum[1]*fEdge[15]-0.2041241452319315*nuVtSqSum[0]*fSkin[14]+0.2041241452319315*nuVtSqSum[0]*fEdge[14]-0.2041241452319315*nuVtSqSum[3]*fSkin[13]+0.2041241452319315*nuVtSqSum[3]*fEdge[13]+0.1767766952966368*nuVtSqSum[1]*fSkin[11]+0.1767766952966368*nuVtSqSum[1]*fEdge[11]-0.2041241452319315*nuVtSqSum[2]*fSkin[10]+0.2041241452319315*nuVtSqSum[2]*fEdge[10]+0.1767766952966368*nuVtSqSum[0]*fSkin[7]+0.1767766952966368*nuVtSqSum[0]*fEdge[7]+0.1767766952966368*nuVtSqSum[3]*fSkin[6]+0.1767766952966368*nuVtSqSum[3]*fEdge[6]+0.1767766952966368*nuVtSqSum[2]*fSkin[3]+0.1767766952966368*nuVtSqSum[2]*fEdge[3]; 
  Gdiff2[7] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[15])+0.2041241452319315*nuVtSqSum[0]*fEdge[15]-0.2041241452319315*nuVtSqSum[1]*fSkin[14]+0.2041241452319315*nuVtSqSum[1]*fEdge[14]-0.2041241452319315*nuVtSqSum[2]*fSkin[13]+0.2041241452319315*nuVtSqSum[2]*fEdge[13]+0.1767766952966368*nuVtSqSum[0]*fSkin[11]+0.1767766952966368*nuVtSqSum[0]*fEdge[11]-0.2041241452319315*nuVtSqSum[3]*fSkin[10]+0.2041241452319315*nuVtSqSum[3]*fEdge[10]+0.1767766952966368*nuVtSqSum[1]*fSkin[7]+0.1767766952966368*nuVtSqSum[1]*fEdge[7]+0.1767766952966368*nuVtSqSum[2]*fSkin[6]+0.1767766952966368*nuVtSqSum[2]*fEdge[6]+0.1767766952966368*fSkin[3]*nuVtSqSum[3]+0.1767766952966368*fEdge[3]*nuVtSqSum[3]; 

  Gdiff[0] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[12])-0.3827327723098713*nuVtSqSum[3]*fEdge[12]-0.3827327723098713*nuVtSqSum[2]*fSkin[9]-0.3827327723098713*nuVtSqSum[2]*fEdge[9]-0.3827327723098713*nuVtSqSum[1]*fSkin[8]-0.3827327723098713*nuVtSqSum[1]*fEdge[8]+0.3977475644174328*nuVtSqSum[3]*fSkin[5]-0.3977475644174328*nuVtSqSum[3]*fEdge[5]-0.3827327723098713*nuVtSqSum[0]*fSkin[4]-0.3827327723098713*nuVtSqSum[0]*fEdge[4]+0.3977475644174328*fSkin[2]*nuVtSqSum[2]-0.3977475644174328*fEdge[2]*nuVtSqSum[2]+0.3977475644174328*fSkin[1]*nuVtSqSum[1]-0.3977475644174328*fEdge[1]*nuVtSqSum[1]+0.3977475644174328*fSkin[0]*nuVtSqSum[0]-0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[12])-0.3827327723098713*nuVtSqSum[2]*fEdge[12]-0.3827327723098713*nuVtSqSum[3]*fSkin[9]-0.3827327723098713*nuVtSqSum[3]*fEdge[9]-0.3827327723098713*nuVtSqSum[0]*fSkin[8]-0.3827327723098713*nuVtSqSum[0]*fEdge[8]+0.3977475644174328*nuVtSqSum[2]*fSkin[5]-0.3977475644174328*nuVtSqSum[2]*fEdge[5]-0.3827327723098713*nuVtSqSum[1]*fSkin[4]-0.3827327723098713*nuVtSqSum[1]*fEdge[4]+0.3977475644174328*fSkin[2]*nuVtSqSum[3]-0.3977475644174328*fEdge[2]*nuVtSqSum[3]+0.3977475644174328*fSkin[0]*nuVtSqSum[1]-0.3977475644174328*fEdge[0]*nuVtSqSum[1]+0.3977475644174328*nuVtSqSum[0]*fSkin[1]-0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[12])-0.3827327723098713*nuVtSqSum[1]*fEdge[12]-0.3827327723098713*nuVtSqSum[0]*fSkin[9]-0.3827327723098713*nuVtSqSum[0]*fEdge[9]-0.3827327723098713*nuVtSqSum[3]*fSkin[8]-0.3827327723098713*nuVtSqSum[3]*fEdge[8]+0.3977475644174328*nuVtSqSum[1]*fSkin[5]-0.3977475644174328*nuVtSqSum[1]*fEdge[5]-0.3827327723098713*nuVtSqSum[2]*fSkin[4]-0.3827327723098713*nuVtSqSum[2]*fEdge[4]+0.3977475644174328*fSkin[1]*nuVtSqSum[3]-0.3977475644174328*fEdge[1]*nuVtSqSum[3]+0.3977475644174328*fSkin[0]*nuVtSqSum[2]-0.3977475644174328*fEdge[0]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[0]*fSkin[2]-0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[15])-0.3827327723098713*nuVtSqSum[3]*fEdge[15]-0.3827327723098713*nuVtSqSum[2]*fSkin[14]-0.3827327723098713*nuVtSqSum[2]*fEdge[14]-0.3827327723098713*nuVtSqSum[1]*fSkin[13]-0.3827327723098713*nuVtSqSum[1]*fEdge[13]+0.3977475644174328*nuVtSqSum[3]*fSkin[11]-0.3977475644174328*nuVtSqSum[3]*fEdge[11]-0.3827327723098713*nuVtSqSum[0]*fSkin[10]-0.3827327723098713*nuVtSqSum[0]*fEdge[10]+0.3977475644174328*nuVtSqSum[2]*fSkin[7]-0.3977475644174328*nuVtSqSum[2]*fEdge[7]+0.3977475644174328*nuVtSqSum[1]*fSkin[6]-0.3977475644174328*nuVtSqSum[1]*fEdge[6]+0.3977475644174328*nuVtSqSum[0]*fSkin[3]-0.3977475644174328*nuVtSqSum[0]*fEdge[3]; 
  Gdiff[4] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[12])-0.3827327723098713*nuVtSqSum[0]*fEdge[12]-0.3827327723098713*nuVtSqSum[1]*fSkin[9]-0.3827327723098713*nuVtSqSum[1]*fEdge[9]-0.3827327723098713*nuVtSqSum[2]*fSkin[8]-0.3827327723098713*nuVtSqSum[2]*fEdge[8]+0.3977475644174328*nuVtSqSum[0]*fSkin[5]-0.3977475644174328*nuVtSqSum[0]*fEdge[5]-0.3827327723098713*nuVtSqSum[3]*fSkin[4]-0.3827327723098713*nuVtSqSum[3]*fEdge[4]+0.3977475644174328*fSkin[0]*nuVtSqSum[3]-0.3977475644174328*fEdge[0]*nuVtSqSum[3]+0.3977475644174328*fSkin[1]*nuVtSqSum[2]-0.3977475644174328*fEdge[1]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[1]*fSkin[2]-0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[5] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[15])-0.3827327723098713*nuVtSqSum[2]*fEdge[15]-0.3827327723098713*nuVtSqSum[3]*fSkin[14]-0.3827327723098713*nuVtSqSum[3]*fEdge[14]-0.3827327723098713*nuVtSqSum[0]*fSkin[13]-0.3827327723098713*nuVtSqSum[0]*fEdge[13]+0.3977475644174328*nuVtSqSum[2]*fSkin[11]-0.3977475644174328*nuVtSqSum[2]*fEdge[11]-0.3827327723098713*nuVtSqSum[1]*fSkin[10]-0.3827327723098713*nuVtSqSum[1]*fEdge[10]+0.3977475644174328*nuVtSqSum[3]*fSkin[7]-0.3977475644174328*nuVtSqSum[3]*fEdge[7]+0.3977475644174328*nuVtSqSum[0]*fSkin[6]-0.3977475644174328*nuVtSqSum[0]*fEdge[6]+0.3977475644174328*nuVtSqSum[1]*fSkin[3]-0.3977475644174328*nuVtSqSum[1]*fEdge[3]; 
  Gdiff[6] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[15])-0.3827327723098713*nuVtSqSum[1]*fEdge[15]-0.3827327723098713*nuVtSqSum[0]*fSkin[14]-0.3827327723098713*nuVtSqSum[0]*fEdge[14]-0.3827327723098713*nuVtSqSum[3]*fSkin[13]-0.3827327723098713*nuVtSqSum[3]*fEdge[13]+0.3977475644174328*nuVtSqSum[1]*fSkin[11]-0.3977475644174328*nuVtSqSum[1]*fEdge[11]-0.3827327723098713*nuVtSqSum[2]*fSkin[10]-0.3827327723098713*nuVtSqSum[2]*fEdge[10]+0.3977475644174328*nuVtSqSum[0]*fSkin[7]-0.3977475644174328*nuVtSqSum[0]*fEdge[7]+0.3977475644174328*nuVtSqSum[3]*fSkin[6]-0.3977475644174328*nuVtSqSum[3]*fEdge[6]+0.3977475644174328*nuVtSqSum[2]*fSkin[3]-0.3977475644174328*nuVtSqSum[2]*fEdge[3]; 
  Gdiff[7] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[15])-0.3827327723098713*nuVtSqSum[0]*fEdge[15]-0.3827327723098713*nuVtSqSum[1]*fSkin[14]-0.3827327723098713*nuVtSqSum[1]*fEdge[14]-0.3827327723098713*nuVtSqSum[2]*fSkin[13]-0.3827327723098713*nuVtSqSum[2]*fEdge[13]+0.3977475644174328*nuVtSqSum[0]*fSkin[11]-0.3977475644174328*nuVtSqSum[0]*fEdge[11]-0.3827327723098713*nuVtSqSum[3]*fSkin[10]-0.3827327723098713*nuVtSqSum[3]*fEdge[10]+0.3977475644174328*nuVtSqSum[1]*fSkin[7]-0.3977475644174328*nuVtSqSum[1]*fEdge[7]+0.3977475644174328*nuVtSqSum[2]*fSkin[6]-0.3977475644174328*nuVtSqSum[2]*fEdge[6]+0.3977475644174328*fSkin[3]*nuVtSqSum[3]-0.3977475644174328*fEdge[3]*nuVtSqSum[3]; 

  Ghat[0] = Gdiff[0]*rdv2+0.3535533905932737*alphaDrSurf[4]*fUpwind[4]+0.3535533905932737*alphaDrSurf[2]*fUpwind[2]+0.3535533905932737*alphaDrSurf[1]*fUpwind[1]+0.3535533905932737*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.3535533905932737*alphaDrSurf[2]*fUpwind[4]+0.3535533905932737*fUpwind[2]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[1]+0.3535533905932737*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2+0.3535533905932737*alphaDrSurf[1]*fUpwind[4]+0.3535533905932737*fUpwind[1]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[0]*fUpwind[2]+0.3535533905932737*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = Gdiff[3]*rdv2+0.3535533905932737*alphaDrSurf[4]*fUpwind[7]+0.3535533905932737*alphaDrSurf[2]*fUpwind[6]+0.3535533905932737*alphaDrSurf[1]*fUpwind[5]+0.3535533905932737*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = Gdiff[4]*rdv2+0.3535533905932737*alphaDrSurf[0]*fUpwind[4]+0.3535533905932737*fUpwind[0]*alphaDrSurf[4]+0.3535533905932737*alphaDrSurf[1]*fUpwind[2]+0.3535533905932737*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[5] = Gdiff[5]*rdv2+0.3535533905932737*alphaDrSurf[2]*fUpwind[7]+0.3535533905932737*alphaDrSurf[4]*fUpwind[6]+0.3535533905932737*alphaDrSurf[0]*fUpwind[5]+0.3535533905932737*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[6] = Gdiff[6]*rdv2+0.3535533905932737*alphaDrSurf[1]*fUpwind[7]+0.3535533905932737*alphaDrSurf[0]*fUpwind[6]+0.3535533905932737*alphaDrSurf[4]*fUpwind[5]+0.3535533905932737*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[7] = Gdiff[7]*rdv2+0.3535533905932737*alphaDrSurf[0]*fUpwind[7]+0.3535533905932737*alphaDrSurf[1]*fUpwind[6]+0.3535533905932737*alphaDrSurf[2]*fUpwind[5]+0.3535533905932737*fUpwind[3]*alphaDrSurf[4]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[4] += 0.75*nuVtSqSum[3]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[1]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[0]*rdvSq4-1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[8] += 0.75*nuVtSqSum[2]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[1]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[9] += 0.75*nuVtSqSum[1]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[4]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[2]*rdvSq4-1.224744871391589*Gdiff2[2]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[10] += 0.75*nuVtSqSum[3]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[6]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[3]*rdvSq4-1.224744871391589*Gdiff2[3]*rdvSq4-1.224744871391589*Ghat[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[12] += 0.75*nuVtSqSum[0]*fSkin[12]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[9]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[8]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[5]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[4]*rdvSq4-1.224744871391589*Gdiff2[4]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[3]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat[4]*rdv2; 
  out[13] += 0.75*nuVtSqSum[2]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2[5]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[3]*rdvSq4-1.224744871391589*Ghat[5]*rdv2; 
  out[14] += 0.75*nuVtSqSum[1]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2[6]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[3]*rdvSq4-1.224744871391589*Ghat[6]*rdv2; 
  out[15] += 0.75*nuVtSqSum[0]*fSkin[15]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[14]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[11]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[7]*rdvSq4-1.224744871391589*Gdiff2[7]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[6]*rdvSq4+0.4330127018922193*fSkin[3]*nuVtSqSum[3]*rdvSq4-1.224744871391589*Ghat[7]*rdv2; 

  } 
} 
