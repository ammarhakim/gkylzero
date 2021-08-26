#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         Cell-center coordinates. 
  // dxv[5]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[16] = {0.0}; 
  double fUpwindQuad[16] = {0.0};
  double fUpwind[16] = {0.0};;
  double Ghat[16] = {0.0}; 
  double Gdiff[16] = {0.0}; 
  double Gdiff2[16] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = (4.0*w[2]+2.0*dxv[2])*nuSum-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = -2.0*sumNuUx[1]; 
  alphaDrSurf[2] = -2.0*sumNuUx[2]; 
  alphaDrSurf[5] = -2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(1, fSkin); 
  } else { 

    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(1, fSkin); 
  } else { 

    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(1, fSkin); 
  } else { 

    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(1, fSkin); 
  } else { 

    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(1, fSkin); 
  } else { 

    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(1, fSkin); 
  } else { 

    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(1, fSkin); 
  } else { 

    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(1, fSkin); 
  } else { 

    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(1, fSkin); 
  } else { 

    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(1, fSkin); 
  } else { 

    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(1, fSkin); 
  } else { 

    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(1, fSkin); 
  } else { 

    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(-1, fEdge); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(1, fSkin); 
  } else { 

    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(1, fSkin); 
  } else { 

    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(-1, fEdge); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(1, fSkin); 
  } else { 

    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(-1, fEdge); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(1, fSkin); 
  } else { 

    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(-1, fEdge); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Gdiff2[0] = 0.2041241452319315*nuVtSqSum[3]*fSkin[16]-0.2041241452319315*nuVtSqSum[3]*fEdge[16]+0.2041241452319315*nuVtSqSum[2]*fSkin[8]-0.2041241452319315*nuVtSqSum[2]*fEdge[8]+0.2041241452319315*nuVtSqSum[1]*fSkin[7]-0.2041241452319315*nuVtSqSum[1]*fEdge[7]+0.1767766952966368*nuVtSqSum[3]*fSkin[6]+0.1767766952966368*nuVtSqSum[3]*fEdge[6]+0.2041241452319315*nuVtSqSum[0]*fSkin[3]-0.2041241452319315*nuVtSqSum[0]*fEdge[3]+0.1767766952966368*fSkin[2]*nuVtSqSum[2]+0.1767766952966368*fEdge[2]*nuVtSqSum[2]+0.1767766952966368*fSkin[1]*nuVtSqSum[1]+0.1767766952966368*fEdge[1]*nuVtSqSum[1]+0.1767766952966368*fSkin[0]*nuVtSqSum[0]+0.1767766952966368*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.2041241452319315*nuVtSqSum[2]*fSkin[16]-0.2041241452319315*nuVtSqSum[2]*fEdge[16]+0.2041241452319315*nuVtSqSum[3]*fSkin[8]-0.2041241452319315*nuVtSqSum[3]*fEdge[8]+0.2041241452319315*nuVtSqSum[0]*fSkin[7]-0.2041241452319315*nuVtSqSum[0]*fEdge[7]+0.1767766952966368*nuVtSqSum[2]*fSkin[6]+0.1767766952966368*nuVtSqSum[2]*fEdge[6]+0.1767766952966368*fSkin[2]*nuVtSqSum[3]+0.1767766952966368*fEdge[2]*nuVtSqSum[3]+0.2041241452319315*nuVtSqSum[1]*fSkin[3]-0.2041241452319315*nuVtSqSum[1]*fEdge[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[1]+0.1767766952966368*fEdge[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fSkin[1]+0.1767766952966368*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.2041241452319315*nuVtSqSum[1]*fSkin[16]-0.2041241452319315*nuVtSqSum[1]*fEdge[16]+0.2041241452319315*nuVtSqSum[0]*fSkin[8]-0.2041241452319315*nuVtSqSum[0]*fEdge[8]+0.2041241452319315*nuVtSqSum[3]*fSkin[7]-0.2041241452319315*nuVtSqSum[3]*fEdge[7]+0.1767766952966368*nuVtSqSum[1]*fSkin[6]+0.1767766952966368*nuVtSqSum[1]*fEdge[6]+0.1767766952966368*fSkin[1]*nuVtSqSum[3]+0.1767766952966368*fEdge[1]*nuVtSqSum[3]+0.2041241452319315*nuVtSqSum[2]*fSkin[3]-0.2041241452319315*nuVtSqSum[2]*fEdge[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[2]+0.1767766952966368*fEdge[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fSkin[2]+0.1767766952966368*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = 0.2041241452319315*nuVtSqSum[3]*fSkin[26]-0.2041241452319315*nuVtSqSum[3]*fEdge[26]+0.2041241452319315*nuVtSqSum[2]*fSkin[19]-0.2041241452319315*nuVtSqSum[2]*fEdge[19]+0.2041241452319315*nuVtSqSum[1]*fSkin[18]-0.2041241452319315*nuVtSqSum[1]*fEdge[18]+0.1767766952966368*nuVtSqSum[3]*fSkin[17]+0.1767766952966368*nuVtSqSum[3]*fEdge[17]+0.2041241452319315*nuVtSqSum[0]*fSkin[11]-0.2041241452319315*nuVtSqSum[0]*fEdge[11]+0.1767766952966368*nuVtSqSum[2]*fSkin[10]+0.1767766952966368*nuVtSqSum[2]*fEdge[10]+0.1767766952966368*nuVtSqSum[1]*fSkin[9]+0.1767766952966368*nuVtSqSum[1]*fEdge[9]+0.1767766952966368*nuVtSqSum[0]*fSkin[4]+0.1767766952966368*nuVtSqSum[0]*fEdge[4]; 
  Gdiff2[4] = 0.2041241452319315*nuVtSqSum[3]*fSkin[27]-0.2041241452319315*nuVtSqSum[3]*fEdge[27]+0.2041241452319315*nuVtSqSum[2]*fSkin[22]-0.2041241452319315*nuVtSqSum[2]*fEdge[22]+0.2041241452319315*nuVtSqSum[1]*fSkin[21]-0.2041241452319315*nuVtSqSum[1]*fEdge[21]+0.1767766952966368*nuVtSqSum[3]*fSkin[20]+0.1767766952966368*nuVtSqSum[3]*fEdge[20]+0.2041241452319315*nuVtSqSum[0]*fSkin[14]-0.2041241452319315*nuVtSqSum[0]*fEdge[14]+0.1767766952966368*nuVtSqSum[2]*fSkin[13]+0.1767766952966368*nuVtSqSum[2]*fEdge[13]+0.1767766952966368*nuVtSqSum[1]*fSkin[12]+0.1767766952966368*nuVtSqSum[1]*fEdge[12]+0.1767766952966368*nuVtSqSum[0]*fSkin[5]+0.1767766952966368*nuVtSqSum[0]*fEdge[5]; 
  Gdiff2[5] = 0.2041241452319315*nuVtSqSum[0]*fSkin[16]-0.2041241452319315*nuVtSqSum[0]*fEdge[16]+0.2041241452319315*nuVtSqSum[1]*fSkin[8]-0.2041241452319315*nuVtSqSum[1]*fEdge[8]+0.2041241452319315*nuVtSqSum[2]*fSkin[7]-0.2041241452319315*nuVtSqSum[2]*fEdge[7]+0.1767766952966368*nuVtSqSum[0]*fSkin[6]+0.1767766952966368*nuVtSqSum[0]*fEdge[6]+0.2041241452319315*fSkin[3]*nuVtSqSum[3]-0.2041241452319315*fEdge[3]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[3]+0.1767766952966368*fEdge[0]*nuVtSqSum[3]+0.1767766952966368*fSkin[1]*nuVtSqSum[2]+0.1767766952966368*fEdge[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fSkin[2]+0.1767766952966368*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[6] = 0.2041241452319315*nuVtSqSum[2]*fSkin[26]-0.2041241452319315*nuVtSqSum[2]*fEdge[26]+0.2041241452319315*nuVtSqSum[3]*fSkin[19]-0.2041241452319315*nuVtSqSum[3]*fEdge[19]+0.2041241452319315*nuVtSqSum[0]*fSkin[18]-0.2041241452319315*nuVtSqSum[0]*fEdge[18]+0.1767766952966368*nuVtSqSum[2]*fSkin[17]+0.1767766952966368*nuVtSqSum[2]*fEdge[17]+0.2041241452319315*nuVtSqSum[1]*fSkin[11]-0.2041241452319315*nuVtSqSum[1]*fEdge[11]+0.1767766952966368*nuVtSqSum[3]*fSkin[10]+0.1767766952966368*nuVtSqSum[3]*fEdge[10]+0.1767766952966368*nuVtSqSum[0]*fSkin[9]+0.1767766952966368*nuVtSqSum[0]*fEdge[9]+0.1767766952966368*nuVtSqSum[1]*fSkin[4]+0.1767766952966368*nuVtSqSum[1]*fEdge[4]; 
  Gdiff2[7] = 0.2041241452319315*nuVtSqSum[1]*fSkin[26]-0.2041241452319315*nuVtSqSum[1]*fEdge[26]+0.2041241452319315*nuVtSqSum[0]*fSkin[19]-0.2041241452319315*nuVtSqSum[0]*fEdge[19]+0.2041241452319315*nuVtSqSum[3]*fSkin[18]-0.2041241452319315*nuVtSqSum[3]*fEdge[18]+0.1767766952966368*nuVtSqSum[1]*fSkin[17]+0.1767766952966368*nuVtSqSum[1]*fEdge[17]+0.2041241452319315*nuVtSqSum[2]*fSkin[11]-0.2041241452319315*nuVtSqSum[2]*fEdge[11]+0.1767766952966368*nuVtSqSum[0]*fSkin[10]+0.1767766952966368*nuVtSqSum[0]*fEdge[10]+0.1767766952966368*nuVtSqSum[3]*fSkin[9]+0.1767766952966368*nuVtSqSum[3]*fEdge[9]+0.1767766952966368*nuVtSqSum[2]*fSkin[4]+0.1767766952966368*nuVtSqSum[2]*fEdge[4]; 
  Gdiff2[8] = 0.2041241452319315*nuVtSqSum[2]*fSkin[27]-0.2041241452319315*nuVtSqSum[2]*fEdge[27]+0.2041241452319315*nuVtSqSum[3]*fSkin[22]-0.2041241452319315*nuVtSqSum[3]*fEdge[22]+0.2041241452319315*nuVtSqSum[0]*fSkin[21]-0.2041241452319315*nuVtSqSum[0]*fEdge[21]+0.1767766952966368*nuVtSqSum[2]*fSkin[20]+0.1767766952966368*nuVtSqSum[2]*fEdge[20]+0.2041241452319315*nuVtSqSum[1]*fSkin[14]-0.2041241452319315*nuVtSqSum[1]*fEdge[14]+0.1767766952966368*nuVtSqSum[3]*fSkin[13]+0.1767766952966368*nuVtSqSum[3]*fEdge[13]+0.1767766952966368*nuVtSqSum[0]*fSkin[12]+0.1767766952966368*nuVtSqSum[0]*fEdge[12]+0.1767766952966368*nuVtSqSum[1]*fSkin[5]+0.1767766952966368*nuVtSqSum[1]*fEdge[5]; 
  Gdiff2[9] = 0.2041241452319315*nuVtSqSum[1]*fSkin[27]-0.2041241452319315*nuVtSqSum[1]*fEdge[27]+0.2041241452319315*nuVtSqSum[0]*fSkin[22]-0.2041241452319315*nuVtSqSum[0]*fEdge[22]+0.2041241452319315*nuVtSqSum[3]*fSkin[21]-0.2041241452319315*nuVtSqSum[3]*fEdge[21]+0.1767766952966368*nuVtSqSum[1]*fSkin[20]+0.1767766952966368*nuVtSqSum[1]*fEdge[20]+0.2041241452319315*nuVtSqSum[2]*fSkin[14]-0.2041241452319315*nuVtSqSum[2]*fEdge[14]+0.1767766952966368*nuVtSqSum[0]*fSkin[13]+0.1767766952966368*nuVtSqSum[0]*fEdge[13]+0.1767766952966368*nuVtSqSum[3]*fSkin[12]+0.1767766952966368*nuVtSqSum[3]*fEdge[12]+0.1767766952966368*nuVtSqSum[2]*fSkin[5]+0.1767766952966368*nuVtSqSum[2]*fEdge[5]; 
  Gdiff2[10] = 0.2041241452319315*nuVtSqSum[3]*fSkin[31]-0.2041241452319315*nuVtSqSum[3]*fEdge[31]+0.2041241452319315*nuVtSqSum[2]*fSkin[30]-0.2041241452319315*nuVtSqSum[2]*fEdge[30]+0.2041241452319315*nuVtSqSum[1]*fSkin[29]-0.2041241452319315*nuVtSqSum[1]*fEdge[29]+0.1767766952966368*nuVtSqSum[3]*fSkin[28]+0.1767766952966368*nuVtSqSum[3]*fEdge[28]+0.2041241452319315*nuVtSqSum[0]*fSkin[25]-0.2041241452319315*nuVtSqSum[0]*fEdge[25]+0.1767766952966368*nuVtSqSum[2]*fSkin[24]+0.1767766952966368*nuVtSqSum[2]*fEdge[24]+0.1767766952966368*nuVtSqSum[1]*fSkin[23]+0.1767766952966368*nuVtSqSum[1]*fEdge[23]+0.1767766952966368*nuVtSqSum[0]*fSkin[15]+0.1767766952966368*nuVtSqSum[0]*fEdge[15]; 
  Gdiff2[11] = 0.2041241452319315*nuVtSqSum[0]*fSkin[26]-0.2041241452319315*nuVtSqSum[0]*fEdge[26]+0.2041241452319315*nuVtSqSum[1]*fSkin[19]-0.2041241452319315*nuVtSqSum[1]*fEdge[19]+0.2041241452319315*nuVtSqSum[2]*fSkin[18]-0.2041241452319315*nuVtSqSum[2]*fEdge[18]+0.1767766952966368*nuVtSqSum[0]*fSkin[17]+0.1767766952966368*nuVtSqSum[0]*fEdge[17]+0.2041241452319315*nuVtSqSum[3]*fSkin[11]-0.2041241452319315*nuVtSqSum[3]*fEdge[11]+0.1767766952966368*nuVtSqSum[1]*fSkin[10]+0.1767766952966368*nuVtSqSum[1]*fEdge[10]+0.1767766952966368*nuVtSqSum[2]*fSkin[9]+0.1767766952966368*nuVtSqSum[2]*fEdge[9]+0.1767766952966368*nuVtSqSum[3]*fSkin[4]+0.1767766952966368*nuVtSqSum[3]*fEdge[4]; 
  Gdiff2[12] = 0.2041241452319315*nuVtSqSum[0]*fSkin[27]-0.2041241452319315*nuVtSqSum[0]*fEdge[27]+0.2041241452319315*nuVtSqSum[1]*fSkin[22]-0.2041241452319315*nuVtSqSum[1]*fEdge[22]+0.2041241452319315*nuVtSqSum[2]*fSkin[21]-0.2041241452319315*nuVtSqSum[2]*fEdge[21]+0.1767766952966368*nuVtSqSum[0]*fSkin[20]+0.1767766952966368*nuVtSqSum[0]*fEdge[20]+0.2041241452319315*nuVtSqSum[3]*fSkin[14]-0.2041241452319315*nuVtSqSum[3]*fEdge[14]+0.1767766952966368*nuVtSqSum[1]*fSkin[13]+0.1767766952966368*nuVtSqSum[1]*fEdge[13]+0.1767766952966368*nuVtSqSum[2]*fSkin[12]+0.1767766952966368*nuVtSqSum[2]*fEdge[12]+0.1767766952966368*nuVtSqSum[3]*fSkin[5]+0.1767766952966368*nuVtSqSum[3]*fEdge[5]; 
  Gdiff2[13] = 0.2041241452319315*nuVtSqSum[2]*fSkin[31]-0.2041241452319315*nuVtSqSum[2]*fEdge[31]+0.2041241452319315*nuVtSqSum[3]*fSkin[30]-0.2041241452319315*nuVtSqSum[3]*fEdge[30]+0.2041241452319315*nuVtSqSum[0]*fSkin[29]-0.2041241452319315*nuVtSqSum[0]*fEdge[29]+0.1767766952966368*nuVtSqSum[2]*fSkin[28]+0.1767766952966368*nuVtSqSum[2]*fEdge[28]+0.2041241452319315*nuVtSqSum[1]*fSkin[25]-0.2041241452319315*nuVtSqSum[1]*fEdge[25]+0.1767766952966368*nuVtSqSum[3]*fSkin[24]+0.1767766952966368*nuVtSqSum[3]*fEdge[24]+0.1767766952966368*nuVtSqSum[0]*fSkin[23]+0.1767766952966368*nuVtSqSum[0]*fEdge[23]+0.1767766952966368*nuVtSqSum[1]*fSkin[15]+0.1767766952966368*nuVtSqSum[1]*fEdge[15]; 
  Gdiff2[14] = 0.2041241452319315*nuVtSqSum[1]*fSkin[31]-0.2041241452319315*nuVtSqSum[1]*fEdge[31]+0.2041241452319315*nuVtSqSum[0]*fSkin[30]-0.2041241452319315*nuVtSqSum[0]*fEdge[30]+0.2041241452319315*nuVtSqSum[3]*fSkin[29]-0.2041241452319315*nuVtSqSum[3]*fEdge[29]+0.1767766952966368*nuVtSqSum[1]*fSkin[28]+0.1767766952966368*nuVtSqSum[1]*fEdge[28]+0.2041241452319315*nuVtSqSum[2]*fSkin[25]-0.2041241452319315*nuVtSqSum[2]*fEdge[25]+0.1767766952966368*nuVtSqSum[0]*fSkin[24]+0.1767766952966368*nuVtSqSum[0]*fEdge[24]+0.1767766952966368*nuVtSqSum[3]*fSkin[23]+0.1767766952966368*nuVtSqSum[3]*fEdge[23]+0.1767766952966368*nuVtSqSum[2]*fSkin[15]+0.1767766952966368*nuVtSqSum[2]*fEdge[15]; 
  Gdiff2[15] = 0.2041241452319315*nuVtSqSum[0]*fSkin[31]-0.2041241452319315*nuVtSqSum[0]*fEdge[31]+0.2041241452319315*nuVtSqSum[1]*fSkin[30]-0.2041241452319315*nuVtSqSum[1]*fEdge[30]+0.2041241452319315*nuVtSqSum[2]*fSkin[29]-0.2041241452319315*nuVtSqSum[2]*fEdge[29]+0.1767766952966368*nuVtSqSum[0]*fSkin[28]+0.1767766952966368*nuVtSqSum[0]*fEdge[28]+0.2041241452319315*nuVtSqSum[3]*fSkin[25]-0.2041241452319315*nuVtSqSum[3]*fEdge[25]+0.1767766952966368*nuVtSqSum[1]*fSkin[24]+0.1767766952966368*nuVtSqSum[1]*fEdge[24]+0.1767766952966368*nuVtSqSum[2]*fSkin[23]+0.1767766952966368*nuVtSqSum[2]*fEdge[23]+0.1767766952966368*nuVtSqSum[3]*fSkin[15]+0.1767766952966368*nuVtSqSum[3]*fEdge[15]; 

  Gdiff[0] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[16])-0.3827327723098713*nuVtSqSum[3]*fEdge[16]-0.3827327723098713*nuVtSqSum[2]*fSkin[8]-0.3827327723098713*nuVtSqSum[2]*fEdge[8]-0.3827327723098713*nuVtSqSum[1]*fSkin[7]-0.3827327723098713*nuVtSqSum[1]*fEdge[7]-0.3977475644174328*nuVtSqSum[3]*fSkin[6]+0.3977475644174328*nuVtSqSum[3]*fEdge[6]-0.3827327723098713*nuVtSqSum[0]*fSkin[3]-0.3827327723098713*nuVtSqSum[0]*fEdge[3]-0.3977475644174328*fSkin[2]*nuVtSqSum[2]+0.3977475644174328*fEdge[2]*nuVtSqSum[2]-0.3977475644174328*fSkin[1]*nuVtSqSum[1]+0.3977475644174328*fEdge[1]*nuVtSqSum[1]-0.3977475644174328*fSkin[0]*nuVtSqSum[0]+0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[16])-0.3827327723098713*nuVtSqSum[2]*fEdge[16]-0.3827327723098713*nuVtSqSum[3]*fSkin[8]-0.3827327723098713*nuVtSqSum[3]*fEdge[8]-0.3827327723098713*nuVtSqSum[0]*fSkin[7]-0.3827327723098713*nuVtSqSum[0]*fEdge[7]-0.3977475644174328*nuVtSqSum[2]*fSkin[6]+0.3977475644174328*nuVtSqSum[2]*fEdge[6]-0.3977475644174328*fSkin[2]*nuVtSqSum[3]+0.3977475644174328*fEdge[2]*nuVtSqSum[3]-0.3827327723098713*nuVtSqSum[1]*fSkin[3]-0.3827327723098713*nuVtSqSum[1]*fEdge[3]-0.3977475644174328*fSkin[0]*nuVtSqSum[1]+0.3977475644174328*fEdge[0]*nuVtSqSum[1]-0.3977475644174328*nuVtSqSum[0]*fSkin[1]+0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[16])-0.3827327723098713*nuVtSqSum[1]*fEdge[16]-0.3827327723098713*nuVtSqSum[0]*fSkin[8]-0.3827327723098713*nuVtSqSum[0]*fEdge[8]-0.3827327723098713*nuVtSqSum[3]*fSkin[7]-0.3827327723098713*nuVtSqSum[3]*fEdge[7]-0.3977475644174328*nuVtSqSum[1]*fSkin[6]+0.3977475644174328*nuVtSqSum[1]*fEdge[6]-0.3977475644174328*fSkin[1]*nuVtSqSum[3]+0.3977475644174328*fEdge[1]*nuVtSqSum[3]-0.3827327723098713*nuVtSqSum[2]*fSkin[3]-0.3827327723098713*nuVtSqSum[2]*fEdge[3]-0.3977475644174328*fSkin[0]*nuVtSqSum[2]+0.3977475644174328*fEdge[0]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[0]*fSkin[2]+0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[26])-0.3827327723098713*nuVtSqSum[3]*fEdge[26]-0.3827327723098713*nuVtSqSum[2]*fSkin[19]-0.3827327723098713*nuVtSqSum[2]*fEdge[19]-0.3827327723098713*nuVtSqSum[1]*fSkin[18]-0.3827327723098713*nuVtSqSum[1]*fEdge[18]-0.3977475644174328*nuVtSqSum[3]*fSkin[17]+0.3977475644174328*nuVtSqSum[3]*fEdge[17]-0.3827327723098713*nuVtSqSum[0]*fSkin[11]-0.3827327723098713*nuVtSqSum[0]*fEdge[11]-0.3977475644174328*nuVtSqSum[2]*fSkin[10]+0.3977475644174328*nuVtSqSum[2]*fEdge[10]-0.3977475644174328*nuVtSqSum[1]*fSkin[9]+0.3977475644174328*nuVtSqSum[1]*fEdge[9]-0.3977475644174328*nuVtSqSum[0]*fSkin[4]+0.3977475644174328*nuVtSqSum[0]*fEdge[4]; 
  Gdiff[4] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[27])-0.3827327723098713*nuVtSqSum[3]*fEdge[27]-0.3827327723098713*nuVtSqSum[2]*fSkin[22]-0.3827327723098713*nuVtSqSum[2]*fEdge[22]-0.3827327723098713*nuVtSqSum[1]*fSkin[21]-0.3827327723098713*nuVtSqSum[1]*fEdge[21]-0.3977475644174328*nuVtSqSum[3]*fSkin[20]+0.3977475644174328*nuVtSqSum[3]*fEdge[20]-0.3827327723098713*nuVtSqSum[0]*fSkin[14]-0.3827327723098713*nuVtSqSum[0]*fEdge[14]-0.3977475644174328*nuVtSqSum[2]*fSkin[13]+0.3977475644174328*nuVtSqSum[2]*fEdge[13]-0.3977475644174328*nuVtSqSum[1]*fSkin[12]+0.3977475644174328*nuVtSqSum[1]*fEdge[12]-0.3977475644174328*nuVtSqSum[0]*fSkin[5]+0.3977475644174328*nuVtSqSum[0]*fEdge[5]; 
  Gdiff[5] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[16])-0.3827327723098713*nuVtSqSum[0]*fEdge[16]-0.3827327723098713*nuVtSqSum[1]*fSkin[8]-0.3827327723098713*nuVtSqSum[1]*fEdge[8]-0.3827327723098713*nuVtSqSum[2]*fSkin[7]-0.3827327723098713*nuVtSqSum[2]*fEdge[7]-0.3977475644174328*nuVtSqSum[0]*fSkin[6]+0.3977475644174328*nuVtSqSum[0]*fEdge[6]-0.3827327723098713*fSkin[3]*nuVtSqSum[3]-0.3827327723098713*fEdge[3]*nuVtSqSum[3]-0.3977475644174328*fSkin[0]*nuVtSqSum[3]+0.3977475644174328*fEdge[0]*nuVtSqSum[3]-0.3977475644174328*fSkin[1]*nuVtSqSum[2]+0.3977475644174328*fEdge[1]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[1]*fSkin[2]+0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[6] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[26])-0.3827327723098713*nuVtSqSum[2]*fEdge[26]-0.3827327723098713*nuVtSqSum[3]*fSkin[19]-0.3827327723098713*nuVtSqSum[3]*fEdge[19]-0.3827327723098713*nuVtSqSum[0]*fSkin[18]-0.3827327723098713*nuVtSqSum[0]*fEdge[18]-0.3977475644174328*nuVtSqSum[2]*fSkin[17]+0.3977475644174328*nuVtSqSum[2]*fEdge[17]-0.3827327723098713*nuVtSqSum[1]*fSkin[11]-0.3827327723098713*nuVtSqSum[1]*fEdge[11]-0.3977475644174328*nuVtSqSum[3]*fSkin[10]+0.3977475644174328*nuVtSqSum[3]*fEdge[10]-0.3977475644174328*nuVtSqSum[0]*fSkin[9]+0.3977475644174328*nuVtSqSum[0]*fEdge[9]-0.3977475644174328*nuVtSqSum[1]*fSkin[4]+0.3977475644174328*nuVtSqSum[1]*fEdge[4]; 
  Gdiff[7] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[26])-0.3827327723098713*nuVtSqSum[1]*fEdge[26]-0.3827327723098713*nuVtSqSum[0]*fSkin[19]-0.3827327723098713*nuVtSqSum[0]*fEdge[19]-0.3827327723098713*nuVtSqSum[3]*fSkin[18]-0.3827327723098713*nuVtSqSum[3]*fEdge[18]-0.3977475644174328*nuVtSqSum[1]*fSkin[17]+0.3977475644174328*nuVtSqSum[1]*fEdge[17]-0.3827327723098713*nuVtSqSum[2]*fSkin[11]-0.3827327723098713*nuVtSqSum[2]*fEdge[11]-0.3977475644174328*nuVtSqSum[0]*fSkin[10]+0.3977475644174328*nuVtSqSum[0]*fEdge[10]-0.3977475644174328*nuVtSqSum[3]*fSkin[9]+0.3977475644174328*nuVtSqSum[3]*fEdge[9]-0.3977475644174328*nuVtSqSum[2]*fSkin[4]+0.3977475644174328*nuVtSqSum[2]*fEdge[4]; 
  Gdiff[8] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[27])-0.3827327723098713*nuVtSqSum[2]*fEdge[27]-0.3827327723098713*nuVtSqSum[3]*fSkin[22]-0.3827327723098713*nuVtSqSum[3]*fEdge[22]-0.3827327723098713*nuVtSqSum[0]*fSkin[21]-0.3827327723098713*nuVtSqSum[0]*fEdge[21]-0.3977475644174328*nuVtSqSum[2]*fSkin[20]+0.3977475644174328*nuVtSqSum[2]*fEdge[20]-0.3827327723098713*nuVtSqSum[1]*fSkin[14]-0.3827327723098713*nuVtSqSum[1]*fEdge[14]-0.3977475644174328*nuVtSqSum[3]*fSkin[13]+0.3977475644174328*nuVtSqSum[3]*fEdge[13]-0.3977475644174328*nuVtSqSum[0]*fSkin[12]+0.3977475644174328*nuVtSqSum[0]*fEdge[12]-0.3977475644174328*nuVtSqSum[1]*fSkin[5]+0.3977475644174328*nuVtSqSum[1]*fEdge[5]; 
  Gdiff[9] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[27])-0.3827327723098713*nuVtSqSum[1]*fEdge[27]-0.3827327723098713*nuVtSqSum[0]*fSkin[22]-0.3827327723098713*nuVtSqSum[0]*fEdge[22]-0.3827327723098713*nuVtSqSum[3]*fSkin[21]-0.3827327723098713*nuVtSqSum[3]*fEdge[21]-0.3977475644174328*nuVtSqSum[1]*fSkin[20]+0.3977475644174328*nuVtSqSum[1]*fEdge[20]-0.3827327723098713*nuVtSqSum[2]*fSkin[14]-0.3827327723098713*nuVtSqSum[2]*fEdge[14]-0.3977475644174328*nuVtSqSum[0]*fSkin[13]+0.3977475644174328*nuVtSqSum[0]*fEdge[13]-0.3977475644174328*nuVtSqSum[3]*fSkin[12]+0.3977475644174328*nuVtSqSum[3]*fEdge[12]-0.3977475644174328*nuVtSqSum[2]*fSkin[5]+0.3977475644174328*nuVtSqSum[2]*fEdge[5]; 
  Gdiff[10] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[31])-0.3827327723098713*nuVtSqSum[3]*fEdge[31]-0.3827327723098713*nuVtSqSum[2]*fSkin[30]-0.3827327723098713*nuVtSqSum[2]*fEdge[30]-0.3827327723098713*nuVtSqSum[1]*fSkin[29]-0.3827327723098713*nuVtSqSum[1]*fEdge[29]-0.3977475644174328*nuVtSqSum[3]*fSkin[28]+0.3977475644174328*nuVtSqSum[3]*fEdge[28]-0.3827327723098713*nuVtSqSum[0]*fSkin[25]-0.3827327723098713*nuVtSqSum[0]*fEdge[25]-0.3977475644174328*nuVtSqSum[2]*fSkin[24]+0.3977475644174328*nuVtSqSum[2]*fEdge[24]-0.3977475644174328*nuVtSqSum[1]*fSkin[23]+0.3977475644174328*nuVtSqSum[1]*fEdge[23]-0.3977475644174328*nuVtSqSum[0]*fSkin[15]+0.3977475644174328*nuVtSqSum[0]*fEdge[15]; 
  Gdiff[11] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[26])-0.3827327723098713*nuVtSqSum[0]*fEdge[26]-0.3827327723098713*nuVtSqSum[1]*fSkin[19]-0.3827327723098713*nuVtSqSum[1]*fEdge[19]-0.3827327723098713*nuVtSqSum[2]*fSkin[18]-0.3827327723098713*nuVtSqSum[2]*fEdge[18]-0.3977475644174328*nuVtSqSum[0]*fSkin[17]+0.3977475644174328*nuVtSqSum[0]*fEdge[17]-0.3827327723098713*nuVtSqSum[3]*fSkin[11]-0.3827327723098713*nuVtSqSum[3]*fEdge[11]-0.3977475644174328*nuVtSqSum[1]*fSkin[10]+0.3977475644174328*nuVtSqSum[1]*fEdge[10]-0.3977475644174328*nuVtSqSum[2]*fSkin[9]+0.3977475644174328*nuVtSqSum[2]*fEdge[9]-0.3977475644174328*nuVtSqSum[3]*fSkin[4]+0.3977475644174328*nuVtSqSum[3]*fEdge[4]; 
  Gdiff[12] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[27])-0.3827327723098713*nuVtSqSum[0]*fEdge[27]-0.3827327723098713*nuVtSqSum[1]*fSkin[22]-0.3827327723098713*nuVtSqSum[1]*fEdge[22]-0.3827327723098713*nuVtSqSum[2]*fSkin[21]-0.3827327723098713*nuVtSqSum[2]*fEdge[21]-0.3977475644174328*nuVtSqSum[0]*fSkin[20]+0.3977475644174328*nuVtSqSum[0]*fEdge[20]-0.3827327723098713*nuVtSqSum[3]*fSkin[14]-0.3827327723098713*nuVtSqSum[3]*fEdge[14]-0.3977475644174328*nuVtSqSum[1]*fSkin[13]+0.3977475644174328*nuVtSqSum[1]*fEdge[13]-0.3977475644174328*nuVtSqSum[2]*fSkin[12]+0.3977475644174328*nuVtSqSum[2]*fEdge[12]-0.3977475644174328*nuVtSqSum[3]*fSkin[5]+0.3977475644174328*nuVtSqSum[3]*fEdge[5]; 
  Gdiff[13] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[31])-0.3827327723098713*nuVtSqSum[2]*fEdge[31]-0.3827327723098713*nuVtSqSum[3]*fSkin[30]-0.3827327723098713*nuVtSqSum[3]*fEdge[30]-0.3827327723098713*nuVtSqSum[0]*fSkin[29]-0.3827327723098713*nuVtSqSum[0]*fEdge[29]-0.3977475644174328*nuVtSqSum[2]*fSkin[28]+0.3977475644174328*nuVtSqSum[2]*fEdge[28]-0.3827327723098713*nuVtSqSum[1]*fSkin[25]-0.3827327723098713*nuVtSqSum[1]*fEdge[25]-0.3977475644174328*nuVtSqSum[3]*fSkin[24]+0.3977475644174328*nuVtSqSum[3]*fEdge[24]-0.3977475644174328*nuVtSqSum[0]*fSkin[23]+0.3977475644174328*nuVtSqSum[0]*fEdge[23]-0.3977475644174328*nuVtSqSum[1]*fSkin[15]+0.3977475644174328*nuVtSqSum[1]*fEdge[15]; 
  Gdiff[14] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[31])-0.3827327723098713*nuVtSqSum[1]*fEdge[31]-0.3827327723098713*nuVtSqSum[0]*fSkin[30]-0.3827327723098713*nuVtSqSum[0]*fEdge[30]-0.3827327723098713*nuVtSqSum[3]*fSkin[29]-0.3827327723098713*nuVtSqSum[3]*fEdge[29]-0.3977475644174328*nuVtSqSum[1]*fSkin[28]+0.3977475644174328*nuVtSqSum[1]*fEdge[28]-0.3827327723098713*nuVtSqSum[2]*fSkin[25]-0.3827327723098713*nuVtSqSum[2]*fEdge[25]-0.3977475644174328*nuVtSqSum[0]*fSkin[24]+0.3977475644174328*nuVtSqSum[0]*fEdge[24]-0.3977475644174328*nuVtSqSum[3]*fSkin[23]+0.3977475644174328*nuVtSqSum[3]*fEdge[23]-0.3977475644174328*nuVtSqSum[2]*fSkin[15]+0.3977475644174328*nuVtSqSum[2]*fEdge[15]; 
  Gdiff[15] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[31])-0.3827327723098713*nuVtSqSum[0]*fEdge[31]-0.3827327723098713*nuVtSqSum[1]*fSkin[30]-0.3827327723098713*nuVtSqSum[1]*fEdge[30]-0.3827327723098713*nuVtSqSum[2]*fSkin[29]-0.3827327723098713*nuVtSqSum[2]*fEdge[29]-0.3977475644174328*nuVtSqSum[0]*fSkin[28]+0.3977475644174328*nuVtSqSum[0]*fEdge[28]-0.3827327723098713*nuVtSqSum[3]*fSkin[25]-0.3827327723098713*nuVtSqSum[3]*fEdge[25]-0.3977475644174328*nuVtSqSum[1]*fSkin[24]+0.3977475644174328*nuVtSqSum[1]*fEdge[24]-0.3977475644174328*nuVtSqSum[2]*fSkin[23]+0.3977475644174328*nuVtSqSum[2]*fEdge[23]-0.3977475644174328*nuVtSqSum[3]*fSkin[15]+0.3977475644174328*nuVtSqSum[3]*fEdge[15]; 

  Ghat[0] = (-1.0*Gdiff[0]*rdv2)+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = (-1.0*Gdiff[1]*rdv2)+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = (-1.0*Gdiff[2]*rdv2)+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = (-1.0*Gdiff[3]*rdv2)+0.25*alphaDrSurf[5]*fUpwind[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = (-1.0*Gdiff[4]*rdv2)+0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[5] = (-1.0*Gdiff[5]*rdv2)+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[6] = (-1.0*Gdiff[6]*rdv2)+0.25*alphaDrSurf[2]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = (-1.0*Gdiff[7]*rdv2)+0.25*alphaDrSurf[1]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[8] = (-1.0*Gdiff[8]*rdv2)+0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  Ghat[9] = (-1.0*Gdiff[9]*rdv2)+0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  Ghat[10] = (-1.0*Gdiff[10]*rdv2)+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  Ghat[11] = (-1.0*Gdiff[11]*rdv2)+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  Ghat[12] = (-1.0*Gdiff[12]*rdv2)+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  Ghat[13] = (-1.0*Gdiff[13]*rdv2)+0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  Ghat[14] = (-1.0*Gdiff[14]*rdv2)+0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  Ghat[15] = (-1.0*Gdiff[15]*rdv2)+0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[10]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.75*nuVtSqSum[3]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[6]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[3]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[1]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[0]*rdvSq4+1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += -0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += -0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.75*nuVtSqSum[2]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[6]*rdvSq4-0.4330127018922193*fSkin[2]*nuVtSqSum[3]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[1]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[8] += 0.75*nuVtSqSum[1]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[6]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[3]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[2]*rdvSq4+1.224744871391589*Gdiff2[2]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[9] += -0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += -0.7071067811865475*Ghat[7]*rdv2; 
  out[11] += 0.75*nuVtSqSum[3]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[18]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[11]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[9]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[4]*rdvSq4+1.224744871391589*Gdiff2[3]*rdvSq4-1.224744871391589*Ghat[3]*rdv2; 
  out[12] += -0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += -0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 0.75*nuVtSqSum[3]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[21]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[14]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[12]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[5]*rdvSq4+1.224744871391589*Gdiff2[4]*rdvSq4-1.224744871391589*Ghat[4]*rdv2; 
  out[15] += -0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 0.75*nuVtSqSum[0]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[7]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[6]*rdvSq4+1.224744871391589*Gdiff2[5]*rdvSq4+0.75*fSkin[3]*nuVtSqSum[3]*rdvSq4-0.4330127018922193*fSkin[0]*nuVtSqSum[3]*rdvSq4-0.4330127018922193*fSkin[1]*nuVtSqSum[2]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat[5]*rdv2; 
  out[17] += -0.7071067811865475*Ghat[11]*rdv2; 
  out[18] += 0.75*nuVtSqSum[2]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[18]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[11]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[9]*rdvSq4+1.224744871391589*Gdiff2[6]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 0.75*nuVtSqSum[1]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[18]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[11]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[9]*rdvSq4+1.224744871391589*Gdiff2[7]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[7]*rdv2; 
  out[20] += -0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += 0.75*nuVtSqSum[2]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[21]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[14]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[12]*rdvSq4+1.224744871391589*Gdiff2[8]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[8]*rdv2; 
  out[22] += 0.75*nuVtSqSum[1]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[21]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[14]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[12]*rdvSq4+1.224744871391589*Gdiff2[9]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[9]*rdv2; 
  out[23] += -0.7071067811865475*Ghat[13]*rdv2; 
  out[24] += -0.7071067811865475*Ghat[14]*rdv2; 
  out[25] += 0.75*nuVtSqSum[3]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[29]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[25]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[24]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[23]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[15]*rdvSq4+1.224744871391589*Gdiff2[10]*rdvSq4-1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 0.75*nuVtSqSum[0]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[18]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[11]*rdvSq4+1.224744871391589*Gdiff2[11]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[10]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[9]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[11]*rdv2; 
  out[27] += 0.75*nuVtSqSum[0]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[21]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[14]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[13]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[12]*rdvSq4+1.224744871391589*Gdiff2[12]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[12]*rdv2; 
  out[28] += -0.7071067811865475*Ghat[15]*rdv2; 
  out[29] += 0.75*nuVtSqSum[2]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[29]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[25]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[24]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[23]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[15]*rdvSq4+1.224744871391589*Gdiff2[13]*rdvSq4-1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 0.75*nuVtSqSum[1]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[29]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[25]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[24]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[23]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[15]*rdvSq4+1.224744871391589*Gdiff2[14]*rdvSq4-1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 0.75*nuVtSqSum[0]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[29]*rdvSq4-0.4330127018922193*nuVtSqSum[0]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[25]*rdvSq4-0.4330127018922193*nuVtSqSum[1]*fSkin[24]*rdvSq4-0.4330127018922193*nuVtSqSum[2]*fSkin[23]*rdvSq4-0.4330127018922193*nuVtSqSum[3]*fSkin[15]*rdvSq4+1.224744871391589*Gdiff2[15]*rdvSq4-1.224744871391589*Ghat[15]*rdv2; 

  } else { 

  alphaDrSurf[0] = (4.0*w[2]+2.0*dxv[2])*nuSum-2.0*sumNuUx[0]; 
  alphaDrSurf[1] = -2.0*sumNuUx[1]; 
  alphaDrSurf[2] = -2.0*sumNuUx[2]; 
  alphaDrSurf[5] = -2.0*sumNuUx[3]; 

  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_2x3v_p1_surfvx_quad_0(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_2x3v_p1_surfvx_quad_1(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_2x3v_p1_surfvx_quad_2(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(1, fEdge); 
  } else { 
    fUpwindQuad[3] = ser_2x3v_p1_surfvx_quad_3(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(1, fEdge); 
  } else { 
    fUpwindQuad[4] = ser_2x3v_p1_surfvx_quad_4(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(1, fEdge); 
  } else { 
    fUpwindQuad[5] = ser_2x3v_p1_surfvx_quad_5(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(1, fEdge); 
  } else { 
    fUpwindQuad[6] = ser_2x3v_p1_surfvx_quad_6(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(1, fEdge); 
  } else { 
    fUpwindQuad[7] = ser_2x3v_p1_surfvx_quad_7(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(1, fEdge); 
  } else { 
    fUpwindQuad[8] = ser_2x3v_p1_surfvx_quad_8(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(1, fEdge); 
  } else { 
    fUpwindQuad[9] = ser_2x3v_p1_surfvx_quad_9(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(1, fEdge); 
  } else { 
    fUpwindQuad[10] = ser_2x3v_p1_surfvx_quad_10(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(1, fEdge); 
  } else { 
    fUpwindQuad[11] = ser_2x3v_p1_surfvx_quad_11(-1, fSkin); 
  } 
  if (alphaDrSurf[5]-alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(1, fEdge); 
  } else { 
    fUpwindQuad[12] = ser_2x3v_p1_surfvx_quad_12(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])-alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(1, fEdge); 
  } else { 
    fUpwindQuad[13] = ser_2x3v_p1_surfvx_quad_13(-1, fSkin); 
  } 
  if ((-alphaDrSurf[5])+alphaDrSurf[2]-alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(1, fEdge); 
  } else { 
    fUpwindQuad[14] = ser_2x3v_p1_surfvx_quad_14(-1, fSkin); 
  } 
  if (alphaDrSurf[5]+alphaDrSurf[2]+alphaDrSurf[1]+alphaDrSurf[0] > 0) { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(1, fEdge); 
  } else { 
    fUpwindQuad[15] = ser_2x3v_p1_surfvx_quad_15(-1, fSkin); 
  } 

  fUpwind[0] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[1] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[2] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[3] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8])+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[4] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[5] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 
  fUpwind[6] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]-1.0*fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[7] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]-1.0*(fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[8] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*fUpwindQuad[12]+fUpwindQuad[11]-1.0*fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]-1.0*fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[9] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12])+fUpwindQuad[11]+fUpwindQuad[10]-1.0*(fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]-1.0*(fUpwindQuad[3]+fUpwindQuad[2])+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[10] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]+fUpwindQuad[13]+fUpwindQuad[12]-1.0*(fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]+fUpwindQuad[8]+fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]+fUpwindQuad[4])+fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]+fUpwindQuad[0]); 
  fUpwind[11] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*fUpwindQuad[8]+fUpwindQuad[7]-1.0*(fUpwindQuad[6]+fUpwindQuad[5])+fUpwindQuad[4]-1.0*fUpwindQuad[3]+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[12] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]+fUpwindQuad[11]-1.0*(fUpwindQuad[10]+fUpwindQuad[9])+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]+fUpwindQuad[5]-1.0*(fUpwindQuad[4]+fUpwindQuad[3])+fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[13] = 0.25*(fUpwindQuad[15]-1.0*fUpwindQuad[14]+fUpwindQuad[13]-1.0*(fUpwindQuad[12]+fUpwindQuad[11])+fUpwindQuad[10]-1.0*fUpwindQuad[9]+fUpwindQuad[8]-1.0*fUpwindQuad[7]+fUpwindQuad[6]-1.0*fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]-1.0*fUpwindQuad[2]+fUpwindQuad[1]-1.0*fUpwindQuad[0]); 
  fUpwind[14] = 0.25*(fUpwindQuad[15]+fUpwindQuad[14]-1.0*(fUpwindQuad[13]+fUpwindQuad[12]+fUpwindQuad[11]+fUpwindQuad[10])+fUpwindQuad[9]+fUpwindQuad[8]-1.0*(fUpwindQuad[7]+fUpwindQuad[6])+fUpwindQuad[5]+fUpwindQuad[4]+fUpwindQuad[3]+fUpwindQuad[2]-1.0*(fUpwindQuad[1]+fUpwindQuad[0])); 
  fUpwind[15] = 0.25*(fUpwindQuad[15]-1.0*(fUpwindQuad[14]+fUpwindQuad[13])+fUpwindQuad[12]-1.0*fUpwindQuad[11]+fUpwindQuad[10]+fUpwindQuad[9]-1.0*(fUpwindQuad[8]+fUpwindQuad[7])+fUpwindQuad[6]+fUpwindQuad[5]-1.0*fUpwindQuad[4]+fUpwindQuad[3]-1.0*(fUpwindQuad[2]+fUpwindQuad[1])+fUpwindQuad[0]); 

  Gdiff2[0] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[16])+0.2041241452319315*nuVtSqSum[3]*fEdge[16]-0.2041241452319315*nuVtSqSum[2]*fSkin[8]+0.2041241452319315*nuVtSqSum[2]*fEdge[8]-0.2041241452319315*nuVtSqSum[1]*fSkin[7]+0.2041241452319315*nuVtSqSum[1]*fEdge[7]+0.1767766952966368*nuVtSqSum[3]*fSkin[6]+0.1767766952966368*nuVtSqSum[3]*fEdge[6]-0.2041241452319315*nuVtSqSum[0]*fSkin[3]+0.2041241452319315*nuVtSqSum[0]*fEdge[3]+0.1767766952966368*fSkin[2]*nuVtSqSum[2]+0.1767766952966368*fEdge[2]*nuVtSqSum[2]+0.1767766952966368*fSkin[1]*nuVtSqSum[1]+0.1767766952966368*fEdge[1]*nuVtSqSum[1]+0.1767766952966368*fSkin[0]*nuVtSqSum[0]+0.1767766952966368*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[16])+0.2041241452319315*nuVtSqSum[2]*fEdge[16]-0.2041241452319315*nuVtSqSum[3]*fSkin[8]+0.2041241452319315*nuVtSqSum[3]*fEdge[8]-0.2041241452319315*nuVtSqSum[0]*fSkin[7]+0.2041241452319315*nuVtSqSum[0]*fEdge[7]+0.1767766952966368*nuVtSqSum[2]*fSkin[6]+0.1767766952966368*nuVtSqSum[2]*fEdge[6]+0.1767766952966368*fSkin[2]*nuVtSqSum[3]+0.1767766952966368*fEdge[2]*nuVtSqSum[3]-0.2041241452319315*nuVtSqSum[1]*fSkin[3]+0.2041241452319315*nuVtSqSum[1]*fEdge[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[1]+0.1767766952966368*fEdge[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fSkin[1]+0.1767766952966368*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[16])+0.2041241452319315*nuVtSqSum[1]*fEdge[16]-0.2041241452319315*nuVtSqSum[0]*fSkin[8]+0.2041241452319315*nuVtSqSum[0]*fEdge[8]-0.2041241452319315*nuVtSqSum[3]*fSkin[7]+0.2041241452319315*nuVtSqSum[3]*fEdge[7]+0.1767766952966368*nuVtSqSum[1]*fSkin[6]+0.1767766952966368*nuVtSqSum[1]*fEdge[6]+0.1767766952966368*fSkin[1]*nuVtSqSum[3]+0.1767766952966368*fEdge[1]*nuVtSqSum[3]-0.2041241452319315*nuVtSqSum[2]*fSkin[3]+0.2041241452319315*nuVtSqSum[2]*fEdge[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[2]+0.1767766952966368*fEdge[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fSkin[2]+0.1767766952966368*nuVtSqSum[0]*fEdge[2]; 
  Gdiff2[3] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[26])+0.2041241452319315*nuVtSqSum[3]*fEdge[26]-0.2041241452319315*nuVtSqSum[2]*fSkin[19]+0.2041241452319315*nuVtSqSum[2]*fEdge[19]-0.2041241452319315*nuVtSqSum[1]*fSkin[18]+0.2041241452319315*nuVtSqSum[1]*fEdge[18]+0.1767766952966368*nuVtSqSum[3]*fSkin[17]+0.1767766952966368*nuVtSqSum[3]*fEdge[17]-0.2041241452319315*nuVtSqSum[0]*fSkin[11]+0.2041241452319315*nuVtSqSum[0]*fEdge[11]+0.1767766952966368*nuVtSqSum[2]*fSkin[10]+0.1767766952966368*nuVtSqSum[2]*fEdge[10]+0.1767766952966368*nuVtSqSum[1]*fSkin[9]+0.1767766952966368*nuVtSqSum[1]*fEdge[9]+0.1767766952966368*nuVtSqSum[0]*fSkin[4]+0.1767766952966368*nuVtSqSum[0]*fEdge[4]; 
  Gdiff2[4] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[27])+0.2041241452319315*nuVtSqSum[3]*fEdge[27]-0.2041241452319315*nuVtSqSum[2]*fSkin[22]+0.2041241452319315*nuVtSqSum[2]*fEdge[22]-0.2041241452319315*nuVtSqSum[1]*fSkin[21]+0.2041241452319315*nuVtSqSum[1]*fEdge[21]+0.1767766952966368*nuVtSqSum[3]*fSkin[20]+0.1767766952966368*nuVtSqSum[3]*fEdge[20]-0.2041241452319315*nuVtSqSum[0]*fSkin[14]+0.2041241452319315*nuVtSqSum[0]*fEdge[14]+0.1767766952966368*nuVtSqSum[2]*fSkin[13]+0.1767766952966368*nuVtSqSum[2]*fEdge[13]+0.1767766952966368*nuVtSqSum[1]*fSkin[12]+0.1767766952966368*nuVtSqSum[1]*fEdge[12]+0.1767766952966368*nuVtSqSum[0]*fSkin[5]+0.1767766952966368*nuVtSqSum[0]*fEdge[5]; 
  Gdiff2[5] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[16])+0.2041241452319315*nuVtSqSum[0]*fEdge[16]-0.2041241452319315*nuVtSqSum[1]*fSkin[8]+0.2041241452319315*nuVtSqSum[1]*fEdge[8]-0.2041241452319315*nuVtSqSum[2]*fSkin[7]+0.2041241452319315*nuVtSqSum[2]*fEdge[7]+0.1767766952966368*nuVtSqSum[0]*fSkin[6]+0.1767766952966368*nuVtSqSum[0]*fEdge[6]-0.2041241452319315*fSkin[3]*nuVtSqSum[3]+0.2041241452319315*fEdge[3]*nuVtSqSum[3]+0.1767766952966368*fSkin[0]*nuVtSqSum[3]+0.1767766952966368*fEdge[0]*nuVtSqSum[3]+0.1767766952966368*fSkin[1]*nuVtSqSum[2]+0.1767766952966368*fEdge[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fSkin[2]+0.1767766952966368*nuVtSqSum[1]*fEdge[2]; 
  Gdiff2[6] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[26])+0.2041241452319315*nuVtSqSum[2]*fEdge[26]-0.2041241452319315*nuVtSqSum[3]*fSkin[19]+0.2041241452319315*nuVtSqSum[3]*fEdge[19]-0.2041241452319315*nuVtSqSum[0]*fSkin[18]+0.2041241452319315*nuVtSqSum[0]*fEdge[18]+0.1767766952966368*nuVtSqSum[2]*fSkin[17]+0.1767766952966368*nuVtSqSum[2]*fEdge[17]-0.2041241452319315*nuVtSqSum[1]*fSkin[11]+0.2041241452319315*nuVtSqSum[1]*fEdge[11]+0.1767766952966368*nuVtSqSum[3]*fSkin[10]+0.1767766952966368*nuVtSqSum[3]*fEdge[10]+0.1767766952966368*nuVtSqSum[0]*fSkin[9]+0.1767766952966368*nuVtSqSum[0]*fEdge[9]+0.1767766952966368*nuVtSqSum[1]*fSkin[4]+0.1767766952966368*nuVtSqSum[1]*fEdge[4]; 
  Gdiff2[7] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[26])+0.2041241452319315*nuVtSqSum[1]*fEdge[26]-0.2041241452319315*nuVtSqSum[0]*fSkin[19]+0.2041241452319315*nuVtSqSum[0]*fEdge[19]-0.2041241452319315*nuVtSqSum[3]*fSkin[18]+0.2041241452319315*nuVtSqSum[3]*fEdge[18]+0.1767766952966368*nuVtSqSum[1]*fSkin[17]+0.1767766952966368*nuVtSqSum[1]*fEdge[17]-0.2041241452319315*nuVtSqSum[2]*fSkin[11]+0.2041241452319315*nuVtSqSum[2]*fEdge[11]+0.1767766952966368*nuVtSqSum[0]*fSkin[10]+0.1767766952966368*nuVtSqSum[0]*fEdge[10]+0.1767766952966368*nuVtSqSum[3]*fSkin[9]+0.1767766952966368*nuVtSqSum[3]*fEdge[9]+0.1767766952966368*nuVtSqSum[2]*fSkin[4]+0.1767766952966368*nuVtSqSum[2]*fEdge[4]; 
  Gdiff2[8] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[27])+0.2041241452319315*nuVtSqSum[2]*fEdge[27]-0.2041241452319315*nuVtSqSum[3]*fSkin[22]+0.2041241452319315*nuVtSqSum[3]*fEdge[22]-0.2041241452319315*nuVtSqSum[0]*fSkin[21]+0.2041241452319315*nuVtSqSum[0]*fEdge[21]+0.1767766952966368*nuVtSqSum[2]*fSkin[20]+0.1767766952966368*nuVtSqSum[2]*fEdge[20]-0.2041241452319315*nuVtSqSum[1]*fSkin[14]+0.2041241452319315*nuVtSqSum[1]*fEdge[14]+0.1767766952966368*nuVtSqSum[3]*fSkin[13]+0.1767766952966368*nuVtSqSum[3]*fEdge[13]+0.1767766952966368*nuVtSqSum[0]*fSkin[12]+0.1767766952966368*nuVtSqSum[0]*fEdge[12]+0.1767766952966368*nuVtSqSum[1]*fSkin[5]+0.1767766952966368*nuVtSqSum[1]*fEdge[5]; 
  Gdiff2[9] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[27])+0.2041241452319315*nuVtSqSum[1]*fEdge[27]-0.2041241452319315*nuVtSqSum[0]*fSkin[22]+0.2041241452319315*nuVtSqSum[0]*fEdge[22]-0.2041241452319315*nuVtSqSum[3]*fSkin[21]+0.2041241452319315*nuVtSqSum[3]*fEdge[21]+0.1767766952966368*nuVtSqSum[1]*fSkin[20]+0.1767766952966368*nuVtSqSum[1]*fEdge[20]-0.2041241452319315*nuVtSqSum[2]*fSkin[14]+0.2041241452319315*nuVtSqSum[2]*fEdge[14]+0.1767766952966368*nuVtSqSum[0]*fSkin[13]+0.1767766952966368*nuVtSqSum[0]*fEdge[13]+0.1767766952966368*nuVtSqSum[3]*fSkin[12]+0.1767766952966368*nuVtSqSum[3]*fEdge[12]+0.1767766952966368*nuVtSqSum[2]*fSkin[5]+0.1767766952966368*nuVtSqSum[2]*fEdge[5]; 
  Gdiff2[10] = (-0.2041241452319315*nuVtSqSum[3]*fSkin[31])+0.2041241452319315*nuVtSqSum[3]*fEdge[31]-0.2041241452319315*nuVtSqSum[2]*fSkin[30]+0.2041241452319315*nuVtSqSum[2]*fEdge[30]-0.2041241452319315*nuVtSqSum[1]*fSkin[29]+0.2041241452319315*nuVtSqSum[1]*fEdge[29]+0.1767766952966368*nuVtSqSum[3]*fSkin[28]+0.1767766952966368*nuVtSqSum[3]*fEdge[28]-0.2041241452319315*nuVtSqSum[0]*fSkin[25]+0.2041241452319315*nuVtSqSum[0]*fEdge[25]+0.1767766952966368*nuVtSqSum[2]*fSkin[24]+0.1767766952966368*nuVtSqSum[2]*fEdge[24]+0.1767766952966368*nuVtSqSum[1]*fSkin[23]+0.1767766952966368*nuVtSqSum[1]*fEdge[23]+0.1767766952966368*nuVtSqSum[0]*fSkin[15]+0.1767766952966368*nuVtSqSum[0]*fEdge[15]; 
  Gdiff2[11] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[26])+0.2041241452319315*nuVtSqSum[0]*fEdge[26]-0.2041241452319315*nuVtSqSum[1]*fSkin[19]+0.2041241452319315*nuVtSqSum[1]*fEdge[19]-0.2041241452319315*nuVtSqSum[2]*fSkin[18]+0.2041241452319315*nuVtSqSum[2]*fEdge[18]+0.1767766952966368*nuVtSqSum[0]*fSkin[17]+0.1767766952966368*nuVtSqSum[0]*fEdge[17]-0.2041241452319315*nuVtSqSum[3]*fSkin[11]+0.2041241452319315*nuVtSqSum[3]*fEdge[11]+0.1767766952966368*nuVtSqSum[1]*fSkin[10]+0.1767766952966368*nuVtSqSum[1]*fEdge[10]+0.1767766952966368*nuVtSqSum[2]*fSkin[9]+0.1767766952966368*nuVtSqSum[2]*fEdge[9]+0.1767766952966368*nuVtSqSum[3]*fSkin[4]+0.1767766952966368*nuVtSqSum[3]*fEdge[4]; 
  Gdiff2[12] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[27])+0.2041241452319315*nuVtSqSum[0]*fEdge[27]-0.2041241452319315*nuVtSqSum[1]*fSkin[22]+0.2041241452319315*nuVtSqSum[1]*fEdge[22]-0.2041241452319315*nuVtSqSum[2]*fSkin[21]+0.2041241452319315*nuVtSqSum[2]*fEdge[21]+0.1767766952966368*nuVtSqSum[0]*fSkin[20]+0.1767766952966368*nuVtSqSum[0]*fEdge[20]-0.2041241452319315*nuVtSqSum[3]*fSkin[14]+0.2041241452319315*nuVtSqSum[3]*fEdge[14]+0.1767766952966368*nuVtSqSum[1]*fSkin[13]+0.1767766952966368*nuVtSqSum[1]*fEdge[13]+0.1767766952966368*nuVtSqSum[2]*fSkin[12]+0.1767766952966368*nuVtSqSum[2]*fEdge[12]+0.1767766952966368*nuVtSqSum[3]*fSkin[5]+0.1767766952966368*nuVtSqSum[3]*fEdge[5]; 
  Gdiff2[13] = (-0.2041241452319315*nuVtSqSum[2]*fSkin[31])+0.2041241452319315*nuVtSqSum[2]*fEdge[31]-0.2041241452319315*nuVtSqSum[3]*fSkin[30]+0.2041241452319315*nuVtSqSum[3]*fEdge[30]-0.2041241452319315*nuVtSqSum[0]*fSkin[29]+0.2041241452319315*nuVtSqSum[0]*fEdge[29]+0.1767766952966368*nuVtSqSum[2]*fSkin[28]+0.1767766952966368*nuVtSqSum[2]*fEdge[28]-0.2041241452319315*nuVtSqSum[1]*fSkin[25]+0.2041241452319315*nuVtSqSum[1]*fEdge[25]+0.1767766952966368*nuVtSqSum[3]*fSkin[24]+0.1767766952966368*nuVtSqSum[3]*fEdge[24]+0.1767766952966368*nuVtSqSum[0]*fSkin[23]+0.1767766952966368*nuVtSqSum[0]*fEdge[23]+0.1767766952966368*nuVtSqSum[1]*fSkin[15]+0.1767766952966368*nuVtSqSum[1]*fEdge[15]; 
  Gdiff2[14] = (-0.2041241452319315*nuVtSqSum[1]*fSkin[31])+0.2041241452319315*nuVtSqSum[1]*fEdge[31]-0.2041241452319315*nuVtSqSum[0]*fSkin[30]+0.2041241452319315*nuVtSqSum[0]*fEdge[30]-0.2041241452319315*nuVtSqSum[3]*fSkin[29]+0.2041241452319315*nuVtSqSum[3]*fEdge[29]+0.1767766952966368*nuVtSqSum[1]*fSkin[28]+0.1767766952966368*nuVtSqSum[1]*fEdge[28]-0.2041241452319315*nuVtSqSum[2]*fSkin[25]+0.2041241452319315*nuVtSqSum[2]*fEdge[25]+0.1767766952966368*nuVtSqSum[0]*fSkin[24]+0.1767766952966368*nuVtSqSum[0]*fEdge[24]+0.1767766952966368*nuVtSqSum[3]*fSkin[23]+0.1767766952966368*nuVtSqSum[3]*fEdge[23]+0.1767766952966368*nuVtSqSum[2]*fSkin[15]+0.1767766952966368*nuVtSqSum[2]*fEdge[15]; 
  Gdiff2[15] = (-0.2041241452319315*nuVtSqSum[0]*fSkin[31])+0.2041241452319315*nuVtSqSum[0]*fEdge[31]-0.2041241452319315*nuVtSqSum[1]*fSkin[30]+0.2041241452319315*nuVtSqSum[1]*fEdge[30]-0.2041241452319315*nuVtSqSum[2]*fSkin[29]+0.2041241452319315*nuVtSqSum[2]*fEdge[29]+0.1767766952966368*nuVtSqSum[0]*fSkin[28]+0.1767766952966368*nuVtSqSum[0]*fEdge[28]-0.2041241452319315*nuVtSqSum[3]*fSkin[25]+0.2041241452319315*nuVtSqSum[3]*fEdge[25]+0.1767766952966368*nuVtSqSum[1]*fSkin[24]+0.1767766952966368*nuVtSqSum[1]*fEdge[24]+0.1767766952966368*nuVtSqSum[2]*fSkin[23]+0.1767766952966368*nuVtSqSum[2]*fEdge[23]+0.1767766952966368*nuVtSqSum[3]*fSkin[15]+0.1767766952966368*nuVtSqSum[3]*fEdge[15]; 

  Gdiff[0] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[16])-0.3827327723098713*nuVtSqSum[3]*fEdge[16]-0.3827327723098713*nuVtSqSum[2]*fSkin[8]-0.3827327723098713*nuVtSqSum[2]*fEdge[8]-0.3827327723098713*nuVtSqSum[1]*fSkin[7]-0.3827327723098713*nuVtSqSum[1]*fEdge[7]+0.3977475644174328*nuVtSqSum[3]*fSkin[6]-0.3977475644174328*nuVtSqSum[3]*fEdge[6]-0.3827327723098713*nuVtSqSum[0]*fSkin[3]-0.3827327723098713*nuVtSqSum[0]*fEdge[3]+0.3977475644174328*fSkin[2]*nuVtSqSum[2]-0.3977475644174328*fEdge[2]*nuVtSqSum[2]+0.3977475644174328*fSkin[1]*nuVtSqSum[1]-0.3977475644174328*fEdge[1]*nuVtSqSum[1]+0.3977475644174328*fSkin[0]*nuVtSqSum[0]-0.3977475644174328*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[16])-0.3827327723098713*nuVtSqSum[2]*fEdge[16]-0.3827327723098713*nuVtSqSum[3]*fSkin[8]-0.3827327723098713*nuVtSqSum[3]*fEdge[8]-0.3827327723098713*nuVtSqSum[0]*fSkin[7]-0.3827327723098713*nuVtSqSum[0]*fEdge[7]+0.3977475644174328*nuVtSqSum[2]*fSkin[6]-0.3977475644174328*nuVtSqSum[2]*fEdge[6]+0.3977475644174328*fSkin[2]*nuVtSqSum[3]-0.3977475644174328*fEdge[2]*nuVtSqSum[3]-0.3827327723098713*nuVtSqSum[1]*fSkin[3]-0.3827327723098713*nuVtSqSum[1]*fEdge[3]+0.3977475644174328*fSkin[0]*nuVtSqSum[1]-0.3977475644174328*fEdge[0]*nuVtSqSum[1]+0.3977475644174328*nuVtSqSum[0]*fSkin[1]-0.3977475644174328*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[16])-0.3827327723098713*nuVtSqSum[1]*fEdge[16]-0.3827327723098713*nuVtSqSum[0]*fSkin[8]-0.3827327723098713*nuVtSqSum[0]*fEdge[8]-0.3827327723098713*nuVtSqSum[3]*fSkin[7]-0.3827327723098713*nuVtSqSum[3]*fEdge[7]+0.3977475644174328*nuVtSqSum[1]*fSkin[6]-0.3977475644174328*nuVtSqSum[1]*fEdge[6]+0.3977475644174328*fSkin[1]*nuVtSqSum[3]-0.3977475644174328*fEdge[1]*nuVtSqSum[3]-0.3827327723098713*nuVtSqSum[2]*fSkin[3]-0.3827327723098713*nuVtSqSum[2]*fEdge[3]+0.3977475644174328*fSkin[0]*nuVtSqSum[2]-0.3977475644174328*fEdge[0]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[0]*fSkin[2]-0.3977475644174328*nuVtSqSum[0]*fEdge[2]; 
  Gdiff[3] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[26])-0.3827327723098713*nuVtSqSum[3]*fEdge[26]-0.3827327723098713*nuVtSqSum[2]*fSkin[19]-0.3827327723098713*nuVtSqSum[2]*fEdge[19]-0.3827327723098713*nuVtSqSum[1]*fSkin[18]-0.3827327723098713*nuVtSqSum[1]*fEdge[18]+0.3977475644174328*nuVtSqSum[3]*fSkin[17]-0.3977475644174328*nuVtSqSum[3]*fEdge[17]-0.3827327723098713*nuVtSqSum[0]*fSkin[11]-0.3827327723098713*nuVtSqSum[0]*fEdge[11]+0.3977475644174328*nuVtSqSum[2]*fSkin[10]-0.3977475644174328*nuVtSqSum[2]*fEdge[10]+0.3977475644174328*nuVtSqSum[1]*fSkin[9]-0.3977475644174328*nuVtSqSum[1]*fEdge[9]+0.3977475644174328*nuVtSqSum[0]*fSkin[4]-0.3977475644174328*nuVtSqSum[0]*fEdge[4]; 
  Gdiff[4] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[27])-0.3827327723098713*nuVtSqSum[3]*fEdge[27]-0.3827327723098713*nuVtSqSum[2]*fSkin[22]-0.3827327723098713*nuVtSqSum[2]*fEdge[22]-0.3827327723098713*nuVtSqSum[1]*fSkin[21]-0.3827327723098713*nuVtSqSum[1]*fEdge[21]+0.3977475644174328*nuVtSqSum[3]*fSkin[20]-0.3977475644174328*nuVtSqSum[3]*fEdge[20]-0.3827327723098713*nuVtSqSum[0]*fSkin[14]-0.3827327723098713*nuVtSqSum[0]*fEdge[14]+0.3977475644174328*nuVtSqSum[2]*fSkin[13]-0.3977475644174328*nuVtSqSum[2]*fEdge[13]+0.3977475644174328*nuVtSqSum[1]*fSkin[12]-0.3977475644174328*nuVtSqSum[1]*fEdge[12]+0.3977475644174328*nuVtSqSum[0]*fSkin[5]-0.3977475644174328*nuVtSqSum[0]*fEdge[5]; 
  Gdiff[5] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[16])-0.3827327723098713*nuVtSqSum[0]*fEdge[16]-0.3827327723098713*nuVtSqSum[1]*fSkin[8]-0.3827327723098713*nuVtSqSum[1]*fEdge[8]-0.3827327723098713*nuVtSqSum[2]*fSkin[7]-0.3827327723098713*nuVtSqSum[2]*fEdge[7]+0.3977475644174328*nuVtSqSum[0]*fSkin[6]-0.3977475644174328*nuVtSqSum[0]*fEdge[6]-0.3827327723098713*fSkin[3]*nuVtSqSum[3]-0.3827327723098713*fEdge[3]*nuVtSqSum[3]+0.3977475644174328*fSkin[0]*nuVtSqSum[3]-0.3977475644174328*fEdge[0]*nuVtSqSum[3]+0.3977475644174328*fSkin[1]*nuVtSqSum[2]-0.3977475644174328*fEdge[1]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[1]*fSkin[2]-0.3977475644174328*nuVtSqSum[1]*fEdge[2]; 
  Gdiff[6] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[26])-0.3827327723098713*nuVtSqSum[2]*fEdge[26]-0.3827327723098713*nuVtSqSum[3]*fSkin[19]-0.3827327723098713*nuVtSqSum[3]*fEdge[19]-0.3827327723098713*nuVtSqSum[0]*fSkin[18]-0.3827327723098713*nuVtSqSum[0]*fEdge[18]+0.3977475644174328*nuVtSqSum[2]*fSkin[17]-0.3977475644174328*nuVtSqSum[2]*fEdge[17]-0.3827327723098713*nuVtSqSum[1]*fSkin[11]-0.3827327723098713*nuVtSqSum[1]*fEdge[11]+0.3977475644174328*nuVtSqSum[3]*fSkin[10]-0.3977475644174328*nuVtSqSum[3]*fEdge[10]+0.3977475644174328*nuVtSqSum[0]*fSkin[9]-0.3977475644174328*nuVtSqSum[0]*fEdge[9]+0.3977475644174328*nuVtSqSum[1]*fSkin[4]-0.3977475644174328*nuVtSqSum[1]*fEdge[4]; 
  Gdiff[7] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[26])-0.3827327723098713*nuVtSqSum[1]*fEdge[26]-0.3827327723098713*nuVtSqSum[0]*fSkin[19]-0.3827327723098713*nuVtSqSum[0]*fEdge[19]-0.3827327723098713*nuVtSqSum[3]*fSkin[18]-0.3827327723098713*nuVtSqSum[3]*fEdge[18]+0.3977475644174328*nuVtSqSum[1]*fSkin[17]-0.3977475644174328*nuVtSqSum[1]*fEdge[17]-0.3827327723098713*nuVtSqSum[2]*fSkin[11]-0.3827327723098713*nuVtSqSum[2]*fEdge[11]+0.3977475644174328*nuVtSqSum[0]*fSkin[10]-0.3977475644174328*nuVtSqSum[0]*fEdge[10]+0.3977475644174328*nuVtSqSum[3]*fSkin[9]-0.3977475644174328*nuVtSqSum[3]*fEdge[9]+0.3977475644174328*nuVtSqSum[2]*fSkin[4]-0.3977475644174328*nuVtSqSum[2]*fEdge[4]; 
  Gdiff[8] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[27])-0.3827327723098713*nuVtSqSum[2]*fEdge[27]-0.3827327723098713*nuVtSqSum[3]*fSkin[22]-0.3827327723098713*nuVtSqSum[3]*fEdge[22]-0.3827327723098713*nuVtSqSum[0]*fSkin[21]-0.3827327723098713*nuVtSqSum[0]*fEdge[21]+0.3977475644174328*nuVtSqSum[2]*fSkin[20]-0.3977475644174328*nuVtSqSum[2]*fEdge[20]-0.3827327723098713*nuVtSqSum[1]*fSkin[14]-0.3827327723098713*nuVtSqSum[1]*fEdge[14]+0.3977475644174328*nuVtSqSum[3]*fSkin[13]-0.3977475644174328*nuVtSqSum[3]*fEdge[13]+0.3977475644174328*nuVtSqSum[0]*fSkin[12]-0.3977475644174328*nuVtSqSum[0]*fEdge[12]+0.3977475644174328*nuVtSqSum[1]*fSkin[5]-0.3977475644174328*nuVtSqSum[1]*fEdge[5]; 
  Gdiff[9] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[27])-0.3827327723098713*nuVtSqSum[1]*fEdge[27]-0.3827327723098713*nuVtSqSum[0]*fSkin[22]-0.3827327723098713*nuVtSqSum[0]*fEdge[22]-0.3827327723098713*nuVtSqSum[3]*fSkin[21]-0.3827327723098713*nuVtSqSum[3]*fEdge[21]+0.3977475644174328*nuVtSqSum[1]*fSkin[20]-0.3977475644174328*nuVtSqSum[1]*fEdge[20]-0.3827327723098713*nuVtSqSum[2]*fSkin[14]-0.3827327723098713*nuVtSqSum[2]*fEdge[14]+0.3977475644174328*nuVtSqSum[0]*fSkin[13]-0.3977475644174328*nuVtSqSum[0]*fEdge[13]+0.3977475644174328*nuVtSqSum[3]*fSkin[12]-0.3977475644174328*nuVtSqSum[3]*fEdge[12]+0.3977475644174328*nuVtSqSum[2]*fSkin[5]-0.3977475644174328*nuVtSqSum[2]*fEdge[5]; 
  Gdiff[10] = (-0.3827327723098713*nuVtSqSum[3]*fSkin[31])-0.3827327723098713*nuVtSqSum[3]*fEdge[31]-0.3827327723098713*nuVtSqSum[2]*fSkin[30]-0.3827327723098713*nuVtSqSum[2]*fEdge[30]-0.3827327723098713*nuVtSqSum[1]*fSkin[29]-0.3827327723098713*nuVtSqSum[1]*fEdge[29]+0.3977475644174328*nuVtSqSum[3]*fSkin[28]-0.3977475644174328*nuVtSqSum[3]*fEdge[28]-0.3827327723098713*nuVtSqSum[0]*fSkin[25]-0.3827327723098713*nuVtSqSum[0]*fEdge[25]+0.3977475644174328*nuVtSqSum[2]*fSkin[24]-0.3977475644174328*nuVtSqSum[2]*fEdge[24]+0.3977475644174328*nuVtSqSum[1]*fSkin[23]-0.3977475644174328*nuVtSqSum[1]*fEdge[23]+0.3977475644174328*nuVtSqSum[0]*fSkin[15]-0.3977475644174328*nuVtSqSum[0]*fEdge[15]; 
  Gdiff[11] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[26])-0.3827327723098713*nuVtSqSum[0]*fEdge[26]-0.3827327723098713*nuVtSqSum[1]*fSkin[19]-0.3827327723098713*nuVtSqSum[1]*fEdge[19]-0.3827327723098713*nuVtSqSum[2]*fSkin[18]-0.3827327723098713*nuVtSqSum[2]*fEdge[18]+0.3977475644174328*nuVtSqSum[0]*fSkin[17]-0.3977475644174328*nuVtSqSum[0]*fEdge[17]-0.3827327723098713*nuVtSqSum[3]*fSkin[11]-0.3827327723098713*nuVtSqSum[3]*fEdge[11]+0.3977475644174328*nuVtSqSum[1]*fSkin[10]-0.3977475644174328*nuVtSqSum[1]*fEdge[10]+0.3977475644174328*nuVtSqSum[2]*fSkin[9]-0.3977475644174328*nuVtSqSum[2]*fEdge[9]+0.3977475644174328*nuVtSqSum[3]*fSkin[4]-0.3977475644174328*nuVtSqSum[3]*fEdge[4]; 
  Gdiff[12] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[27])-0.3827327723098713*nuVtSqSum[0]*fEdge[27]-0.3827327723098713*nuVtSqSum[1]*fSkin[22]-0.3827327723098713*nuVtSqSum[1]*fEdge[22]-0.3827327723098713*nuVtSqSum[2]*fSkin[21]-0.3827327723098713*nuVtSqSum[2]*fEdge[21]+0.3977475644174328*nuVtSqSum[0]*fSkin[20]-0.3977475644174328*nuVtSqSum[0]*fEdge[20]-0.3827327723098713*nuVtSqSum[3]*fSkin[14]-0.3827327723098713*nuVtSqSum[3]*fEdge[14]+0.3977475644174328*nuVtSqSum[1]*fSkin[13]-0.3977475644174328*nuVtSqSum[1]*fEdge[13]+0.3977475644174328*nuVtSqSum[2]*fSkin[12]-0.3977475644174328*nuVtSqSum[2]*fEdge[12]+0.3977475644174328*nuVtSqSum[3]*fSkin[5]-0.3977475644174328*nuVtSqSum[3]*fEdge[5]; 
  Gdiff[13] = (-0.3827327723098713*nuVtSqSum[2]*fSkin[31])-0.3827327723098713*nuVtSqSum[2]*fEdge[31]-0.3827327723098713*nuVtSqSum[3]*fSkin[30]-0.3827327723098713*nuVtSqSum[3]*fEdge[30]-0.3827327723098713*nuVtSqSum[0]*fSkin[29]-0.3827327723098713*nuVtSqSum[0]*fEdge[29]+0.3977475644174328*nuVtSqSum[2]*fSkin[28]-0.3977475644174328*nuVtSqSum[2]*fEdge[28]-0.3827327723098713*nuVtSqSum[1]*fSkin[25]-0.3827327723098713*nuVtSqSum[1]*fEdge[25]+0.3977475644174328*nuVtSqSum[3]*fSkin[24]-0.3977475644174328*nuVtSqSum[3]*fEdge[24]+0.3977475644174328*nuVtSqSum[0]*fSkin[23]-0.3977475644174328*nuVtSqSum[0]*fEdge[23]+0.3977475644174328*nuVtSqSum[1]*fSkin[15]-0.3977475644174328*nuVtSqSum[1]*fEdge[15]; 
  Gdiff[14] = (-0.3827327723098713*nuVtSqSum[1]*fSkin[31])-0.3827327723098713*nuVtSqSum[1]*fEdge[31]-0.3827327723098713*nuVtSqSum[0]*fSkin[30]-0.3827327723098713*nuVtSqSum[0]*fEdge[30]-0.3827327723098713*nuVtSqSum[3]*fSkin[29]-0.3827327723098713*nuVtSqSum[3]*fEdge[29]+0.3977475644174328*nuVtSqSum[1]*fSkin[28]-0.3977475644174328*nuVtSqSum[1]*fEdge[28]-0.3827327723098713*nuVtSqSum[2]*fSkin[25]-0.3827327723098713*nuVtSqSum[2]*fEdge[25]+0.3977475644174328*nuVtSqSum[0]*fSkin[24]-0.3977475644174328*nuVtSqSum[0]*fEdge[24]+0.3977475644174328*nuVtSqSum[3]*fSkin[23]-0.3977475644174328*nuVtSqSum[3]*fEdge[23]+0.3977475644174328*nuVtSqSum[2]*fSkin[15]-0.3977475644174328*nuVtSqSum[2]*fEdge[15]; 
  Gdiff[15] = (-0.3827327723098713*nuVtSqSum[0]*fSkin[31])-0.3827327723098713*nuVtSqSum[0]*fEdge[31]-0.3827327723098713*nuVtSqSum[1]*fSkin[30]-0.3827327723098713*nuVtSqSum[1]*fEdge[30]-0.3827327723098713*nuVtSqSum[2]*fSkin[29]-0.3827327723098713*nuVtSqSum[2]*fEdge[29]+0.3977475644174328*nuVtSqSum[0]*fSkin[28]-0.3977475644174328*nuVtSqSum[0]*fEdge[28]-0.3827327723098713*nuVtSqSum[3]*fSkin[25]-0.3827327723098713*nuVtSqSum[3]*fEdge[25]+0.3977475644174328*nuVtSqSum[1]*fSkin[24]-0.3977475644174328*nuVtSqSum[1]*fEdge[24]+0.3977475644174328*nuVtSqSum[2]*fSkin[23]-0.3977475644174328*nuVtSqSum[2]*fEdge[23]+0.3977475644174328*nuVtSqSum[3]*fSkin[15]-0.3977475644174328*nuVtSqSum[3]*fEdge[15]; 

  Ghat[0] = Gdiff[0]*rdv2+0.25*alphaDrSurf[5]*fUpwind[5]+0.25*alphaDrSurf[2]*fUpwind[2]+0.25*alphaDrSurf[1]*fUpwind[1]+0.25*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.25*alphaDrSurf[2]*fUpwind[5]+0.25*fUpwind[2]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[1]+0.25*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2+0.25*alphaDrSurf[1]*fUpwind[5]+0.25*fUpwind[1]*alphaDrSurf[5]+0.25*alphaDrSurf[0]*fUpwind[2]+0.25*fUpwind[0]*alphaDrSurf[2]; 
  Ghat[3] = Gdiff[3]*rdv2+0.25*alphaDrSurf[5]*fUpwind[11]+0.25*alphaDrSurf[2]*fUpwind[7]+0.25*alphaDrSurf[1]*fUpwind[6]+0.25*alphaDrSurf[0]*fUpwind[3]; 
  Ghat[4] = Gdiff[4]*rdv2+0.25*alphaDrSurf[5]*fUpwind[12]+0.25*alphaDrSurf[2]*fUpwind[9]+0.25*alphaDrSurf[1]*fUpwind[8]+0.25*alphaDrSurf[0]*fUpwind[4]; 
  Ghat[5] = Gdiff[5]*rdv2+0.25*alphaDrSurf[0]*fUpwind[5]+0.25*fUpwind[0]*alphaDrSurf[5]+0.25*alphaDrSurf[1]*fUpwind[2]+0.25*fUpwind[1]*alphaDrSurf[2]; 
  Ghat[6] = Gdiff[6]*rdv2+0.25*alphaDrSurf[2]*fUpwind[11]+0.25*alphaDrSurf[5]*fUpwind[7]+0.25*alphaDrSurf[0]*fUpwind[6]+0.25*alphaDrSurf[1]*fUpwind[3]; 
  Ghat[7] = Gdiff[7]*rdv2+0.25*alphaDrSurf[1]*fUpwind[11]+0.25*alphaDrSurf[0]*fUpwind[7]+0.25*alphaDrSurf[5]*fUpwind[6]+0.25*alphaDrSurf[2]*fUpwind[3]; 
  Ghat[8] = Gdiff[8]*rdv2+0.25*alphaDrSurf[2]*fUpwind[12]+0.25*alphaDrSurf[5]*fUpwind[9]+0.25*alphaDrSurf[0]*fUpwind[8]+0.25*alphaDrSurf[1]*fUpwind[4]; 
  Ghat[9] = Gdiff[9]*rdv2+0.25*alphaDrSurf[1]*fUpwind[12]+0.25*alphaDrSurf[0]*fUpwind[9]+0.25*alphaDrSurf[5]*fUpwind[8]+0.25*alphaDrSurf[2]*fUpwind[4]; 
  Ghat[10] = Gdiff[10]*rdv2+0.25*alphaDrSurf[5]*fUpwind[15]+0.25*alphaDrSurf[2]*fUpwind[14]+0.25*alphaDrSurf[1]*fUpwind[13]+0.25*alphaDrSurf[0]*fUpwind[10]; 
  Ghat[11] = Gdiff[11]*rdv2+0.25*alphaDrSurf[0]*fUpwind[11]+0.25*alphaDrSurf[1]*fUpwind[7]+0.25*alphaDrSurf[2]*fUpwind[6]+0.25*fUpwind[3]*alphaDrSurf[5]; 
  Ghat[12] = Gdiff[12]*rdv2+0.25*alphaDrSurf[0]*fUpwind[12]+0.25*alphaDrSurf[1]*fUpwind[9]+0.25*alphaDrSurf[2]*fUpwind[8]+0.25*fUpwind[4]*alphaDrSurf[5]; 
  Ghat[13] = Gdiff[13]*rdv2+0.25*alphaDrSurf[2]*fUpwind[15]+0.25*alphaDrSurf[5]*fUpwind[14]+0.25*alphaDrSurf[0]*fUpwind[13]+0.25*alphaDrSurf[1]*fUpwind[10]; 
  Ghat[14] = Gdiff[14]*rdv2+0.25*alphaDrSurf[1]*fUpwind[15]+0.25*alphaDrSurf[0]*fUpwind[14]+0.25*alphaDrSurf[5]*fUpwind[13]+0.25*alphaDrSurf[2]*fUpwind[10]; 
  Ghat[15] = Gdiff[15]*rdv2+0.25*alphaDrSurf[0]*fUpwind[15]+0.25*alphaDrSurf[1]*fUpwind[14]+0.25*alphaDrSurf[2]*fUpwind[13]+0.25*alphaDrSurf[5]*fUpwind[10]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[3] += 0.75*nuVtSqSum[3]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[6]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[3]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[1]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[0]*rdvSq4-1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[3]*rdv2; 
  out[5] += 0.7071067811865475*Ghat[4]*rdv2; 
  out[6] += 0.7071067811865475*Ghat[5]*rdv2; 
  out[7] += 0.75*nuVtSqSum[2]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[6]*rdvSq4+0.4330127018922193*fSkin[2]*nuVtSqSum[3]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[1]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[8] += 0.75*nuVtSqSum[1]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[6]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[3]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[2]*rdvSq4-1.224744871391589*Gdiff2[2]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[9] += 0.7071067811865475*Ghat[6]*rdv2; 
  out[10] += 0.7071067811865475*Ghat[7]*rdv2; 
  out[11] += 0.75*nuVtSqSum[3]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[18]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[11]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[9]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[4]*rdvSq4-1.224744871391589*Gdiff2[3]*rdvSq4-1.224744871391589*Ghat[3]*rdv2; 
  out[12] += 0.7071067811865475*Ghat[8]*rdv2; 
  out[13] += 0.7071067811865475*Ghat[9]*rdv2; 
  out[14] += 0.75*nuVtSqSum[3]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[21]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[14]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[12]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[5]*rdvSq4-1.224744871391589*Gdiff2[4]*rdvSq4-1.224744871391589*Ghat[4]*rdv2; 
  out[15] += 0.7071067811865475*Ghat[10]*rdv2; 
  out[16] += 0.75*nuVtSqSum[0]*fSkin[16]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[8]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[7]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[6]*rdvSq4-1.224744871391589*Gdiff2[5]*rdvSq4+0.75*fSkin[3]*nuVtSqSum[3]*rdvSq4+0.4330127018922193*fSkin[0]*nuVtSqSum[3]*rdvSq4+0.4330127018922193*fSkin[1]*nuVtSqSum[2]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[2]*rdvSq4-1.224744871391589*Ghat[5]*rdv2; 
  out[17] += 0.7071067811865475*Ghat[11]*rdv2; 
  out[18] += 0.75*nuVtSqSum[2]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[18]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[11]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[9]*rdvSq4-1.224744871391589*Gdiff2[6]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[6]*rdv2; 
  out[19] += 0.75*nuVtSqSum[1]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[18]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[11]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[9]*rdvSq4-1.224744871391589*Gdiff2[7]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[7]*rdv2; 
  out[20] += 0.7071067811865475*Ghat[12]*rdv2; 
  out[21] += 0.75*nuVtSqSum[2]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[21]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[14]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[12]*rdvSq4-1.224744871391589*Gdiff2[8]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[8]*rdv2; 
  out[22] += 0.75*nuVtSqSum[1]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[21]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[14]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[12]*rdvSq4-1.224744871391589*Gdiff2[9]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[9]*rdv2; 
  out[23] += 0.7071067811865475*Ghat[13]*rdv2; 
  out[24] += 0.7071067811865475*Ghat[14]*rdv2; 
  out[25] += 0.75*nuVtSqSum[3]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[29]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[25]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[24]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[23]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[15]*rdvSq4-1.224744871391589*Gdiff2[10]*rdvSq4-1.224744871391589*Ghat[10]*rdv2; 
  out[26] += 0.75*nuVtSqSum[0]*fSkin[26]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[19]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[18]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[17]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[11]*rdvSq4-1.224744871391589*Gdiff2[11]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[10]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[9]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[4]*rdvSq4-1.224744871391589*Ghat[11]*rdv2; 
  out[27] += 0.75*nuVtSqSum[0]*fSkin[27]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[22]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[21]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[20]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[14]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[13]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[12]*rdvSq4-1.224744871391589*Gdiff2[12]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[5]*rdvSq4-1.224744871391589*Ghat[12]*rdv2; 
  out[28] += 0.7071067811865475*Ghat[15]*rdv2; 
  out[29] += 0.75*nuVtSqSum[2]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[29]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[25]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[24]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[23]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[15]*rdvSq4-1.224744871391589*Gdiff2[13]*rdvSq4-1.224744871391589*Ghat[13]*rdv2; 
  out[30] += 0.75*nuVtSqSum[1]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[0]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[29]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[25]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[24]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[23]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[15]*rdvSq4-1.224744871391589*Gdiff2[14]*rdvSq4-1.224744871391589*Ghat[14]*rdv2; 
  out[31] += 0.75*nuVtSqSum[0]*fSkin[31]*rdvSq4+0.75*nuVtSqSum[1]*fSkin[30]*rdvSq4+0.75*nuVtSqSum[2]*fSkin[29]*rdvSq4+0.4330127018922193*nuVtSqSum[0]*fSkin[28]*rdvSq4+0.75*nuVtSqSum[3]*fSkin[25]*rdvSq4+0.4330127018922193*nuVtSqSum[1]*fSkin[24]*rdvSq4+0.4330127018922193*nuVtSqSum[2]*fSkin[23]*rdvSq4+0.4330127018922193*nuVtSqSum[3]*fSkin[15]*rdvSq4-1.224744871391589*Gdiff2[15]*rdvSq4-1.224744871391589*Ghat[15]*rdv2; 

  } 
} 
