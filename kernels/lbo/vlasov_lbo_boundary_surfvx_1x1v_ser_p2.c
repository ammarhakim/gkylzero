#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x1v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_boundary_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const int edge, const double *fSkin, const double *fEdge, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         Cell-center coordinates. 
  // dxv[2]:       Cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fSkin/Edge:    Distribution function in cells 
  // out:           Incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf[3] = {0.0}; 
  double fUpwindQuad[3] = {0.0};
  double fUpwind[3] = {0.0};;
  double Ghat[3] = {0.0}; 
  double Gdiff[3] = {0.0}; 
  double Gdiff2[3] = {0.0}; 

  if (edge == -1) { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -0.5*(2.0*sumNuUx[2]+((-2.0*w[1])-1.0*dxv[1])*nuSum[2]); 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(1, fSkin); 
  } else { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(-1, fEdge); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(1, fSkin); 
  } else { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(-1, fEdge); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(1, fSkin); 
  } else { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(-1, fEdge); 
  } 

  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

  Gdiff2[0] = 0.2445699350390395*nuVtSqSum[1]*fSkin[7]+0.2445699350390395*nuVtSqSum[1]*fEdge[7]+0.3518228202874282*nuVtSqSum[2]*fSkin[6]-0.3518228202874282*nuVtSqSum[2]*fEdge[6]+0.2445699350390395*nuVtSqSum[0]*fSkin[5]+0.2445699350390395*nuVtSqSum[0]*fEdge[5]+0.25*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[2]*fEdge[4]+0.3518228202874282*nuVtSqSum[1]*fSkin[3]-0.3518228202874282*nuVtSqSum[1]*fEdge[3]+0.3518228202874282*nuVtSqSum[0]*fSkin[2]-0.3518228202874282*nuVtSqSum[0]*fEdge[2]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.21875*nuVtSqSum[2]*fSkin[7]+0.2445699350390395*nuVtSqSum[0]*fSkin[7]+0.21875*nuVtSqSum[2]*fEdge[7]+0.2445699350390395*nuVtSqSum[0]*fEdge[7]+0.3146798968793527*nuVtSqSum[1]*fSkin[6]-0.3146798968793527*nuVtSqSum[1]*fEdge[6]+0.2445699350390395*nuVtSqSum[1]*fSkin[5]+0.2445699350390395*nuVtSqSum[1]*fEdge[5]+0.223606797749979*nuVtSqSum[1]*fSkin[4]+0.223606797749979*nuVtSqSum[1]*fEdge[4]+0.3146798968793526*nuVtSqSum[2]*fSkin[3]+0.3518228202874282*nuVtSqSum[0]*fSkin[3]-0.3146798968793526*nuVtSqSum[2]*fEdge[3]-0.3518228202874282*nuVtSqSum[0]*fEdge[3]+0.223606797749979*fSkin[1]*nuVtSqSum[2]+0.223606797749979*fEdge[1]*nuVtSqSum[2]+0.3518228202874282*nuVtSqSum[1]*fSkin[2]-0.3518228202874282*nuVtSqSum[1]*fEdge[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.21875*nuVtSqSum[1]*fSkin[7]+0.21875*nuVtSqSum[1]*fEdge[7]+0.2247713549138233*nuVtSqSum[2]*fSkin[6]+0.3518228202874282*nuVtSqSum[0]*fSkin[6]-0.2247713549138233*nuVtSqSum[2]*fEdge[6]-0.3518228202874282*nuVtSqSum[0]*fEdge[6]+0.2445699350390395*nuVtSqSum[2]*fSkin[5]+0.2445699350390395*nuVtSqSum[2]*fEdge[5]+0.159719141249985*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[0]*fSkin[4]+0.159719141249985*nuVtSqSum[2]*fEdge[4]+0.25*nuVtSqSum[0]*fEdge[4]+0.3146798968793526*nuVtSqSum[1]*fSkin[3]-0.3146798968793526*nuVtSqSum[1]*fEdge[3]+0.3518228202874282*fSkin[2]*nuVtSqSum[2]-0.3518228202874282*fEdge[2]*nuVtSqSum[2]+0.25*fSkin[0]*nuVtSqSum[2]+0.25*fEdge[0]*nuVtSqSum[2]+0.223606797749979*fSkin[1]*nuVtSqSum[1]+0.223606797749979*fEdge[1]*nuVtSqSum[1]; 

  Gdiff[0] = (-0.6708203932499369*nuVtSqSum[1]*fSkin[7])+0.6708203932499369*nuVtSqSum[1]*fEdge[7]-1.190784930203603*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[2]*fEdge[6]-0.6708203932499369*nuVtSqSum[0]*fSkin[5]+0.6708203932499369*nuVtSqSum[0]*fEdge[5]-0.9375*nuVtSqSum[2]*fSkin[4]+0.9375*nuVtSqSum[2]*fEdge[4]-1.190784930203603*nuVtSqSum[1]*fSkin[3]-1.190784930203603*nuVtSqSum[1]*fEdge[3]-1.190784930203603*nuVtSqSum[0]*fSkin[2]-1.190784930203603*nuVtSqSum[0]*fEdge[2]-0.9375*fSkin[1]*nuVtSqSum[1]+0.9375*fEdge[1]*nuVtSqSum[1]-0.9375*fSkin[0]*nuVtSqSum[0]+0.9375*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = (-0.5999999999999999*nuVtSqSum[2]*fSkin[7])-0.6708203932499369*nuVtSqSum[0]*fSkin[7]+0.5999999999999999*nuVtSqSum[2]*fEdge[7]+0.6708203932499369*nuVtSqSum[0]*fEdge[7]-1.06507042020704*nuVtSqSum[1]*fSkin[6]-1.06507042020704*nuVtSqSum[1]*fEdge[6]-0.6708203932499369*nuVtSqSum[1]*fSkin[5]+0.6708203932499369*nuVtSqSum[1]*fEdge[5]-0.8385254915624212*nuVtSqSum[1]*fSkin[4]+0.8385254915624212*nuVtSqSum[1]*fEdge[4]-1.06507042020704*nuVtSqSum[2]*fSkin[3]-1.190784930203603*nuVtSqSum[0]*fSkin[3]-1.06507042020704*nuVtSqSum[2]*fEdge[3]-1.190784930203603*nuVtSqSum[0]*fEdge[3]-0.8385254915624212*fSkin[1]*nuVtSqSum[2]+0.8385254915624212*fEdge[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fSkin[2]-1.190784930203603*nuVtSqSum[1]*fEdge[2]-0.9375*fSkin[0]*nuVtSqSum[1]+0.9375*fEdge[0]*nuVtSqSum[1]-0.9375*nuVtSqSum[0]*fSkin[1]+0.9375*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = (-0.5999999999999999*nuVtSqSum[1]*fSkin[7])+0.5999999999999999*nuVtSqSum[1]*fEdge[7]-0.7607645858621712*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[0]*fSkin[6]-0.7607645858621712*nuVtSqSum[2]*fEdge[6]-1.190784930203603*nuVtSqSum[0]*fEdge[6]-0.6708203932499369*nuVtSqSum[2]*fSkin[5]+0.6708203932499369*nuVtSqSum[2]*fEdge[5]-0.5989467796874438*nuVtSqSum[2]*fSkin[4]-0.9375*nuVtSqSum[0]*fSkin[4]+0.5989467796874438*nuVtSqSum[2]*fEdge[4]+0.9375*nuVtSqSum[0]*fEdge[4]-1.06507042020704*nuVtSqSum[1]*fSkin[3]-1.06507042020704*nuVtSqSum[1]*fEdge[3]-1.190784930203603*fSkin[2]*nuVtSqSum[2]-1.190784930203603*fEdge[2]*nuVtSqSum[2]-0.9375*fSkin[0]*nuVtSqSum[2]+0.9375*fEdge[0]*nuVtSqSum[2]-0.8385254915624212*fSkin[1]*nuVtSqSum[1]+0.8385254915624212*fEdge[1]*nuVtSqSum[1]; 

  Ghat[0] = (-1.0*Gdiff[0]*rdv2)+0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = (-1.0*Gdiff[1]*rdv2)+0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = (-1.0*Gdiff[2]*rdv2)+0.4517539514526256*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += -0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += -0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += (-1.369306393762915*nuVtSqSum[1]*fSkin[7]*rdvSq4)+1.060660171779821*nuVtSqSum[2]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum[0]*fSkin[5]*rdvSq4-0.6123724356957944*nuVtSqSum[2]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4+1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[3] += (-1.224744871391589*nuVtSqSum[2]*fSkin[7]*rdvSq4)-1.369306393762915*nuVtSqSum[0]*fSkin[7]*rdvSq4+0.9486832980505138*nuVtSqSum[1]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum[1]*fSkin[5]*rdvSq4-0.5477225575051661*nuVtSqSum[1]*fSkin[4]*rdvSq4+0.9486832980505137*nuVtSqSum[2]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum[2]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4+1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[4] += -0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += 5.303300858899106*nuVtSqSum[1]*fSkin[7]*rdvSq4-4.107919181288745*nuVtSqSum[2]*fSkin[6]*rdvSq4+5.303300858899105*nuVtSqSum[0]*fSkin[5]*rdvSq4+2.371708245126284*nuVtSqSum[2]*fSkin[4]*rdvSq4-4.107919181288745*nuVtSqSum[1]*fSkin[3]*rdvSq4-4.107919181288745*nuVtSqSum[0]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[1]*nuVtSqSum[1]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum[0]*rdvSq4+4.743416490252569*Gdiff2[0]*rdvSq4-1.58113883008419*Ghat[0]*rdv2; 
  out[6] += (-1.224744871391589*nuVtSqSum[1]*fSkin[7]*rdvSq4)+0.6776309271789384*nuVtSqSum[2]*fSkin[6]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[6]*rdvSq4-1.369306393762915*nuVtSqSum[2]*fSkin[5]*rdvSq4-0.3912303982179757*nuVtSqSum[2]*fSkin[4]*rdvSq4-0.6123724356957944*nuVtSqSum[0]*fSkin[4]*rdvSq4+0.9486832980505138*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*fSkin[2]*nuVtSqSum[2]*rdvSq4-0.6123724356957944*fSkin[0]*nuVtSqSum[2]*rdvSq4+1.224744871391589*Gdiff2[2]*rdvSq4-0.5477225575051661*fSkin[1]*nuVtSqSum[1]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 4.743416490252569*nuVtSqSum[2]*fSkin[7]*rdvSq4+5.303300858899105*nuVtSqSum[0]*fSkin[7]*rdvSq4-3.674234614174766*nuVtSqSum[1]*fSkin[6]*rdvSq4+5.303300858899106*nuVtSqSum[1]*fSkin[5]*rdvSq4+2.121320343559642*nuVtSqSum[1]*fSkin[4]*rdvSq4-3.674234614174767*nuVtSqSum[2]*fSkin[3]*rdvSq4-4.107919181288746*nuVtSqSum[0]*fSkin[3]*rdvSq4+2.121320343559642*fSkin[1]*nuVtSqSum[2]*rdvSq4-4.107919181288746*nuVtSqSum[1]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum[1]*rdvSq4+2.371708245126284*nuVtSqSum[0]*fSkin[1]*rdvSq4+4.743416490252569*Gdiff2[1]*rdvSq4-1.58113883008419*Ghat[1]*rdv2; 

  } else { 

  alphaDrSurf[0] = 0.5*(nuSum[0]*(2.0*w[1]+dxv[1])-2.0*sumNuUx[0]); 
  alphaDrSurf[1] = 0.5*(2.0*nuSum[1]*w[1]-2.0*sumNuUx[1]+dxv[1]*nuSum[1]); 
  alphaDrSurf[2] = -0.5*(2.0*sumNuUx[2]+((-2.0*w[1])-1.0*dxv[1])*nuSum[2]); 

  if (0.6324555320336759*alphaDrSurf[2]-0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(1, fEdge); 
  } else { 
    fUpwindQuad[0] = ser_1x1v_p2_surfvx_quad_0(-1, fSkin); 
  } 
  if (0.7071067811865475*alphaDrSurf[0]-0.7905694150420947*alphaDrSurf[2] < 0) { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(1, fEdge); 
  } else { 
    fUpwindQuad[1] = ser_1x1v_p2_surfvx_quad_1(-1, fSkin); 
  } 
  if (0.6324555320336759*alphaDrSurf[2]+0.9486832980505137*alphaDrSurf[1]+0.7071067811865475*alphaDrSurf[0] < 0) { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(1, fEdge); 
  } else { 
    fUpwindQuad[2] = ser_1x1v_p2_surfvx_quad_2(-1, fSkin); 
  } 

  fUpwind[0] = 0.392837100659193*fUpwindQuad[2]+0.6285393610547091*fUpwindQuad[1]+0.392837100659193*fUpwindQuad[0]; 
  fUpwind[1] = 0.5270462766947298*fUpwindQuad[2]-0.5270462766947298*fUpwindQuad[0]; 
  fUpwind[2] = 0.3513641844631533*fUpwindQuad[2]-0.7027283689263066*fUpwindQuad[1]+0.3513641844631533*fUpwindQuad[0]; 

  Gdiff2[0] = 0.2445699350390395*nuVtSqSum[1]*fSkin[7]+0.2445699350390395*nuVtSqSum[1]*fEdge[7]-0.3518228202874282*nuVtSqSum[2]*fSkin[6]+0.3518228202874282*nuVtSqSum[2]*fEdge[6]+0.2445699350390395*nuVtSqSum[0]*fSkin[5]+0.2445699350390395*nuVtSqSum[0]*fEdge[5]+0.25*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[2]*fEdge[4]-0.3518228202874282*nuVtSqSum[1]*fSkin[3]+0.3518228202874282*nuVtSqSum[1]*fEdge[3]-0.3518228202874282*nuVtSqSum[0]*fSkin[2]+0.3518228202874282*nuVtSqSum[0]*fEdge[2]+0.25*fSkin[1]*nuVtSqSum[1]+0.25*fEdge[1]*nuVtSqSum[1]+0.25*fSkin[0]*nuVtSqSum[0]+0.25*fEdge[0]*nuVtSqSum[0]; 
  Gdiff2[1] = 0.21875*nuVtSqSum[2]*fSkin[7]+0.2445699350390395*nuVtSqSum[0]*fSkin[7]+0.21875*nuVtSqSum[2]*fEdge[7]+0.2445699350390395*nuVtSqSum[0]*fEdge[7]-0.3146798968793527*nuVtSqSum[1]*fSkin[6]+0.3146798968793527*nuVtSqSum[1]*fEdge[6]+0.2445699350390395*nuVtSqSum[1]*fSkin[5]+0.2445699350390395*nuVtSqSum[1]*fEdge[5]+0.223606797749979*nuVtSqSum[1]*fSkin[4]+0.223606797749979*nuVtSqSum[1]*fEdge[4]-0.3146798968793526*nuVtSqSum[2]*fSkin[3]-0.3518228202874282*nuVtSqSum[0]*fSkin[3]+0.3146798968793526*nuVtSqSum[2]*fEdge[3]+0.3518228202874282*nuVtSqSum[0]*fEdge[3]+0.223606797749979*fSkin[1]*nuVtSqSum[2]+0.223606797749979*fEdge[1]*nuVtSqSum[2]-0.3518228202874282*nuVtSqSum[1]*fSkin[2]+0.3518228202874282*nuVtSqSum[1]*fEdge[2]+0.25*fSkin[0]*nuVtSqSum[1]+0.25*fEdge[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fSkin[1]+0.25*nuVtSqSum[0]*fEdge[1]; 
  Gdiff2[2] = 0.21875*nuVtSqSum[1]*fSkin[7]+0.21875*nuVtSqSum[1]*fEdge[7]-0.2247713549138233*nuVtSqSum[2]*fSkin[6]-0.3518228202874282*nuVtSqSum[0]*fSkin[6]+0.2247713549138233*nuVtSqSum[2]*fEdge[6]+0.3518228202874282*nuVtSqSum[0]*fEdge[6]+0.2445699350390395*nuVtSqSum[2]*fSkin[5]+0.2445699350390395*nuVtSqSum[2]*fEdge[5]+0.159719141249985*nuVtSqSum[2]*fSkin[4]+0.25*nuVtSqSum[0]*fSkin[4]+0.159719141249985*nuVtSqSum[2]*fEdge[4]+0.25*nuVtSqSum[0]*fEdge[4]-0.3146798968793526*nuVtSqSum[1]*fSkin[3]+0.3146798968793526*nuVtSqSum[1]*fEdge[3]-0.3518228202874282*fSkin[2]*nuVtSqSum[2]+0.3518228202874282*fEdge[2]*nuVtSqSum[2]+0.25*fSkin[0]*nuVtSqSum[2]+0.25*fEdge[0]*nuVtSqSum[2]+0.223606797749979*fSkin[1]*nuVtSqSum[1]+0.223606797749979*fEdge[1]*nuVtSqSum[1]; 

  Gdiff[0] = 0.6708203932499369*nuVtSqSum[1]*fSkin[7]-0.6708203932499369*nuVtSqSum[1]*fEdge[7]-1.190784930203603*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[2]*fEdge[6]+0.6708203932499369*nuVtSqSum[0]*fSkin[5]-0.6708203932499369*nuVtSqSum[0]*fEdge[5]+0.9375*nuVtSqSum[2]*fSkin[4]-0.9375*nuVtSqSum[2]*fEdge[4]-1.190784930203603*nuVtSqSum[1]*fSkin[3]-1.190784930203603*nuVtSqSum[1]*fEdge[3]-1.190784930203603*nuVtSqSum[0]*fSkin[2]-1.190784930203603*nuVtSqSum[0]*fEdge[2]+0.9375*fSkin[1]*nuVtSqSum[1]-0.9375*fEdge[1]*nuVtSqSum[1]+0.9375*fSkin[0]*nuVtSqSum[0]-0.9375*fEdge[0]*nuVtSqSum[0]; 
  Gdiff[1] = 0.5999999999999999*nuVtSqSum[2]*fSkin[7]+0.6708203932499369*nuVtSqSum[0]*fSkin[7]-0.5999999999999999*nuVtSqSum[2]*fEdge[7]-0.6708203932499369*nuVtSqSum[0]*fEdge[7]-1.06507042020704*nuVtSqSum[1]*fSkin[6]-1.06507042020704*nuVtSqSum[1]*fEdge[6]+0.6708203932499369*nuVtSqSum[1]*fSkin[5]-0.6708203932499369*nuVtSqSum[1]*fEdge[5]+0.8385254915624212*nuVtSqSum[1]*fSkin[4]-0.8385254915624212*nuVtSqSum[1]*fEdge[4]-1.06507042020704*nuVtSqSum[2]*fSkin[3]-1.190784930203603*nuVtSqSum[0]*fSkin[3]-1.06507042020704*nuVtSqSum[2]*fEdge[3]-1.190784930203603*nuVtSqSum[0]*fEdge[3]+0.8385254915624212*fSkin[1]*nuVtSqSum[2]-0.8385254915624212*fEdge[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fSkin[2]-1.190784930203603*nuVtSqSum[1]*fEdge[2]+0.9375*fSkin[0]*nuVtSqSum[1]-0.9375*fEdge[0]*nuVtSqSum[1]+0.9375*nuVtSqSum[0]*fSkin[1]-0.9375*nuVtSqSum[0]*fEdge[1]; 
  Gdiff[2] = 0.5999999999999999*nuVtSqSum[1]*fSkin[7]-0.5999999999999999*nuVtSqSum[1]*fEdge[7]-0.7607645858621712*nuVtSqSum[2]*fSkin[6]-1.190784930203603*nuVtSqSum[0]*fSkin[6]-0.7607645858621712*nuVtSqSum[2]*fEdge[6]-1.190784930203603*nuVtSqSum[0]*fEdge[6]+0.6708203932499369*nuVtSqSum[2]*fSkin[5]-0.6708203932499369*nuVtSqSum[2]*fEdge[5]+0.5989467796874438*nuVtSqSum[2]*fSkin[4]+0.9375*nuVtSqSum[0]*fSkin[4]-0.5989467796874438*nuVtSqSum[2]*fEdge[4]-0.9375*nuVtSqSum[0]*fEdge[4]-1.06507042020704*nuVtSqSum[1]*fSkin[3]-1.06507042020704*nuVtSqSum[1]*fEdge[3]-1.190784930203603*fSkin[2]*nuVtSqSum[2]-1.190784930203603*fEdge[2]*nuVtSqSum[2]+0.9375*fSkin[0]*nuVtSqSum[2]-0.9375*fEdge[0]*nuVtSqSum[2]+0.8385254915624212*fSkin[1]*nuVtSqSum[1]-0.8385254915624212*fEdge[1]*nuVtSqSum[1]; 

  Ghat[0] = Gdiff[0]*rdv2+0.7071067811865475*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[1]*fUpwind[1]+0.7071067811865475*alphaDrSurf[0]*fUpwind[0]; 
  Ghat[1] = Gdiff[1]*rdv2+0.6324555320336759*alphaDrSurf[1]*fUpwind[2]+0.6324555320336759*fUpwind[1]*alphaDrSurf[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[1]+0.7071067811865475*fUpwind[0]*alphaDrSurf[1]; 
  Ghat[2] = Gdiff[2]*rdv2+0.4517539514526256*alphaDrSurf[2]*fUpwind[2]+0.7071067811865475*alphaDrSurf[0]*fUpwind[2]+0.7071067811865475*fUpwind[0]*alphaDrSurf[2]+0.6324555320336759*alphaDrSurf[1]*fUpwind[1]; 

  out[0] += 0.7071067811865475*Ghat[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat[1]*rdv2; 
  out[2] += 1.369306393762915*nuVtSqSum[1]*fSkin[7]*rdvSq4+1.060660171779821*nuVtSqSum[2]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[5]*rdvSq4+0.6123724356957944*nuVtSqSum[2]*fSkin[4]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[1]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[0]*rdvSq4-1.224744871391589*Gdiff2[0]*rdvSq4-1.224744871391589*Ghat[0]*rdv2; 
  out[3] += 1.224744871391589*nuVtSqSum[2]*fSkin[7]*rdvSq4+1.369306393762915*nuVtSqSum[0]*fSkin[7]*rdvSq4+0.9486832980505138*nuVtSqSum[1]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum[1]*fSkin[5]*rdvSq4+0.5477225575051661*nuVtSqSum[1]*fSkin[4]*rdvSq4+0.9486832980505137*nuVtSqSum[2]*fSkin[3]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[3]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum[2]*rdvSq4+1.060660171779821*nuVtSqSum[1]*fSkin[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[1]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[1]*rdvSq4-1.224744871391589*Gdiff2[1]*rdvSq4-1.224744871391589*Ghat[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat[2]*rdv2; 
  out[5] += 5.303300858899106*nuVtSqSum[1]*fSkin[7]*rdvSq4+4.107919181288745*nuVtSqSum[2]*fSkin[6]*rdvSq4+5.303300858899105*nuVtSqSum[0]*fSkin[5]*rdvSq4+2.371708245126284*nuVtSqSum[2]*fSkin[4]*rdvSq4+4.107919181288745*nuVtSqSum[1]*fSkin[3]*rdvSq4+4.107919181288745*nuVtSqSum[0]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[1]*nuVtSqSum[1]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum[0]*rdvSq4+4.743416490252569*Gdiff2[0]*rdvSq4+1.58113883008419*Ghat[0]*rdv2; 
  out[6] += 1.224744871391589*nuVtSqSum[1]*fSkin[7]*rdvSq4+0.6776309271789384*nuVtSqSum[2]*fSkin[6]*rdvSq4+1.060660171779821*nuVtSqSum[0]*fSkin[6]*rdvSq4+1.369306393762915*nuVtSqSum[2]*fSkin[5]*rdvSq4+0.3912303982179757*nuVtSqSum[2]*fSkin[4]*rdvSq4+0.6123724356957944*nuVtSqSum[0]*fSkin[4]*rdvSq4+0.9486832980505138*nuVtSqSum[1]*fSkin[3]*rdvSq4+1.060660171779821*fSkin[2]*nuVtSqSum[2]*rdvSq4+0.6123724356957944*fSkin[0]*nuVtSqSum[2]*rdvSq4-1.224744871391589*Gdiff2[2]*rdvSq4+0.5477225575051661*fSkin[1]*nuVtSqSum[1]*rdvSq4-1.224744871391589*Ghat[2]*rdv2; 
  out[7] += 4.743416490252569*nuVtSqSum[2]*fSkin[7]*rdvSq4+5.303300858899105*nuVtSqSum[0]*fSkin[7]*rdvSq4+3.674234614174766*nuVtSqSum[1]*fSkin[6]*rdvSq4+5.303300858899106*nuVtSqSum[1]*fSkin[5]*rdvSq4+2.121320343559642*nuVtSqSum[1]*fSkin[4]*rdvSq4+3.674234614174767*nuVtSqSum[2]*fSkin[3]*rdvSq4+4.107919181288746*nuVtSqSum[0]*fSkin[3]*rdvSq4+2.121320343559642*fSkin[1]*nuVtSqSum[2]*rdvSq4+4.107919181288746*nuVtSqSum[1]*fSkin[2]*rdvSq4+2.371708245126284*fSkin[0]*nuVtSqSum[1]*rdvSq4+2.371708245126284*nuVtSqSum[0]*fSkin[1]*rdvSq4+4.743416490252569*Gdiff2[1]*rdvSq4+1.58113883008419*Ghat[1]*rdv2; 

  } 
} 
