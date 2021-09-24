#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x1v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvx_1x1v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[2]:         cell-center coordinates. 
  // dxv[2]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[3]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[3] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[1]-0.7071067811865475*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_l[1] = -1.0*sumNuUx[1]; 
  alphaDrSurf_l[2] = -1.0*sumNuUx[2]; 

  double alphaDrSurf_r[3] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[1]+0.7071067811865475*nuSum[0]*dxv[1]-1.0*sumNuUx[0]; 
  alphaDrSurf_r[1] = -1.0*sumNuUx[1]; 
  alphaDrSurf_r[2] = -1.0*sumNuUx[2]; 

  double fUpwindQuad_l[3] = {0.0};
  double fUpwindQuad_r[3] = {0.0};
  double fUpwind_l[3] = {0.0};;
  double fUpwind_r[3] = {0.0};
  double Ghat_l[3] = {0.0}; 
  double Ghat_r[3] = {0.0}; 
  double Gdiff_l[3] = {0.0}; 
  double Gdiff_r[3] = {0.0}; 
  double Gdiff2_l[3] = {0.0}; 
  double Gdiff2_r[3] = {0.0}; 

  if (0.6324555320336759*alphaDrSurf_l[2]-0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_1x1v_p2_surfvx_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x1v_p2_surfvx_quad_0(-1, fc); 
  } 
  if (0.7071067811865475*alphaDrSurf_l[0]-0.7905694150420947*alphaDrSurf_l[2] > 0) { 
    fUpwindQuad_l[1] = ser_1x1v_p2_surfvx_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x1v_p2_surfvx_quad_1(-1, fc); 
  } 
  if (0.6324555320336759*alphaDrSurf_l[2]+0.9486832980505137*alphaDrSurf_l[1]+0.7071067811865475*alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_1x1v_p2_surfvx_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_1x1v_p2_surfvx_quad_2(-1, fc); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]-0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_1x1v_p2_surfvx_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x1v_p2_surfvx_quad_0(-1, fr); 
  } 
  if (0.7071067811865475*alphaDrSurf_r[0]-0.7905694150420947*alphaDrSurf_r[2] > 0) { 
    fUpwindQuad_r[1] = ser_1x1v_p2_surfvx_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x1v_p2_surfvx_quad_1(-1, fr); 
  } 
  if (0.6324555320336759*alphaDrSurf_r[2]+0.9486832980505137*alphaDrSurf_r[1]+0.7071067811865475*alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_1x1v_p2_surfvx_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_1x1v_p2_surfvx_quad_2(-1, fr); 
  } 
  fUpwind_l[0] = 0.392837100659193*fUpwindQuad_l[2]+0.6285393610547091*fUpwindQuad_l[1]+0.392837100659193*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.5270462766947298*fUpwindQuad_l[2]-0.5270462766947298*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.3513641844631533*fUpwindQuad_l[2]-0.7027283689263066*fUpwindQuad_l[1]+0.3513641844631533*fUpwindQuad_l[0]; 

  fUpwind_r[0] = 0.392837100659193*fUpwindQuad_r[2]+0.6285393610547091*fUpwindQuad_r[1]+0.392837100659193*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.5270462766947298*fUpwindQuad_r[2]-0.5270462766947298*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.3513641844631533*fUpwindQuad_r[2]-0.7027283689263066*fUpwindQuad_r[1]+0.3513641844631533*fUpwindQuad_r[0]; 

  Gdiff2_l[0] = 0.2445699350390395*nuVtSqSum[1]*fl[7]+0.2445699350390395*nuVtSqSum[1]*fc[7]+0.3518228202874282*nuVtSqSum[2]*fl[6]-0.3518228202874282*nuVtSqSum[2]*fc[6]+0.2445699350390395*nuVtSqSum[0]*fl[5]+0.2445699350390395*nuVtSqSum[0]*fc[5]+0.25*nuVtSqSum[2]*fl[4]+0.25*nuVtSqSum[2]*fc[4]+0.3518228202874282*nuVtSqSum[1]*fl[3]-0.3518228202874282*nuVtSqSum[1]*fc[3]+0.3518228202874282*nuVtSqSum[0]*fl[2]-0.3518228202874282*nuVtSqSum[0]*fc[2]+0.25*fl[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fl[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.21875*nuVtSqSum[2]*fl[7]+0.2445699350390395*nuVtSqSum[0]*fl[7]+0.21875*nuVtSqSum[2]*fc[7]+0.2445699350390395*nuVtSqSum[0]*fc[7]+0.3146798968793526*nuVtSqSum[1]*fl[6]-0.3146798968793526*nuVtSqSum[1]*fc[6]+0.2445699350390395*nuVtSqSum[1]*fl[5]+0.2445699350390395*nuVtSqSum[1]*fc[5]+0.223606797749979*nuVtSqSum[1]*fl[4]+0.223606797749979*nuVtSqSum[1]*fc[4]+0.3146798968793526*nuVtSqSum[2]*fl[3]+0.3518228202874282*nuVtSqSum[0]*fl[3]-0.3146798968793526*nuVtSqSum[2]*fc[3]-0.3518228202874282*nuVtSqSum[0]*fc[3]+0.223606797749979*fl[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]+0.3518228202874282*nuVtSqSum[1]*fl[2]-0.3518228202874282*nuVtSqSum[1]*fc[2]+0.25*fl[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fl[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.21875*nuVtSqSum[1]*fl[7]+0.21875*nuVtSqSum[1]*fc[7]+0.2247713549138233*nuVtSqSum[2]*fl[6]+0.3518228202874282*nuVtSqSum[0]*fl[6]-0.2247713549138233*nuVtSqSum[2]*fc[6]-0.3518228202874282*nuVtSqSum[0]*fc[6]+0.2445699350390395*nuVtSqSum[2]*fl[5]+0.2445699350390395*nuVtSqSum[2]*fc[5]+0.159719141249985*nuVtSqSum[2]*fl[4]+0.25*nuVtSqSum[0]*fl[4]+0.159719141249985*nuVtSqSum[2]*fc[4]+0.25*nuVtSqSum[0]*fc[4]+0.3146798968793526*nuVtSqSum[1]*fl[3]-0.3146798968793526*nuVtSqSum[1]*fc[3]+0.3518228202874282*fl[2]*nuVtSqSum[2]-0.3518228202874282*fc[2]*nuVtSqSum[2]+0.25*fl[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fl[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 

  Gdiff2_r[0] = 0.2445699350390395*nuVtSqSum[1]*fr[7]+0.2445699350390395*nuVtSqSum[1]*fc[7]-0.3518228202874282*nuVtSqSum[2]*fr[6]+0.3518228202874282*nuVtSqSum[2]*fc[6]+0.2445699350390395*nuVtSqSum[0]*fr[5]+0.2445699350390395*nuVtSqSum[0]*fc[5]+0.25*nuVtSqSum[2]*fr[4]+0.25*nuVtSqSum[2]*fc[4]-0.3518228202874282*nuVtSqSum[1]*fr[3]+0.3518228202874282*nuVtSqSum[1]*fc[3]-0.3518228202874282*nuVtSqSum[0]*fr[2]+0.3518228202874282*nuVtSqSum[0]*fc[2]+0.25*fr[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fr[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = 0.21875*nuVtSqSum[2]*fr[7]+0.2445699350390395*nuVtSqSum[0]*fr[7]+0.21875*nuVtSqSum[2]*fc[7]+0.2445699350390395*nuVtSqSum[0]*fc[7]-0.3146798968793526*nuVtSqSum[1]*fr[6]+0.3146798968793526*nuVtSqSum[1]*fc[6]+0.2445699350390395*nuVtSqSum[1]*fr[5]+0.2445699350390395*nuVtSqSum[1]*fc[5]+0.223606797749979*nuVtSqSum[1]*fr[4]+0.223606797749979*nuVtSqSum[1]*fc[4]-0.3146798968793526*nuVtSqSum[2]*fr[3]-0.3518228202874282*nuVtSqSum[0]*fr[3]+0.3146798968793526*nuVtSqSum[2]*fc[3]+0.3518228202874282*nuVtSqSum[0]*fc[3]+0.223606797749979*fr[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]-0.3518228202874282*nuVtSqSum[1]*fr[2]+0.3518228202874282*nuVtSqSum[1]*fc[2]+0.25*fr[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fr[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = 0.21875*nuVtSqSum[1]*fr[7]+0.21875*nuVtSqSum[1]*fc[7]-0.2247713549138233*nuVtSqSum[2]*fr[6]-0.3518228202874282*nuVtSqSum[0]*fr[6]+0.2247713549138233*nuVtSqSum[2]*fc[6]+0.3518228202874282*nuVtSqSum[0]*fc[6]+0.2445699350390395*nuVtSqSum[2]*fr[5]+0.2445699350390395*nuVtSqSum[2]*fc[5]+0.159719141249985*nuVtSqSum[2]*fr[4]+0.25*nuVtSqSum[0]*fr[4]+0.159719141249985*nuVtSqSum[2]*fc[4]+0.25*nuVtSqSum[0]*fc[4]-0.3146798968793526*nuVtSqSum[1]*fr[3]+0.3146798968793526*nuVtSqSum[1]*fc[3]-0.3518228202874282*fr[2]*nuVtSqSum[2]+0.3518228202874282*fc[2]*nuVtSqSum[2]+0.25*fr[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fr[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 

  Gdiff_l[0] = (-0.6708203932499369*nuVtSqSum[1]*fl[7])+0.6708203932499369*nuVtSqSum[1]*fc[7]-1.190784930203603*nuVtSqSum[2]*fl[6]-1.190784930203603*nuVtSqSum[2]*fc[6]-0.6708203932499369*nuVtSqSum[0]*fl[5]+0.6708203932499369*nuVtSqSum[0]*fc[5]-0.9375*nuVtSqSum[2]*fl[4]+0.9375*nuVtSqSum[2]*fc[4]-1.190784930203603*nuVtSqSum[1]*fl[3]-1.190784930203603*nuVtSqSum[1]*fc[3]-1.190784930203603*nuVtSqSum[0]*fl[2]-1.190784930203603*nuVtSqSum[0]*fc[2]-0.9375*fl[1]*nuVtSqSum[1]+0.9375*fc[1]*nuVtSqSum[1]-0.9375*fl[0]*nuVtSqSum[0]+0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.5999999999999999*nuVtSqSum[2]*fl[7])-0.6708203932499369*nuVtSqSum[0]*fl[7]+0.5999999999999999*nuVtSqSum[2]*fc[7]+0.6708203932499369*nuVtSqSum[0]*fc[7]-1.06507042020704*nuVtSqSum[1]*fl[6]-1.06507042020704*nuVtSqSum[1]*fc[6]-0.6708203932499369*nuVtSqSum[1]*fl[5]+0.6708203932499369*nuVtSqSum[1]*fc[5]-0.8385254915624212*nuVtSqSum[1]*fl[4]+0.8385254915624212*nuVtSqSum[1]*fc[4]-1.06507042020704*nuVtSqSum[2]*fl[3]-1.190784930203603*nuVtSqSum[0]*fl[3]-1.06507042020704*nuVtSqSum[2]*fc[3]-1.190784930203603*nuVtSqSum[0]*fc[3]-0.8385254915624212*fl[1]*nuVtSqSum[2]+0.8385254915624212*fc[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fl[2]-1.190784930203603*nuVtSqSum[1]*fc[2]-0.9375*fl[0]*nuVtSqSum[1]+0.9375*fc[0]*nuVtSqSum[1]-0.9375*nuVtSqSum[0]*fl[1]+0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.5999999999999999*nuVtSqSum[1]*fl[7])+0.5999999999999999*nuVtSqSum[1]*fc[7]-0.7607645858621712*nuVtSqSum[2]*fl[6]-1.190784930203603*nuVtSqSum[0]*fl[6]-0.7607645858621712*nuVtSqSum[2]*fc[6]-1.190784930203603*nuVtSqSum[0]*fc[6]-0.6708203932499369*nuVtSqSum[2]*fl[5]+0.6708203932499369*nuVtSqSum[2]*fc[5]-0.5989467796874438*nuVtSqSum[2]*fl[4]-0.9375*nuVtSqSum[0]*fl[4]+0.5989467796874438*nuVtSqSum[2]*fc[4]+0.9375*nuVtSqSum[0]*fc[4]-1.06507042020704*nuVtSqSum[1]*fl[3]-1.06507042020704*nuVtSqSum[1]*fc[3]-1.190784930203603*fl[2]*nuVtSqSum[2]-1.190784930203603*fc[2]*nuVtSqSum[2]-0.9375*fl[0]*nuVtSqSum[2]+0.9375*fc[0]*nuVtSqSum[2]-0.8385254915624212*fl[1]*nuVtSqSum[1]+0.8385254915624212*fc[1]*nuVtSqSum[1]; 

  Gdiff_r[0] = 0.6708203932499369*nuVtSqSum[1]*fr[7]-0.6708203932499369*nuVtSqSum[1]*fc[7]-1.190784930203603*nuVtSqSum[2]*fr[6]-1.190784930203603*nuVtSqSum[2]*fc[6]+0.6708203932499369*nuVtSqSum[0]*fr[5]-0.6708203932499369*nuVtSqSum[0]*fc[5]+0.9375*nuVtSqSum[2]*fr[4]-0.9375*nuVtSqSum[2]*fc[4]-1.190784930203603*nuVtSqSum[1]*fr[3]-1.190784930203603*nuVtSqSum[1]*fc[3]-1.190784930203603*nuVtSqSum[0]*fr[2]-1.190784930203603*nuVtSqSum[0]*fc[2]+0.9375*fr[1]*nuVtSqSum[1]-0.9375*fc[1]*nuVtSqSum[1]+0.9375*fr[0]*nuVtSqSum[0]-0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = 0.5999999999999999*nuVtSqSum[2]*fr[7]+0.6708203932499369*nuVtSqSum[0]*fr[7]-0.5999999999999999*nuVtSqSum[2]*fc[7]-0.6708203932499369*nuVtSqSum[0]*fc[7]-1.06507042020704*nuVtSqSum[1]*fr[6]-1.06507042020704*nuVtSqSum[1]*fc[6]+0.6708203932499369*nuVtSqSum[1]*fr[5]-0.6708203932499369*nuVtSqSum[1]*fc[5]+0.8385254915624212*nuVtSqSum[1]*fr[4]-0.8385254915624212*nuVtSqSum[1]*fc[4]-1.06507042020704*nuVtSqSum[2]*fr[3]-1.190784930203603*nuVtSqSum[0]*fr[3]-1.06507042020704*nuVtSqSum[2]*fc[3]-1.190784930203603*nuVtSqSum[0]*fc[3]+0.8385254915624212*fr[1]*nuVtSqSum[2]-0.8385254915624212*fc[1]*nuVtSqSum[2]-1.190784930203603*nuVtSqSum[1]*fr[2]-1.190784930203603*nuVtSqSum[1]*fc[2]+0.9375*fr[0]*nuVtSqSum[1]-0.9375*fc[0]*nuVtSqSum[1]+0.9375*nuVtSqSum[0]*fr[1]-0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = 0.5999999999999999*nuVtSqSum[1]*fr[7]-0.5999999999999999*nuVtSqSum[1]*fc[7]-0.7607645858621712*nuVtSqSum[2]*fr[6]-1.190784930203603*nuVtSqSum[0]*fr[6]-0.7607645858621712*nuVtSqSum[2]*fc[6]-1.190784930203603*nuVtSqSum[0]*fc[6]+0.6708203932499369*nuVtSqSum[2]*fr[5]-0.6708203932499369*nuVtSqSum[2]*fc[5]+0.5989467796874438*nuVtSqSum[2]*fr[4]+0.9375*nuVtSqSum[0]*fr[4]-0.5989467796874438*nuVtSqSum[2]*fc[4]-0.9375*nuVtSqSum[0]*fc[4]-1.06507042020704*nuVtSqSum[1]*fr[3]-1.06507042020704*nuVtSqSum[1]*fc[3]-1.190784930203603*fr[2]*nuVtSqSum[2]-1.190784930203603*fc[2]*nuVtSqSum[2]+0.9375*fr[0]*nuVtSqSum[2]-0.9375*fc[0]*nuVtSqSum[2]+0.8385254915624212*fr[1]*nuVtSqSum[1]-0.8385254915624212*fc[1]*nuVtSqSum[1]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.7071067811865475*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[1]*fUpwind_l[1]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[2]+0.6324555320336759*fUpwind_l[1]*alphaDrSurf_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[1]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.4517539514526256*alphaDrSurf_l[2]*fUpwind_l[2]+0.7071067811865475*alphaDrSurf_l[0]*fUpwind_l[2]+0.7071067811865475*fUpwind_l[0]*alphaDrSurf_l[2]+0.6324555320336759*alphaDrSurf_l[1]*fUpwind_l[1]; 

  Ghat_r[0] = Gdiff_r[0]*rdv2+0.7071067811865475*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[1]*fUpwind_r[1]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = Gdiff_r[1]*rdv2+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[2]+0.6324555320336759*fUpwind_r[1]*alphaDrSurf_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[1]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = Gdiff_r[2]*rdv2+0.4517539514526256*alphaDrSurf_r[2]*fUpwind_r[2]+0.7071067811865475*alphaDrSurf_r[0]*fUpwind_r[2]+0.7071067811865475*fUpwind_r[0]*alphaDrSurf_r[2]+0.6324555320336759*alphaDrSurf_r[1]*fUpwind_r[1]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[3] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[5] += 4.743416490252569*Gdiff2_r[0]*rdvSq4+4.743416490252569*Gdiff2_l[0]*rdvSq4-1.58113883008419*Ghat_r[0]*rdv2+1.58113883008419*Ghat_l[0]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 4.743416490252569*Gdiff2_r[1]*rdvSq4+4.743416490252569*Gdiff2_l[1]*rdvSq4-1.58113883008419*Ghat_r[1]*rdv2+1.58113883008419*Ghat_l[1]*rdv2; 
} 
