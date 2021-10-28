#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_tensor_1x2v_p2_surfvy_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvy_1x2v_tensor_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[2]; 
  double rdvSq4 = 4.0/(dxv[2]*dxv[2]); 

  const double *sumNuUy = &nuUSum[3]; 

  double alphaDrSurf_l[9] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[2]-0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[2]-0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_l[4] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*sumNuUy[2]-0.7071067811865475*dxv[2]*nuSum[2]; 

  double alphaDrSurf_r[9] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[2]+0.7071067811865475*nuSum[0]*dxv[2]-1.414213562373095*sumNuUy[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[2]+0.7071067811865475*nuSum[1]*dxv[2]-1.414213562373095*sumNuUy[1]; 
  alphaDrSurf_r[4] = 1.414213562373095*nuSum[2]*w[2]-1.414213562373095*sumNuUy[2]+0.7071067811865475*dxv[2]*nuSum[2]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[9] = {0.0};;
  double fUpwind_r[9] = {0.0};
  double Ghat_l[9] = {0.0}; 
  double Ghat_r[9] = {0.0}; 
  double Gdiff_l[9] = {0.0}; 
  double Gdiff_r[9] = {0.0}; 
  double Gdiff2_l[9] = {0.0}; 
  double Gdiff2_r[9] = {0.0}; 

  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = tensor_1x2v_p2_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = tensor_1x2v_p2_surfvy_quad_0(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[1] = tensor_1x2v_p2_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = tensor_1x2v_p2_surfvy_quad_1(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = tensor_1x2v_p2_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = tensor_1x2v_p2_surfvy_quad_2(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = tensor_1x2v_p2_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = tensor_1x2v_p2_surfvy_quad_3(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[4] = tensor_1x2v_p2_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = tensor_1x2v_p2_surfvy_quad_4(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = tensor_1x2v_p2_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = tensor_1x2v_p2_surfvy_quad_5(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = tensor_1x2v_p2_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = tensor_1x2v_p2_surfvy_quad_6(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[7] = tensor_1x2v_p2_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = tensor_1x2v_p2_surfvy_quad_7(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = tensor_1x2v_p2_surfvy_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = tensor_1x2v_p2_surfvy_quad_8(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = tensor_1x2v_p2_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = tensor_1x2v_p2_surfvy_quad_0(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[1] = tensor_1x2v_p2_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = tensor_1x2v_p2_surfvy_quad_1(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = tensor_1x2v_p2_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = tensor_1x2v_p2_surfvy_quad_2(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = tensor_1x2v_p2_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = tensor_1x2v_p2_surfvy_quad_3(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[4] = tensor_1x2v_p2_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = tensor_1x2v_p2_surfvy_quad_4(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = tensor_1x2v_p2_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = tensor_1x2v_p2_surfvy_quad_5(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = tensor_1x2v_p2_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = tensor_1x2v_p2_surfvy_quad_6(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[7] = tensor_1x2v_p2_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = tensor_1x2v_p2_surfvy_quad_7(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = tensor_1x2v_p2_surfvy_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = tensor_1x2v_p2_surfvy_quad_8(-1, fr); 
  } 
  fUpwind_l[0] = 0.154320987654321*fUpwindQuad_l[8]+0.2469135802469136*fUpwindQuad_l[7]+0.154320987654321*fUpwindQuad_l[6]+0.2469135802469136*fUpwindQuad_l[5]+0.3950617283950617*fUpwindQuad_l[4]+0.2469135802469136*fUpwindQuad_l[3]+0.154320987654321*fUpwindQuad_l[2]+0.2469135802469136*fUpwindQuad_l[1]+0.154320987654321*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.2070433312499806*fUpwindQuad_l[8]-0.2070433312499806*fUpwindQuad_l[6]+0.3312693299999688*fUpwindQuad_l[5]-0.3312693299999688*fUpwindQuad_l[3]+0.2070433312499806*fUpwindQuad_l[2]-0.2070433312499806*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.2070433312499806*fUpwindQuad_l[8]+0.3312693299999688*fUpwindQuad_l[7]+0.2070433312499806*fUpwindQuad_l[6]-0.2070433312499806*fUpwindQuad_l[2]-0.3312693299999688*fUpwindQuad_l[1]-0.2070433312499806*fUpwindQuad_l[0]; 
  fUpwind_l[3] = 0.2777777777777778*fUpwindQuad_l[8]-0.2777777777777778*fUpwindQuad_l[6]-0.2777777777777778*fUpwindQuad_l[2]+0.2777777777777778*fUpwindQuad_l[0]; 
  fUpwind_l[4] = 0.138028887499987*fUpwindQuad_l[8]-0.2760577749999741*fUpwindQuad_l[7]+0.138028887499987*fUpwindQuad_l[6]+0.2208462199999792*fUpwindQuad_l[5]-0.4416924399999584*fUpwindQuad_l[4]+0.2208462199999792*fUpwindQuad_l[3]+0.138028887499987*fUpwindQuad_l[2]-0.2760577749999741*fUpwindQuad_l[1]+0.138028887499987*fUpwindQuad_l[0]; 
  fUpwind_l[5] = 0.138028887499987*fUpwindQuad_l[8]+0.2208462199999792*fUpwindQuad_l[7]+0.138028887499987*fUpwindQuad_l[6]-0.2760577749999741*fUpwindQuad_l[5]-0.4416924399999584*fUpwindQuad_l[4]-0.2760577749999741*fUpwindQuad_l[3]+0.138028887499987*fUpwindQuad_l[2]+0.2208462199999792*fUpwindQuad_l[1]+0.138028887499987*fUpwindQuad_l[0]; 
  fUpwind_l[6] = 0.1851851851851853*fUpwindQuad_l[8]-0.3703703703703705*fUpwindQuad_l[7]+0.1851851851851853*fUpwindQuad_l[6]-0.1851851851851853*fUpwindQuad_l[2]+0.3703703703703705*fUpwindQuad_l[1]-0.1851851851851853*fUpwindQuad_l[0]; 
  fUpwind_l[7] = 0.1851851851851853*fUpwindQuad_l[8]-0.1851851851851853*fUpwindQuad_l[6]-0.3703703703703705*fUpwindQuad_l[5]+0.3703703703703705*fUpwindQuad_l[3]+0.1851851851851853*fUpwindQuad_l[2]-0.1851851851851853*fUpwindQuad_l[0]; 
  fUpwind_l[8] = 0.1234567901234568*fUpwindQuad_l[8]-0.2469135802469136*fUpwindQuad_l[7]+0.1234567901234568*fUpwindQuad_l[6]-0.2469135802469136*fUpwindQuad_l[5]+0.4938271604938271*fUpwindQuad_l[4]-0.2469135802469136*fUpwindQuad_l[3]+0.1234567901234568*fUpwindQuad_l[2]-0.2469135802469136*fUpwindQuad_l[1]+0.1234567901234568*fUpwindQuad_l[0]; 

  fUpwind_r[0] = 0.154320987654321*fUpwindQuad_r[8]+0.2469135802469136*fUpwindQuad_r[7]+0.154320987654321*fUpwindQuad_r[6]+0.2469135802469136*fUpwindQuad_r[5]+0.3950617283950617*fUpwindQuad_r[4]+0.2469135802469136*fUpwindQuad_r[3]+0.154320987654321*fUpwindQuad_r[2]+0.2469135802469136*fUpwindQuad_r[1]+0.154320987654321*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.2070433312499806*fUpwindQuad_r[8]-0.2070433312499806*fUpwindQuad_r[6]+0.3312693299999688*fUpwindQuad_r[5]-0.3312693299999688*fUpwindQuad_r[3]+0.2070433312499806*fUpwindQuad_r[2]-0.2070433312499806*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.2070433312499806*fUpwindQuad_r[8]+0.3312693299999688*fUpwindQuad_r[7]+0.2070433312499806*fUpwindQuad_r[6]-0.2070433312499806*fUpwindQuad_r[2]-0.3312693299999688*fUpwindQuad_r[1]-0.2070433312499806*fUpwindQuad_r[0]; 
  fUpwind_r[3] = 0.2777777777777778*fUpwindQuad_r[8]-0.2777777777777778*fUpwindQuad_r[6]-0.2777777777777778*fUpwindQuad_r[2]+0.2777777777777778*fUpwindQuad_r[0]; 
  fUpwind_r[4] = 0.138028887499987*fUpwindQuad_r[8]-0.2760577749999741*fUpwindQuad_r[7]+0.138028887499987*fUpwindQuad_r[6]+0.2208462199999792*fUpwindQuad_r[5]-0.4416924399999584*fUpwindQuad_r[4]+0.2208462199999792*fUpwindQuad_r[3]+0.138028887499987*fUpwindQuad_r[2]-0.2760577749999741*fUpwindQuad_r[1]+0.138028887499987*fUpwindQuad_r[0]; 
  fUpwind_r[5] = 0.138028887499987*fUpwindQuad_r[8]+0.2208462199999792*fUpwindQuad_r[7]+0.138028887499987*fUpwindQuad_r[6]-0.2760577749999741*fUpwindQuad_r[5]-0.4416924399999584*fUpwindQuad_r[4]-0.2760577749999741*fUpwindQuad_r[3]+0.138028887499987*fUpwindQuad_r[2]+0.2208462199999792*fUpwindQuad_r[1]+0.138028887499987*fUpwindQuad_r[0]; 
  fUpwind_r[6] = 0.1851851851851853*fUpwindQuad_r[8]-0.3703703703703705*fUpwindQuad_r[7]+0.1851851851851853*fUpwindQuad_r[6]-0.1851851851851853*fUpwindQuad_r[2]+0.3703703703703705*fUpwindQuad_r[1]-0.1851851851851853*fUpwindQuad_r[0]; 
  fUpwind_r[7] = 0.1851851851851853*fUpwindQuad_r[8]-0.1851851851851853*fUpwindQuad_r[6]-0.3703703703703705*fUpwindQuad_r[5]+0.3703703703703705*fUpwindQuad_r[3]+0.1851851851851853*fUpwindQuad_r[2]-0.1851851851851853*fUpwindQuad_r[0]; 
  fUpwind_r[8] = 0.1234567901234568*fUpwindQuad_r[8]-0.2469135802469136*fUpwindQuad_r[7]+0.1234567901234568*fUpwindQuad_r[6]-0.2469135802469136*fUpwindQuad_r[5]+0.4938271604938271*fUpwindQuad_r[4]-0.2469135802469136*fUpwindQuad_r[3]+0.1234567901234568*fUpwindQuad_r[2]-0.2469135802469136*fUpwindQuad_r[1]+0.1234567901234568*fUpwindQuad_r[0]; 

  Gdiff2_l[0] = 0.2445699350390395*nuVtSqSum[2]*fl[21]+0.2445699350390395*nuVtSqSum[2]*fc[21]+0.2445699350390395*nuVtSqSum[1]*fl[15]+0.2445699350390395*nuVtSqSum[1]*fc[15]+0.3518228202874282*nuVtSqSum[2]*fl[13]-0.3518228202874282*nuVtSqSum[2]*fc[13]+0.2445699350390395*nuVtSqSum[0]*fl[9]+0.2445699350390395*nuVtSqSum[0]*fc[9]+0.25*nuVtSqSum[2]*fl[7]+0.25*nuVtSqSum[2]*fc[7]+0.3518228202874282*nuVtSqSum[1]*fl[5]-0.3518228202874282*nuVtSqSum[1]*fc[5]+0.3518228202874282*nuVtSqSum[0]*fl[3]-0.3518228202874282*nuVtSqSum[0]*fc[3]+0.25*fl[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fl[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.21875*nuVtSqSum[1]*fl[21]+0.21875*nuVtSqSum[1]*fc[21]+0.21875*nuVtSqSum[2]*fl[15]+0.2445699350390395*nuVtSqSum[0]*fl[15]+0.21875*nuVtSqSum[2]*fc[15]+0.2445699350390395*nuVtSqSum[0]*fc[15]+0.3146798968793526*nuVtSqSum[1]*fl[13]-0.3146798968793526*nuVtSqSum[1]*fc[13]+0.2445699350390395*nuVtSqSum[1]*fl[9]+0.2445699350390395*nuVtSqSum[1]*fc[9]+0.223606797749979*nuVtSqSum[1]*fl[7]+0.223606797749979*nuVtSqSum[1]*fc[7]+0.3146798968793526*nuVtSqSum[2]*fl[5]+0.3518228202874282*nuVtSqSum[0]*fl[5]-0.3146798968793526*nuVtSqSum[2]*fc[5]-0.3518228202874282*nuVtSqSum[0]*fc[5]+0.3518228202874282*nuVtSqSum[1]*fl[3]-0.3518228202874282*nuVtSqSum[1]*fc[3]+0.223606797749979*fl[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]+0.25*fl[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fl[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.2445699350390395*nuVtSqSum[2]*fl[24]+0.2445699350390395*nuVtSqSum[2]*fc[24]+0.2445699350390395*nuVtSqSum[1]*fl[19]+0.2445699350390395*nuVtSqSum[1]*fc[19]+0.3518228202874282*nuVtSqSum[2]*fl[17]-0.3518228202874282*nuVtSqSum[2]*fc[17]+0.2445699350390395*nuVtSqSum[0]*fl[16]+0.2445699350390395*nuVtSqSum[0]*fc[16]+0.25*nuVtSqSum[2]*fl[11]+0.25*nuVtSqSum[2]*fc[11]+0.3518228202874282*nuVtSqSum[1]*fl[10]-0.3518228202874282*nuVtSqSum[1]*fc[10]+0.3518228202874282*nuVtSqSum[0]*fl[6]-0.3518228202874282*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fl[4]+0.25*nuVtSqSum[1]*fc[4]+0.25*nuVtSqSum[0]*fl[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_l[3] = 0.21875*nuVtSqSum[1]*fl[24]+0.21875*nuVtSqSum[1]*fc[24]+0.21875*nuVtSqSum[2]*fl[19]+0.2445699350390395*nuVtSqSum[0]*fl[19]+0.21875*nuVtSqSum[2]*fc[19]+0.2445699350390395*nuVtSqSum[0]*fc[19]+0.3146798968793526*nuVtSqSum[1]*fl[17]-0.3146798968793526*nuVtSqSum[1]*fc[17]+0.2445699350390395*nuVtSqSum[1]*fl[16]+0.2445699350390395*nuVtSqSum[1]*fc[16]+0.223606797749979*nuVtSqSum[1]*fl[11]+0.223606797749979*nuVtSqSum[1]*fc[11]+0.3146798968793526*nuVtSqSum[2]*fl[10]+0.3518228202874282*nuVtSqSum[0]*fl[10]-0.3146798968793526*nuVtSqSum[2]*fc[10]-0.3518228202874282*nuVtSqSum[0]*fc[10]+0.3518228202874282*nuVtSqSum[1]*fl[6]-0.3518228202874282*nuVtSqSum[1]*fc[6]+0.223606797749979*nuVtSqSum[2]*fl[4]+0.25*nuVtSqSum[0]*fl[4]+0.223606797749979*nuVtSqSum[2]*fc[4]+0.25*nuVtSqSum[0]*fc[4]+0.25*nuVtSqSum[1]*fl[2]+0.25*nuVtSqSum[1]*fc[2]; 
  Gdiff2_l[4] = 0.15625*nuVtSqSum[2]*fl[21]+0.2445699350390395*nuVtSqSum[0]*fl[21]+0.15625*nuVtSqSum[2]*fc[21]+0.2445699350390395*nuVtSqSum[0]*fc[21]+0.21875*nuVtSqSum[1]*fl[15]+0.21875*nuVtSqSum[1]*fc[15]+0.2247713549138233*nuVtSqSum[2]*fl[13]+0.3518228202874282*nuVtSqSum[0]*fl[13]-0.2247713549138233*nuVtSqSum[2]*fc[13]-0.3518228202874282*nuVtSqSum[0]*fc[13]+0.2445699350390395*nuVtSqSum[2]*fl[9]+0.2445699350390395*nuVtSqSum[2]*fc[9]+0.159719141249985*nuVtSqSum[2]*fl[7]+0.25*nuVtSqSum[0]*fl[7]+0.159719141249985*nuVtSqSum[2]*fc[7]+0.25*nuVtSqSum[0]*fc[7]+0.3146798968793526*nuVtSqSum[1]*fl[5]-0.3146798968793526*nuVtSqSum[1]*fc[5]+0.3518228202874282*nuVtSqSum[2]*fl[3]-0.3518228202874282*nuVtSqSum[2]*fc[3]+0.25*fl[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fl[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 
  Gdiff2_l[5] = 0.2445699350390395*nuVtSqSum[2]*fl[26]+0.2445699350390395*nuVtSqSum[2]*fc[26]+0.2445699350390395*nuVtSqSum[1]*fl[25]+0.2445699350390395*nuVtSqSum[1]*fc[25]+0.3518228202874282*nuVtSqSum[2]*fl[23]-0.3518228202874282*nuVtSqSum[2]*fc[23]+0.2445699350390395*nuVtSqSum[0]*fl[22]+0.2445699350390395*nuVtSqSum[0]*fc[22]+0.25*nuVtSqSum[2]*fl[20]+0.25*nuVtSqSum[2]*fc[20]+0.3518228202874282*nuVtSqSum[1]*fl[18]-0.3518228202874282*nuVtSqSum[1]*fc[18]+0.3518228202874282*nuVtSqSum[0]*fl[14]-0.3518228202874282*nuVtSqSum[0]*fc[14]+0.25*nuVtSqSum[1]*fl[12]+0.25*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[0]*fl[8]+0.25*nuVtSqSum[0]*fc[8]; 
  Gdiff2_l[6] = 0.15625*nuVtSqSum[2]*fl[24]+0.2445699350390395*nuVtSqSum[0]*fl[24]+0.15625*nuVtSqSum[2]*fc[24]+0.2445699350390395*nuVtSqSum[0]*fc[24]+0.21875*nuVtSqSum[1]*fl[19]+0.21875*nuVtSqSum[1]*fc[19]+0.2247713549138233*nuVtSqSum[2]*fl[17]+0.3518228202874282*nuVtSqSum[0]*fl[17]-0.2247713549138233*nuVtSqSum[2]*fc[17]-0.3518228202874282*nuVtSqSum[0]*fc[17]+0.2445699350390395*nuVtSqSum[2]*fl[16]+0.2445699350390395*nuVtSqSum[2]*fc[16]+0.159719141249985*nuVtSqSum[2]*fl[11]+0.25*nuVtSqSum[0]*fl[11]+0.159719141249985*nuVtSqSum[2]*fc[11]+0.25*nuVtSqSum[0]*fc[11]+0.3146798968793526*nuVtSqSum[1]*fl[10]-0.3146798968793526*nuVtSqSum[1]*fc[10]+0.3518228202874282*nuVtSqSum[2]*fl[6]-0.3518228202874282*nuVtSqSum[2]*fc[6]+0.223606797749979*nuVtSqSum[1]*fl[4]+0.223606797749979*nuVtSqSum[1]*fc[4]+0.25*fl[2]*nuVtSqSum[2]+0.25*fc[2]*nuVtSqSum[2]; 
  Gdiff2_l[7] = 0.21875*nuVtSqSum[1]*fl[26]+0.21875*nuVtSqSum[1]*fc[26]+0.21875*nuVtSqSum[2]*fl[25]+0.2445699350390395*nuVtSqSum[0]*fl[25]+0.21875*nuVtSqSum[2]*fc[25]+0.2445699350390395*nuVtSqSum[0]*fc[25]+0.3146798968793526*nuVtSqSum[1]*fl[23]-0.3146798968793526*nuVtSqSum[1]*fc[23]+0.2445699350390395*nuVtSqSum[1]*fl[22]+0.2445699350390395*nuVtSqSum[1]*fc[22]+0.223606797749979*nuVtSqSum[1]*fl[20]+0.223606797749979*nuVtSqSum[1]*fc[20]+0.3146798968793526*nuVtSqSum[2]*fl[18]+0.3518228202874282*nuVtSqSum[0]*fl[18]-0.3146798968793526*nuVtSqSum[2]*fc[18]-0.3518228202874282*nuVtSqSum[0]*fc[18]+0.3518228202874282*nuVtSqSum[1]*fl[14]-0.3518228202874282*nuVtSqSum[1]*fc[14]+0.223606797749979*nuVtSqSum[2]*fl[12]+0.25*nuVtSqSum[0]*fl[12]+0.223606797749979*nuVtSqSum[2]*fc[12]+0.25*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[1]*fl[8]+0.25*nuVtSqSum[1]*fc[8]; 
  Gdiff2_l[8] = 0.15625*nuVtSqSum[2]*fl[26]+0.2445699350390395*nuVtSqSum[0]*fl[26]+0.15625*nuVtSqSum[2]*fc[26]+0.2445699350390395*nuVtSqSum[0]*fc[26]+0.21875*nuVtSqSum[1]*fl[25]+0.21875*nuVtSqSum[1]*fc[25]+0.2247713549138233*nuVtSqSum[2]*fl[23]+0.3518228202874282*nuVtSqSum[0]*fl[23]-0.2247713549138233*nuVtSqSum[2]*fc[23]-0.3518228202874282*nuVtSqSum[0]*fc[23]+0.2445699350390395*nuVtSqSum[2]*fl[22]+0.2445699350390395*nuVtSqSum[2]*fc[22]+0.159719141249985*nuVtSqSum[2]*fl[20]+0.25*nuVtSqSum[0]*fl[20]+0.159719141249985*nuVtSqSum[2]*fc[20]+0.25*nuVtSqSum[0]*fc[20]+0.3146798968793526*nuVtSqSum[1]*fl[18]-0.3146798968793526*nuVtSqSum[1]*fc[18]+0.3518228202874282*nuVtSqSum[2]*fl[14]-0.3518228202874282*nuVtSqSum[2]*fc[14]+0.223606797749979*nuVtSqSum[1]*fl[12]+0.223606797749979*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[2]*fl[8]+0.25*nuVtSqSum[2]*fc[8]; 

  Gdiff2_r[0] = 0.2445699350390395*nuVtSqSum[2]*fr[21]+0.2445699350390395*nuVtSqSum[2]*fc[21]+0.2445699350390395*nuVtSqSum[1]*fr[15]+0.2445699350390395*nuVtSqSum[1]*fc[15]-0.3518228202874282*nuVtSqSum[2]*fr[13]+0.3518228202874282*nuVtSqSum[2]*fc[13]+0.2445699350390395*nuVtSqSum[0]*fr[9]+0.2445699350390395*nuVtSqSum[0]*fc[9]+0.25*nuVtSqSum[2]*fr[7]+0.25*nuVtSqSum[2]*fc[7]-0.3518228202874282*nuVtSqSum[1]*fr[5]+0.3518228202874282*nuVtSqSum[1]*fc[5]-0.3518228202874282*nuVtSqSum[0]*fr[3]+0.3518228202874282*nuVtSqSum[0]*fc[3]+0.25*fr[1]*nuVtSqSum[1]+0.25*fc[1]*nuVtSqSum[1]+0.25*fr[0]*nuVtSqSum[0]+0.25*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = 0.21875*nuVtSqSum[1]*fr[21]+0.21875*nuVtSqSum[1]*fc[21]+0.21875*nuVtSqSum[2]*fr[15]+0.2445699350390395*nuVtSqSum[0]*fr[15]+0.21875*nuVtSqSum[2]*fc[15]+0.2445699350390395*nuVtSqSum[0]*fc[15]-0.3146798968793526*nuVtSqSum[1]*fr[13]+0.3146798968793526*nuVtSqSum[1]*fc[13]+0.2445699350390395*nuVtSqSum[1]*fr[9]+0.2445699350390395*nuVtSqSum[1]*fc[9]+0.223606797749979*nuVtSqSum[1]*fr[7]+0.223606797749979*nuVtSqSum[1]*fc[7]-0.3146798968793526*nuVtSqSum[2]*fr[5]-0.3518228202874282*nuVtSqSum[0]*fr[5]+0.3146798968793526*nuVtSqSum[2]*fc[5]+0.3518228202874282*nuVtSqSum[0]*fc[5]-0.3518228202874282*nuVtSqSum[1]*fr[3]+0.3518228202874282*nuVtSqSum[1]*fc[3]+0.223606797749979*fr[1]*nuVtSqSum[2]+0.223606797749979*fc[1]*nuVtSqSum[2]+0.25*fr[0]*nuVtSqSum[1]+0.25*fc[0]*nuVtSqSum[1]+0.25*nuVtSqSum[0]*fr[1]+0.25*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = 0.2445699350390395*nuVtSqSum[2]*fr[24]+0.2445699350390395*nuVtSqSum[2]*fc[24]+0.2445699350390395*nuVtSqSum[1]*fr[19]+0.2445699350390395*nuVtSqSum[1]*fc[19]-0.3518228202874282*nuVtSqSum[2]*fr[17]+0.3518228202874282*nuVtSqSum[2]*fc[17]+0.2445699350390395*nuVtSqSum[0]*fr[16]+0.2445699350390395*nuVtSqSum[0]*fc[16]+0.25*nuVtSqSum[2]*fr[11]+0.25*nuVtSqSum[2]*fc[11]-0.3518228202874282*nuVtSqSum[1]*fr[10]+0.3518228202874282*nuVtSqSum[1]*fc[10]-0.3518228202874282*nuVtSqSum[0]*fr[6]+0.3518228202874282*nuVtSqSum[0]*fc[6]+0.25*nuVtSqSum[1]*fr[4]+0.25*nuVtSqSum[1]*fc[4]+0.25*nuVtSqSum[0]*fr[2]+0.25*nuVtSqSum[0]*fc[2]; 
  Gdiff2_r[3] = 0.21875*nuVtSqSum[1]*fr[24]+0.21875*nuVtSqSum[1]*fc[24]+0.21875*nuVtSqSum[2]*fr[19]+0.2445699350390395*nuVtSqSum[0]*fr[19]+0.21875*nuVtSqSum[2]*fc[19]+0.2445699350390395*nuVtSqSum[0]*fc[19]-0.3146798968793526*nuVtSqSum[1]*fr[17]+0.3146798968793526*nuVtSqSum[1]*fc[17]+0.2445699350390395*nuVtSqSum[1]*fr[16]+0.2445699350390395*nuVtSqSum[1]*fc[16]+0.223606797749979*nuVtSqSum[1]*fr[11]+0.223606797749979*nuVtSqSum[1]*fc[11]-0.3146798968793526*nuVtSqSum[2]*fr[10]-0.3518228202874282*nuVtSqSum[0]*fr[10]+0.3146798968793526*nuVtSqSum[2]*fc[10]+0.3518228202874282*nuVtSqSum[0]*fc[10]-0.3518228202874282*nuVtSqSum[1]*fr[6]+0.3518228202874282*nuVtSqSum[1]*fc[6]+0.223606797749979*nuVtSqSum[2]*fr[4]+0.25*nuVtSqSum[0]*fr[4]+0.223606797749979*nuVtSqSum[2]*fc[4]+0.25*nuVtSqSum[0]*fc[4]+0.25*nuVtSqSum[1]*fr[2]+0.25*nuVtSqSum[1]*fc[2]; 
  Gdiff2_r[4] = 0.15625*nuVtSqSum[2]*fr[21]+0.2445699350390395*nuVtSqSum[0]*fr[21]+0.15625*nuVtSqSum[2]*fc[21]+0.2445699350390395*nuVtSqSum[0]*fc[21]+0.21875*nuVtSqSum[1]*fr[15]+0.21875*nuVtSqSum[1]*fc[15]-0.2247713549138233*nuVtSqSum[2]*fr[13]-0.3518228202874282*nuVtSqSum[0]*fr[13]+0.2247713549138233*nuVtSqSum[2]*fc[13]+0.3518228202874282*nuVtSqSum[0]*fc[13]+0.2445699350390395*nuVtSqSum[2]*fr[9]+0.2445699350390395*nuVtSqSum[2]*fc[9]+0.159719141249985*nuVtSqSum[2]*fr[7]+0.25*nuVtSqSum[0]*fr[7]+0.159719141249985*nuVtSqSum[2]*fc[7]+0.25*nuVtSqSum[0]*fc[7]-0.3146798968793526*nuVtSqSum[1]*fr[5]+0.3146798968793526*nuVtSqSum[1]*fc[5]-0.3518228202874282*nuVtSqSum[2]*fr[3]+0.3518228202874282*nuVtSqSum[2]*fc[3]+0.25*fr[0]*nuVtSqSum[2]+0.25*fc[0]*nuVtSqSum[2]+0.223606797749979*fr[1]*nuVtSqSum[1]+0.223606797749979*fc[1]*nuVtSqSum[1]; 
  Gdiff2_r[5] = 0.2445699350390395*nuVtSqSum[2]*fr[26]+0.2445699350390395*nuVtSqSum[2]*fc[26]+0.2445699350390395*nuVtSqSum[1]*fr[25]+0.2445699350390395*nuVtSqSum[1]*fc[25]-0.3518228202874282*nuVtSqSum[2]*fr[23]+0.3518228202874282*nuVtSqSum[2]*fc[23]+0.2445699350390395*nuVtSqSum[0]*fr[22]+0.2445699350390395*nuVtSqSum[0]*fc[22]+0.25*nuVtSqSum[2]*fr[20]+0.25*nuVtSqSum[2]*fc[20]-0.3518228202874282*nuVtSqSum[1]*fr[18]+0.3518228202874282*nuVtSqSum[1]*fc[18]-0.3518228202874282*nuVtSqSum[0]*fr[14]+0.3518228202874282*nuVtSqSum[0]*fc[14]+0.25*nuVtSqSum[1]*fr[12]+0.25*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[0]*fr[8]+0.25*nuVtSqSum[0]*fc[8]; 
  Gdiff2_r[6] = 0.15625*nuVtSqSum[2]*fr[24]+0.2445699350390395*nuVtSqSum[0]*fr[24]+0.15625*nuVtSqSum[2]*fc[24]+0.2445699350390395*nuVtSqSum[0]*fc[24]+0.21875*nuVtSqSum[1]*fr[19]+0.21875*nuVtSqSum[1]*fc[19]-0.2247713549138233*nuVtSqSum[2]*fr[17]-0.3518228202874282*nuVtSqSum[0]*fr[17]+0.2247713549138233*nuVtSqSum[2]*fc[17]+0.3518228202874282*nuVtSqSum[0]*fc[17]+0.2445699350390395*nuVtSqSum[2]*fr[16]+0.2445699350390395*nuVtSqSum[2]*fc[16]+0.159719141249985*nuVtSqSum[2]*fr[11]+0.25*nuVtSqSum[0]*fr[11]+0.159719141249985*nuVtSqSum[2]*fc[11]+0.25*nuVtSqSum[0]*fc[11]-0.3146798968793526*nuVtSqSum[1]*fr[10]+0.3146798968793526*nuVtSqSum[1]*fc[10]-0.3518228202874282*nuVtSqSum[2]*fr[6]+0.3518228202874282*nuVtSqSum[2]*fc[6]+0.223606797749979*nuVtSqSum[1]*fr[4]+0.223606797749979*nuVtSqSum[1]*fc[4]+0.25*fr[2]*nuVtSqSum[2]+0.25*fc[2]*nuVtSqSum[2]; 
  Gdiff2_r[7] = 0.21875*nuVtSqSum[1]*fr[26]+0.21875*nuVtSqSum[1]*fc[26]+0.21875*nuVtSqSum[2]*fr[25]+0.2445699350390395*nuVtSqSum[0]*fr[25]+0.21875*nuVtSqSum[2]*fc[25]+0.2445699350390395*nuVtSqSum[0]*fc[25]-0.3146798968793526*nuVtSqSum[1]*fr[23]+0.3146798968793526*nuVtSqSum[1]*fc[23]+0.2445699350390395*nuVtSqSum[1]*fr[22]+0.2445699350390395*nuVtSqSum[1]*fc[22]+0.223606797749979*nuVtSqSum[1]*fr[20]+0.223606797749979*nuVtSqSum[1]*fc[20]-0.3146798968793526*nuVtSqSum[2]*fr[18]-0.3518228202874282*nuVtSqSum[0]*fr[18]+0.3146798968793526*nuVtSqSum[2]*fc[18]+0.3518228202874282*nuVtSqSum[0]*fc[18]-0.3518228202874282*nuVtSqSum[1]*fr[14]+0.3518228202874282*nuVtSqSum[1]*fc[14]+0.223606797749979*nuVtSqSum[2]*fr[12]+0.25*nuVtSqSum[0]*fr[12]+0.223606797749979*nuVtSqSum[2]*fc[12]+0.25*nuVtSqSum[0]*fc[12]+0.25*nuVtSqSum[1]*fr[8]+0.25*nuVtSqSum[1]*fc[8]; 
  Gdiff2_r[8] = 0.15625*nuVtSqSum[2]*fr[26]+0.2445699350390395*nuVtSqSum[0]*fr[26]+0.15625*nuVtSqSum[2]*fc[26]+0.2445699350390395*nuVtSqSum[0]*fc[26]+0.21875*nuVtSqSum[1]*fr[25]+0.21875*nuVtSqSum[1]*fc[25]-0.2247713549138233*nuVtSqSum[2]*fr[23]-0.3518228202874282*nuVtSqSum[0]*fr[23]+0.2247713549138233*nuVtSqSum[2]*fc[23]+0.3518228202874282*nuVtSqSum[0]*fc[23]+0.2445699350390395*nuVtSqSum[2]*fr[22]+0.2445699350390395*nuVtSqSum[2]*fc[22]+0.159719141249985*nuVtSqSum[2]*fr[20]+0.25*nuVtSqSum[0]*fr[20]+0.159719141249985*nuVtSqSum[2]*fc[20]+0.25*nuVtSqSum[0]*fc[20]-0.3146798968793526*nuVtSqSum[1]*fr[18]+0.3146798968793526*nuVtSqSum[1]*fc[18]-0.3518228202874282*nuVtSqSum[2]*fr[14]+0.3518228202874282*nuVtSqSum[2]*fc[14]+0.223606797749979*nuVtSqSum[1]*fr[12]+0.223606797749979*nuVtSqSum[1]*fc[12]+0.25*nuVtSqSum[2]*fr[8]+0.25*nuVtSqSum[2]*fc[8]; 

  Gdiff_l[0] = (-0.6708203932499369*nuVtSqSum[2]*fl[21])+0.6708203932499369*nuVtSqSum[2]*fc[21]-0.6708203932499369*nuVtSqSum[1]*fl[15]+0.6708203932499369*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[2]*fl[13]-1.190784930203603*nuVtSqSum[2]*fc[13]-0.6708203932499369*nuVtSqSum[0]*fl[9]+0.6708203932499369*nuVtSqSum[0]*fc[9]-0.9375*nuVtSqSum[2]*fl[7]+0.9375*nuVtSqSum[2]*fc[7]-1.190784930203603*nuVtSqSum[1]*fl[5]-1.190784930203603*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[0]*fl[3]-1.190784930203603*nuVtSqSum[0]*fc[3]-0.9375*fl[1]*nuVtSqSum[1]+0.9375*fc[1]*nuVtSqSum[1]-0.9375*fl[0]*nuVtSqSum[0]+0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.6*nuVtSqSum[1]*fl[21])+0.6*nuVtSqSum[1]*fc[21]-0.5999999999999999*nuVtSqSum[2]*fl[15]-0.6708203932499369*nuVtSqSum[0]*fl[15]+0.5999999999999999*nuVtSqSum[2]*fc[15]+0.6708203932499369*nuVtSqSum[0]*fc[15]-1.06507042020704*nuVtSqSum[1]*fl[13]-1.06507042020704*nuVtSqSum[1]*fc[13]-0.6708203932499369*nuVtSqSum[1]*fl[9]+0.6708203932499369*nuVtSqSum[1]*fc[9]-0.8385254915624212*nuVtSqSum[1]*fl[7]+0.8385254915624212*nuVtSqSum[1]*fc[7]-1.06507042020704*nuVtSqSum[2]*fl[5]-1.190784930203603*nuVtSqSum[0]*fl[5]-1.06507042020704*nuVtSqSum[2]*fc[5]-1.190784930203603*nuVtSqSum[0]*fc[5]-1.190784930203603*nuVtSqSum[1]*fl[3]-1.190784930203603*nuVtSqSum[1]*fc[3]-0.8385254915624212*fl[1]*nuVtSqSum[2]+0.8385254915624212*fc[1]*nuVtSqSum[2]-0.9375*fl[0]*nuVtSqSum[1]+0.9375*fc[0]*nuVtSqSum[1]-0.9375*nuVtSqSum[0]*fl[1]+0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.6708203932499369*nuVtSqSum[2]*fl[24])+0.6708203932499369*nuVtSqSum[2]*fc[24]-0.6708203932499369*nuVtSqSum[1]*fl[19]+0.6708203932499369*nuVtSqSum[1]*fc[19]-1.190784930203603*nuVtSqSum[2]*fl[17]-1.190784930203603*nuVtSqSum[2]*fc[17]-0.6708203932499369*nuVtSqSum[0]*fl[16]+0.6708203932499369*nuVtSqSum[0]*fc[16]-0.9375000000000001*nuVtSqSum[2]*fl[11]+0.9375000000000001*nuVtSqSum[2]*fc[11]-1.190784930203603*nuVtSqSum[1]*fl[10]-1.190784930203603*nuVtSqSum[1]*fc[10]-1.190784930203603*nuVtSqSum[0]*fl[6]-1.190784930203603*nuVtSqSum[0]*fc[6]-0.9375*nuVtSqSum[1]*fl[4]+0.9375*nuVtSqSum[1]*fc[4]-0.9375*nuVtSqSum[0]*fl[2]+0.9375*nuVtSqSum[0]*fc[2]; 
  Gdiff_l[3] = (-0.6*nuVtSqSum[1]*fl[24])+0.6*nuVtSqSum[1]*fc[24]-0.6*nuVtSqSum[2]*fl[19]-0.6708203932499369*nuVtSqSum[0]*fl[19]+0.6*nuVtSqSum[2]*fc[19]+0.6708203932499369*nuVtSqSum[0]*fc[19]-1.06507042020704*nuVtSqSum[1]*fl[17]-1.06507042020704*nuVtSqSum[1]*fc[17]-0.6708203932499369*nuVtSqSum[1]*fl[16]+0.6708203932499369*nuVtSqSum[1]*fc[16]-0.8385254915624211*nuVtSqSum[1]*fl[11]+0.8385254915624211*nuVtSqSum[1]*fc[11]-1.06507042020704*nuVtSqSum[2]*fl[10]-1.190784930203603*nuVtSqSum[0]*fl[10]-1.06507042020704*nuVtSqSum[2]*fc[10]-1.190784930203603*nuVtSqSum[0]*fc[10]-1.190784930203603*nuVtSqSum[1]*fl[6]-1.190784930203603*nuVtSqSum[1]*fc[6]-0.8385254915624212*nuVtSqSum[2]*fl[4]-0.9375*nuVtSqSum[0]*fl[4]+0.8385254915624212*nuVtSqSum[2]*fc[4]+0.9375*nuVtSqSum[0]*fc[4]-0.9375*nuVtSqSum[1]*fl[2]+0.9375*nuVtSqSum[1]*fc[2]; 
  Gdiff_l[4] = (-0.4285714285714285*nuVtSqSum[2]*fl[21])-0.6708203932499369*nuVtSqSum[0]*fl[21]+0.4285714285714285*nuVtSqSum[2]*fc[21]+0.6708203932499369*nuVtSqSum[0]*fc[21]-0.5999999999999999*nuVtSqSum[1]*fl[15]+0.5999999999999999*nuVtSqSum[1]*fc[15]-0.7607645858621712*nuVtSqSum[2]*fl[13]-1.190784930203603*nuVtSqSum[0]*fl[13]-0.7607645858621712*nuVtSqSum[2]*fc[13]-1.190784930203603*nuVtSqSum[0]*fc[13]-0.6708203932499369*nuVtSqSum[2]*fl[9]+0.6708203932499369*nuVtSqSum[2]*fc[9]-0.5989467796874438*nuVtSqSum[2]*fl[7]-0.9375*nuVtSqSum[0]*fl[7]+0.5989467796874438*nuVtSqSum[2]*fc[7]+0.9375*nuVtSqSum[0]*fc[7]-1.06507042020704*nuVtSqSum[1]*fl[5]-1.06507042020704*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[2]*fl[3]-1.190784930203603*nuVtSqSum[2]*fc[3]-0.9375*fl[0]*nuVtSqSum[2]+0.9375*fc[0]*nuVtSqSum[2]-0.8385254915624212*fl[1]*nuVtSqSum[1]+0.8385254915624212*fc[1]*nuVtSqSum[1]; 
  Gdiff_l[5] = (-0.6708203932499369*nuVtSqSum[2]*fl[26])+0.6708203932499369*nuVtSqSum[2]*fc[26]-0.6708203932499369*nuVtSqSum[1]*fl[25]+0.6708203932499369*nuVtSqSum[1]*fc[25]-1.190784930203603*nuVtSqSum[2]*fl[23]-1.190784930203603*nuVtSqSum[2]*fc[23]-0.6708203932499369*nuVtSqSum[0]*fl[22]+0.6708203932499369*nuVtSqSum[0]*fc[22]-0.9375*nuVtSqSum[2]*fl[20]+0.9375*nuVtSqSum[2]*fc[20]-1.190784930203603*nuVtSqSum[1]*fl[18]-1.190784930203603*nuVtSqSum[1]*fc[18]-1.190784930203603*nuVtSqSum[0]*fl[14]-1.190784930203603*nuVtSqSum[0]*fc[14]-0.9375000000000001*nuVtSqSum[1]*fl[12]+0.9375000000000001*nuVtSqSum[1]*fc[12]-0.9375*nuVtSqSum[0]*fl[8]+0.9375*nuVtSqSum[0]*fc[8]; 
  Gdiff_l[6] = (-0.4285714285714285*nuVtSqSum[2]*fl[24])-0.6708203932499369*nuVtSqSum[0]*fl[24]+0.4285714285714285*nuVtSqSum[2]*fc[24]+0.6708203932499369*nuVtSqSum[0]*fc[24]-0.5999999999999999*nuVtSqSum[1]*fl[19]+0.5999999999999999*nuVtSqSum[1]*fc[19]-0.7607645858621712*nuVtSqSum[2]*fl[17]-1.190784930203603*nuVtSqSum[0]*fl[17]-0.7607645858621712*nuVtSqSum[2]*fc[17]-1.190784930203603*nuVtSqSum[0]*fc[17]-0.6708203932499369*nuVtSqSum[2]*fl[16]+0.6708203932499369*nuVtSqSum[2]*fc[16]-0.5989467796874438*nuVtSqSum[2]*fl[11]-0.9375*nuVtSqSum[0]*fl[11]+0.5989467796874438*nuVtSqSum[2]*fc[11]+0.9375*nuVtSqSum[0]*fc[11]-1.06507042020704*nuVtSqSum[1]*fl[10]-1.06507042020704*nuVtSqSum[1]*fc[10]-1.190784930203603*nuVtSqSum[2]*fl[6]-1.190784930203603*nuVtSqSum[2]*fc[6]-0.8385254915624211*nuVtSqSum[1]*fl[4]+0.8385254915624211*nuVtSqSum[1]*fc[4]-0.9375000000000001*fl[2]*nuVtSqSum[2]+0.9375000000000001*fc[2]*nuVtSqSum[2]; 
  Gdiff_l[7] = (-0.5999999999999999*nuVtSqSum[1]*fl[26])+0.5999999999999999*nuVtSqSum[1]*fc[26]-0.5999999999999999*nuVtSqSum[2]*fl[25]-0.6708203932499369*nuVtSqSum[0]*fl[25]+0.5999999999999999*nuVtSqSum[2]*fc[25]+0.6708203932499369*nuVtSqSum[0]*fc[25]-1.06507042020704*nuVtSqSum[1]*fl[23]-1.06507042020704*nuVtSqSum[1]*fc[23]-0.6708203932499369*nuVtSqSum[1]*fl[22]+0.6708203932499369*nuVtSqSum[1]*fc[22]-0.8385254915624211*nuVtSqSum[1]*fl[20]+0.8385254915624211*nuVtSqSum[1]*fc[20]-1.06507042020704*nuVtSqSum[2]*fl[18]-1.190784930203603*nuVtSqSum[0]*fl[18]-1.06507042020704*nuVtSqSum[2]*fc[18]-1.190784930203603*nuVtSqSum[0]*fc[18]-1.190784930203603*nuVtSqSum[1]*fl[14]-1.190784930203603*nuVtSqSum[1]*fc[14]-0.8385254915624212*nuVtSqSum[2]*fl[12]-0.9375*nuVtSqSum[0]*fl[12]+0.8385254915624212*nuVtSqSum[2]*fc[12]+0.9375*nuVtSqSum[0]*fc[12]-0.9375000000000001*nuVtSqSum[1]*fl[8]+0.9375000000000001*nuVtSqSum[1]*fc[8]; 
  Gdiff_l[8] = (-0.4285714285714285*nuVtSqSum[2]*fl[26])-0.6708203932499369*nuVtSqSum[0]*fl[26]+0.4285714285714285*nuVtSqSum[2]*fc[26]+0.6708203932499369*nuVtSqSum[0]*fc[26]-0.6*nuVtSqSum[1]*fl[25]+0.6*nuVtSqSum[1]*fc[25]-0.7607645858621711*nuVtSqSum[2]*fl[23]-1.190784930203603*nuVtSqSum[0]*fl[23]-0.7607645858621711*nuVtSqSum[2]*fc[23]-1.190784930203603*nuVtSqSum[0]*fc[23]-0.6708203932499369*nuVtSqSum[2]*fl[22]+0.6708203932499369*nuVtSqSum[2]*fc[22]-0.5989467796874438*nuVtSqSum[2]*fl[20]-0.9375*nuVtSqSum[0]*fl[20]+0.5989467796874438*nuVtSqSum[2]*fc[20]+0.9375*nuVtSqSum[0]*fc[20]-1.06507042020704*nuVtSqSum[1]*fl[18]-1.06507042020704*nuVtSqSum[1]*fc[18]-1.190784930203603*nuVtSqSum[2]*fl[14]-1.190784930203603*nuVtSqSum[2]*fc[14]-0.8385254915624211*nuVtSqSum[1]*fl[12]+0.8385254915624211*nuVtSqSum[1]*fc[12]-0.9375*nuVtSqSum[2]*fl[8]+0.9375*nuVtSqSum[2]*fc[8]; 

  Gdiff_r[0] = 0.6708203932499369*nuVtSqSum[2]*fr[21]-0.6708203932499369*nuVtSqSum[2]*fc[21]+0.6708203932499369*nuVtSqSum[1]*fr[15]-0.6708203932499369*nuVtSqSum[1]*fc[15]-1.190784930203603*nuVtSqSum[2]*fr[13]-1.190784930203603*nuVtSqSum[2]*fc[13]+0.6708203932499369*nuVtSqSum[0]*fr[9]-0.6708203932499369*nuVtSqSum[0]*fc[9]+0.9375*nuVtSqSum[2]*fr[7]-0.9375*nuVtSqSum[2]*fc[7]-1.190784930203603*nuVtSqSum[1]*fr[5]-1.190784930203603*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[0]*fr[3]-1.190784930203603*nuVtSqSum[0]*fc[3]+0.9375*fr[1]*nuVtSqSum[1]-0.9375*fc[1]*nuVtSqSum[1]+0.9375*fr[0]*nuVtSqSum[0]-0.9375*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = 0.6*nuVtSqSum[1]*fr[21]-0.6*nuVtSqSum[1]*fc[21]+0.5999999999999999*nuVtSqSum[2]*fr[15]+0.6708203932499369*nuVtSqSum[0]*fr[15]-0.5999999999999999*nuVtSqSum[2]*fc[15]-0.6708203932499369*nuVtSqSum[0]*fc[15]-1.06507042020704*nuVtSqSum[1]*fr[13]-1.06507042020704*nuVtSqSum[1]*fc[13]+0.6708203932499369*nuVtSqSum[1]*fr[9]-0.6708203932499369*nuVtSqSum[1]*fc[9]+0.8385254915624212*nuVtSqSum[1]*fr[7]-0.8385254915624212*nuVtSqSum[1]*fc[7]-1.06507042020704*nuVtSqSum[2]*fr[5]-1.190784930203603*nuVtSqSum[0]*fr[5]-1.06507042020704*nuVtSqSum[2]*fc[5]-1.190784930203603*nuVtSqSum[0]*fc[5]-1.190784930203603*nuVtSqSum[1]*fr[3]-1.190784930203603*nuVtSqSum[1]*fc[3]+0.8385254915624212*fr[1]*nuVtSqSum[2]-0.8385254915624212*fc[1]*nuVtSqSum[2]+0.9375*fr[0]*nuVtSqSum[1]-0.9375*fc[0]*nuVtSqSum[1]+0.9375*nuVtSqSum[0]*fr[1]-0.9375*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = 0.6708203932499369*nuVtSqSum[2]*fr[24]-0.6708203932499369*nuVtSqSum[2]*fc[24]+0.6708203932499369*nuVtSqSum[1]*fr[19]-0.6708203932499369*nuVtSqSum[1]*fc[19]-1.190784930203603*nuVtSqSum[2]*fr[17]-1.190784930203603*nuVtSqSum[2]*fc[17]+0.6708203932499369*nuVtSqSum[0]*fr[16]-0.6708203932499369*nuVtSqSum[0]*fc[16]+0.9375000000000001*nuVtSqSum[2]*fr[11]-0.9375000000000001*nuVtSqSum[2]*fc[11]-1.190784930203603*nuVtSqSum[1]*fr[10]-1.190784930203603*nuVtSqSum[1]*fc[10]-1.190784930203603*nuVtSqSum[0]*fr[6]-1.190784930203603*nuVtSqSum[0]*fc[6]+0.9375*nuVtSqSum[1]*fr[4]-0.9375*nuVtSqSum[1]*fc[4]+0.9375*nuVtSqSum[0]*fr[2]-0.9375*nuVtSqSum[0]*fc[2]; 
  Gdiff_r[3] = 0.6*nuVtSqSum[1]*fr[24]-0.6*nuVtSqSum[1]*fc[24]+0.6*nuVtSqSum[2]*fr[19]+0.6708203932499369*nuVtSqSum[0]*fr[19]-0.6*nuVtSqSum[2]*fc[19]-0.6708203932499369*nuVtSqSum[0]*fc[19]-1.06507042020704*nuVtSqSum[1]*fr[17]-1.06507042020704*nuVtSqSum[1]*fc[17]+0.6708203932499369*nuVtSqSum[1]*fr[16]-0.6708203932499369*nuVtSqSum[1]*fc[16]+0.8385254915624211*nuVtSqSum[1]*fr[11]-0.8385254915624211*nuVtSqSum[1]*fc[11]-1.06507042020704*nuVtSqSum[2]*fr[10]-1.190784930203603*nuVtSqSum[0]*fr[10]-1.06507042020704*nuVtSqSum[2]*fc[10]-1.190784930203603*nuVtSqSum[0]*fc[10]-1.190784930203603*nuVtSqSum[1]*fr[6]-1.190784930203603*nuVtSqSum[1]*fc[6]+0.8385254915624212*nuVtSqSum[2]*fr[4]+0.9375*nuVtSqSum[0]*fr[4]-0.8385254915624212*nuVtSqSum[2]*fc[4]-0.9375*nuVtSqSum[0]*fc[4]+0.9375*nuVtSqSum[1]*fr[2]-0.9375*nuVtSqSum[1]*fc[2]; 
  Gdiff_r[4] = 0.4285714285714285*nuVtSqSum[2]*fr[21]+0.6708203932499369*nuVtSqSum[0]*fr[21]-0.4285714285714285*nuVtSqSum[2]*fc[21]-0.6708203932499369*nuVtSqSum[0]*fc[21]+0.5999999999999999*nuVtSqSum[1]*fr[15]-0.5999999999999999*nuVtSqSum[1]*fc[15]-0.7607645858621712*nuVtSqSum[2]*fr[13]-1.190784930203603*nuVtSqSum[0]*fr[13]-0.7607645858621712*nuVtSqSum[2]*fc[13]-1.190784930203603*nuVtSqSum[0]*fc[13]+0.6708203932499369*nuVtSqSum[2]*fr[9]-0.6708203932499369*nuVtSqSum[2]*fc[9]+0.5989467796874438*nuVtSqSum[2]*fr[7]+0.9375*nuVtSqSum[0]*fr[7]-0.5989467796874438*nuVtSqSum[2]*fc[7]-0.9375*nuVtSqSum[0]*fc[7]-1.06507042020704*nuVtSqSum[1]*fr[5]-1.06507042020704*nuVtSqSum[1]*fc[5]-1.190784930203603*nuVtSqSum[2]*fr[3]-1.190784930203603*nuVtSqSum[2]*fc[3]+0.9375*fr[0]*nuVtSqSum[2]-0.9375*fc[0]*nuVtSqSum[2]+0.8385254915624212*fr[1]*nuVtSqSum[1]-0.8385254915624212*fc[1]*nuVtSqSum[1]; 
  Gdiff_r[5] = 0.6708203932499369*nuVtSqSum[2]*fr[26]-0.6708203932499369*nuVtSqSum[2]*fc[26]+0.6708203932499369*nuVtSqSum[1]*fr[25]-0.6708203932499369*nuVtSqSum[1]*fc[25]-1.190784930203603*nuVtSqSum[2]*fr[23]-1.190784930203603*nuVtSqSum[2]*fc[23]+0.6708203932499369*nuVtSqSum[0]*fr[22]-0.6708203932499369*nuVtSqSum[0]*fc[22]+0.9375*nuVtSqSum[2]*fr[20]-0.9375*nuVtSqSum[2]*fc[20]-1.190784930203603*nuVtSqSum[1]*fr[18]-1.190784930203603*nuVtSqSum[1]*fc[18]-1.190784930203603*nuVtSqSum[0]*fr[14]-1.190784930203603*nuVtSqSum[0]*fc[14]+0.9375000000000001*nuVtSqSum[1]*fr[12]-0.9375000000000001*nuVtSqSum[1]*fc[12]+0.9375*nuVtSqSum[0]*fr[8]-0.9375*nuVtSqSum[0]*fc[8]; 
  Gdiff_r[6] = 0.4285714285714285*nuVtSqSum[2]*fr[24]+0.6708203932499369*nuVtSqSum[0]*fr[24]-0.4285714285714285*nuVtSqSum[2]*fc[24]-0.6708203932499369*nuVtSqSum[0]*fc[24]+0.5999999999999999*nuVtSqSum[1]*fr[19]-0.5999999999999999*nuVtSqSum[1]*fc[19]-0.7607645858621712*nuVtSqSum[2]*fr[17]-1.190784930203603*nuVtSqSum[0]*fr[17]-0.7607645858621712*nuVtSqSum[2]*fc[17]-1.190784930203603*nuVtSqSum[0]*fc[17]+0.6708203932499369*nuVtSqSum[2]*fr[16]-0.6708203932499369*nuVtSqSum[2]*fc[16]+0.5989467796874438*nuVtSqSum[2]*fr[11]+0.9375*nuVtSqSum[0]*fr[11]-0.5989467796874438*nuVtSqSum[2]*fc[11]-0.9375*nuVtSqSum[0]*fc[11]-1.06507042020704*nuVtSqSum[1]*fr[10]-1.06507042020704*nuVtSqSum[1]*fc[10]-1.190784930203603*nuVtSqSum[2]*fr[6]-1.190784930203603*nuVtSqSum[2]*fc[6]+0.8385254915624211*nuVtSqSum[1]*fr[4]-0.8385254915624211*nuVtSqSum[1]*fc[4]+0.9375000000000001*fr[2]*nuVtSqSum[2]-0.9375000000000001*fc[2]*nuVtSqSum[2]; 
  Gdiff_r[7] = 0.5999999999999999*nuVtSqSum[1]*fr[26]-0.5999999999999999*nuVtSqSum[1]*fc[26]+0.5999999999999999*nuVtSqSum[2]*fr[25]+0.6708203932499369*nuVtSqSum[0]*fr[25]-0.5999999999999999*nuVtSqSum[2]*fc[25]-0.6708203932499369*nuVtSqSum[0]*fc[25]-1.06507042020704*nuVtSqSum[1]*fr[23]-1.06507042020704*nuVtSqSum[1]*fc[23]+0.6708203932499369*nuVtSqSum[1]*fr[22]-0.6708203932499369*nuVtSqSum[1]*fc[22]+0.8385254915624211*nuVtSqSum[1]*fr[20]-0.8385254915624211*nuVtSqSum[1]*fc[20]-1.06507042020704*nuVtSqSum[2]*fr[18]-1.190784930203603*nuVtSqSum[0]*fr[18]-1.06507042020704*nuVtSqSum[2]*fc[18]-1.190784930203603*nuVtSqSum[0]*fc[18]-1.190784930203603*nuVtSqSum[1]*fr[14]-1.190784930203603*nuVtSqSum[1]*fc[14]+0.8385254915624212*nuVtSqSum[2]*fr[12]+0.9375*nuVtSqSum[0]*fr[12]-0.8385254915624212*nuVtSqSum[2]*fc[12]-0.9375*nuVtSqSum[0]*fc[12]+0.9375000000000001*nuVtSqSum[1]*fr[8]-0.9375000000000001*nuVtSqSum[1]*fc[8]; 
  Gdiff_r[8] = 0.4285714285714285*nuVtSqSum[2]*fr[26]+0.6708203932499369*nuVtSqSum[0]*fr[26]-0.4285714285714285*nuVtSqSum[2]*fc[26]-0.6708203932499369*nuVtSqSum[0]*fc[26]+0.6*nuVtSqSum[1]*fr[25]-0.6*nuVtSqSum[1]*fc[25]-0.7607645858621711*nuVtSqSum[2]*fr[23]-1.190784930203603*nuVtSqSum[0]*fr[23]-0.7607645858621711*nuVtSqSum[2]*fc[23]-1.190784930203603*nuVtSqSum[0]*fc[23]+0.6708203932499369*nuVtSqSum[2]*fr[22]-0.6708203932499369*nuVtSqSum[2]*fc[22]+0.5989467796874438*nuVtSqSum[2]*fr[20]+0.9375*nuVtSqSum[0]*fr[20]-0.5989467796874438*nuVtSqSum[2]*fc[20]-0.9375*nuVtSqSum[0]*fc[20]-1.06507042020704*nuVtSqSum[1]*fr[18]-1.06507042020704*nuVtSqSum[1]*fc[18]-1.190784930203603*nuVtSqSum[2]*fr[14]-1.190784930203603*nuVtSqSum[2]*fc[14]+0.8385254915624211*nuVtSqSum[1]*fr[12]-0.8385254915624211*nuVtSqSum[1]*fc[12]+0.9375*nuVtSqSum[2]*fr[8]-0.9375*nuVtSqSum[2]*fc[8]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.5*alphaDrSurf_l[4]*fUpwind_l[8]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[8]+0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]; 
  Ghat_l[8] = Gdiff_l[8]*rdv2+0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[8]+0.5*alphaDrSurf_l[0]*fUpwind_l[8]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[4]*fUpwind_l[5]; 

  Ghat_r[0] = Gdiff_r[0]*rdv2+0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = Gdiff_r[1]*rdv2+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = Gdiff_r[2]*rdv2+0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Ghat_r[3] = Gdiff_r[3]*rdv2+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Ghat_r[4] = Gdiff_r[4]*rdv2+0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Ghat_r[5] = Gdiff_r[5]*rdv2+0.5*alphaDrSurf_r[4]*fUpwind_r[8]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Ghat_r[6] = Gdiff_r[6]*rdv2+0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[7] = Gdiff_r[7]*rdv2+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[8]+0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]; 
  Ghat_r[8] = Gdiff_r[8]*rdv2+0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[8]+0.5*alphaDrSurf_r[0]*fUpwind_r[8]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[4]*fUpwind_r[5]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[5] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[6] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_l[4]*rdv2-0.7071067811865475*Ghat_r[4]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_l[5]*rdv2-0.7071067811865475*Ghat_r[5]*rdv2; 
  out[9] += 4.743416490252569*Gdiff2_r[0]*rdvSq4+4.743416490252569*Gdiff2_l[0]*rdvSq4-1.58113883008419*Ghat_r[0]*rdv2+1.58113883008419*Ghat_l[0]*rdv2; 
  out[10] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_l[6]*rdv2-0.7071067811865475*Ghat_r[6]*rdv2; 
  out[12] += 0.7071067811865475*Ghat_l[7]*rdv2-0.7071067811865475*Ghat_r[7]*rdv2; 
  out[13] += 1.224744871391589*Gdiff2_r[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2-1.224744871391589*Ghat_l[4]*rdv2; 
  out[14] += 1.224744871391589*Gdiff2_r[5]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2-1.224744871391589*Ghat_l[5]*rdv2; 
  out[15] += 4.743416490252569*Gdiff2_r[1]*rdvSq4+4.743416490252569*Gdiff2_l[1]*rdvSq4-1.58113883008419*Ghat_r[1]*rdv2+1.58113883008419*Ghat_l[1]*rdv2; 
  out[16] += 4.743416490252569*Gdiff2_r[2]*rdvSq4+4.743416490252569*Gdiff2_l[2]*rdvSq4-1.58113883008419*Ghat_r[2]*rdv2+1.58113883008419*Ghat_l[2]*rdv2; 
  out[17] += 1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2-1.224744871391589*Ghat_l[6]*rdv2; 
  out[18] += 1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2-1.224744871391589*Ghat_l[7]*rdv2; 
  out[19] += 4.743416490252569*Gdiff2_r[3]*rdvSq4+4.743416490252569*Gdiff2_l[3]*rdvSq4-1.58113883008419*Ghat_r[3]*rdv2+1.58113883008419*Ghat_l[3]*rdv2; 
  out[20] += 0.7071067811865475*Ghat_l[8]*rdv2-0.7071067811865475*Ghat_r[8]*rdv2; 
  out[21] += 4.743416490252569*Gdiff2_r[4]*rdvSq4+4.743416490252569*Gdiff2_l[4]*rdvSq4-1.58113883008419*Ghat_r[4]*rdv2+1.58113883008419*Ghat_l[4]*rdv2; 
  out[22] += 4.743416490252569*Gdiff2_r[5]*rdvSq4+4.743416490252569*Gdiff2_l[5]*rdvSq4-1.58113883008419*Ghat_r[5]*rdv2+1.58113883008419*Ghat_l[5]*rdv2; 
  out[23] += 1.224744871391589*Gdiff2_r[8]*rdvSq4-1.224744871391589*Gdiff2_l[8]*rdvSq4-1.224744871391589*Ghat_r[8]*rdv2-1.224744871391589*Ghat_l[8]*rdv2; 
  out[24] += 4.743416490252569*Gdiff2_r[6]*rdvSq4+4.743416490252569*Gdiff2_l[6]*rdvSq4-1.58113883008419*Ghat_r[6]*rdv2+1.58113883008419*Ghat_l[6]*rdv2; 
  out[25] += 4.743416490252569*Gdiff2_r[7]*rdvSq4+4.743416490252569*Gdiff2_l[7]*rdvSq4-1.58113883008419*Ghat_r[7]*rdv2+1.58113883008419*Ghat_l[7]*rdv2; 
  out[26] += 4.743416490252569*Gdiff2_r[8]*rdvSq4+4.743416490252569*Gdiff2_l[8]*rdvSq4-1.58113883008419*Ghat_r[8]*rdv2+1.58113883008419*Ghat_l[8]*rdv2; 
} 
