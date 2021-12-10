#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_1x2v_p2_surfvx_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvx_1x2v_ser_p2(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[3]:         cell-center coordinates. 
  // dxv[3]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[6]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[3]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[1]; 
  double rdvSq4 = 4.0/(dxv[1]*dxv[1]); 

  const double *sumNuUx = &nuUSum[0]; 

  double alphaDrSurf_l[8] = {0.0}; 
  alphaDrSurf_l[0] = 1.414213562373095*nuSum[0]*w[1]-0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_l[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*sumNuUx[1]-0.7071067811865475*dxv[1]*nuSum[1]; 
  alphaDrSurf_l[4] = (-1.414213562373095*sumNuUx[2])+1.414213562373095*w[1]*nuSum[2]-0.7071067811865475*dxv[1]*nuSum[2]; 

  double alphaDrSurf_r[8] = {0.0}; 
  alphaDrSurf_r[0] = 1.414213562373095*nuSum[0]*w[1]+0.7071067811865475*nuSum[0]*dxv[1]-1.414213562373095*sumNuUx[0]; 
  alphaDrSurf_r[1] = 1.414213562373095*nuSum[1]*w[1]-1.414213562373095*sumNuUx[1]+0.7071067811865475*dxv[1]*nuSum[1]; 
  alphaDrSurf_r[4] = (-1.414213562373095*sumNuUx[2])+1.414213562373095*w[1]*nuSum[2]+0.7071067811865475*dxv[1]*nuSum[2]; 

  double fUpwindQuad_l[9] = {0.0};
  double fUpwindQuad_r[9] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Gdrag_l[8] = {0.0}; 
  double Gdrag_r[8] = {0.0}; 

  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_1x2v_p2_surfvx_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_1x2v_p2_surfvx_quad_0(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[1] = ser_1x2v_p2_surfvx_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_1x2v_p2_surfvx_quad_1(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_1x2v_p2_surfvx_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_1x2v_p2_surfvx_quad_2(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_1x2v_p2_surfvx_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_1x2v_p2_surfvx_quad_3(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[4] = ser_1x2v_p2_surfvx_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_1x2v_p2_surfvx_quad_4(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_1x2v_p2_surfvx_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_1x2v_p2_surfvx_quad_5(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]-0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_1x2v_p2_surfvx_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_1x2v_p2_surfvx_quad_6(-1, fc); 
  } 
  if (0.5*alphaDrSurf_l[0]-0.5590169943749475*alphaDrSurf_l[4] < 0) { 
    fUpwindQuad_l[7] = ser_1x2v_p2_surfvx_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_1x2v_p2_surfvx_quad_7(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_l[4]+0.6708203932499369*alphaDrSurf_l[1]+0.5*alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_1x2v_p2_surfvx_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_1x2v_p2_surfvx_quad_8(-1, fc); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_1x2v_p2_surfvx_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_1x2v_p2_surfvx_quad_0(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[1] = ser_1x2v_p2_surfvx_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_1x2v_p2_surfvx_quad_1(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_1x2v_p2_surfvx_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_1x2v_p2_surfvx_quad_2(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_1x2v_p2_surfvx_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_1x2v_p2_surfvx_quad_3(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[4] = ser_1x2v_p2_surfvx_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_1x2v_p2_surfvx_quad_4(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_1x2v_p2_surfvx_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_1x2v_p2_surfvx_quad_5(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]-0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_1x2v_p2_surfvx_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_1x2v_p2_surfvx_quad_6(-1, fr); 
  } 
  if (0.5*alphaDrSurf_r[0]-0.5590169943749475*alphaDrSurf_r[4] < 0) { 
    fUpwindQuad_r[7] = ser_1x2v_p2_surfvx_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_1x2v_p2_surfvx_quad_7(-1, fr); 
  } 
  if (0.4472135954999579*alphaDrSurf_r[4]+0.6708203932499369*alphaDrSurf_r[1]+0.5*alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_1x2v_p2_surfvx_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_1x2v_p2_surfvx_quad_8(-1, fr); 
  } 
  fUpwind_l[0] = 0.154320987654321*fUpwindQuad_l[8]+0.2469135802469136*fUpwindQuad_l[7]+0.154320987654321*fUpwindQuad_l[6]+0.2469135802469136*fUpwindQuad_l[5]+0.3950617283950617*fUpwindQuad_l[4]+0.2469135802469136*fUpwindQuad_l[3]+0.154320987654321*fUpwindQuad_l[2]+0.2469135802469136*fUpwindQuad_l[1]+0.154320987654321*fUpwindQuad_l[0]; 
  fUpwind_l[1] = 0.2070433312499806*fUpwindQuad_l[8]-0.2070433312499806*fUpwindQuad_l[6]+0.3312693299999688*fUpwindQuad_l[5]-0.3312693299999688*fUpwindQuad_l[3]+0.2070433312499806*fUpwindQuad_l[2]-0.2070433312499806*fUpwindQuad_l[0]; 
  fUpwind_l[2] = 0.2070433312499806*fUpwindQuad_l[8]+0.3312693299999688*fUpwindQuad_l[7]+0.2070433312499806*fUpwindQuad_l[6]-0.2070433312499806*fUpwindQuad_l[2]-0.3312693299999688*fUpwindQuad_l[1]-0.2070433312499806*fUpwindQuad_l[0]; 
  fUpwind_l[3] = 0.2777777777777778*fUpwindQuad_l[8]-0.2777777777777778*fUpwindQuad_l[6]-0.2777777777777778*fUpwindQuad_l[2]+0.2777777777777778*fUpwindQuad_l[0]; 
  fUpwind_l[4] = 0.138028887499987*fUpwindQuad_l[8]-0.2760577749999741*fUpwindQuad_l[7]+0.138028887499987*fUpwindQuad_l[6]+0.2208462199999792*fUpwindQuad_l[5]-0.4416924399999584*fUpwindQuad_l[4]+0.2208462199999792*fUpwindQuad_l[3]+0.138028887499987*fUpwindQuad_l[2]-0.2760577749999741*fUpwindQuad_l[1]+0.138028887499987*fUpwindQuad_l[0]; 
  fUpwind_l[5] = 0.138028887499987*fUpwindQuad_l[8]+0.2208462199999792*fUpwindQuad_l[7]+0.138028887499987*fUpwindQuad_l[6]-0.2760577749999741*fUpwindQuad_l[5]-0.4416924399999584*fUpwindQuad_l[4]-0.2760577749999741*fUpwindQuad_l[3]+0.138028887499987*fUpwindQuad_l[2]+0.2208462199999792*fUpwindQuad_l[1]+0.138028887499987*fUpwindQuad_l[0]; 
  fUpwind_l[6] = 0.1851851851851853*fUpwindQuad_l[8]-0.3703703703703705*fUpwindQuad_l[7]+0.1851851851851853*fUpwindQuad_l[6]-0.1851851851851853*fUpwindQuad_l[2]+0.3703703703703705*fUpwindQuad_l[1]-0.1851851851851853*fUpwindQuad_l[0]; 
  fUpwind_l[7] = 0.1851851851851853*fUpwindQuad_l[8]-0.1851851851851853*fUpwindQuad_l[6]-0.3703703703703705*fUpwindQuad_l[5]+0.3703703703703705*fUpwindQuad_l[3]+0.1851851851851853*fUpwindQuad_l[2]-0.1851851851851853*fUpwindQuad_l[0]; 

  fUpwind_r[0] = 0.154320987654321*fUpwindQuad_r[8]+0.2469135802469136*fUpwindQuad_r[7]+0.154320987654321*fUpwindQuad_r[6]+0.2469135802469136*fUpwindQuad_r[5]+0.3950617283950617*fUpwindQuad_r[4]+0.2469135802469136*fUpwindQuad_r[3]+0.154320987654321*fUpwindQuad_r[2]+0.2469135802469136*fUpwindQuad_r[1]+0.154320987654321*fUpwindQuad_r[0]; 
  fUpwind_r[1] = 0.2070433312499806*fUpwindQuad_r[8]-0.2070433312499806*fUpwindQuad_r[6]+0.3312693299999688*fUpwindQuad_r[5]-0.3312693299999688*fUpwindQuad_r[3]+0.2070433312499806*fUpwindQuad_r[2]-0.2070433312499806*fUpwindQuad_r[0]; 
  fUpwind_r[2] = 0.2070433312499806*fUpwindQuad_r[8]+0.3312693299999688*fUpwindQuad_r[7]+0.2070433312499806*fUpwindQuad_r[6]-0.2070433312499806*fUpwindQuad_r[2]-0.3312693299999688*fUpwindQuad_r[1]-0.2070433312499806*fUpwindQuad_r[0]; 
  fUpwind_r[3] = 0.2777777777777778*fUpwindQuad_r[8]-0.2777777777777778*fUpwindQuad_r[6]-0.2777777777777778*fUpwindQuad_r[2]+0.2777777777777778*fUpwindQuad_r[0]; 
  fUpwind_r[4] = 0.138028887499987*fUpwindQuad_r[8]-0.2760577749999741*fUpwindQuad_r[7]+0.138028887499987*fUpwindQuad_r[6]+0.2208462199999792*fUpwindQuad_r[5]-0.4416924399999584*fUpwindQuad_r[4]+0.2208462199999792*fUpwindQuad_r[3]+0.138028887499987*fUpwindQuad_r[2]-0.2760577749999741*fUpwindQuad_r[1]+0.138028887499987*fUpwindQuad_r[0]; 
  fUpwind_r[5] = 0.138028887499987*fUpwindQuad_r[8]+0.2208462199999792*fUpwindQuad_r[7]+0.138028887499987*fUpwindQuad_r[6]-0.2760577749999741*fUpwindQuad_r[5]-0.4416924399999584*fUpwindQuad_r[4]-0.2760577749999741*fUpwindQuad_r[3]+0.138028887499987*fUpwindQuad_r[2]+0.2208462199999792*fUpwindQuad_r[1]+0.138028887499987*fUpwindQuad_r[0]; 
  fUpwind_r[6] = 0.1851851851851853*fUpwindQuad_r[8]-0.3703703703703705*fUpwindQuad_r[7]+0.1851851851851853*fUpwindQuad_r[6]-0.1851851851851853*fUpwindQuad_r[2]+0.3703703703703705*fUpwindQuad_r[1]-0.1851851851851853*fUpwindQuad_r[0]; 
  fUpwind_r[7] = 0.1851851851851853*fUpwindQuad_r[8]-0.1851851851851853*fUpwindQuad_r[6]-0.3703703703703705*fUpwindQuad_r[5]+0.3703703703703705*fUpwindQuad_r[3]+0.1851851851851853*fUpwindQuad_r[2]-0.1851851851851853*fUpwindQuad_r[0]; 

  Gdrag_l[0] = 0.5*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[1]*fUpwind_l[1]+0.5*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[4]+0.4472135954999579*fUpwind_l[1]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[1]+0.5*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.5000000000000001*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[1]*fUpwind_l[3]+0.5*alphaDrSurf_l[0]*fUpwind_l[2]; 
  Gdrag_l[3] = 0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[6]+0.4472135954999579*fUpwind_l[3]*alphaDrSurf_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[3]+0.5*alphaDrSurf_l[1]*fUpwind_l[2]; 
  Gdrag_l[4] = 0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[4]+0.5*alphaDrSurf_l[0]*fUpwind_l[4]+0.5*fUpwind_l[0]*alphaDrSurf_l[4]+0.4472135954999579*alphaDrSurf_l[1]*fUpwind_l[1]; 
  Gdrag_l[5] = 0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[5]; 
  Gdrag_l[6] = 0.31943828249997*alphaDrSurf_l[4]*fUpwind_l[6]+0.5*alphaDrSurf_l[0]*fUpwind_l[6]+0.5000000000000001*fUpwind_l[2]*alphaDrSurf_l[4]+0.447213595499958*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.4472135954999579*alphaDrSurf_l[4]*fUpwind_l[7]+0.5*alphaDrSurf_l[0]*fUpwind_l[7]+0.5000000000000001*alphaDrSurf_l[1]*fUpwind_l[5]; 

  Gdrag_r[0] = 0.5*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[1]*fUpwind_r[1]+0.5*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[4]+0.4472135954999579*fUpwind_r[1]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[1]+0.5*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.5000000000000001*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[1]*fUpwind_r[3]+0.5*alphaDrSurf_r[0]*fUpwind_r[2]; 
  Gdrag_r[3] = 0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[6]+0.4472135954999579*fUpwind_r[3]*alphaDrSurf_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[3]+0.5*alphaDrSurf_r[1]*fUpwind_r[2]; 
  Gdrag_r[4] = 0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[4]+0.5*alphaDrSurf_r[0]*fUpwind_r[4]+0.5*fUpwind_r[0]*alphaDrSurf_r[4]+0.4472135954999579*alphaDrSurf_r[1]*fUpwind_r[1]; 
  Gdrag_r[5] = 0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[5]; 
  Gdrag_r[6] = 0.31943828249997*alphaDrSurf_r[4]*fUpwind_r[6]+0.5*alphaDrSurf_r[0]*fUpwind_r[6]+0.5000000000000001*fUpwind_r[2]*alphaDrSurf_r[4]+0.447213595499958*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.4472135954999579*alphaDrSurf_r[4]*fUpwind_r[7]+0.5*alphaDrSurf_r[0]*fUpwind_r[7]+0.5000000000000001*alphaDrSurf_r[1]*fUpwind_r[5]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[6] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[8] += 1.58113883008419*Gdrag_r[0]*rdv2-1.58113883008419*Gdrag_l[0]*rdv2; 
  out[9] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[12] += 1.58113883008419*Gdrag_r[1]*rdv2-1.58113883008419*Gdrag_l[1]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[14] += 1.58113883008419*Gdrag_r[2]*rdv2-1.58113883008419*Gdrag_l[2]*rdv2; 
  out[15] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[16] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[18] += 1.58113883008419*Gdrag_r[3]*rdv2-1.58113883008419*Gdrag_l[3]*rdv2; 
  out[19] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
} 
