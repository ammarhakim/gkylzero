#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_4x_p1_surfx3_quad.h> 
#include <gkyl_basis_ser_4x_p1_upwind.h> 
GKYL_CU_DH void vlasov_poisson_extem_surfvy_1x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double dv3 = dxv[3], wv3 = w[3]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[2]; 
  const double *A2 = &vecA[4]; 
  double alpha[8] = {0.0}; 

  alpha[0] = -3.464101615137754*A1[1]*dx10*wv1; 
  alpha[2] = -1.0*A1[1]*dv1*dx10; 

  double fUpwindQuad_l[8] = {0.0};
  double fUpwindQuad_r[8] = {0.0};
  double fUpwind_l[8] = {0.0};;
  double fUpwind_r[8] = {0.0};
  double Ghat_l[8] = {0.0}; 
  double Ghat_r[8] = {0.0}; 

  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx3_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx3_quad_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_4x_p1_surfx3_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_4x_p1_surfx3_quad_0_l(fr); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx3_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx3_quad_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_4x_p1_surfx3_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_4x_p1_surfx3_quad_1_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx3_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx3_quad_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_4x_p1_surfx3_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_4x_p1_surfx3_quad_2_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx3_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx3_quad_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_4x_p1_surfx3_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_4x_p1_surfx3_quad_3_l(fr); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx3_quad_4_r(fl); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx3_quad_4_r(fc); 
  } else { 
    fUpwindQuad_l[4] = ser_4x_p1_surfx3_quad_4_l(fc); 
    fUpwindQuad_r[4] = ser_4x_p1_surfx3_quad_4_l(fr); 
  } 
  if (alpha[0]-alpha[2] > 0) { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx3_quad_5_r(fl); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx3_quad_5_r(fc); 
  } else { 
    fUpwindQuad_l[5] = ser_4x_p1_surfx3_quad_5_l(fc); 
    fUpwindQuad_r[5] = ser_4x_p1_surfx3_quad_5_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx3_quad_6_r(fl); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx3_quad_6_r(fc); 
  } else { 
    fUpwindQuad_l[6] = ser_4x_p1_surfx3_quad_6_l(fc); 
    fUpwindQuad_r[6] = ser_4x_p1_surfx3_quad_6_l(fr); 
  } 
  if (alpha[2]+alpha[0] > 0) { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx3_quad_7_r(fl); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx3_quad_7_r(fc); 
  } else { 
    fUpwindQuad_l[7] = ser_4x_p1_surfx3_quad_7_l(fc); 
    fUpwindQuad_r[7] = ser_4x_p1_surfx3_quad_7_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_4x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_4x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.3535533905932737*(alpha[2]*fUpwind_l[2]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.3535533905932737*(alpha[2]*fUpwind_l[4]+alpha[0]*fUpwind_l[1]); 
  Ghat_l[2] += 0.3535533905932737*(alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.3535533905932737*(alpha[2]*fUpwind_l[6]+alpha[0]*fUpwind_l[3]); 
  Ghat_l[4] += 0.3535533905932737*(alpha[0]*fUpwind_l[4]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[5] += 0.3535533905932737*(alpha[2]*fUpwind_l[7]+alpha[0]*fUpwind_l[5]); 
  Ghat_l[6] += 0.3535533905932737*(alpha[0]*fUpwind_l[6]+alpha[2]*fUpwind_l[3]); 
  Ghat_l[7] += 0.3535533905932737*(alpha[0]*fUpwind_l[7]+alpha[2]*fUpwind_l[5]); 

  Ghat_r[0] += 0.3535533905932737*(alpha[2]*fUpwind_r[2]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.3535533905932737*(alpha[2]*fUpwind_r[4]+alpha[0]*fUpwind_r[1]); 
  Ghat_r[2] += 0.3535533905932737*(alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.3535533905932737*(alpha[2]*fUpwind_r[6]+alpha[0]*fUpwind_r[3]); 
  Ghat_r[4] += 0.3535533905932737*(alpha[0]*fUpwind_r[4]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[5] += 0.3535533905932737*(alpha[2]*fUpwind_r[7]+alpha[0]*fUpwind_r[5]); 
  Ghat_r[6] += 0.3535533905932737*(alpha[0]*fUpwind_r[6]+alpha[2]*fUpwind_r[3]); 
  Ghat_r[7] += 0.3535533905932737*(alpha[0]*fUpwind_r[7]+alpha[2]*fUpwind_r[5]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv11; 
  out[6] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[7] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[8] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv11; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv11; 
  out[10] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 
  out[11] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv11; 
  out[12] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv11; 
  out[13] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv11; 
  out[14] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv11; 
  out[15] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv11; 

} 
