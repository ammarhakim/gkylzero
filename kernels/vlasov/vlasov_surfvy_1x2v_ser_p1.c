#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_3x_p1_surfx3_quad.h> 
#include <gkyl_basis_ser_3x_p1_upwind.h> 
GKYL_CU_DH void vlasov_surfvy_1x2v_ser_p1(const double *w, const double *dxv, const double *qmem, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // qmem:      q/m*EM fields.
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv11 = 2/dxv[2]; 
  const double dv1 = dxv[1], wv1 = w[1]; 
  const double dv2 = dxv[2], wv2 = w[2]; 
  const double *E1 = &qmem[2]; 
  const double *B2 = &qmem[10]; 

  double alpha[4] = {0.0}; 

  alpha[0] = 1.414213562373095*E1[0]-1.414213562373095*B2[0]*wv1; 
  alpha[1] = 1.414213562373095*E1[1]-1.414213562373095*B2[1]*wv1; 
  alpha[2] = -0.408248290463863*B2[0]*dv1; 
  alpha[3] = -0.408248290463863*B2[1]*dv1; 

  double fUpwindQuad_l[4] = {0.0};
  double fUpwindQuad_r[4] = {0.0};
  double fUpwind_l[4] = {0.0};;
  double fUpwind_r[4] = {0.0};
  double Ghat_l[4] = {0.0}; 
  double Ghat_r[4] = {0.0}; 

  if (alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx3_quad_0_r(fl); 
    fUpwindQuad_r[0] = ser_3x_p1_surfx3_quad_0_r(fc); 
  } else { 
    fUpwindQuad_l[0] = ser_3x_p1_surfx3_quad_0_l(fc); 
    fUpwindQuad_r[0] = ser_3x_p1_surfx3_quad_0_l(fr); 
  } 
  if ((-alpha[3])-alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx3_quad_1_r(fl); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx3_quad_1_r(fc); 
  } else { 
    fUpwindQuad_l[1] = ser_3x_p1_surfx3_quad_1_l(fc); 
    fUpwindQuad_r[1] = ser_3x_p1_surfx3_quad_1_l(fr); 
  } 
  if ((-alpha[3])+alpha[2]-alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx3_quad_2_r(fl); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx3_quad_2_r(fc); 
  } else { 
    fUpwindQuad_l[2] = ser_3x_p1_surfx3_quad_2_l(fc); 
    fUpwindQuad_r[2] = ser_3x_p1_surfx3_quad_2_l(fr); 
  } 
  if (alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 
    fUpwindQuad_l[3] = ser_3x_p1_surfx3_quad_3_r(fl); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx3_quad_3_r(fc); 
  } else { 
    fUpwindQuad_l[3] = ser_3x_p1_surfx3_quad_3_l(fc); 
    fUpwindQuad_r[3] = ser_3x_p1_surfx3_quad_3_l(fr); 
  } 

  // Project nodal basis back onto modal basis. 
  ser_3x_p1_upwind(fUpwindQuad_l, fUpwind_l); 
  ser_3x_p1_upwind(fUpwindQuad_r, fUpwind_r); 

  Ghat_l[0] += 0.5*(alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.5*(alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.5*(alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.5*(alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 

  Ghat_r[0] += 0.5*(alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.5*(alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.5*(alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.5*(alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv11; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv11; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv11; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv11; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv11; 
  out[5] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv11; 
  out[6] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv11; 
  out[7] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv11; 

} 
