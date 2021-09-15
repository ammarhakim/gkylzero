#include <gkyl_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvx_quad.h> 
GKYL_CU_DH void vlasov_poisson_extem_surfvx_2x3v_ser_p1(const double *w, const double *dxv, const double *fac_phi, const double *vecA, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w:         Cell-center coordinates.
  // dxv[NDIM]: Cell spacing.
  // fac_phi:   potential (scaled by appropriate factors).
  // vecA:      vector potential (scaled by appropriate factors). Unused in pure Vlasov-Poisson. 
  // fl/fc/fr:  Input Distribution function in left/center/right cells 
  // out:       Output distribution function in center cell 
  const double dv10 = 2/dxv[2]; 
  const double dv1 = dxv[2], wv1 = w[2]; 
  const double dv2 = dxv[3], wv2 = w[3]; 
  const double dv3 = dxv[4], wv3 = w[4]; 
  const double *phi = &fac_phi[0]; 
  const double dx10 = 2/dxv[0]; 
  const double dx11 = 2/dxv[1]; 
  const double *A0 = &vecA[0]; 
  const double *A1 = &vecA[4]; 
  const double *A2 = &vecA[8]; 
  double alpha[16] = {0.0}; 

  alpha[0] = 3.464101615137754*A2[1]*dx10*wv3-3.464101615137754*A0[2]*dx11*wv2+3.464101615137754*A1[1]*dx10*wv2-3.464101615137754*phi[1]*dx10; 
  alpha[1] = -3.464101615137754*A0[3]*dx11*wv2; 
  alpha[2] = 3.464101615137754*A2[3]*dx10*wv3+3.464101615137754*A1[3]*dx10*wv2-3.464101615137754*phi[3]*dx10; 
  alpha[3] = A1[1]*dv2*dx10-1.0*A0[2]*dv2*dx11; 
  alpha[4] = A2[1]*dv3*dx10; 
  alpha[6] = -1.0*A0[3]*dv2*dx11; 
  alpha[7] = A1[3]*dv2*dx10; 
  alpha[9] = A2[3]*dv3*dx10; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 

  if (alpha[9]+alpha[7]+alpha[6]-alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[0] = ser_2x3v_p1_surfvx_quad_0(1, fl); 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvx_quad_0(1, fc); 
  } else { 

    fUpwindQuad_l[0] = ser_2x3v_p1_surfvx_quad_0(-1, fc); 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvx_quad_0(-1, fr); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]-alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[1] = ser_2x3v_p1_surfvx_quad_1(1, fl); 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvx_quad_1(1, fc); 
  } else { 

    fUpwindQuad_l[1] = ser_2x3v_p1_surfvx_quad_1(-1, fc); 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvx_quad_1(-1, fr); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]-alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[2] = ser_2x3v_p1_surfvx_quad_2(1, fl); 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvx_quad_2(1, fc); 
  } else { 

    fUpwindQuad_l[2] = ser_2x3v_p1_surfvx_quad_2(-1, fc); 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvx_quad_2(-1, fr); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]-alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[3] = ser_2x3v_p1_surfvx_quad_3(1, fl); 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvx_quad_3(1, fc); 
  } else { 

    fUpwindQuad_l[3] = ser_2x3v_p1_surfvx_quad_3(-1, fc); 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvx_quad_3(-1, fr); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]-alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[4] = ser_2x3v_p1_surfvx_quad_4(1, fl); 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvx_quad_4(1, fc); 
  } else { 

    fUpwindQuad_l[4] = ser_2x3v_p1_surfvx_quad_4(-1, fc); 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvx_quad_4(-1, fr); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]-alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[5] = ser_2x3v_p1_surfvx_quad_5(1, fl); 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvx_quad_5(1, fc); 
  } else { 

    fUpwindQuad_l[5] = ser_2x3v_p1_surfvx_quad_5(-1, fc); 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvx_quad_5(-1, fr); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]-alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[6] = ser_2x3v_p1_surfvx_quad_6(1, fl); 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvx_quad_6(1, fc); 
  } else { 

    fUpwindQuad_l[6] = ser_2x3v_p1_surfvx_quad_6(-1, fc); 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvx_quad_6(-1, fr); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]-alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[7] = ser_2x3v_p1_surfvx_quad_7(1, fl); 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvx_quad_7(1, fc); 
  } else { 

    fUpwindQuad_l[7] = ser_2x3v_p1_surfvx_quad_7(-1, fc); 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvx_quad_7(-1, fr); 
  } 
  if ((-alpha[9])+alpha[7]+alpha[6]+alpha[4]-alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[8] = ser_2x3v_p1_surfvx_quad_8(1, fl); 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvx_quad_8(1, fc); 
  } else { 

    fUpwindQuad_l[8] = ser_2x3v_p1_surfvx_quad_8(-1, fc); 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvx_quad_8(-1, fr); 
  } 
  if ((-alpha[9])+alpha[7]-alpha[6]+alpha[4]-alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[9] = ser_2x3v_p1_surfvx_quad_9(1, fl); 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvx_quad_9(1, fc); 
  } else { 

    fUpwindQuad_l[9] = ser_2x3v_p1_surfvx_quad_9(-1, fc); 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvx_quad_9(-1, fr); 
  } 
  if (alpha[9]-alpha[7]+alpha[6]+alpha[4]-alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[10] = ser_2x3v_p1_surfvx_quad_10(1, fl); 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvx_quad_10(1, fc); 
  } else { 

    fUpwindQuad_l[10] = ser_2x3v_p1_surfvx_quad_10(-1, fc); 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvx_quad_10(-1, fr); 
  } 
  if (alpha[9]-alpha[7]-alpha[6]+alpha[4]-alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[11] = ser_2x3v_p1_surfvx_quad_11(1, fl); 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvx_quad_11(1, fc); 
  } else { 

    fUpwindQuad_l[11] = ser_2x3v_p1_surfvx_quad_11(-1, fc); 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvx_quad_11(-1, fr); 
  } 
  if ((-alpha[9])-alpha[7]-alpha[6]+alpha[4]+alpha[3]-alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[12] = ser_2x3v_p1_surfvx_quad_12(1, fl); 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvx_quad_12(1, fc); 
  } else { 

    fUpwindQuad_l[12] = ser_2x3v_p1_surfvx_quad_12(-1, fc); 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvx_quad_12(-1, fr); 
  } 
  if ((-alpha[9])-alpha[7]+alpha[6]+alpha[4]+alpha[3]-alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[13] = ser_2x3v_p1_surfvx_quad_13(1, fl); 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvx_quad_13(1, fc); 
  } else { 

    fUpwindQuad_l[13] = ser_2x3v_p1_surfvx_quad_13(-1, fc); 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvx_quad_13(-1, fr); 
  } 
  if (alpha[9]+alpha[7]-alpha[6]+alpha[4]+alpha[3]+alpha[2]-alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[14] = ser_2x3v_p1_surfvx_quad_14(1, fl); 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvx_quad_14(1, fc); 
  } else { 

    fUpwindQuad_l[14] = ser_2x3v_p1_surfvx_quad_14(-1, fc); 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvx_quad_14(-1, fr); 
  } 
  if (alpha[9]+alpha[7]+alpha[6]+alpha[4]+alpha[3]+alpha[2]+alpha[1]+alpha[0] > 0) { 

    fUpwindQuad_l[15] = ser_2x3v_p1_surfvx_quad_15(1, fl); 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvx_quad_15(1, fc); 
  } else { 

    fUpwindQuad_l[15] = ser_2x3v_p1_surfvx_quad_15(-1, fc); 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvx_quad_15(-1, fr); 
  } 
  fUpwind_l[0] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[1] = 0.25*(fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[2] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[3] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8])+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[4] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[5] = 0.25*(fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 
  fUpwind_l[6] = 0.25*(fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[7] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*(fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[8] = 0.25*(fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[9] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12])+fUpwindQuad_l[11]+fUpwindQuad_l[10]-1.0*(fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]-1.0*(fUpwindQuad_l[3]+fUpwindQuad_l[2])+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[10] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]+fUpwindQuad_l[13]+fUpwindQuad_l[12]-1.0*(fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]+fUpwindQuad_l[8]+fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]+fUpwindQuad_l[4])+fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]+fUpwindQuad_l[0]); 
  fUpwind_l[11] = 0.25*(fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*fUpwindQuad_l[8]+fUpwindQuad_l[7]-1.0*(fUpwindQuad_l[6]+fUpwindQuad_l[5])+fUpwindQuad_l[4]-1.0*fUpwindQuad_l[3]+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[12] = 0.25*(fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]+fUpwindQuad_l[11]-1.0*(fUpwindQuad_l[10]+fUpwindQuad_l[9])+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*(fUpwindQuad_l[4]+fUpwindQuad_l[3])+fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[13] = 0.25*(fUpwindQuad_l[15]-1.0*fUpwindQuad_l[14]+fUpwindQuad_l[13]-1.0*(fUpwindQuad_l[12]+fUpwindQuad_l[11])+fUpwindQuad_l[10]-1.0*fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*fUpwindQuad_l[7]+fUpwindQuad_l[6]-1.0*fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*fUpwindQuad_l[2]+fUpwindQuad_l[1]-1.0*fUpwindQuad_l[0]); 
  fUpwind_l[14] = 0.25*(fUpwindQuad_l[15]+fUpwindQuad_l[14]-1.0*(fUpwindQuad_l[13]+fUpwindQuad_l[12]+fUpwindQuad_l[11]+fUpwindQuad_l[10])+fUpwindQuad_l[9]+fUpwindQuad_l[8]-1.0*(fUpwindQuad_l[7]+fUpwindQuad_l[6])+fUpwindQuad_l[5]+fUpwindQuad_l[4]+fUpwindQuad_l[3]+fUpwindQuad_l[2]-1.0*(fUpwindQuad_l[1]+fUpwindQuad_l[0])); 
  fUpwind_l[15] = 0.25*(fUpwindQuad_l[15]-1.0*(fUpwindQuad_l[14]+fUpwindQuad_l[13])+fUpwindQuad_l[12]-1.0*fUpwindQuad_l[11]+fUpwindQuad_l[10]+fUpwindQuad_l[9]-1.0*(fUpwindQuad_l[8]+fUpwindQuad_l[7])+fUpwindQuad_l[6]+fUpwindQuad_l[5]-1.0*fUpwindQuad_l[4]+fUpwindQuad_l[3]-1.0*(fUpwindQuad_l[2]+fUpwindQuad_l[1])+fUpwindQuad_l[0]); 

  fUpwind_r[0] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[1] = 0.25*(fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[2] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[3] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8])+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[4] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[5] = 0.25*(fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 
  fUpwind_r[6] = 0.25*(fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[7] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*(fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[8] = 0.25*(fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[9] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12])+fUpwindQuad_r[11]+fUpwindQuad_r[10]-1.0*(fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]-1.0*(fUpwindQuad_r[3]+fUpwindQuad_r[2])+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[10] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]+fUpwindQuad_r[13]+fUpwindQuad_r[12]-1.0*(fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]+fUpwindQuad_r[8]+fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]+fUpwindQuad_r[4])+fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]+fUpwindQuad_r[0]); 
  fUpwind_r[11] = 0.25*(fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*fUpwindQuad_r[8]+fUpwindQuad_r[7]-1.0*(fUpwindQuad_r[6]+fUpwindQuad_r[5])+fUpwindQuad_r[4]-1.0*fUpwindQuad_r[3]+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[12] = 0.25*(fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]+fUpwindQuad_r[11]-1.0*(fUpwindQuad_r[10]+fUpwindQuad_r[9])+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*(fUpwindQuad_r[4]+fUpwindQuad_r[3])+fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[13] = 0.25*(fUpwindQuad_r[15]-1.0*fUpwindQuad_r[14]+fUpwindQuad_r[13]-1.0*(fUpwindQuad_r[12]+fUpwindQuad_r[11])+fUpwindQuad_r[10]-1.0*fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*fUpwindQuad_r[7]+fUpwindQuad_r[6]-1.0*fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*fUpwindQuad_r[2]+fUpwindQuad_r[1]-1.0*fUpwindQuad_r[0]); 
  fUpwind_r[14] = 0.25*(fUpwindQuad_r[15]+fUpwindQuad_r[14]-1.0*(fUpwindQuad_r[13]+fUpwindQuad_r[12]+fUpwindQuad_r[11]+fUpwindQuad_r[10])+fUpwindQuad_r[9]+fUpwindQuad_r[8]-1.0*(fUpwindQuad_r[7]+fUpwindQuad_r[6])+fUpwindQuad_r[5]+fUpwindQuad_r[4]+fUpwindQuad_r[3]+fUpwindQuad_r[2]-1.0*(fUpwindQuad_r[1]+fUpwindQuad_r[0])); 
  fUpwind_r[15] = 0.25*(fUpwindQuad_r[15]-1.0*(fUpwindQuad_r[14]+fUpwindQuad_r[13])+fUpwindQuad_r[12]-1.0*fUpwindQuad_r[11]+fUpwindQuad_r[10]+fUpwindQuad_r[9]-1.0*(fUpwindQuad_r[8]+fUpwindQuad_r[7])+fUpwindQuad_r[6]+fUpwindQuad_r[5]-1.0*fUpwindQuad_r[4]+fUpwindQuad_r[3]-1.0*(fUpwindQuad_r[2]+fUpwindQuad_r[1])+fUpwindQuad_r[0]); 

  Ghat_l[0] += 0.25*(alpha[9]*fUpwind_l[9]+alpha[7]*fUpwind_l[7]+alpha[6]*fUpwind_l[6]+alpha[4]*fUpwind_l[4]+alpha[3]*fUpwind_l[3]+alpha[2]*fUpwind_l[2]+alpha[1]*fUpwind_l[1]+alpha[0]*fUpwind_l[0]); 
  Ghat_l[1] += 0.25*(alpha[9]*fUpwind_l[12]+alpha[7]*fUpwind_l[11]+alpha[4]*fUpwind_l[8]+alpha[3]*fUpwind_l[6]+fUpwind_l[3]*alpha[6]+alpha[2]*fUpwind_l[5]+alpha[0]*fUpwind_l[1]+fUpwind_l[0]*alpha[1]); 
  Ghat_l[2] += 0.25*(alpha[6]*fUpwind_l[11]+alpha[4]*fUpwind_l[9]+fUpwind_l[4]*alpha[9]+alpha[3]*fUpwind_l[7]+fUpwind_l[3]*alpha[7]+alpha[1]*fUpwind_l[5]+alpha[0]*fUpwind_l[2]+fUpwind_l[0]*alpha[2]); 
  Ghat_l[3] += 0.25*(alpha[9]*fUpwind_l[14]+alpha[4]*fUpwind_l[10]+alpha[2]*fUpwind_l[7]+fUpwind_l[2]*alpha[7]+alpha[1]*fUpwind_l[6]+fUpwind_l[1]*alpha[6]+alpha[0]*fUpwind_l[3]+fUpwind_l[0]*alpha[3]); 
  Ghat_l[4] += 0.25*(alpha[7]*fUpwind_l[14]+alpha[6]*fUpwind_l[13]+alpha[3]*fUpwind_l[10]+alpha[2]*fUpwind_l[9]+fUpwind_l[2]*alpha[9]+alpha[1]*fUpwind_l[8]+alpha[0]*fUpwind_l[4]+fUpwind_l[0]*alpha[4]); 
  Ghat_l[5] += 0.25*(alpha[4]*fUpwind_l[12]+alpha[3]*fUpwind_l[11]+fUpwind_l[8]*alpha[9]+alpha[6]*fUpwind_l[7]+fUpwind_l[6]*alpha[7]+alpha[0]*fUpwind_l[5]+alpha[1]*fUpwind_l[2]+fUpwind_l[1]*alpha[2]); 
  Ghat_l[6] += 0.25*(alpha[9]*fUpwind_l[15]+alpha[4]*fUpwind_l[13]+alpha[2]*fUpwind_l[11]+fUpwind_l[5]*alpha[7]+alpha[0]*fUpwind_l[6]+fUpwind_l[0]*alpha[6]+alpha[1]*fUpwind_l[3]+fUpwind_l[1]*alpha[3]); 
  Ghat_l[7] += 0.25*(alpha[4]*fUpwind_l[14]+alpha[1]*fUpwind_l[11]+alpha[9]*fUpwind_l[10]+alpha[0]*fUpwind_l[7]+fUpwind_l[0]*alpha[7]+fUpwind_l[5]*alpha[6]+alpha[2]*fUpwind_l[3]+fUpwind_l[2]*alpha[3]); 
  Ghat_l[8] += 0.25*(alpha[7]*fUpwind_l[15]+alpha[3]*fUpwind_l[13]+alpha[2]*fUpwind_l[12]+alpha[6]*fUpwind_l[10]+fUpwind_l[5]*alpha[9]+alpha[0]*fUpwind_l[8]+alpha[1]*fUpwind_l[4]+fUpwind_l[1]*alpha[4]); 
  Ghat_l[9] += 0.25*(alpha[6]*fUpwind_l[15]+alpha[3]*fUpwind_l[14]+alpha[1]*fUpwind_l[12]+alpha[7]*fUpwind_l[10]+alpha[0]*fUpwind_l[9]+fUpwind_l[0]*alpha[9]+alpha[2]*fUpwind_l[4]+fUpwind_l[2]*alpha[4]); 
  Ghat_l[10] += 0.25*(alpha[2]*fUpwind_l[14]+alpha[1]*fUpwind_l[13]+alpha[0]*fUpwind_l[10]+alpha[7]*fUpwind_l[9]+fUpwind_l[7]*alpha[9]+alpha[6]*fUpwind_l[8]+alpha[3]*fUpwind_l[4]+fUpwind_l[3]*alpha[4]); 
  Ghat_l[11] += 0.25*(alpha[4]*fUpwind_l[15]+alpha[9]*fUpwind_l[13]+alpha[0]*fUpwind_l[11]+alpha[1]*fUpwind_l[7]+fUpwind_l[1]*alpha[7]+alpha[2]*fUpwind_l[6]+fUpwind_l[2]*alpha[6]+alpha[3]*fUpwind_l[5]); 
  Ghat_l[12] += 0.25*(alpha[3]*fUpwind_l[15]+alpha[6]*fUpwind_l[14]+alpha[7]*fUpwind_l[13]+alpha[0]*fUpwind_l[12]+alpha[1]*fUpwind_l[9]+fUpwind_l[1]*alpha[9]+alpha[2]*fUpwind_l[8]+alpha[4]*fUpwind_l[5]); 
  Ghat_l[13] += 0.25*(alpha[2]*fUpwind_l[15]+alpha[0]*fUpwind_l[13]+alpha[7]*fUpwind_l[12]+alpha[9]*fUpwind_l[11]+alpha[1]*fUpwind_l[10]+alpha[3]*fUpwind_l[8]+alpha[4]*fUpwind_l[6]+fUpwind_l[4]*alpha[6]); 
  Ghat_l[14] += 0.25*(alpha[1]*fUpwind_l[15]+alpha[0]*fUpwind_l[14]+alpha[6]*fUpwind_l[12]+alpha[2]*fUpwind_l[10]+alpha[3]*fUpwind_l[9]+fUpwind_l[3]*alpha[9]+alpha[4]*fUpwind_l[7]+fUpwind_l[4]*alpha[7]); 
  Ghat_l[15] += 0.25*(alpha[0]*fUpwind_l[15]+alpha[1]*fUpwind_l[14]+alpha[2]*fUpwind_l[13]+alpha[3]*fUpwind_l[12]+alpha[4]*fUpwind_l[11]+alpha[6]*fUpwind_l[9]+fUpwind_l[6]*alpha[9]+alpha[7]*fUpwind_l[8]); 

  Ghat_r[0] += 0.25*(alpha[9]*fUpwind_r[9]+alpha[7]*fUpwind_r[7]+alpha[6]*fUpwind_r[6]+alpha[4]*fUpwind_r[4]+alpha[3]*fUpwind_r[3]+alpha[2]*fUpwind_r[2]+alpha[1]*fUpwind_r[1]+alpha[0]*fUpwind_r[0]); 
  Ghat_r[1] += 0.25*(alpha[9]*fUpwind_r[12]+alpha[7]*fUpwind_r[11]+alpha[4]*fUpwind_r[8]+alpha[3]*fUpwind_r[6]+fUpwind_r[3]*alpha[6]+alpha[2]*fUpwind_r[5]+alpha[0]*fUpwind_r[1]+fUpwind_r[0]*alpha[1]); 
  Ghat_r[2] += 0.25*(alpha[6]*fUpwind_r[11]+alpha[4]*fUpwind_r[9]+fUpwind_r[4]*alpha[9]+alpha[3]*fUpwind_r[7]+fUpwind_r[3]*alpha[7]+alpha[1]*fUpwind_r[5]+alpha[0]*fUpwind_r[2]+fUpwind_r[0]*alpha[2]); 
  Ghat_r[3] += 0.25*(alpha[9]*fUpwind_r[14]+alpha[4]*fUpwind_r[10]+alpha[2]*fUpwind_r[7]+fUpwind_r[2]*alpha[7]+alpha[1]*fUpwind_r[6]+fUpwind_r[1]*alpha[6]+alpha[0]*fUpwind_r[3]+fUpwind_r[0]*alpha[3]); 
  Ghat_r[4] += 0.25*(alpha[7]*fUpwind_r[14]+alpha[6]*fUpwind_r[13]+alpha[3]*fUpwind_r[10]+alpha[2]*fUpwind_r[9]+fUpwind_r[2]*alpha[9]+alpha[1]*fUpwind_r[8]+alpha[0]*fUpwind_r[4]+fUpwind_r[0]*alpha[4]); 
  Ghat_r[5] += 0.25*(alpha[4]*fUpwind_r[12]+alpha[3]*fUpwind_r[11]+fUpwind_r[8]*alpha[9]+alpha[6]*fUpwind_r[7]+fUpwind_r[6]*alpha[7]+alpha[0]*fUpwind_r[5]+alpha[1]*fUpwind_r[2]+fUpwind_r[1]*alpha[2]); 
  Ghat_r[6] += 0.25*(alpha[9]*fUpwind_r[15]+alpha[4]*fUpwind_r[13]+alpha[2]*fUpwind_r[11]+fUpwind_r[5]*alpha[7]+alpha[0]*fUpwind_r[6]+fUpwind_r[0]*alpha[6]+alpha[1]*fUpwind_r[3]+fUpwind_r[1]*alpha[3]); 
  Ghat_r[7] += 0.25*(alpha[4]*fUpwind_r[14]+alpha[1]*fUpwind_r[11]+alpha[9]*fUpwind_r[10]+alpha[0]*fUpwind_r[7]+fUpwind_r[0]*alpha[7]+fUpwind_r[5]*alpha[6]+alpha[2]*fUpwind_r[3]+fUpwind_r[2]*alpha[3]); 
  Ghat_r[8] += 0.25*(alpha[7]*fUpwind_r[15]+alpha[3]*fUpwind_r[13]+alpha[2]*fUpwind_r[12]+alpha[6]*fUpwind_r[10]+fUpwind_r[5]*alpha[9]+alpha[0]*fUpwind_r[8]+alpha[1]*fUpwind_r[4]+fUpwind_r[1]*alpha[4]); 
  Ghat_r[9] += 0.25*(alpha[6]*fUpwind_r[15]+alpha[3]*fUpwind_r[14]+alpha[1]*fUpwind_r[12]+alpha[7]*fUpwind_r[10]+alpha[0]*fUpwind_r[9]+fUpwind_r[0]*alpha[9]+alpha[2]*fUpwind_r[4]+fUpwind_r[2]*alpha[4]); 
  Ghat_r[10] += 0.25*(alpha[2]*fUpwind_r[14]+alpha[1]*fUpwind_r[13]+alpha[0]*fUpwind_r[10]+alpha[7]*fUpwind_r[9]+fUpwind_r[7]*alpha[9]+alpha[6]*fUpwind_r[8]+alpha[3]*fUpwind_r[4]+fUpwind_r[3]*alpha[4]); 
  Ghat_r[11] += 0.25*(alpha[4]*fUpwind_r[15]+alpha[9]*fUpwind_r[13]+alpha[0]*fUpwind_r[11]+alpha[1]*fUpwind_r[7]+fUpwind_r[1]*alpha[7]+alpha[2]*fUpwind_r[6]+fUpwind_r[2]*alpha[6]+alpha[3]*fUpwind_r[5]); 
  Ghat_r[12] += 0.25*(alpha[3]*fUpwind_r[15]+alpha[6]*fUpwind_r[14]+alpha[7]*fUpwind_r[13]+alpha[0]*fUpwind_r[12]+alpha[1]*fUpwind_r[9]+fUpwind_r[1]*alpha[9]+alpha[2]*fUpwind_r[8]+alpha[4]*fUpwind_r[5]); 
  Ghat_r[13] += 0.25*(alpha[2]*fUpwind_r[15]+alpha[0]*fUpwind_r[13]+alpha[7]*fUpwind_r[12]+alpha[9]*fUpwind_r[11]+alpha[1]*fUpwind_r[10]+alpha[3]*fUpwind_r[8]+alpha[4]*fUpwind_r[6]+fUpwind_r[4]*alpha[6]); 
  Ghat_r[14] += 0.25*(alpha[1]*fUpwind_r[15]+alpha[0]*fUpwind_r[14]+alpha[6]*fUpwind_r[12]+alpha[2]*fUpwind_r[10]+alpha[3]*fUpwind_r[9]+fUpwind_r[3]*alpha[9]+alpha[4]*fUpwind_r[7]+fUpwind_r[4]*alpha[7]); 
  Ghat_r[15] += 0.25*(alpha[0]*fUpwind_r[15]+alpha[1]*fUpwind_r[14]+alpha[2]*fUpwind_r[13]+alpha[3]*fUpwind_r[12]+alpha[4]*fUpwind_r[11]+alpha[6]*fUpwind_r[9]+fUpwind_r[6]*alpha[9]+alpha[7]*fUpwind_r[8]); 

  out[0] += (0.7071067811865475*Ghat_l[0]-0.7071067811865475*Ghat_r[0])*dv10; 
  out[1] += (0.7071067811865475*Ghat_l[1]-0.7071067811865475*Ghat_r[1])*dv10; 
  out[2] += (0.7071067811865475*Ghat_l[2]-0.7071067811865475*Ghat_r[2])*dv10; 
  out[3] += -1.224744871391589*(Ghat_r[0]+Ghat_l[0])*dv10; 
  out[4] += (0.7071067811865475*Ghat_l[3]-0.7071067811865475*Ghat_r[3])*dv10; 
  out[5] += (0.7071067811865475*Ghat_l[4]-0.7071067811865475*Ghat_r[4])*dv10; 
  out[6] += (0.7071067811865475*Ghat_l[5]-0.7071067811865475*Ghat_r[5])*dv10; 
  out[7] += -1.224744871391589*(Ghat_r[1]+Ghat_l[1])*dv10; 
  out[8] += -1.224744871391589*(Ghat_r[2]+Ghat_l[2])*dv10; 
  out[9] += (0.7071067811865475*Ghat_l[6]-0.7071067811865475*Ghat_r[6])*dv10; 
  out[10] += (0.7071067811865475*Ghat_l[7]-0.7071067811865475*Ghat_r[7])*dv10; 
  out[11] += -1.224744871391589*(Ghat_r[3]+Ghat_l[3])*dv10; 
  out[12] += (0.7071067811865475*Ghat_l[8]-0.7071067811865475*Ghat_r[8])*dv10; 
  out[13] += (0.7071067811865475*Ghat_l[9]-0.7071067811865475*Ghat_r[9])*dv10; 
  out[14] += -1.224744871391589*(Ghat_r[4]+Ghat_l[4])*dv10; 
  out[15] += (0.7071067811865475*Ghat_l[10]-0.7071067811865475*Ghat_r[10])*dv10; 
  out[16] += -1.224744871391589*(Ghat_r[5]+Ghat_l[5])*dv10; 
  out[17] += (0.7071067811865475*Ghat_l[11]-0.7071067811865475*Ghat_r[11])*dv10; 
  out[18] += -1.224744871391589*(Ghat_r[6]+Ghat_l[6])*dv10; 
  out[19] += -1.224744871391589*(Ghat_r[7]+Ghat_l[7])*dv10; 
  out[20] += (0.7071067811865475*Ghat_l[12]-0.7071067811865475*Ghat_r[12])*dv10; 
  out[21] += -1.224744871391589*(Ghat_r[8]+Ghat_l[8])*dv10; 
  out[22] += -1.224744871391589*(Ghat_r[9]+Ghat_l[9])*dv10; 
  out[23] += (0.7071067811865475*Ghat_l[13]-0.7071067811865475*Ghat_r[13])*dv10; 
  out[24] += (0.7071067811865475*Ghat_l[14]-0.7071067811865475*Ghat_r[14])*dv10; 
  out[25] += -1.224744871391589*(Ghat_r[10]+Ghat_l[10])*dv10; 
  out[26] += -1.224744871391589*(Ghat_r[11]+Ghat_l[11])*dv10; 
  out[27] += -1.224744871391589*(Ghat_r[12]+Ghat_l[12])*dv10; 
  out[28] += (0.7071067811865475*Ghat_l[15]-0.7071067811865475*Ghat_r[15])*dv10; 
  out[29] += -1.224744871391589*(Ghat_r[13]+Ghat_l[13])*dv10; 
  out[30] += -1.224744871391589*(Ghat_r[14]+Ghat_l[14])*dv10; 
  out[31] += -1.224744871391589*(Ghat_r[15]+Ghat_l[15])*dv10; 

} 
