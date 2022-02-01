#include <gkyl_lbo_vlasov_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvy_quad.h> 
GKYL_CU_DH void lbo_vlasov_drag_surfvy_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[3]; 
  double rdvSq4 = 4.0/(dxv[3]*dxv[3]); 

  const double *sumNuUy = &nuUSum[4]; 

  double alphaDrSurf_l[16] = {0.0}; 
  alphaDrSurf_l[0] = 2.0*nuSum[0]*w[3]-1.0*nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_l[1] = 2.0*nuSum[1]*w[3]-1.0*nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_l[2] = 2.0*nuSum[2]*w[3]-1.0*nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_l[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]-1.0*dxv[3]*nuSum[3]; 

  double alphaDrSurf_r[16] = {0.0}; 
  alphaDrSurf_r[0] = 2.0*nuSum[0]*w[3]+nuSum[0]*dxv[3]-2.0*sumNuUy[0]; 
  alphaDrSurf_r[1] = 2.0*nuSum[1]*w[3]+nuSum[1]*dxv[3]-2.0*sumNuUy[1]; 
  alphaDrSurf_r[2] = 2.0*nuSum[2]*w[3]+nuSum[2]*dxv[3]-2.0*sumNuUy[2]; 
  alphaDrSurf_r[5] = 2.0*nuSum[3]*w[3]-2.0*sumNuUy[3]+dxv[3]*nuSum[3]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Gdrag_l[16] = {0.0}; 
  double Gdrag_r[16] = {0.0}; 

  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[0] = ser_2x3v_p1_surfvy_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x3v_p1_surfvy_quad_0(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[1] = ser_2x3v_p1_surfvy_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x3v_p1_surfvy_quad_1(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[2] = ser_2x3v_p1_surfvy_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x3v_p1_surfvy_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[3] = ser_2x3v_p1_surfvy_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x3v_p1_surfvy_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[4] = ser_2x3v_p1_surfvy_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_2x3v_p1_surfvy_quad_4(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[5] = ser_2x3v_p1_surfvy_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_2x3v_p1_surfvy_quad_5(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[6] = ser_2x3v_p1_surfvy_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_2x3v_p1_surfvy_quad_6(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[7] = ser_2x3v_p1_surfvy_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_2x3v_p1_surfvy_quad_7(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[8] = ser_2x3v_p1_surfvy_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_2x3v_p1_surfvy_quad_8(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[9] = ser_2x3v_p1_surfvy_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = ser_2x3v_p1_surfvy_quad_9(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[10] = ser_2x3v_p1_surfvy_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = ser_2x3v_p1_surfvy_quad_10(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[11] = ser_2x3v_p1_surfvy_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = ser_2x3v_p1_surfvy_quad_11(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[12] = ser_2x3v_p1_surfvy_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = ser_2x3v_p1_surfvy_quad_12(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[13] = ser_2x3v_p1_surfvy_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = ser_2x3v_p1_surfvy_quad_13(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[14] = ser_2x3v_p1_surfvy_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = ser_2x3v_p1_surfvy_quad_14(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] < 0) { 
    fUpwindQuad_l[15] = ser_2x3v_p1_surfvy_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = ser_2x3v_p1_surfvy_quad_15(-1, fc); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvy_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvy_quad_0(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvy_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvy_quad_1(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvy_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvy_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvy_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvy_quad_3(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvy_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvy_quad_4(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvy_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvy_quad_5(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvy_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvy_quad_6(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvy_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvy_quad_7(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvy_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvy_quad_8(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvy_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvy_quad_9(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvy_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvy_quad_10(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvy_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvy_quad_11(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvy_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvy_quad_12(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvy_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvy_quad_13(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvy_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvy_quad_14(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] < 0) { 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvy_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvy_quad_15(-1, fr); 
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

  Gdrag_l[0] = 0.25*alphaDrSurf_l[5]*fUpwind_l[5]+0.25*alphaDrSurf_l[2]*fUpwind_l[2]+0.25*alphaDrSurf_l[1]*fUpwind_l[1]+0.25*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Gdrag_l[1] = 0.25*alphaDrSurf_l[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Gdrag_l[2] = 0.25*alphaDrSurf_l[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Gdrag_l[3] = 0.25*alphaDrSurf_l[5]*fUpwind_l[11]+0.25*alphaDrSurf_l[2]*fUpwind_l[7]+0.25*alphaDrSurf_l[1]*fUpwind_l[6]+0.25*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Gdrag_l[4] = 0.25*alphaDrSurf_l[5]*fUpwind_l[12]+0.25*alphaDrSurf_l[2]*fUpwind_l[9]+0.25*alphaDrSurf_l[1]*fUpwind_l[8]+0.25*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Gdrag_l[5] = 0.25*alphaDrSurf_l[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Gdrag_l[6] = 0.25*alphaDrSurf_l[2]*fUpwind_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[7]+0.25*alphaDrSurf_l[0]*fUpwind_l[6]+0.25*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Gdrag_l[7] = 0.25*alphaDrSurf_l[1]*fUpwind_l[11]+0.25*alphaDrSurf_l[0]*fUpwind_l[7]+0.25*alphaDrSurf_l[5]*fUpwind_l[6]+0.25*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Gdrag_l[8] = 0.25*alphaDrSurf_l[2]*fUpwind_l[12]+0.25*alphaDrSurf_l[5]*fUpwind_l[9]+0.25*alphaDrSurf_l[0]*fUpwind_l[8]+0.25*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Gdrag_l[9] = 0.25*alphaDrSurf_l[1]*fUpwind_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[9]+0.25*alphaDrSurf_l[5]*fUpwind_l[8]+0.25*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Gdrag_l[10] = 0.25*alphaDrSurf_l[5]*fUpwind_l[15]+0.25*alphaDrSurf_l[2]*fUpwind_l[14]+0.25*alphaDrSurf_l[1]*fUpwind_l[13]+0.25*alphaDrSurf_l[0]*fUpwind_l[10]; 
  Gdrag_l[11] = 0.25*alphaDrSurf_l[0]*fUpwind_l[11]+0.25*alphaDrSurf_l[1]*fUpwind_l[7]+0.25*alphaDrSurf_l[2]*fUpwind_l[6]+0.25*fUpwind_l[3]*alphaDrSurf_l[5]; 
  Gdrag_l[12] = 0.25*alphaDrSurf_l[0]*fUpwind_l[12]+0.25*alphaDrSurf_l[1]*fUpwind_l[9]+0.25*alphaDrSurf_l[2]*fUpwind_l[8]+0.25*fUpwind_l[4]*alphaDrSurf_l[5]; 
  Gdrag_l[13] = 0.25*alphaDrSurf_l[2]*fUpwind_l[15]+0.25*alphaDrSurf_l[5]*fUpwind_l[14]+0.25*alphaDrSurf_l[0]*fUpwind_l[13]+0.25*alphaDrSurf_l[1]*fUpwind_l[10]; 
  Gdrag_l[14] = 0.25*alphaDrSurf_l[1]*fUpwind_l[15]+0.25*alphaDrSurf_l[0]*fUpwind_l[14]+0.25*alphaDrSurf_l[5]*fUpwind_l[13]+0.25*alphaDrSurf_l[2]*fUpwind_l[10]; 
  Gdrag_l[15] = 0.25*alphaDrSurf_l[0]*fUpwind_l[15]+0.25*alphaDrSurf_l[1]*fUpwind_l[14]+0.25*alphaDrSurf_l[2]*fUpwind_l[13]+0.25*alphaDrSurf_l[5]*fUpwind_l[10]; 

  Gdrag_r[0] = 0.25*alphaDrSurf_r[5]*fUpwind_r[5]+0.25*alphaDrSurf_r[2]*fUpwind_r[2]+0.25*alphaDrSurf_r[1]*fUpwind_r[1]+0.25*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Gdrag_r[1] = 0.25*alphaDrSurf_r[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Gdrag_r[2] = 0.25*alphaDrSurf_r[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Gdrag_r[3] = 0.25*alphaDrSurf_r[5]*fUpwind_r[11]+0.25*alphaDrSurf_r[2]*fUpwind_r[7]+0.25*alphaDrSurf_r[1]*fUpwind_r[6]+0.25*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Gdrag_r[4] = 0.25*alphaDrSurf_r[5]*fUpwind_r[12]+0.25*alphaDrSurf_r[2]*fUpwind_r[9]+0.25*alphaDrSurf_r[1]*fUpwind_r[8]+0.25*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Gdrag_r[5] = 0.25*alphaDrSurf_r[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Gdrag_r[6] = 0.25*alphaDrSurf_r[2]*fUpwind_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[7]+0.25*alphaDrSurf_r[0]*fUpwind_r[6]+0.25*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Gdrag_r[7] = 0.25*alphaDrSurf_r[1]*fUpwind_r[11]+0.25*alphaDrSurf_r[0]*fUpwind_r[7]+0.25*alphaDrSurf_r[5]*fUpwind_r[6]+0.25*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Gdrag_r[8] = 0.25*alphaDrSurf_r[2]*fUpwind_r[12]+0.25*alphaDrSurf_r[5]*fUpwind_r[9]+0.25*alphaDrSurf_r[0]*fUpwind_r[8]+0.25*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Gdrag_r[9] = 0.25*alphaDrSurf_r[1]*fUpwind_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[9]+0.25*alphaDrSurf_r[5]*fUpwind_r[8]+0.25*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Gdrag_r[10] = 0.25*alphaDrSurf_r[5]*fUpwind_r[15]+0.25*alphaDrSurf_r[2]*fUpwind_r[14]+0.25*alphaDrSurf_r[1]*fUpwind_r[13]+0.25*alphaDrSurf_r[0]*fUpwind_r[10]; 
  Gdrag_r[11] = 0.25*alphaDrSurf_r[0]*fUpwind_r[11]+0.25*alphaDrSurf_r[1]*fUpwind_r[7]+0.25*alphaDrSurf_r[2]*fUpwind_r[6]+0.25*fUpwind_r[3]*alphaDrSurf_r[5]; 
  Gdrag_r[12] = 0.25*alphaDrSurf_r[0]*fUpwind_r[12]+0.25*alphaDrSurf_r[1]*fUpwind_r[9]+0.25*alphaDrSurf_r[2]*fUpwind_r[8]+0.25*fUpwind_r[4]*alphaDrSurf_r[5]; 
  Gdrag_r[13] = 0.25*alphaDrSurf_r[2]*fUpwind_r[15]+0.25*alphaDrSurf_r[5]*fUpwind_r[14]+0.25*alphaDrSurf_r[0]*fUpwind_r[13]+0.25*alphaDrSurf_r[1]*fUpwind_r[10]; 
  Gdrag_r[14] = 0.25*alphaDrSurf_r[1]*fUpwind_r[15]+0.25*alphaDrSurf_r[0]*fUpwind_r[14]+0.25*alphaDrSurf_r[5]*fUpwind_r[13]+0.25*alphaDrSurf_r[2]*fUpwind_r[10]; 
  Gdrag_r[15] = 0.25*alphaDrSurf_r[0]*fUpwind_r[15]+0.25*alphaDrSurf_r[1]*fUpwind_r[14]+0.25*alphaDrSurf_r[2]*fUpwind_r[13]+0.25*alphaDrSurf_r[5]*fUpwind_r[10]; 

  out[0] += 0.7071067811865475*Gdrag_r[0]*rdv2-0.7071067811865475*Gdrag_l[0]*rdv2; 
  out[1] += 0.7071067811865475*Gdrag_r[1]*rdv2-0.7071067811865475*Gdrag_l[1]*rdv2; 
  out[2] += 0.7071067811865475*Gdrag_r[2]*rdv2-0.7071067811865475*Gdrag_l[2]*rdv2; 
  out[3] += 0.7071067811865475*Gdrag_r[3]*rdv2-0.7071067811865475*Gdrag_l[3]*rdv2; 
  out[4] += 1.224744871391589*Gdrag_r[0]*rdv2+1.224744871391589*Gdrag_l[0]*rdv2; 
  out[5] += 0.7071067811865475*Gdrag_r[4]*rdv2-0.7071067811865475*Gdrag_l[4]*rdv2; 
  out[6] += 0.7071067811865475*Gdrag_r[5]*rdv2-0.7071067811865475*Gdrag_l[5]*rdv2; 
  out[7] += 0.7071067811865475*Gdrag_r[6]*rdv2-0.7071067811865475*Gdrag_l[6]*rdv2; 
  out[8] += 0.7071067811865475*Gdrag_r[7]*rdv2-0.7071067811865475*Gdrag_l[7]*rdv2; 
  out[9] += 1.224744871391589*Gdrag_r[1]*rdv2+1.224744871391589*Gdrag_l[1]*rdv2; 
  out[10] += 1.224744871391589*Gdrag_r[2]*rdv2+1.224744871391589*Gdrag_l[2]*rdv2; 
  out[11] += 1.224744871391589*Gdrag_r[3]*rdv2+1.224744871391589*Gdrag_l[3]*rdv2; 
  out[12] += 0.7071067811865475*Gdrag_r[8]*rdv2-0.7071067811865475*Gdrag_l[8]*rdv2; 
  out[13] += 0.7071067811865475*Gdrag_r[9]*rdv2-0.7071067811865475*Gdrag_l[9]*rdv2; 
  out[14] += 0.7071067811865475*Gdrag_r[10]*rdv2-0.7071067811865475*Gdrag_l[10]*rdv2; 
  out[15] += 1.224744871391589*Gdrag_r[4]*rdv2+1.224744871391589*Gdrag_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Gdrag_r[11]*rdv2-0.7071067811865475*Gdrag_l[11]*rdv2; 
  out[17] += 1.224744871391589*Gdrag_r[5]*rdv2+1.224744871391589*Gdrag_l[5]*rdv2; 
  out[18] += 1.224744871391589*Gdrag_r[6]*rdv2+1.224744871391589*Gdrag_l[6]*rdv2; 
  out[19] += 1.224744871391589*Gdrag_r[7]*rdv2+1.224744871391589*Gdrag_l[7]*rdv2; 
  out[20] += 0.7071067811865475*Gdrag_r[12]*rdv2-0.7071067811865475*Gdrag_l[12]*rdv2; 
  out[21] += 0.7071067811865475*Gdrag_r[13]*rdv2-0.7071067811865475*Gdrag_l[13]*rdv2; 
  out[22] += 0.7071067811865475*Gdrag_r[14]*rdv2-0.7071067811865475*Gdrag_l[14]*rdv2; 
  out[23] += 1.224744871391589*Gdrag_r[8]*rdv2+1.224744871391589*Gdrag_l[8]*rdv2; 
  out[24] += 1.224744871391589*Gdrag_r[9]*rdv2+1.224744871391589*Gdrag_l[9]*rdv2; 
  out[25] += 1.224744871391589*Gdrag_r[10]*rdv2+1.224744871391589*Gdrag_l[10]*rdv2; 
  out[26] += 1.224744871391589*Gdrag_r[11]*rdv2+1.224744871391589*Gdrag_l[11]*rdv2; 
  out[27] += 0.7071067811865475*Gdrag_r[15]*rdv2-0.7071067811865475*Gdrag_l[15]*rdv2; 
  out[28] += 1.224744871391589*Gdrag_r[12]*rdv2+1.224744871391589*Gdrag_l[12]*rdv2; 
  out[29] += 1.224744871391589*Gdrag_r[13]*rdv2+1.224744871391589*Gdrag_l[13]*rdv2; 
  out[30] += 1.224744871391589*Gdrag_r[14]*rdv2+1.224744871391589*Gdrag_l[14]*rdv2; 
  out[31] += 1.224744871391589*Gdrag_r[15]*rdv2+1.224744871391589*Gdrag_l[15]*rdv2; 
} 
