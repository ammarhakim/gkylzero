#include <gkyl_vlasov_lbo_kernels.h> 
#include <gkyl_basis_ser_2x3v_p1_surfvz_quad.h> 
GKYL_CU_DH void vlasov_lbo_surfvz_2x3v_ser_p1(const double *w, const double *dxv, const double *nuSum, const double *nuUSum, const double *nuVtSqSum, const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[5]:         cell-center coordinates. 
  // dxv[5]:       cell spacing. 
  // nuSum:         collisionalities added (self and cross species collisionalities). 
  // nuUSum[12]:    sum of bulk velocities times their respective collisionalities. 
  // nuVtSqSum[4]: sum of thermal speeds squared time their respective collisionalities. 
  // fl/fc/fr:      distribution function in cells 
  // out:           incremented distribution function in cell 
  double rdv2 = 2.0/dxv[4]; 
  double rdvSq4 = 4.0/(dxv[4]*dxv[4]); 

  const double *sumNuUz = &nuUSum[8]; 

  double alphaDrSurf_l[16] = {0.0}; 
  alphaDrSurf_l[0] = 4.0*nuSum[0]*w[4]-2.0*nuSum[0]*dxv[4]-2.0*sumNuUz[0]; 
  alphaDrSurf_l[1] = -2.0*sumNuUz[1]; 
  alphaDrSurf_l[2] = -2.0*sumNuUz[2]; 
  alphaDrSurf_l[5] = -2.0*sumNuUz[3]; 

  double alphaDrSurf_r[16] = {0.0}; 
  alphaDrSurf_r[0] = 4.0*nuSum[0]*w[4]+2.0*nuSum[0]*dxv[4]-2.0*sumNuUz[0]; 
  alphaDrSurf_r[1] = -2.0*sumNuUz[1]; 
  alphaDrSurf_r[2] = -2.0*sumNuUz[2]; 
  alphaDrSurf_r[5] = -2.0*sumNuUz[3]; 

  double fUpwindQuad_l[16] = {0.0};
  double fUpwindQuad_r[16] = {0.0};
  double fUpwind_l[16] = {0.0};;
  double fUpwind_r[16] = {0.0};
  double Ghat_l[16] = {0.0}; 
  double Ghat_r[16] = {0.0}; 
  double Gdiff_l[16] = {0.0}; 
  double Gdiff_r[16] = {0.0}; 
  double Gdiff2_l[16] = {0.0}; 
  double Gdiff2_r[16] = {0.0}; 

  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[0] = ser_2x3v_p1_surfvz_quad_0(1, fl); 
  } else { 
    fUpwindQuad_l[0] = ser_2x3v_p1_surfvz_quad_0(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[1] = ser_2x3v_p1_surfvz_quad_1(1, fl); 
  } else { 
    fUpwindQuad_l[1] = ser_2x3v_p1_surfvz_quad_1(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[2] = ser_2x3v_p1_surfvz_quad_2(1, fl); 
  } else { 
    fUpwindQuad_l[2] = ser_2x3v_p1_surfvz_quad_2(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[3] = ser_2x3v_p1_surfvz_quad_3(1, fl); 
  } else { 
    fUpwindQuad_l[3] = ser_2x3v_p1_surfvz_quad_3(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[4] = ser_2x3v_p1_surfvz_quad_4(1, fl); 
  } else { 
    fUpwindQuad_l[4] = ser_2x3v_p1_surfvz_quad_4(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[5] = ser_2x3v_p1_surfvz_quad_5(1, fl); 
  } else { 
    fUpwindQuad_l[5] = ser_2x3v_p1_surfvz_quad_5(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[6] = ser_2x3v_p1_surfvz_quad_6(1, fl); 
  } else { 
    fUpwindQuad_l[6] = ser_2x3v_p1_surfvz_quad_6(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[7] = ser_2x3v_p1_surfvz_quad_7(1, fl); 
  } else { 
    fUpwindQuad_l[7] = ser_2x3v_p1_surfvz_quad_7(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[8] = ser_2x3v_p1_surfvz_quad_8(1, fl); 
  } else { 
    fUpwindQuad_l[8] = ser_2x3v_p1_surfvz_quad_8(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[9] = ser_2x3v_p1_surfvz_quad_9(1, fl); 
  } else { 
    fUpwindQuad_l[9] = ser_2x3v_p1_surfvz_quad_9(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[10] = ser_2x3v_p1_surfvz_quad_10(1, fl); 
  } else { 
    fUpwindQuad_l[10] = ser_2x3v_p1_surfvz_quad_10(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[11] = ser_2x3v_p1_surfvz_quad_11(1, fl); 
  } else { 
    fUpwindQuad_l[11] = ser_2x3v_p1_surfvz_quad_11(-1, fc); 
  } 
  if (alphaDrSurf_l[5]-alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[12] = ser_2x3v_p1_surfvz_quad_12(1, fl); 
  } else { 
    fUpwindQuad_l[12] = ser_2x3v_p1_surfvz_quad_12(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])-alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[13] = ser_2x3v_p1_surfvz_quad_13(1, fl); 
  } else { 
    fUpwindQuad_l[13] = ser_2x3v_p1_surfvz_quad_13(-1, fc); 
  } 
  if ((-alphaDrSurf_l[5])+alphaDrSurf_l[2]-alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[14] = ser_2x3v_p1_surfvz_quad_14(1, fl); 
  } else { 
    fUpwindQuad_l[14] = ser_2x3v_p1_surfvz_quad_14(-1, fc); 
  } 
  if (alphaDrSurf_l[5]+alphaDrSurf_l[2]+alphaDrSurf_l[1]+alphaDrSurf_l[0] > 0) { 
    fUpwindQuad_l[15] = ser_2x3v_p1_surfvz_quad_15(1, fl); 
  } else { 
    fUpwindQuad_l[15] = ser_2x3v_p1_surfvz_quad_15(-1, fc); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvz_quad_0(1, fc); 
  } else { 
    fUpwindQuad_r[0] = ser_2x3v_p1_surfvz_quad_0(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvz_quad_1(1, fc); 
  } else { 
    fUpwindQuad_r[1] = ser_2x3v_p1_surfvz_quad_1(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvz_quad_2(1, fc); 
  } else { 
    fUpwindQuad_r[2] = ser_2x3v_p1_surfvz_quad_2(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvz_quad_3(1, fc); 
  } else { 
    fUpwindQuad_r[3] = ser_2x3v_p1_surfvz_quad_3(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvz_quad_4(1, fc); 
  } else { 
    fUpwindQuad_r[4] = ser_2x3v_p1_surfvz_quad_4(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvz_quad_5(1, fc); 
  } else { 
    fUpwindQuad_r[5] = ser_2x3v_p1_surfvz_quad_5(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvz_quad_6(1, fc); 
  } else { 
    fUpwindQuad_r[6] = ser_2x3v_p1_surfvz_quad_6(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvz_quad_7(1, fc); 
  } else { 
    fUpwindQuad_r[7] = ser_2x3v_p1_surfvz_quad_7(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvz_quad_8(1, fc); 
  } else { 
    fUpwindQuad_r[8] = ser_2x3v_p1_surfvz_quad_8(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvz_quad_9(1, fc); 
  } else { 
    fUpwindQuad_r[9] = ser_2x3v_p1_surfvz_quad_9(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvz_quad_10(1, fc); 
  } else { 
    fUpwindQuad_r[10] = ser_2x3v_p1_surfvz_quad_10(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvz_quad_11(1, fc); 
  } else { 
    fUpwindQuad_r[11] = ser_2x3v_p1_surfvz_quad_11(-1, fr); 
  } 
  if (alphaDrSurf_r[5]-alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvz_quad_12(1, fc); 
  } else { 
    fUpwindQuad_r[12] = ser_2x3v_p1_surfvz_quad_12(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])-alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvz_quad_13(1, fc); 
  } else { 
    fUpwindQuad_r[13] = ser_2x3v_p1_surfvz_quad_13(-1, fr); 
  } 
  if ((-alphaDrSurf_r[5])+alphaDrSurf_r[2]-alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvz_quad_14(1, fc); 
  } else { 
    fUpwindQuad_r[14] = ser_2x3v_p1_surfvz_quad_14(-1, fr); 
  } 
  if (alphaDrSurf_r[5]+alphaDrSurf_r[2]+alphaDrSurf_r[1]+alphaDrSurf_r[0] > 0) { 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvz_quad_15(1, fc); 
  } else { 
    fUpwindQuad_r[15] = ser_2x3v_p1_surfvz_quad_15(-1, fr); 
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

  Gdiff2_l[0] = 0.2041241452319315*nuVtSqSum[3]*fl[20]-0.2041241452319315*nuVtSqSum[3]*fc[20]+0.2041241452319315*nuVtSqSum[2]*fl[13]-0.2041241452319315*nuVtSqSum[2]*fc[13]+0.2041241452319315*nuVtSqSum[1]*fl[12]-0.2041241452319315*nuVtSqSum[1]*fc[12]+0.1767766952966368*nuVtSqSum[3]*fl[6]+0.1767766952966368*nuVtSqSum[3]*fc[6]+0.2041241452319315*nuVtSqSum[0]*fl[5]-0.2041241452319315*nuVtSqSum[0]*fc[5]+0.1767766952966368*fl[2]*nuVtSqSum[2]+0.1767766952966368*fc[2]*nuVtSqSum[2]+0.1767766952966368*fl[1]*nuVtSqSum[1]+0.1767766952966368*fc[1]*nuVtSqSum[1]+0.1767766952966368*fl[0]*nuVtSqSum[0]+0.1767766952966368*fc[0]*nuVtSqSum[0]; 
  Gdiff2_l[1] = 0.2041241452319315*nuVtSqSum[2]*fl[20]-0.2041241452319315*nuVtSqSum[2]*fc[20]+0.2041241452319315*nuVtSqSum[3]*fl[13]-0.2041241452319315*nuVtSqSum[3]*fc[13]+0.2041241452319315*nuVtSqSum[0]*fl[12]-0.2041241452319315*nuVtSqSum[0]*fc[12]+0.1767766952966368*nuVtSqSum[2]*fl[6]+0.1767766952966368*nuVtSqSum[2]*fc[6]+0.2041241452319315*nuVtSqSum[1]*fl[5]-0.2041241452319315*nuVtSqSum[1]*fc[5]+0.1767766952966368*fl[2]*nuVtSqSum[3]+0.1767766952966368*fc[2]*nuVtSqSum[3]+0.1767766952966368*fl[0]*nuVtSqSum[1]+0.1767766952966368*fc[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fl[1]+0.1767766952966368*nuVtSqSum[0]*fc[1]; 
  Gdiff2_l[2] = 0.2041241452319315*nuVtSqSum[1]*fl[20]-0.2041241452319315*nuVtSqSum[1]*fc[20]+0.2041241452319315*nuVtSqSum[0]*fl[13]-0.2041241452319315*nuVtSqSum[0]*fc[13]+0.2041241452319315*nuVtSqSum[3]*fl[12]-0.2041241452319315*nuVtSqSum[3]*fc[12]+0.1767766952966368*nuVtSqSum[1]*fl[6]+0.1767766952966368*nuVtSqSum[1]*fc[6]+0.2041241452319315*nuVtSqSum[2]*fl[5]-0.2041241452319315*nuVtSqSum[2]*fc[5]+0.1767766952966368*fl[1]*nuVtSqSum[3]+0.1767766952966368*fc[1]*nuVtSqSum[3]+0.1767766952966368*fl[0]*nuVtSqSum[2]+0.1767766952966368*fc[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fl[2]+0.1767766952966368*nuVtSqSum[0]*fc[2]; 
  Gdiff2_l[3] = 0.2041241452319315*nuVtSqSum[3]*fl[27]-0.2041241452319315*nuVtSqSum[3]*fc[27]+0.2041241452319315*nuVtSqSum[2]*fl[22]-0.2041241452319315*nuVtSqSum[2]*fc[22]+0.2041241452319315*nuVtSqSum[1]*fl[21]-0.2041241452319315*nuVtSqSum[1]*fc[21]+0.1767766952966368*nuVtSqSum[3]*fl[16]+0.1767766952966368*nuVtSqSum[3]*fc[16]+0.2041241452319315*nuVtSqSum[0]*fl[14]-0.2041241452319315*nuVtSqSum[0]*fc[14]+0.1767766952966368*nuVtSqSum[2]*fl[8]+0.1767766952966368*nuVtSqSum[2]*fc[8]+0.1767766952966368*nuVtSqSum[1]*fl[7]+0.1767766952966368*nuVtSqSum[1]*fc[7]+0.1767766952966368*nuVtSqSum[0]*fl[3]+0.1767766952966368*nuVtSqSum[0]*fc[3]; 
  Gdiff2_l[4] = 0.2041241452319315*nuVtSqSum[3]*fl[28]-0.2041241452319315*nuVtSqSum[3]*fc[28]+0.2041241452319315*nuVtSqSum[2]*fl[24]-0.2041241452319315*nuVtSqSum[2]*fc[24]+0.2041241452319315*nuVtSqSum[1]*fl[23]-0.2041241452319315*nuVtSqSum[1]*fc[23]+0.1767766952966368*nuVtSqSum[3]*fl[17]+0.1767766952966368*nuVtSqSum[3]*fc[17]+0.2041241452319315*nuVtSqSum[0]*fl[15]-0.2041241452319315*nuVtSqSum[0]*fc[15]+0.1767766952966368*nuVtSqSum[2]*fl[10]+0.1767766952966368*nuVtSqSum[2]*fc[10]+0.1767766952966368*nuVtSqSum[1]*fl[9]+0.1767766952966368*nuVtSqSum[1]*fc[9]+0.1767766952966368*nuVtSqSum[0]*fl[4]+0.1767766952966368*nuVtSqSum[0]*fc[4]; 
  Gdiff2_l[5] = 0.2041241452319315*nuVtSqSum[0]*fl[20]-0.2041241452319315*nuVtSqSum[0]*fc[20]+0.2041241452319315*nuVtSqSum[1]*fl[13]-0.2041241452319315*nuVtSqSum[1]*fc[13]+0.2041241452319315*nuVtSqSum[2]*fl[12]-0.2041241452319315*nuVtSqSum[2]*fc[12]+0.1767766952966368*nuVtSqSum[0]*fl[6]+0.1767766952966368*nuVtSqSum[0]*fc[6]+0.2041241452319315*nuVtSqSum[3]*fl[5]-0.2041241452319315*nuVtSqSum[3]*fc[5]+0.1767766952966368*fl[0]*nuVtSqSum[3]+0.1767766952966368*fc[0]*nuVtSqSum[3]+0.1767766952966368*fl[1]*nuVtSqSum[2]+0.1767766952966368*fc[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fl[2]+0.1767766952966368*nuVtSqSum[1]*fc[2]; 
  Gdiff2_l[6] = 0.2041241452319315*nuVtSqSum[2]*fl[27]-0.2041241452319315*nuVtSqSum[2]*fc[27]+0.2041241452319315*nuVtSqSum[3]*fl[22]-0.2041241452319315*nuVtSqSum[3]*fc[22]+0.2041241452319315*nuVtSqSum[0]*fl[21]-0.2041241452319315*nuVtSqSum[0]*fc[21]+0.1767766952966368*nuVtSqSum[2]*fl[16]+0.1767766952966368*nuVtSqSum[2]*fc[16]+0.2041241452319315*nuVtSqSum[1]*fl[14]-0.2041241452319315*nuVtSqSum[1]*fc[14]+0.1767766952966368*nuVtSqSum[3]*fl[8]+0.1767766952966368*nuVtSqSum[3]*fc[8]+0.1767766952966368*nuVtSqSum[0]*fl[7]+0.1767766952966368*nuVtSqSum[0]*fc[7]+0.1767766952966368*nuVtSqSum[1]*fl[3]+0.1767766952966368*nuVtSqSum[1]*fc[3]; 
  Gdiff2_l[7] = 0.2041241452319315*nuVtSqSum[1]*fl[27]-0.2041241452319315*nuVtSqSum[1]*fc[27]+0.2041241452319315*nuVtSqSum[0]*fl[22]-0.2041241452319315*nuVtSqSum[0]*fc[22]+0.2041241452319315*nuVtSqSum[3]*fl[21]-0.2041241452319315*nuVtSqSum[3]*fc[21]+0.1767766952966368*nuVtSqSum[1]*fl[16]+0.1767766952966368*nuVtSqSum[1]*fc[16]+0.2041241452319315*nuVtSqSum[2]*fl[14]-0.2041241452319315*nuVtSqSum[2]*fc[14]+0.1767766952966368*nuVtSqSum[0]*fl[8]+0.1767766952966368*nuVtSqSum[0]*fc[8]+0.1767766952966368*nuVtSqSum[3]*fl[7]+0.1767766952966368*nuVtSqSum[3]*fc[7]+0.1767766952966368*nuVtSqSum[2]*fl[3]+0.1767766952966368*nuVtSqSum[2]*fc[3]; 
  Gdiff2_l[8] = 0.2041241452319315*nuVtSqSum[2]*fl[28]-0.2041241452319315*nuVtSqSum[2]*fc[28]+0.2041241452319315*nuVtSqSum[3]*fl[24]-0.2041241452319315*nuVtSqSum[3]*fc[24]+0.2041241452319315*nuVtSqSum[0]*fl[23]-0.2041241452319315*nuVtSqSum[0]*fc[23]+0.1767766952966368*nuVtSqSum[2]*fl[17]+0.1767766952966368*nuVtSqSum[2]*fc[17]+0.2041241452319315*nuVtSqSum[1]*fl[15]-0.2041241452319315*nuVtSqSum[1]*fc[15]+0.1767766952966368*nuVtSqSum[3]*fl[10]+0.1767766952966368*nuVtSqSum[3]*fc[10]+0.1767766952966368*nuVtSqSum[0]*fl[9]+0.1767766952966368*nuVtSqSum[0]*fc[9]+0.1767766952966368*nuVtSqSum[1]*fl[4]+0.1767766952966368*nuVtSqSum[1]*fc[4]; 
  Gdiff2_l[9] = 0.2041241452319315*nuVtSqSum[1]*fl[28]-0.2041241452319315*nuVtSqSum[1]*fc[28]+0.2041241452319315*nuVtSqSum[0]*fl[24]-0.2041241452319315*nuVtSqSum[0]*fc[24]+0.2041241452319315*nuVtSqSum[3]*fl[23]-0.2041241452319315*nuVtSqSum[3]*fc[23]+0.1767766952966368*nuVtSqSum[1]*fl[17]+0.1767766952966368*nuVtSqSum[1]*fc[17]+0.2041241452319315*nuVtSqSum[2]*fl[15]-0.2041241452319315*nuVtSqSum[2]*fc[15]+0.1767766952966368*nuVtSqSum[0]*fl[10]+0.1767766952966368*nuVtSqSum[0]*fc[10]+0.1767766952966368*nuVtSqSum[3]*fl[9]+0.1767766952966368*nuVtSqSum[3]*fc[9]+0.1767766952966368*nuVtSqSum[2]*fl[4]+0.1767766952966368*nuVtSqSum[2]*fc[4]; 
  Gdiff2_l[10] = 0.2041241452319315*nuVtSqSum[3]*fl[31]-0.2041241452319315*nuVtSqSum[3]*fc[31]+0.2041241452319315*nuVtSqSum[2]*fl[30]-0.2041241452319315*nuVtSqSum[2]*fc[30]+0.2041241452319315*nuVtSqSum[1]*fl[29]-0.2041241452319315*nuVtSqSum[1]*fc[29]+0.1767766952966368*nuVtSqSum[3]*fl[26]+0.1767766952966368*nuVtSqSum[3]*fc[26]+0.2041241452319315*nuVtSqSum[0]*fl[25]-0.2041241452319315*nuVtSqSum[0]*fc[25]+0.1767766952966368*nuVtSqSum[2]*fl[19]+0.1767766952966368*nuVtSqSum[2]*fc[19]+0.1767766952966368*nuVtSqSum[1]*fl[18]+0.1767766952966368*nuVtSqSum[1]*fc[18]+0.1767766952966368*nuVtSqSum[0]*fl[11]+0.1767766952966368*nuVtSqSum[0]*fc[11]; 
  Gdiff2_l[11] = 0.2041241452319315*nuVtSqSum[0]*fl[27]-0.2041241452319315*nuVtSqSum[0]*fc[27]+0.2041241452319315*nuVtSqSum[1]*fl[22]-0.2041241452319315*nuVtSqSum[1]*fc[22]+0.2041241452319315*nuVtSqSum[2]*fl[21]-0.2041241452319315*nuVtSqSum[2]*fc[21]+0.1767766952966368*nuVtSqSum[0]*fl[16]+0.1767766952966368*nuVtSqSum[0]*fc[16]+0.2041241452319315*nuVtSqSum[3]*fl[14]-0.2041241452319315*nuVtSqSum[3]*fc[14]+0.1767766952966368*nuVtSqSum[1]*fl[8]+0.1767766952966368*nuVtSqSum[1]*fc[8]+0.1767766952966368*nuVtSqSum[2]*fl[7]+0.1767766952966368*nuVtSqSum[2]*fc[7]+0.1767766952966368*fl[3]*nuVtSqSum[3]+0.1767766952966368*fc[3]*nuVtSqSum[3]; 
  Gdiff2_l[12] = 0.2041241452319315*nuVtSqSum[0]*fl[28]-0.2041241452319315*nuVtSqSum[0]*fc[28]+0.2041241452319315*nuVtSqSum[1]*fl[24]-0.2041241452319315*nuVtSqSum[1]*fc[24]+0.2041241452319315*nuVtSqSum[2]*fl[23]-0.2041241452319315*nuVtSqSum[2]*fc[23]+0.1767766952966368*nuVtSqSum[0]*fl[17]+0.1767766952966368*nuVtSqSum[0]*fc[17]+0.2041241452319315*nuVtSqSum[3]*fl[15]-0.2041241452319315*nuVtSqSum[3]*fc[15]+0.1767766952966368*nuVtSqSum[1]*fl[10]+0.1767766952966368*nuVtSqSum[1]*fc[10]+0.1767766952966368*nuVtSqSum[2]*fl[9]+0.1767766952966368*nuVtSqSum[2]*fc[9]+0.1767766952966368*nuVtSqSum[3]*fl[4]+0.1767766952966368*nuVtSqSum[3]*fc[4]; 
  Gdiff2_l[13] = 0.2041241452319315*nuVtSqSum[2]*fl[31]-0.2041241452319315*nuVtSqSum[2]*fc[31]+0.2041241452319315*nuVtSqSum[3]*fl[30]-0.2041241452319315*nuVtSqSum[3]*fc[30]+0.2041241452319315*nuVtSqSum[0]*fl[29]-0.2041241452319315*nuVtSqSum[0]*fc[29]+0.1767766952966368*nuVtSqSum[2]*fl[26]+0.1767766952966368*nuVtSqSum[2]*fc[26]+0.2041241452319315*nuVtSqSum[1]*fl[25]-0.2041241452319315*nuVtSqSum[1]*fc[25]+0.1767766952966368*nuVtSqSum[3]*fl[19]+0.1767766952966368*nuVtSqSum[3]*fc[19]+0.1767766952966368*nuVtSqSum[0]*fl[18]+0.1767766952966368*nuVtSqSum[0]*fc[18]+0.1767766952966368*nuVtSqSum[1]*fl[11]+0.1767766952966368*nuVtSqSum[1]*fc[11]; 
  Gdiff2_l[14] = 0.2041241452319315*nuVtSqSum[1]*fl[31]-0.2041241452319315*nuVtSqSum[1]*fc[31]+0.2041241452319315*nuVtSqSum[0]*fl[30]-0.2041241452319315*nuVtSqSum[0]*fc[30]+0.2041241452319315*nuVtSqSum[3]*fl[29]-0.2041241452319315*nuVtSqSum[3]*fc[29]+0.1767766952966368*nuVtSqSum[1]*fl[26]+0.1767766952966368*nuVtSqSum[1]*fc[26]+0.2041241452319315*nuVtSqSum[2]*fl[25]-0.2041241452319315*nuVtSqSum[2]*fc[25]+0.1767766952966368*nuVtSqSum[0]*fl[19]+0.1767766952966368*nuVtSqSum[0]*fc[19]+0.1767766952966368*nuVtSqSum[3]*fl[18]+0.1767766952966368*nuVtSqSum[3]*fc[18]+0.1767766952966368*nuVtSqSum[2]*fl[11]+0.1767766952966368*nuVtSqSum[2]*fc[11]; 
  Gdiff2_l[15] = 0.2041241452319315*nuVtSqSum[0]*fl[31]-0.2041241452319315*nuVtSqSum[0]*fc[31]+0.2041241452319315*nuVtSqSum[1]*fl[30]-0.2041241452319315*nuVtSqSum[1]*fc[30]+0.2041241452319315*nuVtSqSum[2]*fl[29]-0.2041241452319315*nuVtSqSum[2]*fc[29]+0.1767766952966368*nuVtSqSum[0]*fl[26]+0.1767766952966368*nuVtSqSum[0]*fc[26]+0.2041241452319315*nuVtSqSum[3]*fl[25]-0.2041241452319315*nuVtSqSum[3]*fc[25]+0.1767766952966368*nuVtSqSum[1]*fl[19]+0.1767766952966368*nuVtSqSum[1]*fc[19]+0.1767766952966368*nuVtSqSum[2]*fl[18]+0.1767766952966368*nuVtSqSum[2]*fc[18]+0.1767766952966368*nuVtSqSum[3]*fl[11]+0.1767766952966368*nuVtSqSum[3]*fc[11]; 

  Gdiff2_r[0] = (-0.2041241452319315*nuVtSqSum[3]*fr[20])+0.2041241452319315*nuVtSqSum[3]*fc[20]-0.2041241452319315*nuVtSqSum[2]*fr[13]+0.2041241452319315*nuVtSqSum[2]*fc[13]-0.2041241452319315*nuVtSqSum[1]*fr[12]+0.2041241452319315*nuVtSqSum[1]*fc[12]+0.1767766952966368*nuVtSqSum[3]*fr[6]+0.1767766952966368*nuVtSqSum[3]*fc[6]-0.2041241452319315*nuVtSqSum[0]*fr[5]+0.2041241452319315*nuVtSqSum[0]*fc[5]+0.1767766952966368*fr[2]*nuVtSqSum[2]+0.1767766952966368*fc[2]*nuVtSqSum[2]+0.1767766952966368*fr[1]*nuVtSqSum[1]+0.1767766952966368*fc[1]*nuVtSqSum[1]+0.1767766952966368*fr[0]*nuVtSqSum[0]+0.1767766952966368*fc[0]*nuVtSqSum[0]; 
  Gdiff2_r[1] = (-0.2041241452319315*nuVtSqSum[2]*fr[20])+0.2041241452319315*nuVtSqSum[2]*fc[20]-0.2041241452319315*nuVtSqSum[3]*fr[13]+0.2041241452319315*nuVtSqSum[3]*fc[13]-0.2041241452319315*nuVtSqSum[0]*fr[12]+0.2041241452319315*nuVtSqSum[0]*fc[12]+0.1767766952966368*nuVtSqSum[2]*fr[6]+0.1767766952966368*nuVtSqSum[2]*fc[6]-0.2041241452319315*nuVtSqSum[1]*fr[5]+0.2041241452319315*nuVtSqSum[1]*fc[5]+0.1767766952966368*fr[2]*nuVtSqSum[3]+0.1767766952966368*fc[2]*nuVtSqSum[3]+0.1767766952966368*fr[0]*nuVtSqSum[1]+0.1767766952966368*fc[0]*nuVtSqSum[1]+0.1767766952966368*nuVtSqSum[0]*fr[1]+0.1767766952966368*nuVtSqSum[0]*fc[1]; 
  Gdiff2_r[2] = (-0.2041241452319315*nuVtSqSum[1]*fr[20])+0.2041241452319315*nuVtSqSum[1]*fc[20]-0.2041241452319315*nuVtSqSum[0]*fr[13]+0.2041241452319315*nuVtSqSum[0]*fc[13]-0.2041241452319315*nuVtSqSum[3]*fr[12]+0.2041241452319315*nuVtSqSum[3]*fc[12]+0.1767766952966368*nuVtSqSum[1]*fr[6]+0.1767766952966368*nuVtSqSum[1]*fc[6]-0.2041241452319315*nuVtSqSum[2]*fr[5]+0.2041241452319315*nuVtSqSum[2]*fc[5]+0.1767766952966368*fr[1]*nuVtSqSum[3]+0.1767766952966368*fc[1]*nuVtSqSum[3]+0.1767766952966368*fr[0]*nuVtSqSum[2]+0.1767766952966368*fc[0]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[0]*fr[2]+0.1767766952966368*nuVtSqSum[0]*fc[2]; 
  Gdiff2_r[3] = (-0.2041241452319315*nuVtSqSum[3]*fr[27])+0.2041241452319315*nuVtSqSum[3]*fc[27]-0.2041241452319315*nuVtSqSum[2]*fr[22]+0.2041241452319315*nuVtSqSum[2]*fc[22]-0.2041241452319315*nuVtSqSum[1]*fr[21]+0.2041241452319315*nuVtSqSum[1]*fc[21]+0.1767766952966368*nuVtSqSum[3]*fr[16]+0.1767766952966368*nuVtSqSum[3]*fc[16]-0.2041241452319315*nuVtSqSum[0]*fr[14]+0.2041241452319315*nuVtSqSum[0]*fc[14]+0.1767766952966368*nuVtSqSum[2]*fr[8]+0.1767766952966368*nuVtSqSum[2]*fc[8]+0.1767766952966368*nuVtSqSum[1]*fr[7]+0.1767766952966368*nuVtSqSum[1]*fc[7]+0.1767766952966368*nuVtSqSum[0]*fr[3]+0.1767766952966368*nuVtSqSum[0]*fc[3]; 
  Gdiff2_r[4] = (-0.2041241452319315*nuVtSqSum[3]*fr[28])+0.2041241452319315*nuVtSqSum[3]*fc[28]-0.2041241452319315*nuVtSqSum[2]*fr[24]+0.2041241452319315*nuVtSqSum[2]*fc[24]-0.2041241452319315*nuVtSqSum[1]*fr[23]+0.2041241452319315*nuVtSqSum[1]*fc[23]+0.1767766952966368*nuVtSqSum[3]*fr[17]+0.1767766952966368*nuVtSqSum[3]*fc[17]-0.2041241452319315*nuVtSqSum[0]*fr[15]+0.2041241452319315*nuVtSqSum[0]*fc[15]+0.1767766952966368*nuVtSqSum[2]*fr[10]+0.1767766952966368*nuVtSqSum[2]*fc[10]+0.1767766952966368*nuVtSqSum[1]*fr[9]+0.1767766952966368*nuVtSqSum[1]*fc[9]+0.1767766952966368*nuVtSqSum[0]*fr[4]+0.1767766952966368*nuVtSqSum[0]*fc[4]; 
  Gdiff2_r[5] = (-0.2041241452319315*nuVtSqSum[0]*fr[20])+0.2041241452319315*nuVtSqSum[0]*fc[20]-0.2041241452319315*nuVtSqSum[1]*fr[13]+0.2041241452319315*nuVtSqSum[1]*fc[13]-0.2041241452319315*nuVtSqSum[2]*fr[12]+0.2041241452319315*nuVtSqSum[2]*fc[12]+0.1767766952966368*nuVtSqSum[0]*fr[6]+0.1767766952966368*nuVtSqSum[0]*fc[6]-0.2041241452319315*nuVtSqSum[3]*fr[5]+0.2041241452319315*nuVtSqSum[3]*fc[5]+0.1767766952966368*fr[0]*nuVtSqSum[3]+0.1767766952966368*fc[0]*nuVtSqSum[3]+0.1767766952966368*fr[1]*nuVtSqSum[2]+0.1767766952966368*fc[1]*nuVtSqSum[2]+0.1767766952966368*nuVtSqSum[1]*fr[2]+0.1767766952966368*nuVtSqSum[1]*fc[2]; 
  Gdiff2_r[6] = (-0.2041241452319315*nuVtSqSum[2]*fr[27])+0.2041241452319315*nuVtSqSum[2]*fc[27]-0.2041241452319315*nuVtSqSum[3]*fr[22]+0.2041241452319315*nuVtSqSum[3]*fc[22]-0.2041241452319315*nuVtSqSum[0]*fr[21]+0.2041241452319315*nuVtSqSum[0]*fc[21]+0.1767766952966368*nuVtSqSum[2]*fr[16]+0.1767766952966368*nuVtSqSum[2]*fc[16]-0.2041241452319315*nuVtSqSum[1]*fr[14]+0.2041241452319315*nuVtSqSum[1]*fc[14]+0.1767766952966368*nuVtSqSum[3]*fr[8]+0.1767766952966368*nuVtSqSum[3]*fc[8]+0.1767766952966368*nuVtSqSum[0]*fr[7]+0.1767766952966368*nuVtSqSum[0]*fc[7]+0.1767766952966368*nuVtSqSum[1]*fr[3]+0.1767766952966368*nuVtSqSum[1]*fc[3]; 
  Gdiff2_r[7] = (-0.2041241452319315*nuVtSqSum[1]*fr[27])+0.2041241452319315*nuVtSqSum[1]*fc[27]-0.2041241452319315*nuVtSqSum[0]*fr[22]+0.2041241452319315*nuVtSqSum[0]*fc[22]-0.2041241452319315*nuVtSqSum[3]*fr[21]+0.2041241452319315*nuVtSqSum[3]*fc[21]+0.1767766952966368*nuVtSqSum[1]*fr[16]+0.1767766952966368*nuVtSqSum[1]*fc[16]-0.2041241452319315*nuVtSqSum[2]*fr[14]+0.2041241452319315*nuVtSqSum[2]*fc[14]+0.1767766952966368*nuVtSqSum[0]*fr[8]+0.1767766952966368*nuVtSqSum[0]*fc[8]+0.1767766952966368*nuVtSqSum[3]*fr[7]+0.1767766952966368*nuVtSqSum[3]*fc[7]+0.1767766952966368*nuVtSqSum[2]*fr[3]+0.1767766952966368*nuVtSqSum[2]*fc[3]; 
  Gdiff2_r[8] = (-0.2041241452319315*nuVtSqSum[2]*fr[28])+0.2041241452319315*nuVtSqSum[2]*fc[28]-0.2041241452319315*nuVtSqSum[3]*fr[24]+0.2041241452319315*nuVtSqSum[3]*fc[24]-0.2041241452319315*nuVtSqSum[0]*fr[23]+0.2041241452319315*nuVtSqSum[0]*fc[23]+0.1767766952966368*nuVtSqSum[2]*fr[17]+0.1767766952966368*nuVtSqSum[2]*fc[17]-0.2041241452319315*nuVtSqSum[1]*fr[15]+0.2041241452319315*nuVtSqSum[1]*fc[15]+0.1767766952966368*nuVtSqSum[3]*fr[10]+0.1767766952966368*nuVtSqSum[3]*fc[10]+0.1767766952966368*nuVtSqSum[0]*fr[9]+0.1767766952966368*nuVtSqSum[0]*fc[9]+0.1767766952966368*nuVtSqSum[1]*fr[4]+0.1767766952966368*nuVtSqSum[1]*fc[4]; 
  Gdiff2_r[9] = (-0.2041241452319315*nuVtSqSum[1]*fr[28])+0.2041241452319315*nuVtSqSum[1]*fc[28]-0.2041241452319315*nuVtSqSum[0]*fr[24]+0.2041241452319315*nuVtSqSum[0]*fc[24]-0.2041241452319315*nuVtSqSum[3]*fr[23]+0.2041241452319315*nuVtSqSum[3]*fc[23]+0.1767766952966368*nuVtSqSum[1]*fr[17]+0.1767766952966368*nuVtSqSum[1]*fc[17]-0.2041241452319315*nuVtSqSum[2]*fr[15]+0.2041241452319315*nuVtSqSum[2]*fc[15]+0.1767766952966368*nuVtSqSum[0]*fr[10]+0.1767766952966368*nuVtSqSum[0]*fc[10]+0.1767766952966368*nuVtSqSum[3]*fr[9]+0.1767766952966368*nuVtSqSum[3]*fc[9]+0.1767766952966368*nuVtSqSum[2]*fr[4]+0.1767766952966368*nuVtSqSum[2]*fc[4]; 
  Gdiff2_r[10] = (-0.2041241452319315*nuVtSqSum[3]*fr[31])+0.2041241452319315*nuVtSqSum[3]*fc[31]-0.2041241452319315*nuVtSqSum[2]*fr[30]+0.2041241452319315*nuVtSqSum[2]*fc[30]-0.2041241452319315*nuVtSqSum[1]*fr[29]+0.2041241452319315*nuVtSqSum[1]*fc[29]+0.1767766952966368*nuVtSqSum[3]*fr[26]+0.1767766952966368*nuVtSqSum[3]*fc[26]-0.2041241452319315*nuVtSqSum[0]*fr[25]+0.2041241452319315*nuVtSqSum[0]*fc[25]+0.1767766952966368*nuVtSqSum[2]*fr[19]+0.1767766952966368*nuVtSqSum[2]*fc[19]+0.1767766952966368*nuVtSqSum[1]*fr[18]+0.1767766952966368*nuVtSqSum[1]*fc[18]+0.1767766952966368*nuVtSqSum[0]*fr[11]+0.1767766952966368*nuVtSqSum[0]*fc[11]; 
  Gdiff2_r[11] = (-0.2041241452319315*nuVtSqSum[0]*fr[27])+0.2041241452319315*nuVtSqSum[0]*fc[27]-0.2041241452319315*nuVtSqSum[1]*fr[22]+0.2041241452319315*nuVtSqSum[1]*fc[22]-0.2041241452319315*nuVtSqSum[2]*fr[21]+0.2041241452319315*nuVtSqSum[2]*fc[21]+0.1767766952966368*nuVtSqSum[0]*fr[16]+0.1767766952966368*nuVtSqSum[0]*fc[16]-0.2041241452319315*nuVtSqSum[3]*fr[14]+0.2041241452319315*nuVtSqSum[3]*fc[14]+0.1767766952966368*nuVtSqSum[1]*fr[8]+0.1767766952966368*nuVtSqSum[1]*fc[8]+0.1767766952966368*nuVtSqSum[2]*fr[7]+0.1767766952966368*nuVtSqSum[2]*fc[7]+0.1767766952966368*fr[3]*nuVtSqSum[3]+0.1767766952966368*fc[3]*nuVtSqSum[3]; 
  Gdiff2_r[12] = (-0.2041241452319315*nuVtSqSum[0]*fr[28])+0.2041241452319315*nuVtSqSum[0]*fc[28]-0.2041241452319315*nuVtSqSum[1]*fr[24]+0.2041241452319315*nuVtSqSum[1]*fc[24]-0.2041241452319315*nuVtSqSum[2]*fr[23]+0.2041241452319315*nuVtSqSum[2]*fc[23]+0.1767766952966368*nuVtSqSum[0]*fr[17]+0.1767766952966368*nuVtSqSum[0]*fc[17]-0.2041241452319315*nuVtSqSum[3]*fr[15]+0.2041241452319315*nuVtSqSum[3]*fc[15]+0.1767766952966368*nuVtSqSum[1]*fr[10]+0.1767766952966368*nuVtSqSum[1]*fc[10]+0.1767766952966368*nuVtSqSum[2]*fr[9]+0.1767766952966368*nuVtSqSum[2]*fc[9]+0.1767766952966368*nuVtSqSum[3]*fr[4]+0.1767766952966368*nuVtSqSum[3]*fc[4]; 
  Gdiff2_r[13] = (-0.2041241452319315*nuVtSqSum[2]*fr[31])+0.2041241452319315*nuVtSqSum[2]*fc[31]-0.2041241452319315*nuVtSqSum[3]*fr[30]+0.2041241452319315*nuVtSqSum[3]*fc[30]-0.2041241452319315*nuVtSqSum[0]*fr[29]+0.2041241452319315*nuVtSqSum[0]*fc[29]+0.1767766952966368*nuVtSqSum[2]*fr[26]+0.1767766952966368*nuVtSqSum[2]*fc[26]-0.2041241452319315*nuVtSqSum[1]*fr[25]+0.2041241452319315*nuVtSqSum[1]*fc[25]+0.1767766952966368*nuVtSqSum[3]*fr[19]+0.1767766952966368*nuVtSqSum[3]*fc[19]+0.1767766952966368*nuVtSqSum[0]*fr[18]+0.1767766952966368*nuVtSqSum[0]*fc[18]+0.1767766952966368*nuVtSqSum[1]*fr[11]+0.1767766952966368*nuVtSqSum[1]*fc[11]; 
  Gdiff2_r[14] = (-0.2041241452319315*nuVtSqSum[1]*fr[31])+0.2041241452319315*nuVtSqSum[1]*fc[31]-0.2041241452319315*nuVtSqSum[0]*fr[30]+0.2041241452319315*nuVtSqSum[0]*fc[30]-0.2041241452319315*nuVtSqSum[3]*fr[29]+0.2041241452319315*nuVtSqSum[3]*fc[29]+0.1767766952966368*nuVtSqSum[1]*fr[26]+0.1767766952966368*nuVtSqSum[1]*fc[26]-0.2041241452319315*nuVtSqSum[2]*fr[25]+0.2041241452319315*nuVtSqSum[2]*fc[25]+0.1767766952966368*nuVtSqSum[0]*fr[19]+0.1767766952966368*nuVtSqSum[0]*fc[19]+0.1767766952966368*nuVtSqSum[3]*fr[18]+0.1767766952966368*nuVtSqSum[3]*fc[18]+0.1767766952966368*nuVtSqSum[2]*fr[11]+0.1767766952966368*nuVtSqSum[2]*fc[11]; 
  Gdiff2_r[15] = (-0.2041241452319315*nuVtSqSum[0]*fr[31])+0.2041241452319315*nuVtSqSum[0]*fc[31]-0.2041241452319315*nuVtSqSum[1]*fr[30]+0.2041241452319315*nuVtSqSum[1]*fc[30]-0.2041241452319315*nuVtSqSum[2]*fr[29]+0.2041241452319315*nuVtSqSum[2]*fc[29]+0.1767766952966368*nuVtSqSum[0]*fr[26]+0.1767766952966368*nuVtSqSum[0]*fc[26]-0.2041241452319315*nuVtSqSum[3]*fr[25]+0.2041241452319315*nuVtSqSum[3]*fc[25]+0.1767766952966368*nuVtSqSum[1]*fr[19]+0.1767766952966368*nuVtSqSum[1]*fc[19]+0.1767766952966368*nuVtSqSum[2]*fr[18]+0.1767766952966368*nuVtSqSum[2]*fc[18]+0.1767766952966368*nuVtSqSum[3]*fr[11]+0.1767766952966368*nuVtSqSum[3]*fc[11]; 

  Gdiff_l[0] = (-0.3827327723098713*nuVtSqSum[3]*fl[20])-0.3827327723098713*nuVtSqSum[3]*fc[20]-0.3827327723098713*nuVtSqSum[2]*fl[13]-0.3827327723098713*nuVtSqSum[2]*fc[13]-0.3827327723098713*nuVtSqSum[1]*fl[12]-0.3827327723098713*nuVtSqSum[1]*fc[12]-0.3977475644174328*nuVtSqSum[3]*fl[6]+0.3977475644174328*nuVtSqSum[3]*fc[6]-0.3827327723098713*nuVtSqSum[0]*fl[5]-0.3827327723098713*nuVtSqSum[0]*fc[5]-0.3977475644174328*fl[2]*nuVtSqSum[2]+0.3977475644174328*fc[2]*nuVtSqSum[2]-0.3977475644174328*fl[1]*nuVtSqSum[1]+0.3977475644174328*fc[1]*nuVtSqSum[1]-0.3977475644174328*fl[0]*nuVtSqSum[0]+0.3977475644174328*fc[0]*nuVtSqSum[0]; 
  Gdiff_l[1] = (-0.3827327723098713*nuVtSqSum[2]*fl[20])-0.3827327723098713*nuVtSqSum[2]*fc[20]-0.3827327723098713*nuVtSqSum[3]*fl[13]-0.3827327723098713*nuVtSqSum[3]*fc[13]-0.3827327723098713*nuVtSqSum[0]*fl[12]-0.3827327723098713*nuVtSqSum[0]*fc[12]-0.3977475644174328*nuVtSqSum[2]*fl[6]+0.3977475644174328*nuVtSqSum[2]*fc[6]-0.3827327723098713*nuVtSqSum[1]*fl[5]-0.3827327723098713*nuVtSqSum[1]*fc[5]-0.3977475644174328*fl[2]*nuVtSqSum[3]+0.3977475644174328*fc[2]*nuVtSqSum[3]-0.3977475644174328*fl[0]*nuVtSqSum[1]+0.3977475644174328*fc[0]*nuVtSqSum[1]-0.3977475644174328*nuVtSqSum[0]*fl[1]+0.3977475644174328*nuVtSqSum[0]*fc[1]; 
  Gdiff_l[2] = (-0.3827327723098713*nuVtSqSum[1]*fl[20])-0.3827327723098713*nuVtSqSum[1]*fc[20]-0.3827327723098713*nuVtSqSum[0]*fl[13]-0.3827327723098713*nuVtSqSum[0]*fc[13]-0.3827327723098713*nuVtSqSum[3]*fl[12]-0.3827327723098713*nuVtSqSum[3]*fc[12]-0.3977475644174328*nuVtSqSum[1]*fl[6]+0.3977475644174328*nuVtSqSum[1]*fc[6]-0.3827327723098713*nuVtSqSum[2]*fl[5]-0.3827327723098713*nuVtSqSum[2]*fc[5]-0.3977475644174328*fl[1]*nuVtSqSum[3]+0.3977475644174328*fc[1]*nuVtSqSum[3]-0.3977475644174328*fl[0]*nuVtSqSum[2]+0.3977475644174328*fc[0]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[0]*fl[2]+0.3977475644174328*nuVtSqSum[0]*fc[2]; 
  Gdiff_l[3] = (-0.3827327723098713*nuVtSqSum[3]*fl[27])-0.3827327723098713*nuVtSqSum[3]*fc[27]-0.3827327723098713*nuVtSqSum[2]*fl[22]-0.3827327723098713*nuVtSqSum[2]*fc[22]-0.3827327723098713*nuVtSqSum[1]*fl[21]-0.3827327723098713*nuVtSqSum[1]*fc[21]-0.3977475644174328*nuVtSqSum[3]*fl[16]+0.3977475644174328*nuVtSqSum[3]*fc[16]-0.3827327723098713*nuVtSqSum[0]*fl[14]-0.3827327723098713*nuVtSqSum[0]*fc[14]-0.3977475644174328*nuVtSqSum[2]*fl[8]+0.3977475644174328*nuVtSqSum[2]*fc[8]-0.3977475644174328*nuVtSqSum[1]*fl[7]+0.3977475644174328*nuVtSqSum[1]*fc[7]-0.3977475644174328*nuVtSqSum[0]*fl[3]+0.3977475644174328*nuVtSqSum[0]*fc[3]; 
  Gdiff_l[4] = (-0.3827327723098713*nuVtSqSum[3]*fl[28])-0.3827327723098713*nuVtSqSum[3]*fc[28]-0.3827327723098713*nuVtSqSum[2]*fl[24]-0.3827327723098713*nuVtSqSum[2]*fc[24]-0.3827327723098713*nuVtSqSum[1]*fl[23]-0.3827327723098713*nuVtSqSum[1]*fc[23]-0.3977475644174328*nuVtSqSum[3]*fl[17]+0.3977475644174328*nuVtSqSum[3]*fc[17]-0.3827327723098713*nuVtSqSum[0]*fl[15]-0.3827327723098713*nuVtSqSum[0]*fc[15]-0.3977475644174328*nuVtSqSum[2]*fl[10]+0.3977475644174328*nuVtSqSum[2]*fc[10]-0.3977475644174328*nuVtSqSum[1]*fl[9]+0.3977475644174328*nuVtSqSum[1]*fc[9]-0.3977475644174328*nuVtSqSum[0]*fl[4]+0.3977475644174328*nuVtSqSum[0]*fc[4]; 
  Gdiff_l[5] = (-0.3827327723098713*nuVtSqSum[0]*fl[20])-0.3827327723098713*nuVtSqSum[0]*fc[20]-0.3827327723098713*nuVtSqSum[1]*fl[13]-0.3827327723098713*nuVtSqSum[1]*fc[13]-0.3827327723098713*nuVtSqSum[2]*fl[12]-0.3827327723098713*nuVtSqSum[2]*fc[12]-0.3977475644174328*nuVtSqSum[0]*fl[6]+0.3977475644174328*nuVtSqSum[0]*fc[6]-0.3827327723098713*nuVtSqSum[3]*fl[5]-0.3827327723098713*nuVtSqSum[3]*fc[5]-0.3977475644174328*fl[0]*nuVtSqSum[3]+0.3977475644174328*fc[0]*nuVtSqSum[3]-0.3977475644174328*fl[1]*nuVtSqSum[2]+0.3977475644174328*fc[1]*nuVtSqSum[2]-0.3977475644174328*nuVtSqSum[1]*fl[2]+0.3977475644174328*nuVtSqSum[1]*fc[2]; 
  Gdiff_l[6] = (-0.3827327723098713*nuVtSqSum[2]*fl[27])-0.3827327723098713*nuVtSqSum[2]*fc[27]-0.3827327723098713*nuVtSqSum[3]*fl[22]-0.3827327723098713*nuVtSqSum[3]*fc[22]-0.3827327723098713*nuVtSqSum[0]*fl[21]-0.3827327723098713*nuVtSqSum[0]*fc[21]-0.3977475644174328*nuVtSqSum[2]*fl[16]+0.3977475644174328*nuVtSqSum[2]*fc[16]-0.3827327723098713*nuVtSqSum[1]*fl[14]-0.3827327723098713*nuVtSqSum[1]*fc[14]-0.3977475644174328*nuVtSqSum[3]*fl[8]+0.3977475644174328*nuVtSqSum[3]*fc[8]-0.3977475644174328*nuVtSqSum[0]*fl[7]+0.3977475644174328*nuVtSqSum[0]*fc[7]-0.3977475644174328*nuVtSqSum[1]*fl[3]+0.3977475644174328*nuVtSqSum[1]*fc[3]; 
  Gdiff_l[7] = (-0.3827327723098713*nuVtSqSum[1]*fl[27])-0.3827327723098713*nuVtSqSum[1]*fc[27]-0.3827327723098713*nuVtSqSum[0]*fl[22]-0.3827327723098713*nuVtSqSum[0]*fc[22]-0.3827327723098713*nuVtSqSum[3]*fl[21]-0.3827327723098713*nuVtSqSum[3]*fc[21]-0.3977475644174328*nuVtSqSum[1]*fl[16]+0.3977475644174328*nuVtSqSum[1]*fc[16]-0.3827327723098713*nuVtSqSum[2]*fl[14]-0.3827327723098713*nuVtSqSum[2]*fc[14]-0.3977475644174328*nuVtSqSum[0]*fl[8]+0.3977475644174328*nuVtSqSum[0]*fc[8]-0.3977475644174328*nuVtSqSum[3]*fl[7]+0.3977475644174328*nuVtSqSum[3]*fc[7]-0.3977475644174328*nuVtSqSum[2]*fl[3]+0.3977475644174328*nuVtSqSum[2]*fc[3]; 
  Gdiff_l[8] = (-0.3827327723098713*nuVtSqSum[2]*fl[28])-0.3827327723098713*nuVtSqSum[2]*fc[28]-0.3827327723098713*nuVtSqSum[3]*fl[24]-0.3827327723098713*nuVtSqSum[3]*fc[24]-0.3827327723098713*nuVtSqSum[0]*fl[23]-0.3827327723098713*nuVtSqSum[0]*fc[23]-0.3977475644174328*nuVtSqSum[2]*fl[17]+0.3977475644174328*nuVtSqSum[2]*fc[17]-0.3827327723098713*nuVtSqSum[1]*fl[15]-0.3827327723098713*nuVtSqSum[1]*fc[15]-0.3977475644174328*nuVtSqSum[3]*fl[10]+0.3977475644174328*nuVtSqSum[3]*fc[10]-0.3977475644174328*nuVtSqSum[0]*fl[9]+0.3977475644174328*nuVtSqSum[0]*fc[9]-0.3977475644174328*nuVtSqSum[1]*fl[4]+0.3977475644174328*nuVtSqSum[1]*fc[4]; 
  Gdiff_l[9] = (-0.3827327723098713*nuVtSqSum[1]*fl[28])-0.3827327723098713*nuVtSqSum[1]*fc[28]-0.3827327723098713*nuVtSqSum[0]*fl[24]-0.3827327723098713*nuVtSqSum[0]*fc[24]-0.3827327723098713*nuVtSqSum[3]*fl[23]-0.3827327723098713*nuVtSqSum[3]*fc[23]-0.3977475644174328*nuVtSqSum[1]*fl[17]+0.3977475644174328*nuVtSqSum[1]*fc[17]-0.3827327723098713*nuVtSqSum[2]*fl[15]-0.3827327723098713*nuVtSqSum[2]*fc[15]-0.3977475644174328*nuVtSqSum[0]*fl[10]+0.3977475644174328*nuVtSqSum[0]*fc[10]-0.3977475644174328*nuVtSqSum[3]*fl[9]+0.3977475644174328*nuVtSqSum[3]*fc[9]-0.3977475644174328*nuVtSqSum[2]*fl[4]+0.3977475644174328*nuVtSqSum[2]*fc[4]; 
  Gdiff_l[10] = (-0.3827327723098713*nuVtSqSum[3]*fl[31])-0.3827327723098713*nuVtSqSum[3]*fc[31]-0.3827327723098713*nuVtSqSum[2]*fl[30]-0.3827327723098713*nuVtSqSum[2]*fc[30]-0.3827327723098713*nuVtSqSum[1]*fl[29]-0.3827327723098713*nuVtSqSum[1]*fc[29]-0.3977475644174328*nuVtSqSum[3]*fl[26]+0.3977475644174328*nuVtSqSum[3]*fc[26]-0.3827327723098713*nuVtSqSum[0]*fl[25]-0.3827327723098713*nuVtSqSum[0]*fc[25]-0.3977475644174328*nuVtSqSum[2]*fl[19]+0.3977475644174328*nuVtSqSum[2]*fc[19]-0.3977475644174328*nuVtSqSum[1]*fl[18]+0.3977475644174328*nuVtSqSum[1]*fc[18]-0.3977475644174328*nuVtSqSum[0]*fl[11]+0.3977475644174328*nuVtSqSum[0]*fc[11]; 
  Gdiff_l[11] = (-0.3827327723098713*nuVtSqSum[0]*fl[27])-0.3827327723098713*nuVtSqSum[0]*fc[27]-0.3827327723098713*nuVtSqSum[1]*fl[22]-0.3827327723098713*nuVtSqSum[1]*fc[22]-0.3827327723098713*nuVtSqSum[2]*fl[21]-0.3827327723098713*nuVtSqSum[2]*fc[21]-0.3977475644174328*nuVtSqSum[0]*fl[16]+0.3977475644174328*nuVtSqSum[0]*fc[16]-0.3827327723098713*nuVtSqSum[3]*fl[14]-0.3827327723098713*nuVtSqSum[3]*fc[14]-0.3977475644174328*nuVtSqSum[1]*fl[8]+0.3977475644174328*nuVtSqSum[1]*fc[8]-0.3977475644174328*nuVtSqSum[2]*fl[7]+0.3977475644174328*nuVtSqSum[2]*fc[7]-0.3977475644174328*fl[3]*nuVtSqSum[3]+0.3977475644174328*fc[3]*nuVtSqSum[3]; 
  Gdiff_l[12] = (-0.3827327723098713*nuVtSqSum[0]*fl[28])-0.3827327723098713*nuVtSqSum[0]*fc[28]-0.3827327723098713*nuVtSqSum[1]*fl[24]-0.3827327723098713*nuVtSqSum[1]*fc[24]-0.3827327723098713*nuVtSqSum[2]*fl[23]-0.3827327723098713*nuVtSqSum[2]*fc[23]-0.3977475644174328*nuVtSqSum[0]*fl[17]+0.3977475644174328*nuVtSqSum[0]*fc[17]-0.3827327723098713*nuVtSqSum[3]*fl[15]-0.3827327723098713*nuVtSqSum[3]*fc[15]-0.3977475644174328*nuVtSqSum[1]*fl[10]+0.3977475644174328*nuVtSqSum[1]*fc[10]-0.3977475644174328*nuVtSqSum[2]*fl[9]+0.3977475644174328*nuVtSqSum[2]*fc[9]-0.3977475644174328*nuVtSqSum[3]*fl[4]+0.3977475644174328*nuVtSqSum[3]*fc[4]; 
  Gdiff_l[13] = (-0.3827327723098713*nuVtSqSum[2]*fl[31])-0.3827327723098713*nuVtSqSum[2]*fc[31]-0.3827327723098713*nuVtSqSum[3]*fl[30]-0.3827327723098713*nuVtSqSum[3]*fc[30]-0.3827327723098713*nuVtSqSum[0]*fl[29]-0.3827327723098713*nuVtSqSum[0]*fc[29]-0.3977475644174328*nuVtSqSum[2]*fl[26]+0.3977475644174328*nuVtSqSum[2]*fc[26]-0.3827327723098713*nuVtSqSum[1]*fl[25]-0.3827327723098713*nuVtSqSum[1]*fc[25]-0.3977475644174328*nuVtSqSum[3]*fl[19]+0.3977475644174328*nuVtSqSum[3]*fc[19]-0.3977475644174328*nuVtSqSum[0]*fl[18]+0.3977475644174328*nuVtSqSum[0]*fc[18]-0.3977475644174328*nuVtSqSum[1]*fl[11]+0.3977475644174328*nuVtSqSum[1]*fc[11]; 
  Gdiff_l[14] = (-0.3827327723098713*nuVtSqSum[1]*fl[31])-0.3827327723098713*nuVtSqSum[1]*fc[31]-0.3827327723098713*nuVtSqSum[0]*fl[30]-0.3827327723098713*nuVtSqSum[0]*fc[30]-0.3827327723098713*nuVtSqSum[3]*fl[29]-0.3827327723098713*nuVtSqSum[3]*fc[29]-0.3977475644174328*nuVtSqSum[1]*fl[26]+0.3977475644174328*nuVtSqSum[1]*fc[26]-0.3827327723098713*nuVtSqSum[2]*fl[25]-0.3827327723098713*nuVtSqSum[2]*fc[25]-0.3977475644174328*nuVtSqSum[0]*fl[19]+0.3977475644174328*nuVtSqSum[0]*fc[19]-0.3977475644174328*nuVtSqSum[3]*fl[18]+0.3977475644174328*nuVtSqSum[3]*fc[18]-0.3977475644174328*nuVtSqSum[2]*fl[11]+0.3977475644174328*nuVtSqSum[2]*fc[11]; 
  Gdiff_l[15] = (-0.3827327723098713*nuVtSqSum[0]*fl[31])-0.3827327723098713*nuVtSqSum[0]*fc[31]-0.3827327723098713*nuVtSqSum[1]*fl[30]-0.3827327723098713*nuVtSqSum[1]*fc[30]-0.3827327723098713*nuVtSqSum[2]*fl[29]-0.3827327723098713*nuVtSqSum[2]*fc[29]-0.3977475644174328*nuVtSqSum[0]*fl[26]+0.3977475644174328*nuVtSqSum[0]*fc[26]-0.3827327723098713*nuVtSqSum[3]*fl[25]-0.3827327723098713*nuVtSqSum[3]*fc[25]-0.3977475644174328*nuVtSqSum[1]*fl[19]+0.3977475644174328*nuVtSqSum[1]*fc[19]-0.3977475644174328*nuVtSqSum[2]*fl[18]+0.3977475644174328*nuVtSqSum[2]*fc[18]-0.3977475644174328*nuVtSqSum[3]*fl[11]+0.3977475644174328*nuVtSqSum[3]*fc[11]; 

  Gdiff_r[0] = (-0.3827327723098713*nuVtSqSum[3]*fr[20])-0.3827327723098713*nuVtSqSum[3]*fc[20]-0.3827327723098713*nuVtSqSum[2]*fr[13]-0.3827327723098713*nuVtSqSum[2]*fc[13]-0.3827327723098713*nuVtSqSum[1]*fr[12]-0.3827327723098713*nuVtSqSum[1]*fc[12]+0.3977475644174328*nuVtSqSum[3]*fr[6]-0.3977475644174328*nuVtSqSum[3]*fc[6]-0.3827327723098713*nuVtSqSum[0]*fr[5]-0.3827327723098713*nuVtSqSum[0]*fc[5]+0.3977475644174328*fr[2]*nuVtSqSum[2]-0.3977475644174328*fc[2]*nuVtSqSum[2]+0.3977475644174328*fr[1]*nuVtSqSum[1]-0.3977475644174328*fc[1]*nuVtSqSum[1]+0.3977475644174328*fr[0]*nuVtSqSum[0]-0.3977475644174328*fc[0]*nuVtSqSum[0]; 
  Gdiff_r[1] = (-0.3827327723098713*nuVtSqSum[2]*fr[20])-0.3827327723098713*nuVtSqSum[2]*fc[20]-0.3827327723098713*nuVtSqSum[3]*fr[13]-0.3827327723098713*nuVtSqSum[3]*fc[13]-0.3827327723098713*nuVtSqSum[0]*fr[12]-0.3827327723098713*nuVtSqSum[0]*fc[12]+0.3977475644174328*nuVtSqSum[2]*fr[6]-0.3977475644174328*nuVtSqSum[2]*fc[6]-0.3827327723098713*nuVtSqSum[1]*fr[5]-0.3827327723098713*nuVtSqSum[1]*fc[5]+0.3977475644174328*fr[2]*nuVtSqSum[3]-0.3977475644174328*fc[2]*nuVtSqSum[3]+0.3977475644174328*fr[0]*nuVtSqSum[1]-0.3977475644174328*fc[0]*nuVtSqSum[1]+0.3977475644174328*nuVtSqSum[0]*fr[1]-0.3977475644174328*nuVtSqSum[0]*fc[1]; 
  Gdiff_r[2] = (-0.3827327723098713*nuVtSqSum[1]*fr[20])-0.3827327723098713*nuVtSqSum[1]*fc[20]-0.3827327723098713*nuVtSqSum[0]*fr[13]-0.3827327723098713*nuVtSqSum[0]*fc[13]-0.3827327723098713*nuVtSqSum[3]*fr[12]-0.3827327723098713*nuVtSqSum[3]*fc[12]+0.3977475644174328*nuVtSqSum[1]*fr[6]-0.3977475644174328*nuVtSqSum[1]*fc[6]-0.3827327723098713*nuVtSqSum[2]*fr[5]-0.3827327723098713*nuVtSqSum[2]*fc[5]+0.3977475644174328*fr[1]*nuVtSqSum[3]-0.3977475644174328*fc[1]*nuVtSqSum[3]+0.3977475644174328*fr[0]*nuVtSqSum[2]-0.3977475644174328*fc[0]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[0]*fr[2]-0.3977475644174328*nuVtSqSum[0]*fc[2]; 
  Gdiff_r[3] = (-0.3827327723098713*nuVtSqSum[3]*fr[27])-0.3827327723098713*nuVtSqSum[3]*fc[27]-0.3827327723098713*nuVtSqSum[2]*fr[22]-0.3827327723098713*nuVtSqSum[2]*fc[22]-0.3827327723098713*nuVtSqSum[1]*fr[21]-0.3827327723098713*nuVtSqSum[1]*fc[21]+0.3977475644174328*nuVtSqSum[3]*fr[16]-0.3977475644174328*nuVtSqSum[3]*fc[16]-0.3827327723098713*nuVtSqSum[0]*fr[14]-0.3827327723098713*nuVtSqSum[0]*fc[14]+0.3977475644174328*nuVtSqSum[2]*fr[8]-0.3977475644174328*nuVtSqSum[2]*fc[8]+0.3977475644174328*nuVtSqSum[1]*fr[7]-0.3977475644174328*nuVtSqSum[1]*fc[7]+0.3977475644174328*nuVtSqSum[0]*fr[3]-0.3977475644174328*nuVtSqSum[0]*fc[3]; 
  Gdiff_r[4] = (-0.3827327723098713*nuVtSqSum[3]*fr[28])-0.3827327723098713*nuVtSqSum[3]*fc[28]-0.3827327723098713*nuVtSqSum[2]*fr[24]-0.3827327723098713*nuVtSqSum[2]*fc[24]-0.3827327723098713*nuVtSqSum[1]*fr[23]-0.3827327723098713*nuVtSqSum[1]*fc[23]+0.3977475644174328*nuVtSqSum[3]*fr[17]-0.3977475644174328*nuVtSqSum[3]*fc[17]-0.3827327723098713*nuVtSqSum[0]*fr[15]-0.3827327723098713*nuVtSqSum[0]*fc[15]+0.3977475644174328*nuVtSqSum[2]*fr[10]-0.3977475644174328*nuVtSqSum[2]*fc[10]+0.3977475644174328*nuVtSqSum[1]*fr[9]-0.3977475644174328*nuVtSqSum[1]*fc[9]+0.3977475644174328*nuVtSqSum[0]*fr[4]-0.3977475644174328*nuVtSqSum[0]*fc[4]; 
  Gdiff_r[5] = (-0.3827327723098713*nuVtSqSum[0]*fr[20])-0.3827327723098713*nuVtSqSum[0]*fc[20]-0.3827327723098713*nuVtSqSum[1]*fr[13]-0.3827327723098713*nuVtSqSum[1]*fc[13]-0.3827327723098713*nuVtSqSum[2]*fr[12]-0.3827327723098713*nuVtSqSum[2]*fc[12]+0.3977475644174328*nuVtSqSum[0]*fr[6]-0.3977475644174328*nuVtSqSum[0]*fc[6]-0.3827327723098713*nuVtSqSum[3]*fr[5]-0.3827327723098713*nuVtSqSum[3]*fc[5]+0.3977475644174328*fr[0]*nuVtSqSum[3]-0.3977475644174328*fc[0]*nuVtSqSum[3]+0.3977475644174328*fr[1]*nuVtSqSum[2]-0.3977475644174328*fc[1]*nuVtSqSum[2]+0.3977475644174328*nuVtSqSum[1]*fr[2]-0.3977475644174328*nuVtSqSum[1]*fc[2]; 
  Gdiff_r[6] = (-0.3827327723098713*nuVtSqSum[2]*fr[27])-0.3827327723098713*nuVtSqSum[2]*fc[27]-0.3827327723098713*nuVtSqSum[3]*fr[22]-0.3827327723098713*nuVtSqSum[3]*fc[22]-0.3827327723098713*nuVtSqSum[0]*fr[21]-0.3827327723098713*nuVtSqSum[0]*fc[21]+0.3977475644174328*nuVtSqSum[2]*fr[16]-0.3977475644174328*nuVtSqSum[2]*fc[16]-0.3827327723098713*nuVtSqSum[1]*fr[14]-0.3827327723098713*nuVtSqSum[1]*fc[14]+0.3977475644174328*nuVtSqSum[3]*fr[8]-0.3977475644174328*nuVtSqSum[3]*fc[8]+0.3977475644174328*nuVtSqSum[0]*fr[7]-0.3977475644174328*nuVtSqSum[0]*fc[7]+0.3977475644174328*nuVtSqSum[1]*fr[3]-0.3977475644174328*nuVtSqSum[1]*fc[3]; 
  Gdiff_r[7] = (-0.3827327723098713*nuVtSqSum[1]*fr[27])-0.3827327723098713*nuVtSqSum[1]*fc[27]-0.3827327723098713*nuVtSqSum[0]*fr[22]-0.3827327723098713*nuVtSqSum[0]*fc[22]-0.3827327723098713*nuVtSqSum[3]*fr[21]-0.3827327723098713*nuVtSqSum[3]*fc[21]+0.3977475644174328*nuVtSqSum[1]*fr[16]-0.3977475644174328*nuVtSqSum[1]*fc[16]-0.3827327723098713*nuVtSqSum[2]*fr[14]-0.3827327723098713*nuVtSqSum[2]*fc[14]+0.3977475644174328*nuVtSqSum[0]*fr[8]-0.3977475644174328*nuVtSqSum[0]*fc[8]+0.3977475644174328*nuVtSqSum[3]*fr[7]-0.3977475644174328*nuVtSqSum[3]*fc[7]+0.3977475644174328*nuVtSqSum[2]*fr[3]-0.3977475644174328*nuVtSqSum[2]*fc[3]; 
  Gdiff_r[8] = (-0.3827327723098713*nuVtSqSum[2]*fr[28])-0.3827327723098713*nuVtSqSum[2]*fc[28]-0.3827327723098713*nuVtSqSum[3]*fr[24]-0.3827327723098713*nuVtSqSum[3]*fc[24]-0.3827327723098713*nuVtSqSum[0]*fr[23]-0.3827327723098713*nuVtSqSum[0]*fc[23]+0.3977475644174328*nuVtSqSum[2]*fr[17]-0.3977475644174328*nuVtSqSum[2]*fc[17]-0.3827327723098713*nuVtSqSum[1]*fr[15]-0.3827327723098713*nuVtSqSum[1]*fc[15]+0.3977475644174328*nuVtSqSum[3]*fr[10]-0.3977475644174328*nuVtSqSum[3]*fc[10]+0.3977475644174328*nuVtSqSum[0]*fr[9]-0.3977475644174328*nuVtSqSum[0]*fc[9]+0.3977475644174328*nuVtSqSum[1]*fr[4]-0.3977475644174328*nuVtSqSum[1]*fc[4]; 
  Gdiff_r[9] = (-0.3827327723098713*nuVtSqSum[1]*fr[28])-0.3827327723098713*nuVtSqSum[1]*fc[28]-0.3827327723098713*nuVtSqSum[0]*fr[24]-0.3827327723098713*nuVtSqSum[0]*fc[24]-0.3827327723098713*nuVtSqSum[3]*fr[23]-0.3827327723098713*nuVtSqSum[3]*fc[23]+0.3977475644174328*nuVtSqSum[1]*fr[17]-0.3977475644174328*nuVtSqSum[1]*fc[17]-0.3827327723098713*nuVtSqSum[2]*fr[15]-0.3827327723098713*nuVtSqSum[2]*fc[15]+0.3977475644174328*nuVtSqSum[0]*fr[10]-0.3977475644174328*nuVtSqSum[0]*fc[10]+0.3977475644174328*nuVtSqSum[3]*fr[9]-0.3977475644174328*nuVtSqSum[3]*fc[9]+0.3977475644174328*nuVtSqSum[2]*fr[4]-0.3977475644174328*nuVtSqSum[2]*fc[4]; 
  Gdiff_r[10] = (-0.3827327723098713*nuVtSqSum[3]*fr[31])-0.3827327723098713*nuVtSqSum[3]*fc[31]-0.3827327723098713*nuVtSqSum[2]*fr[30]-0.3827327723098713*nuVtSqSum[2]*fc[30]-0.3827327723098713*nuVtSqSum[1]*fr[29]-0.3827327723098713*nuVtSqSum[1]*fc[29]+0.3977475644174328*nuVtSqSum[3]*fr[26]-0.3977475644174328*nuVtSqSum[3]*fc[26]-0.3827327723098713*nuVtSqSum[0]*fr[25]-0.3827327723098713*nuVtSqSum[0]*fc[25]+0.3977475644174328*nuVtSqSum[2]*fr[19]-0.3977475644174328*nuVtSqSum[2]*fc[19]+0.3977475644174328*nuVtSqSum[1]*fr[18]-0.3977475644174328*nuVtSqSum[1]*fc[18]+0.3977475644174328*nuVtSqSum[0]*fr[11]-0.3977475644174328*nuVtSqSum[0]*fc[11]; 
  Gdiff_r[11] = (-0.3827327723098713*nuVtSqSum[0]*fr[27])-0.3827327723098713*nuVtSqSum[0]*fc[27]-0.3827327723098713*nuVtSqSum[1]*fr[22]-0.3827327723098713*nuVtSqSum[1]*fc[22]-0.3827327723098713*nuVtSqSum[2]*fr[21]-0.3827327723098713*nuVtSqSum[2]*fc[21]+0.3977475644174328*nuVtSqSum[0]*fr[16]-0.3977475644174328*nuVtSqSum[0]*fc[16]-0.3827327723098713*nuVtSqSum[3]*fr[14]-0.3827327723098713*nuVtSqSum[3]*fc[14]+0.3977475644174328*nuVtSqSum[1]*fr[8]-0.3977475644174328*nuVtSqSum[1]*fc[8]+0.3977475644174328*nuVtSqSum[2]*fr[7]-0.3977475644174328*nuVtSqSum[2]*fc[7]+0.3977475644174328*fr[3]*nuVtSqSum[3]-0.3977475644174328*fc[3]*nuVtSqSum[3]; 
  Gdiff_r[12] = (-0.3827327723098713*nuVtSqSum[0]*fr[28])-0.3827327723098713*nuVtSqSum[0]*fc[28]-0.3827327723098713*nuVtSqSum[1]*fr[24]-0.3827327723098713*nuVtSqSum[1]*fc[24]-0.3827327723098713*nuVtSqSum[2]*fr[23]-0.3827327723098713*nuVtSqSum[2]*fc[23]+0.3977475644174328*nuVtSqSum[0]*fr[17]-0.3977475644174328*nuVtSqSum[0]*fc[17]-0.3827327723098713*nuVtSqSum[3]*fr[15]-0.3827327723098713*nuVtSqSum[3]*fc[15]+0.3977475644174328*nuVtSqSum[1]*fr[10]-0.3977475644174328*nuVtSqSum[1]*fc[10]+0.3977475644174328*nuVtSqSum[2]*fr[9]-0.3977475644174328*nuVtSqSum[2]*fc[9]+0.3977475644174328*nuVtSqSum[3]*fr[4]-0.3977475644174328*nuVtSqSum[3]*fc[4]; 
  Gdiff_r[13] = (-0.3827327723098713*nuVtSqSum[2]*fr[31])-0.3827327723098713*nuVtSqSum[2]*fc[31]-0.3827327723098713*nuVtSqSum[3]*fr[30]-0.3827327723098713*nuVtSqSum[3]*fc[30]-0.3827327723098713*nuVtSqSum[0]*fr[29]-0.3827327723098713*nuVtSqSum[0]*fc[29]+0.3977475644174328*nuVtSqSum[2]*fr[26]-0.3977475644174328*nuVtSqSum[2]*fc[26]-0.3827327723098713*nuVtSqSum[1]*fr[25]-0.3827327723098713*nuVtSqSum[1]*fc[25]+0.3977475644174328*nuVtSqSum[3]*fr[19]-0.3977475644174328*nuVtSqSum[3]*fc[19]+0.3977475644174328*nuVtSqSum[0]*fr[18]-0.3977475644174328*nuVtSqSum[0]*fc[18]+0.3977475644174328*nuVtSqSum[1]*fr[11]-0.3977475644174328*nuVtSqSum[1]*fc[11]; 
  Gdiff_r[14] = (-0.3827327723098713*nuVtSqSum[1]*fr[31])-0.3827327723098713*nuVtSqSum[1]*fc[31]-0.3827327723098713*nuVtSqSum[0]*fr[30]-0.3827327723098713*nuVtSqSum[0]*fc[30]-0.3827327723098713*nuVtSqSum[3]*fr[29]-0.3827327723098713*nuVtSqSum[3]*fc[29]+0.3977475644174328*nuVtSqSum[1]*fr[26]-0.3977475644174328*nuVtSqSum[1]*fc[26]-0.3827327723098713*nuVtSqSum[2]*fr[25]-0.3827327723098713*nuVtSqSum[2]*fc[25]+0.3977475644174328*nuVtSqSum[0]*fr[19]-0.3977475644174328*nuVtSqSum[0]*fc[19]+0.3977475644174328*nuVtSqSum[3]*fr[18]-0.3977475644174328*nuVtSqSum[3]*fc[18]+0.3977475644174328*nuVtSqSum[2]*fr[11]-0.3977475644174328*nuVtSqSum[2]*fc[11]; 
  Gdiff_r[15] = (-0.3827327723098713*nuVtSqSum[0]*fr[31])-0.3827327723098713*nuVtSqSum[0]*fc[31]-0.3827327723098713*nuVtSqSum[1]*fr[30]-0.3827327723098713*nuVtSqSum[1]*fc[30]-0.3827327723098713*nuVtSqSum[2]*fr[29]-0.3827327723098713*nuVtSqSum[2]*fc[29]+0.3977475644174328*nuVtSqSum[0]*fr[26]-0.3977475644174328*nuVtSqSum[0]*fc[26]-0.3827327723098713*nuVtSqSum[3]*fr[25]-0.3827327723098713*nuVtSqSum[3]*fc[25]+0.3977475644174328*nuVtSqSum[1]*fr[19]-0.3977475644174328*nuVtSqSum[1]*fc[19]+0.3977475644174328*nuVtSqSum[2]*fr[18]-0.3977475644174328*nuVtSqSum[2]*fc[18]+0.3977475644174328*nuVtSqSum[3]*fr[11]-0.3977475644174328*nuVtSqSum[3]*fc[11]; 

  Ghat_l[0] = Gdiff_l[0]*rdv2+0.25*alphaDrSurf_l[5]*fUpwind_l[5]+0.25*alphaDrSurf_l[2]*fUpwind_l[2]+0.25*alphaDrSurf_l[1]*fUpwind_l[1]+0.25*alphaDrSurf_l[0]*fUpwind_l[0]; 
  Ghat_l[1] = Gdiff_l[1]*rdv2+0.25*alphaDrSurf_l[2]*fUpwind_l[5]+0.25*fUpwind_l[2]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[1]+0.25*fUpwind_l[0]*alphaDrSurf_l[1]; 
  Ghat_l[2] = Gdiff_l[2]*rdv2+0.25*alphaDrSurf_l[1]*fUpwind_l[5]+0.25*fUpwind_l[1]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[0]*fUpwind_l[2]+0.25*fUpwind_l[0]*alphaDrSurf_l[2]; 
  Ghat_l[3] = Gdiff_l[3]*rdv2+0.25*alphaDrSurf_l[5]*fUpwind_l[11]+0.25*alphaDrSurf_l[2]*fUpwind_l[7]+0.25*alphaDrSurf_l[1]*fUpwind_l[6]+0.25*alphaDrSurf_l[0]*fUpwind_l[3]; 
  Ghat_l[4] = Gdiff_l[4]*rdv2+0.25*alphaDrSurf_l[5]*fUpwind_l[12]+0.25*alphaDrSurf_l[2]*fUpwind_l[9]+0.25*alphaDrSurf_l[1]*fUpwind_l[8]+0.25*alphaDrSurf_l[0]*fUpwind_l[4]; 
  Ghat_l[5] = Gdiff_l[5]*rdv2+0.25*alphaDrSurf_l[0]*fUpwind_l[5]+0.25*fUpwind_l[0]*alphaDrSurf_l[5]+0.25*alphaDrSurf_l[1]*fUpwind_l[2]+0.25*fUpwind_l[1]*alphaDrSurf_l[2]; 
  Ghat_l[6] = Gdiff_l[6]*rdv2+0.25*alphaDrSurf_l[2]*fUpwind_l[11]+0.25*alphaDrSurf_l[5]*fUpwind_l[7]+0.25*alphaDrSurf_l[0]*fUpwind_l[6]+0.25*alphaDrSurf_l[1]*fUpwind_l[3]; 
  Ghat_l[7] = Gdiff_l[7]*rdv2+0.25*alphaDrSurf_l[1]*fUpwind_l[11]+0.25*alphaDrSurf_l[0]*fUpwind_l[7]+0.25*alphaDrSurf_l[5]*fUpwind_l[6]+0.25*alphaDrSurf_l[2]*fUpwind_l[3]; 
  Ghat_l[8] = Gdiff_l[8]*rdv2+0.25*alphaDrSurf_l[2]*fUpwind_l[12]+0.25*alphaDrSurf_l[5]*fUpwind_l[9]+0.25*alphaDrSurf_l[0]*fUpwind_l[8]+0.25*alphaDrSurf_l[1]*fUpwind_l[4]; 
  Ghat_l[9] = Gdiff_l[9]*rdv2+0.25*alphaDrSurf_l[1]*fUpwind_l[12]+0.25*alphaDrSurf_l[0]*fUpwind_l[9]+0.25*alphaDrSurf_l[5]*fUpwind_l[8]+0.25*alphaDrSurf_l[2]*fUpwind_l[4]; 
  Ghat_l[10] = Gdiff_l[10]*rdv2+0.25*alphaDrSurf_l[5]*fUpwind_l[15]+0.25*alphaDrSurf_l[2]*fUpwind_l[14]+0.25*alphaDrSurf_l[1]*fUpwind_l[13]+0.25*alphaDrSurf_l[0]*fUpwind_l[10]; 
  Ghat_l[11] = Gdiff_l[11]*rdv2+0.25*alphaDrSurf_l[0]*fUpwind_l[11]+0.25*alphaDrSurf_l[1]*fUpwind_l[7]+0.25*alphaDrSurf_l[2]*fUpwind_l[6]+0.25*fUpwind_l[3]*alphaDrSurf_l[5]; 
  Ghat_l[12] = Gdiff_l[12]*rdv2+0.25*alphaDrSurf_l[0]*fUpwind_l[12]+0.25*alphaDrSurf_l[1]*fUpwind_l[9]+0.25*alphaDrSurf_l[2]*fUpwind_l[8]+0.25*fUpwind_l[4]*alphaDrSurf_l[5]; 
  Ghat_l[13] = Gdiff_l[13]*rdv2+0.25*alphaDrSurf_l[2]*fUpwind_l[15]+0.25*alphaDrSurf_l[5]*fUpwind_l[14]+0.25*alphaDrSurf_l[0]*fUpwind_l[13]+0.25*alphaDrSurf_l[1]*fUpwind_l[10]; 
  Ghat_l[14] = Gdiff_l[14]*rdv2+0.25*alphaDrSurf_l[1]*fUpwind_l[15]+0.25*alphaDrSurf_l[0]*fUpwind_l[14]+0.25*alphaDrSurf_l[5]*fUpwind_l[13]+0.25*alphaDrSurf_l[2]*fUpwind_l[10]; 
  Ghat_l[15] = Gdiff_l[15]*rdv2+0.25*alphaDrSurf_l[0]*fUpwind_l[15]+0.25*alphaDrSurf_l[1]*fUpwind_l[14]+0.25*alphaDrSurf_l[2]*fUpwind_l[13]+0.25*alphaDrSurf_l[5]*fUpwind_l[10]; 

  Ghat_r[0] = Gdiff_r[0]*rdv2+0.25*alphaDrSurf_r[5]*fUpwind_r[5]+0.25*alphaDrSurf_r[2]*fUpwind_r[2]+0.25*alphaDrSurf_r[1]*fUpwind_r[1]+0.25*alphaDrSurf_r[0]*fUpwind_r[0]; 
  Ghat_r[1] = Gdiff_r[1]*rdv2+0.25*alphaDrSurf_r[2]*fUpwind_r[5]+0.25*fUpwind_r[2]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[1]+0.25*fUpwind_r[0]*alphaDrSurf_r[1]; 
  Ghat_r[2] = Gdiff_r[2]*rdv2+0.25*alphaDrSurf_r[1]*fUpwind_r[5]+0.25*fUpwind_r[1]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[0]*fUpwind_r[2]+0.25*fUpwind_r[0]*alphaDrSurf_r[2]; 
  Ghat_r[3] = Gdiff_r[3]*rdv2+0.25*alphaDrSurf_r[5]*fUpwind_r[11]+0.25*alphaDrSurf_r[2]*fUpwind_r[7]+0.25*alphaDrSurf_r[1]*fUpwind_r[6]+0.25*alphaDrSurf_r[0]*fUpwind_r[3]; 
  Ghat_r[4] = Gdiff_r[4]*rdv2+0.25*alphaDrSurf_r[5]*fUpwind_r[12]+0.25*alphaDrSurf_r[2]*fUpwind_r[9]+0.25*alphaDrSurf_r[1]*fUpwind_r[8]+0.25*alphaDrSurf_r[0]*fUpwind_r[4]; 
  Ghat_r[5] = Gdiff_r[5]*rdv2+0.25*alphaDrSurf_r[0]*fUpwind_r[5]+0.25*fUpwind_r[0]*alphaDrSurf_r[5]+0.25*alphaDrSurf_r[1]*fUpwind_r[2]+0.25*fUpwind_r[1]*alphaDrSurf_r[2]; 
  Ghat_r[6] = Gdiff_r[6]*rdv2+0.25*alphaDrSurf_r[2]*fUpwind_r[11]+0.25*alphaDrSurf_r[5]*fUpwind_r[7]+0.25*alphaDrSurf_r[0]*fUpwind_r[6]+0.25*alphaDrSurf_r[1]*fUpwind_r[3]; 
  Ghat_r[7] = Gdiff_r[7]*rdv2+0.25*alphaDrSurf_r[1]*fUpwind_r[11]+0.25*alphaDrSurf_r[0]*fUpwind_r[7]+0.25*alphaDrSurf_r[5]*fUpwind_r[6]+0.25*alphaDrSurf_r[2]*fUpwind_r[3]; 
  Ghat_r[8] = Gdiff_r[8]*rdv2+0.25*alphaDrSurf_r[2]*fUpwind_r[12]+0.25*alphaDrSurf_r[5]*fUpwind_r[9]+0.25*alphaDrSurf_r[0]*fUpwind_r[8]+0.25*alphaDrSurf_r[1]*fUpwind_r[4]; 
  Ghat_r[9] = Gdiff_r[9]*rdv2+0.25*alphaDrSurf_r[1]*fUpwind_r[12]+0.25*alphaDrSurf_r[0]*fUpwind_r[9]+0.25*alphaDrSurf_r[5]*fUpwind_r[8]+0.25*alphaDrSurf_r[2]*fUpwind_r[4]; 
  Ghat_r[10] = Gdiff_r[10]*rdv2+0.25*alphaDrSurf_r[5]*fUpwind_r[15]+0.25*alphaDrSurf_r[2]*fUpwind_r[14]+0.25*alphaDrSurf_r[1]*fUpwind_r[13]+0.25*alphaDrSurf_r[0]*fUpwind_r[10]; 
  Ghat_r[11] = Gdiff_r[11]*rdv2+0.25*alphaDrSurf_r[0]*fUpwind_r[11]+0.25*alphaDrSurf_r[1]*fUpwind_r[7]+0.25*alphaDrSurf_r[2]*fUpwind_r[6]+0.25*fUpwind_r[3]*alphaDrSurf_r[5]; 
  Ghat_r[12] = Gdiff_r[12]*rdv2+0.25*alphaDrSurf_r[0]*fUpwind_r[12]+0.25*alphaDrSurf_r[1]*fUpwind_r[9]+0.25*alphaDrSurf_r[2]*fUpwind_r[8]+0.25*fUpwind_r[4]*alphaDrSurf_r[5]; 
  Ghat_r[13] = Gdiff_r[13]*rdv2+0.25*alphaDrSurf_r[2]*fUpwind_r[15]+0.25*alphaDrSurf_r[5]*fUpwind_r[14]+0.25*alphaDrSurf_r[0]*fUpwind_r[13]+0.25*alphaDrSurf_r[1]*fUpwind_r[10]; 
  Ghat_r[14] = Gdiff_r[14]*rdv2+0.25*alphaDrSurf_r[1]*fUpwind_r[15]+0.25*alphaDrSurf_r[0]*fUpwind_r[14]+0.25*alphaDrSurf_r[5]*fUpwind_r[13]+0.25*alphaDrSurf_r[2]*fUpwind_r[10]; 
  Ghat_r[15] = Gdiff_r[15]*rdv2+0.25*alphaDrSurf_r[0]*fUpwind_r[15]+0.25*alphaDrSurf_r[1]*fUpwind_r[14]+0.25*alphaDrSurf_r[2]*fUpwind_r[13]+0.25*alphaDrSurf_r[5]*fUpwind_r[10]; 

  out[0] += 0.7071067811865475*Ghat_l[0]*rdv2-0.7071067811865475*Ghat_r[0]*rdv2; 
  out[1] += 0.7071067811865475*Ghat_l[1]*rdv2-0.7071067811865475*Ghat_r[1]*rdv2; 
  out[2] += 0.7071067811865475*Ghat_l[2]*rdv2-0.7071067811865475*Ghat_r[2]*rdv2; 
  out[3] += 0.7071067811865475*Ghat_l[3]*rdv2-0.7071067811865475*Ghat_r[3]*rdv2; 
  out[4] += 0.7071067811865475*Ghat_l[4]*rdv2-0.7071067811865475*Ghat_r[4]*rdv2; 
  out[5] += 1.224744871391589*Gdiff2_r[0]*rdvSq4-1.224744871391589*Gdiff2_l[0]*rdvSq4-1.224744871391589*Ghat_r[0]*rdv2-1.224744871391589*Ghat_l[0]*rdv2; 
  out[6] += 0.7071067811865475*Ghat_l[5]*rdv2-0.7071067811865475*Ghat_r[5]*rdv2; 
  out[7] += 0.7071067811865475*Ghat_l[6]*rdv2-0.7071067811865475*Ghat_r[6]*rdv2; 
  out[8] += 0.7071067811865475*Ghat_l[7]*rdv2-0.7071067811865475*Ghat_r[7]*rdv2; 
  out[9] += 0.7071067811865475*Ghat_l[8]*rdv2-0.7071067811865475*Ghat_r[8]*rdv2; 
  out[10] += 0.7071067811865475*Ghat_l[9]*rdv2-0.7071067811865475*Ghat_r[9]*rdv2; 
  out[11] += 0.7071067811865475*Ghat_l[10]*rdv2-0.7071067811865475*Ghat_r[10]*rdv2; 
  out[12] += 1.224744871391589*Gdiff2_r[1]*rdvSq4-1.224744871391589*Gdiff2_l[1]*rdvSq4-1.224744871391589*Ghat_r[1]*rdv2-1.224744871391589*Ghat_l[1]*rdv2; 
  out[13] += 1.224744871391589*Gdiff2_r[2]*rdvSq4-1.224744871391589*Gdiff2_l[2]*rdvSq4-1.224744871391589*Ghat_r[2]*rdv2-1.224744871391589*Ghat_l[2]*rdv2; 
  out[14] += 1.224744871391589*Gdiff2_r[3]*rdvSq4-1.224744871391589*Gdiff2_l[3]*rdvSq4-1.224744871391589*Ghat_r[3]*rdv2-1.224744871391589*Ghat_l[3]*rdv2; 
  out[15] += 1.224744871391589*Gdiff2_r[4]*rdvSq4-1.224744871391589*Gdiff2_l[4]*rdvSq4-1.224744871391589*Ghat_r[4]*rdv2-1.224744871391589*Ghat_l[4]*rdv2; 
  out[16] += 0.7071067811865475*Ghat_l[11]*rdv2-0.7071067811865475*Ghat_r[11]*rdv2; 
  out[17] += 0.7071067811865475*Ghat_l[12]*rdv2-0.7071067811865475*Ghat_r[12]*rdv2; 
  out[18] += 0.7071067811865475*Ghat_l[13]*rdv2-0.7071067811865475*Ghat_r[13]*rdv2; 
  out[19] += 0.7071067811865475*Ghat_l[14]*rdv2-0.7071067811865475*Ghat_r[14]*rdv2; 
  out[20] += 1.224744871391589*Gdiff2_r[5]*rdvSq4-1.224744871391589*Gdiff2_l[5]*rdvSq4-1.224744871391589*Ghat_r[5]*rdv2-1.224744871391589*Ghat_l[5]*rdv2; 
  out[21] += 1.224744871391589*Gdiff2_r[6]*rdvSq4-1.224744871391589*Gdiff2_l[6]*rdvSq4-1.224744871391589*Ghat_r[6]*rdv2-1.224744871391589*Ghat_l[6]*rdv2; 
  out[22] += 1.224744871391589*Gdiff2_r[7]*rdvSq4-1.224744871391589*Gdiff2_l[7]*rdvSq4-1.224744871391589*Ghat_r[7]*rdv2-1.224744871391589*Ghat_l[7]*rdv2; 
  out[23] += 1.224744871391589*Gdiff2_r[8]*rdvSq4-1.224744871391589*Gdiff2_l[8]*rdvSq4-1.224744871391589*Ghat_r[8]*rdv2-1.224744871391589*Ghat_l[8]*rdv2; 
  out[24] += 1.224744871391589*Gdiff2_r[9]*rdvSq4-1.224744871391589*Gdiff2_l[9]*rdvSq4-1.224744871391589*Ghat_r[9]*rdv2-1.224744871391589*Ghat_l[9]*rdv2; 
  out[25] += 1.224744871391589*Gdiff2_r[10]*rdvSq4-1.224744871391589*Gdiff2_l[10]*rdvSq4-1.224744871391589*Ghat_r[10]*rdv2-1.224744871391589*Ghat_l[10]*rdv2; 
  out[26] += 0.7071067811865475*Ghat_l[15]*rdv2-0.7071067811865475*Ghat_r[15]*rdv2; 
  out[27] += 1.224744871391589*Gdiff2_r[11]*rdvSq4-1.224744871391589*Gdiff2_l[11]*rdvSq4-1.224744871391589*Ghat_r[11]*rdv2-1.224744871391589*Ghat_l[11]*rdv2; 
  out[28] += 1.224744871391589*Gdiff2_r[12]*rdvSq4-1.224744871391589*Gdiff2_l[12]*rdvSq4-1.224744871391589*Ghat_r[12]*rdv2-1.224744871391589*Ghat_l[12]*rdv2; 
  out[29] += 1.224744871391589*Gdiff2_r[13]*rdvSq4-1.224744871391589*Gdiff2_l[13]*rdvSq4-1.224744871391589*Ghat_r[13]*rdv2-1.224744871391589*Ghat_l[13]*rdv2; 
  out[30] += 1.224744871391589*Gdiff2_r[14]*rdvSq4-1.224744871391589*Gdiff2_l[14]*rdvSq4-1.224744871391589*Ghat_r[14]*rdv2-1.224744871391589*Ghat_l[14]*rdv2; 
  out[31] += 1.224744871391589*Gdiff2_r[15]*rdvSq4-1.224744871391589*Gdiff2_l[15]*rdvSq4-1.224744871391589*Ghat_r[15]*rdv2-1.224744871391589*Ghat_l[15]*rdv2; 
} 
