#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vx_ser_p1_invx(const double *dxv, const double gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  double dv1_sq = 4.0/dxv[1]/dxv[1]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  const double* G_R = fpo_g_stencil[2]; 
  
  const double *fpo_d2g_surf_C_vxvx = &fpo_d2gdv2_surf[0]; 
  const double *fpo_d2g_surf_C_vyvy = &fpo_d2gdv2_surf[32]; 
  const double *fpo_d2g_surf_C_vzvz = &fpo_d2gdv2_surf[64]; 
  
  const double* d2G_surf_C = fpo_d2g_surf_C_vxvx; 
  
  double *diff_coeff_vxvx = &diff_coeff[0]; 
  double *diff_coeff_vyvy = &diff_coeff[160]; 
  double *diff_coeff_vzvz = &diff_coeff[320]; 
  
  double *out = diff_coeff_vxvx; 
  
  out[0] = -(0.5412658773652741*G_R[2]*dv1_sq*gamma)+0.5412658773652741*G_L[2]*dv1_sq*gamma+0.5625*G_R[0]*dv1_sq*gamma+0.5625*G_L[0]*dv1_sq*gamma-1.125*G_C[0]*dv1_sq*gamma; 
  out[1] = -(0.5412658773652741*G_R[5]*dv1_sq*gamma)+0.5412658773652741*G_L[5]*dv1_sq*gamma+0.5625*G_R[1]*dv1_sq*gamma+0.5625*G_L[1]*dv1_sq*gamma-1.125*G_C[1]*dv1_sq*gamma; 
  out[2] = -(0.4375*G_R[2]*dv1_sq*gamma)-0.4375*G_L[2]*dv1_sq*gamma-2.875*G_C[2]*dv1_sq*gamma+0.5412658773652741*G_R[0]*dv1_sq*gamma-0.5412658773652741*G_L[0]*dv1_sq*gamma; 
  out[3] = -(0.5412658773652741*G_R[7]*dv1_sq*gamma)+0.5412658773652741*G_L[7]*dv1_sq*gamma+0.5625*G_R[3]*dv1_sq*gamma+0.5625*G_L[3]*dv1_sq*gamma-1.125*G_C[3]*dv1_sq*gamma; 
  out[4] = -(0.5412658773652741*G_R[9]*dv1_sq*gamma)+0.5412658773652741*G_L[9]*dv1_sq*gamma+0.5625*G_R[4]*dv1_sq*gamma+0.5625*G_L[4]*dv1_sq*gamma-1.125*G_C[4]*dv1_sq*gamma; 
  out[5] = -(0.4375*G_R[5]*dv1_sq*gamma)-0.4375*G_L[5]*dv1_sq*gamma-2.875*G_C[5]*dv1_sq*gamma+0.5412658773652741*G_R[1]*dv1_sq*gamma-0.5412658773652741*G_L[1]*dv1_sq*gamma; 
  out[6] = -(0.5412658773652741*G_R[11]*dv1_sq*gamma)+0.5412658773652741*G_L[11]*dv1_sq*gamma+0.5625*G_R[6]*dv1_sq*gamma+0.5625*G_L[6]*dv1_sq*gamma-1.125*G_C[6]*dv1_sq*gamma; 
  out[7] = -(0.4375*G_R[7]*dv1_sq*gamma)-0.4375*G_L[7]*dv1_sq*gamma-2.875*G_C[7]*dv1_sq*gamma+0.5412658773652741*G_R[3]*dv1_sq*gamma-0.5412658773652741*G_L[3]*dv1_sq*gamma; 
  out[8] = -(0.5412658773652741*G_R[12]*dv1_sq*gamma)+0.5412658773652741*G_L[12]*dv1_sq*gamma+0.5625*G_R[8]*dv1_sq*gamma+0.5625*G_L[8]*dv1_sq*gamma-1.125*G_C[8]*dv1_sq*gamma; 
  out[9] = -(0.4375*G_R[9]*dv1_sq*gamma)-0.4375*G_L[9]*dv1_sq*gamma-2.875*G_C[9]*dv1_sq*gamma+0.5412658773652741*G_R[4]*dv1_sq*gamma-0.5412658773652741*G_L[4]*dv1_sq*gamma; 
  out[10] = -(0.5412658773652741*G_R[14]*dv1_sq*gamma)+0.5412658773652741*G_L[14]*dv1_sq*gamma+0.5625*G_R[10]*dv1_sq*gamma+0.5625*G_L[10]*dv1_sq*gamma-1.125*G_C[10]*dv1_sq*gamma; 
  out[11] = -(0.4375*G_R[11]*dv1_sq*gamma)-0.4375*G_L[11]*dv1_sq*gamma-2.875*G_C[11]*dv1_sq*gamma+0.5412658773652741*G_R[6]*dv1_sq*gamma-0.5412658773652741*G_L[6]*dv1_sq*gamma; 
  out[12] = -(0.4375*G_R[12]*dv1_sq*gamma)-0.4375*G_L[12]*dv1_sq*gamma-2.875*G_C[12]*dv1_sq*gamma+0.5412658773652741*G_R[8]*dv1_sq*gamma-0.5412658773652741*G_L[8]*dv1_sq*gamma; 
  out[13] = -(0.5412658773652741*G_R[15]*dv1_sq*gamma)+0.5412658773652741*G_L[15]*dv1_sq*gamma+0.5625*G_R[13]*dv1_sq*gamma+0.5625*G_L[13]*dv1_sq*gamma-1.125*G_C[13]*dv1_sq*gamma; 
  out[14] = -(0.4375*G_R[14]*dv1_sq*gamma)-0.4375*G_L[14]*dv1_sq*gamma-2.875*G_C[14]*dv1_sq*gamma+0.5412658773652741*G_R[10]*dv1_sq*gamma-0.5412658773652741*G_L[10]*dv1_sq*gamma; 
  out[15] = -(0.4375*G_R[15]*dv1_sq*gamma)-0.4375*G_L[15]*dv1_sq*gamma-2.875*G_C[15]*dv1_sq*gamma+0.5412658773652741*G_R[13]*dv1_sq*gamma-0.5412658773652741*G_L[13]*dv1_sq*gamma; 
  out[16] = 0.7261843774138906*G_R[2]*dv1_sq*gamma-0.7261843774138906*G_L[2]*dv1_sq*gamma-0.4192627457812106*G_R[0]*dv1_sq*gamma-0.4192627457812106*G_L[0]*dv1_sq*gamma+0.8385254915624212*G_C[0]*dv1_sq*gamma; 
  out[17] = 0.7261843774138907*G_R[5]*dv1_sq*gamma-0.7261843774138907*G_L[5]*dv1_sq*gamma-0.41926274578121053*G_R[1]*dv1_sq*gamma-0.41926274578121053*G_L[1]*dv1_sq*gamma+0.8385254915624211*G_C[1]*dv1_sq*gamma; 
  out[18] = 0.7261843774138907*G_R[7]*dv1_sq*gamma-0.7261843774138907*G_L[7]*dv1_sq*gamma-0.41926274578121053*G_R[3]*dv1_sq*gamma-0.41926274578121053*G_L[3]*dv1_sq*gamma+0.8385254915624211*G_C[3]*dv1_sq*gamma; 
  out[19] = 0.7261843774138907*G_R[9]*dv1_sq*gamma-0.7261843774138907*G_L[9]*dv1_sq*gamma-0.41926274578121053*G_R[4]*dv1_sq*gamma-0.41926274578121053*G_L[4]*dv1_sq*gamma+0.8385254915624211*G_C[4]*dv1_sq*gamma; 
  out[20] = 0.7261843774138906*G_R[11]*dv1_sq*gamma-0.7261843774138906*G_L[11]*dv1_sq*gamma-0.4192627457812106*G_R[6]*dv1_sq*gamma-0.4192627457812106*G_L[6]*dv1_sq*gamma+0.8385254915624212*G_C[6]*dv1_sq*gamma; 
  out[21] = 0.7261843774138906*G_R[12]*dv1_sq*gamma-0.7261843774138906*G_L[12]*dv1_sq*gamma-0.4192627457812106*G_R[8]*dv1_sq*gamma-0.4192627457812106*G_L[8]*dv1_sq*gamma+0.8385254915624212*G_C[8]*dv1_sq*gamma; 
  out[22] = 0.7261843774138906*G_R[14]*dv1_sq*gamma-0.7261843774138906*G_L[14]*dv1_sq*gamma-0.4192627457812106*G_R[10]*dv1_sq*gamma-0.4192627457812106*G_L[10]*dv1_sq*gamma+0.8385254915624212*G_C[10]*dv1_sq*gamma; 
  out[23] = 0.7261843774138907*G_R[15]*dv1_sq*gamma-0.7261843774138907*G_L[15]*dv1_sq*gamma-0.41926274578121053*G_R[13]*dv1_sq*gamma-0.41926274578121053*G_L[13]*dv1_sq*gamma+0.8385254915624211*G_C[13]*dv1_sq*gamma; 
} 
