#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_diff_coeff_diag_1x3v_vz_ser_p1_upvz(const double *dxv, const double *gamma, const double* fpo_g_stencil[3], const double* fpo_d2gdv2_surf, double *diff_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_g_stencil[3]: 3 cell stencil of Rosenbluth potential G. 
  // fpo_d2gdv2_surf: Surface projection of d2G/dv2 in center cell. 
  // diff_coeff: Output array for diffusion tensor. 

  // Use cell-average value for gamma. 
 double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1_sq = 4.0/dxv[3]/dxv[3]; 

  const double* G_L = fpo_g_stencil[0]; 
  const double* G_C = fpo_g_stencil[1]; 
  
  const double *fpo_d2g_surf_C_vxvx = &fpo_d2gdv2_surf[0]; 
  const double *fpo_d2g_surf_C_vyvy = &fpo_d2gdv2_surf[32]; 
  const double *fpo_d2g_surf_C_vzvz = &fpo_d2gdv2_surf[64]; 
  
  const double* d2G_surf_C = fpo_d2g_surf_C_vzvz; 
  
  double *diff_coeff_vxvx = &diff_coeff[0]; 
  double *diff_coeff_vyvy = &diff_coeff[160]; 
  double *diff_coeff_vzvz = &diff_coeff[320]; 
  
  double *out = diff_coeff_vzvz; 
  
  out[0] = 1.4433756729740645*G_L[4]*dv1_sq*gamma_avg+2.886751345948129*G_C[4]*dv1_sq*gamma_avg+1.25*G_L[0]*dv1_sq*gamma_avg-1.25*G_C[0]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[0]*gamma_avg; 
  out[1] = 1.4433756729740645*G_L[8]*dv1_sq*gamma_avg+2.886751345948129*G_C[8]*dv1_sq*gamma_avg+1.25*G_L[1]*dv1_sq*gamma_avg-1.25*G_C[1]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[1]*gamma_avg; 
  out[2] = 1.4433756729740645*G_L[9]*dv1_sq*gamma_avg+2.886751345948129*G_C[9]*dv1_sq*gamma_avg+1.25*G_L[2]*dv1_sq*gamma_avg-1.25*G_C[2]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[2]*gamma_avg; 
  out[3] = 1.4433756729740645*G_L[10]*dv1_sq*gamma_avg+2.886751345948129*G_C[10]*dv1_sq*gamma_avg+1.25*G_L[3]*dv1_sq*gamma_avg-1.25*G_C[3]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[3]*gamma_avg; 
  out[4] = 0.5208333333333334*G_L[4]*dv1_sq*gamma_avg+0.10416666666666667*G_C[4]*dv1_sq*gamma_avg+0.18042195912175807*G_L[0]*dv1_sq*gamma_avg-0.18042195912175807*G_C[0]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[0]*gamma_avg; 
  out[5] = 1.4433756729740645*G_L[12]*dv1_sq*gamma_avg+2.886751345948129*G_C[12]*dv1_sq*gamma_avg+1.25*G_L[5]*dv1_sq*gamma_avg-1.25*G_C[5]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[4]*gamma_avg; 
  out[6] = 1.4433756729740645*G_L[13]*dv1_sq*gamma_avg+2.886751345948129*G_C[13]*dv1_sq*gamma_avg+1.25*G_L[6]*dv1_sq*gamma_avg-1.25*G_C[6]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[5]*gamma_avg; 
  out[7] = 1.4433756729740645*G_L[14]*dv1_sq*gamma_avg+2.886751345948129*G_C[14]*dv1_sq*gamma_avg+1.25*G_L[7]*dv1_sq*gamma_avg-1.25*G_C[7]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[6]*gamma_avg; 
  out[8] = 0.5208333333333334*G_L[8]*dv1_sq*gamma_avg+0.10416666666666667*G_C[8]*dv1_sq*gamma_avg+0.18042195912175807*G_L[1]*dv1_sq*gamma_avg-0.18042195912175807*G_C[1]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[1]*gamma_avg; 
  out[9] = 0.5208333333333334*G_L[9]*dv1_sq*gamma_avg+0.10416666666666667*G_C[9]*dv1_sq*gamma_avg+0.18042195912175807*G_L[2]*dv1_sq*gamma_avg-0.18042195912175807*G_C[2]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[2]*gamma_avg; 
  out[10] = 0.5208333333333334*G_L[10]*dv1_sq*gamma_avg+0.10416666666666667*G_C[10]*dv1_sq*gamma_avg+0.18042195912175807*G_L[3]*dv1_sq*gamma_avg-0.18042195912175807*G_C[3]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[3]*gamma_avg; 
  out[11] = 1.4433756729740645*G_L[15]*dv1_sq*gamma_avg+2.886751345948129*G_C[15]*dv1_sq*gamma_avg+1.25*G_L[11]*dv1_sq*gamma_avg-1.25*G_C[11]*dv1_sq*gamma_avg+0.2357022603955158*d2G_surf_C[7]*gamma_avg; 
  out[12] = 0.5208333333333334*G_L[12]*dv1_sq*gamma_avg+0.10416666666666667*G_C[12]*dv1_sq*gamma_avg+0.18042195912175807*G_L[5]*dv1_sq*gamma_avg-0.18042195912175807*G_C[5]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[4]*gamma_avg; 
  out[13] = 0.5208333333333334*G_L[13]*dv1_sq*gamma_avg+0.10416666666666667*G_C[13]*dv1_sq*gamma_avg+0.18042195912175807*G_L[6]*dv1_sq*gamma_avg-0.18042195912175807*G_C[6]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[5]*gamma_avg; 
  out[14] = 0.5208333333333334*G_L[14]*dv1_sq*gamma_avg+0.10416666666666667*G_C[14]*dv1_sq*gamma_avg+0.18042195912175807*G_L[7]*dv1_sq*gamma_avg-0.18042195912175807*G_C[7]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[6]*gamma_avg; 
  out[15] = 0.5208333333333334*G_L[15]*dv1_sq*gamma_avg+0.10416666666666667*G_C[15]*dv1_sq*gamma_avg+0.18042195912175807*G_L[11]*dv1_sq*gamma_avg-0.18042195912175807*G_C[11]*dv1_sq*gamma_avg+0.34020690871988585*d2G_surf_C[7]*gamma_avg; 
  out[32] = -(1.0489329895978423*G_L[4]*dv1_sq*gamma_avg)-1.3716816017817939*G_C[4]*dv1_sq*gamma_avg-0.6987712429686844*G_L[0]*dv1_sq*gamma_avg+0.6987712429686844*G_C[0]*dv1_sq*gamma_avg+0.2635231383473649*d2G_surf_C[0]*gamma_avg; 
  out[33] = -(1.048932989597842*G_L[8]*dv1_sq*gamma_avg)-1.3716816017817937*G_C[8]*dv1_sq*gamma_avg-0.6987712429686844*G_L[1]*dv1_sq*gamma_avg+0.6987712429686844*G_C[1]*dv1_sq*gamma_avg+0.26352313834736496*d2G_surf_C[1]*gamma_avg; 
  out[34] = -(1.048932989597842*G_L[9]*dv1_sq*gamma_avg)-1.3716816017817937*G_C[9]*dv1_sq*gamma_avg-0.6987712429686844*G_L[2]*dv1_sq*gamma_avg+0.6987712429686844*G_C[2]*dv1_sq*gamma_avg+0.26352313834736496*d2G_surf_C[2]*gamma_avg; 
  out[35] = -(1.048932989597842*G_L[10]*dv1_sq*gamma_avg)-1.3716816017817937*G_C[10]*dv1_sq*gamma_avg-0.6987712429686844*G_L[3]*dv1_sq*gamma_avg+0.6987712429686844*G_C[3]*dv1_sq*gamma_avg+0.26352313834736496*d2G_surf_C[3]*gamma_avg; 
  out[36] = -(1.0489329895978423*G_L[12]*dv1_sq*gamma_avg)-1.3716816017817939*G_C[12]*dv1_sq*gamma_avg-0.6987712429686844*G_L[5]*dv1_sq*gamma_avg+0.6987712429686844*G_C[5]*dv1_sq*gamma_avg+0.2635231383473649*d2G_surf_C[4]*gamma_avg; 
  out[37] = -(1.0489329895978423*G_L[13]*dv1_sq*gamma_avg)-1.3716816017817939*G_C[13]*dv1_sq*gamma_avg-0.6987712429686844*G_L[6]*dv1_sq*gamma_avg+0.6987712429686844*G_C[6]*dv1_sq*gamma_avg+0.2635231383473649*d2G_surf_C[5]*gamma_avg; 
  out[38] = -(1.0489329895978423*G_L[14]*dv1_sq*gamma_avg)-1.3716816017817939*G_C[14]*dv1_sq*gamma_avg-0.6987712429686844*G_L[7]*dv1_sq*gamma_avg+0.6987712429686844*G_C[7]*dv1_sq*gamma_avg+0.2635231383473649*d2G_surf_C[6]*gamma_avg; 
  out[39] = -(1.048932989597842*G_L[15]*dv1_sq*gamma_avg)-1.3716816017817937*G_C[15]*dv1_sq*gamma_avg-0.6987712429686844*G_L[11]*dv1_sq*gamma_avg+0.6987712429686844*G_C[11]*dv1_sq*gamma_avg+0.26352313834736496*d2G_surf_C[7]*gamma_avg; 
} 
