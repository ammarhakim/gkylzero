#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  // Use cell-average value for gamma. 
 double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[2]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[8]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[16]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vy; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[40]; 
  double *drag_coeff_vz = &drag_coeff[80]; 
  
  double *out = drag_coeff_vy; 
  
  out[0] = -(0.180421959121758*H_L[3]*dv1*gamma_avg)+2.70632938682637*H_C[3]*dv1*gamma_avg-0.1875*H_L[0]*dv1*gamma_avg+0.1875*H_C[0]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[0]*gamma_avg; 
  out[1] = -(0.180421959121758*H_L[6]*dv1*gamma_avg)+2.70632938682637*H_C[6]*dv1*gamma_avg-0.1875*H_L[1]*dv1*gamma_avg+0.1875*H_C[1]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[1]*gamma_avg; 
  out[2] = -(0.180421959121758*H_L[7]*dv1*gamma_avg)+2.70632938682637*H_C[7]*dv1*gamma_avg-0.1875*H_L[2]*dv1*gamma_avg+0.1875*H_C[2]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[2]*gamma_avg; 
  out[3] = 1.6875*H_L[3]*dv1*gamma_avg+2.6875*H_C[3]*dv1*gamma_avg+1.4072912811497125*H_L[0]*dv1*gamma_avg-1.4072912811497125*H_C[0]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[0]*gamma_avg; 
  out[4] = -(0.180421959121758*H_L[10]*dv1*gamma_avg)+2.70632938682637*H_C[10]*dv1*gamma_avg-0.1875*H_L[4]*dv1*gamma_avg+0.1875*H_C[4]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[3]*gamma_avg; 
  out[5] = -(0.180421959121758*H_L[11]*dv1*gamma_avg)+2.70632938682637*H_C[11]*dv1*gamma_avg-0.1875*H_L[5]*dv1*gamma_avg+0.1875*H_C[5]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[4]*gamma_avg; 
  out[6] = 1.6875*H_L[6]*dv1*gamma_avg+2.6875*H_C[6]*dv1*gamma_avg+1.4072912811497125*H_L[1]*dv1*gamma_avg-1.4072912811497125*H_C[1]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[1]*gamma_avg; 
  out[7] = 1.6875*H_L[7]*dv1*gamma_avg+2.6875*H_C[7]*dv1*gamma_avg+1.4072912811497125*H_L[2]*dv1*gamma_avg-1.4072912811497125*H_C[2]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[2]*gamma_avg; 
  out[8] = -(0.180421959121758*H_L[13]*dv1*gamma_avg)+2.70632938682637*H_C[13]*dv1*gamma_avg-0.1875*H_L[8]*dv1*gamma_avg+0.1875*H_C[8]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[5]*gamma_avg; 
  out[9] = -(0.180421959121758*H_L[14]*dv1*gamma_avg)+2.70632938682637*H_C[14]*dv1*gamma_avg-0.1875*H_L[9]*dv1*gamma_avg+0.1875*H_C[9]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[6]*gamma_avg; 
  out[10] = 1.6875*H_L[10]*dv1*gamma_avg+2.6875*H_C[10]*dv1*gamma_avg+1.4072912811497125*H_L[4]*dv1*gamma_avg-1.4072912811497125*H_C[4]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[3]*gamma_avg; 
  out[11] = 1.6875*H_L[11]*dv1*gamma_avg+2.6875*H_C[11]*dv1*gamma_avg+1.4072912811497125*H_L[5]*dv1*gamma_avg-1.4072912811497125*H_C[5]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[4]*gamma_avg; 
  out[12] = -(0.180421959121758*H_L[15]*dv1*gamma_avg)+2.70632938682637*H_C[15]*dv1*gamma_avg-0.1875*H_L[12]*dv1*gamma_avg+0.1875*H_C[12]*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[7]*gamma_avg; 
  out[13] = 1.6875*H_L[13]*dv1*gamma_avg+2.6875*H_C[13]*dv1*gamma_avg+1.4072912811497125*H_L[8]*dv1*gamma_avg-1.4072912811497125*H_C[8]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[5]*gamma_avg; 
  out[14] = 1.6875*H_L[14]*dv1*gamma_avg+2.6875*H_C[14]*dv1*gamma_avg+1.4072912811497125*H_L[9]*dv1*gamma_avg-1.4072912811497125*H_C[9]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[6]*gamma_avg; 
  out[15] = 1.6875*H_L[15]*dv1*gamma_avg+2.6875*H_C[15]*dv1*gamma_avg+1.4072912811497125*H_L[12]*dv1*gamma_avg-1.4072912811497125*H_C[12]*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[7]*gamma_avg; 
  out[24] = -(0.4034357652299393*H_L[3]*dv1*gamma_avg)-1.6944302139657443*H_C[3]*dv1*gamma_avg-0.4192627457812106*H_L[0]*dv1*gamma_avg+0.4192627457812106*H_C[0]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[0]*gamma_avg; 
  out[25] = -(0.4034357652299393*H_L[6]*dv1*gamma_avg)-1.694430213965745*H_C[6]*dv1*gamma_avg-0.4192627457812105*H_L[1]*dv1*gamma_avg+0.4192627457812105*H_C[1]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[1]*gamma_avg; 
  out[26] = -(0.4034357652299393*H_L[7]*dv1*gamma_avg)-1.694430213965745*H_C[7]*dv1*gamma_avg-0.4192627457812105*H_L[2]*dv1*gamma_avg+0.4192627457812105*H_C[2]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[2]*gamma_avg; 
  out[27] = -(0.4034357652299393*H_L[10]*dv1*gamma_avg)-1.694430213965745*H_C[10]*dv1*gamma_avg-0.4192627457812105*H_L[4]*dv1*gamma_avg+0.4192627457812105*H_C[4]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[3]*gamma_avg; 
  out[28] = -(0.4034357652299393*H_L[11]*dv1*gamma_avg)-1.6944302139657443*H_C[11]*dv1*gamma_avg-0.4192627457812106*H_L[5]*dv1*gamma_avg+0.4192627457812106*H_C[5]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[4]*gamma_avg; 
  out[29] = -(0.4034357652299393*H_L[13]*dv1*gamma_avg)-1.6944302139657443*H_C[13]*dv1*gamma_avg-0.4192627457812106*H_L[8]*dv1*gamma_avg+0.4192627457812106*H_C[8]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[5]*gamma_avg; 
  out[30] = -(0.4034357652299393*H_L[14]*dv1*gamma_avg)-1.6944302139657443*H_C[14]*dv1*gamma_avg-0.4192627457812106*H_L[9]*dv1*gamma_avg+0.4192627457812106*H_C[9]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[6]*gamma_avg; 
  out[31] = -(0.4034357652299393*H_L[15]*dv1*gamma_avg)-1.694430213965745*H_C[15]*dv1*gamma_avg-0.4192627457812105*H_L[12]*dv1*gamma_avg+0.4192627457812105*H_C[12]*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[7]*gamma_avg; 
} 
