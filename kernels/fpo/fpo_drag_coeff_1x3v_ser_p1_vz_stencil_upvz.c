#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vz_ser_p1_upvz(const double *dxv, const double gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  double dv1 = 2.0/dxv[3]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[8]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[16]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vz; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[40]; 
  double *drag_coeff_vz = &drag_coeff[80]; 
  
  double *out = drag_coeff_vz; 
  
  out[0] = -(0.180421959121758*H_L[4]*dv1*gamma)+2.70632938682637*H_C[4]*dv1*gamma-0.1875*H_L[0]*dv1*gamma+0.1875*H_C[0]*dv1*gamma+0.23570226039551595*dHdv_surf_C[0]*gamma; 
  out[1] = -(0.180421959121758*H_L[8]*dv1*gamma)+2.70632938682637*H_C[8]*dv1*gamma-0.1875*H_L[1]*dv1*gamma+0.1875*H_C[1]*dv1*gamma+0.23570226039551595*dHdv_surf_C[1]*gamma; 
  out[2] = -(0.180421959121758*H_L[9]*dv1*gamma)+2.70632938682637*H_C[9]*dv1*gamma-0.1875*H_L[2]*dv1*gamma+0.1875*H_C[2]*dv1*gamma+0.23570226039551595*dHdv_surf_C[2]*gamma; 
  out[3] = -(0.180421959121758*H_L[10]*dv1*gamma)+2.70632938682637*H_C[10]*dv1*gamma-0.1875*H_L[3]*dv1*gamma+0.1875*H_C[3]*dv1*gamma+0.23570226039551595*dHdv_surf_C[3]*gamma; 
  out[4] = 1.6875*H_L[4]*dv1*gamma+2.6875*H_C[4]*dv1*gamma+1.4072912811497125*H_L[0]*dv1*gamma-1.4072912811497125*H_C[0]*dv1*gamma+0.4082482904638632*dHdv_surf_C[0]*gamma; 
  out[5] = -(0.180421959121758*H_L[12]*dv1*gamma)+2.70632938682637*H_C[12]*dv1*gamma-0.1875*H_L[5]*dv1*gamma+0.1875*H_C[5]*dv1*gamma+0.23570226039551595*dHdv_surf_C[4]*gamma; 
  out[6] = -(0.180421959121758*H_L[13]*dv1*gamma)+2.70632938682637*H_C[13]*dv1*gamma-0.1875*H_L[6]*dv1*gamma+0.1875*H_C[6]*dv1*gamma+0.23570226039551595*dHdv_surf_C[5]*gamma; 
  out[7] = -(0.180421959121758*H_L[14]*dv1*gamma)+2.70632938682637*H_C[14]*dv1*gamma-0.1875*H_L[7]*dv1*gamma+0.1875*H_C[7]*dv1*gamma+0.23570226039551595*dHdv_surf_C[6]*gamma; 
  out[8] = 1.6875*H_L[8]*dv1*gamma+2.6875*H_C[8]*dv1*gamma+1.4072912811497125*H_L[1]*dv1*gamma-1.4072912811497125*H_C[1]*dv1*gamma+0.4082482904638632*dHdv_surf_C[1]*gamma; 
  out[9] = 1.6875*H_L[9]*dv1*gamma+2.6875*H_C[9]*dv1*gamma+1.4072912811497125*H_L[2]*dv1*gamma-1.4072912811497125*H_C[2]*dv1*gamma+0.4082482904638632*dHdv_surf_C[2]*gamma; 
  out[10] = 1.6875*H_L[10]*dv1*gamma+2.6875*H_C[10]*dv1*gamma+1.4072912811497125*H_L[3]*dv1*gamma-1.4072912811497125*H_C[3]*dv1*gamma+0.4082482904638632*dHdv_surf_C[3]*gamma; 
  out[11] = -(0.180421959121758*H_L[15]*dv1*gamma)+2.70632938682637*H_C[15]*dv1*gamma-0.1875*H_L[11]*dv1*gamma+0.1875*H_C[11]*dv1*gamma+0.23570226039551595*dHdv_surf_C[7]*gamma; 
  out[12] = 1.6875*H_L[12]*dv1*gamma+2.6875*H_C[12]*dv1*gamma+1.4072912811497125*H_L[5]*dv1*gamma-1.4072912811497125*H_C[5]*dv1*gamma+0.4082482904638632*dHdv_surf_C[4]*gamma; 
  out[13] = 1.6875*H_L[13]*dv1*gamma+2.6875*H_C[13]*dv1*gamma+1.4072912811497125*H_L[6]*dv1*gamma-1.4072912811497125*H_C[6]*dv1*gamma+0.4082482904638632*dHdv_surf_C[5]*gamma; 
  out[14] = 1.6875*H_L[14]*dv1*gamma+2.6875*H_C[14]*dv1*gamma+1.4072912811497125*H_L[7]*dv1*gamma-1.4072912811497125*H_C[7]*dv1*gamma+0.4082482904638632*dHdv_surf_C[6]*gamma; 
  out[15] = 1.6875*H_L[15]*dv1*gamma+2.6875*H_C[15]*dv1*gamma+1.4072912811497125*H_L[11]*dv1*gamma-1.4072912811497125*H_C[11]*dv1*gamma+0.4082482904638632*dHdv_surf_C[7]*gamma; 
  out[32] = -(0.4034357652299393*H_L[4]*dv1*gamma)-1.6944302139657443*H_C[4]*dv1*gamma-0.4192627457812106*H_L[0]*dv1*gamma+0.4192627457812106*H_C[0]*dv1*gamma+0.5270462766947301*dHdv_surf_C[0]*gamma; 
  out[33] = -(0.4034357652299393*H_L[8]*dv1*gamma)-1.694430213965745*H_C[8]*dv1*gamma-0.4192627457812105*H_L[1]*dv1*gamma+0.4192627457812105*H_C[1]*dv1*gamma+0.5270462766947301*dHdv_surf_C[1]*gamma; 
  out[34] = -(0.4034357652299393*H_L[9]*dv1*gamma)-1.694430213965745*H_C[9]*dv1*gamma-0.4192627457812105*H_L[2]*dv1*gamma+0.4192627457812105*H_C[2]*dv1*gamma+0.5270462766947301*dHdv_surf_C[2]*gamma; 
  out[35] = -(0.4034357652299393*H_L[10]*dv1*gamma)-1.694430213965745*H_C[10]*dv1*gamma-0.4192627457812105*H_L[3]*dv1*gamma+0.4192627457812105*H_C[3]*dv1*gamma+0.5270462766947301*dHdv_surf_C[3]*gamma; 
  out[36] = -(0.4034357652299393*H_L[12]*dv1*gamma)-1.6944302139657443*H_C[12]*dv1*gamma-0.4192627457812106*H_L[5]*dv1*gamma+0.4192627457812106*H_C[5]*dv1*gamma+0.5270462766947301*dHdv_surf_C[4]*gamma; 
  out[37] = -(0.4034357652299393*H_L[13]*dv1*gamma)-1.6944302139657443*H_C[13]*dv1*gamma-0.4192627457812106*H_L[6]*dv1*gamma+0.4192627457812106*H_C[6]*dv1*gamma+0.5270462766947301*dHdv_surf_C[5]*gamma; 
  out[38] = -(0.4034357652299393*H_L[14]*dv1*gamma)-1.6944302139657443*H_C[14]*dv1*gamma-0.4192627457812106*H_L[7]*dv1*gamma+0.4192627457812106*H_C[7]*dv1*gamma+0.5270462766947301*dHdv_surf_C[6]*gamma; 
  out[39] = -(0.4034357652299393*H_L[15]*dv1*gamma)-1.694430213965745*H_C[15]*dv1*gamma-0.4192627457812105*H_L[11]*dv1*gamma+0.4192627457812105*H_C[11]*dv1*gamma+0.5270462766947301*dHdv_surf_C[7]*gamma; 
} 
