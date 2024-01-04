#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_lovx(const double *dxv, const double gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  double dv1 = 2.0/dxv[1]; 

  const double* H_C = fpo_h_stencil[0]; 
  const double* H_R = fpo_h_stencil[1]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[8]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[16]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vx; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[40]; 
  double *drag_coeff_vz = &drag_coeff[80]; 
  
  double *out = drag_coeff_vx; 
  
  out[0] = -(0.180421959121758*H_R[2]*dv1*gamma)+2.70632938682637*H_C[2]*dv1*gamma+0.1875*H_R[0]*dv1*gamma-0.1875*H_C[0]*dv1*gamma+0.23570226039551595*dHdv_surf_C[0]*gamma; 
  out[1] = -(0.180421959121758*H_R[5]*dv1*gamma)+2.70632938682637*H_C[5]*dv1*gamma+0.1875*H_R[1]*dv1*gamma-0.1875*H_C[1]*dv1*gamma+0.23570226039551595*dHdv_surf_C[1]*gamma; 
  out[2] = -(1.6875*H_R[2]*dv1*gamma)-2.6875*H_C[2]*dv1*gamma+1.4072912811497125*H_R[0]*dv1*gamma-1.4072912811497125*H_C[0]*dv1*gamma-0.4082482904638632*dHdv_surf_C[0]*gamma; 
  out[3] = -(0.180421959121758*H_R[7]*dv1*gamma)+2.70632938682637*H_C[7]*dv1*gamma+0.1875*H_R[3]*dv1*gamma-0.1875*H_C[3]*dv1*gamma+0.23570226039551595*dHdv_surf_C[2]*gamma; 
  out[4] = -(0.180421959121758*H_R[9]*dv1*gamma)+2.70632938682637*H_C[9]*dv1*gamma+0.1875*H_R[4]*dv1*gamma-0.1875*H_C[4]*dv1*gamma+0.23570226039551595*dHdv_surf_C[3]*gamma; 
  out[5] = -(1.6875*H_R[5]*dv1*gamma)-2.6875*H_C[5]*dv1*gamma+1.4072912811497125*H_R[1]*dv1*gamma-1.4072912811497125*H_C[1]*dv1*gamma-0.4082482904638632*dHdv_surf_C[1]*gamma; 
  out[6] = -(0.180421959121758*H_R[11]*dv1*gamma)+2.70632938682637*H_C[11]*dv1*gamma+0.1875*H_R[6]*dv1*gamma-0.1875*H_C[6]*dv1*gamma+0.23570226039551595*dHdv_surf_C[4]*gamma; 
  out[7] = -(1.6875*H_R[7]*dv1*gamma)-2.6875*H_C[7]*dv1*gamma+1.4072912811497125*H_R[3]*dv1*gamma-1.4072912811497125*H_C[3]*dv1*gamma-0.4082482904638632*dHdv_surf_C[2]*gamma; 
  out[8] = -(0.180421959121758*H_R[12]*dv1*gamma)+2.70632938682637*H_C[12]*dv1*gamma+0.1875*H_R[8]*dv1*gamma-0.1875*H_C[8]*dv1*gamma+0.23570226039551595*dHdv_surf_C[5]*gamma; 
  out[9] = -(1.6875*H_R[9]*dv1*gamma)-2.6875*H_C[9]*dv1*gamma+1.4072912811497125*H_R[4]*dv1*gamma-1.4072912811497125*H_C[4]*dv1*gamma-0.4082482904638632*dHdv_surf_C[3]*gamma; 
  out[10] = -(0.180421959121758*H_R[14]*dv1*gamma)+2.70632938682637*H_C[14]*dv1*gamma+0.1875*H_R[10]*dv1*gamma-0.1875*H_C[10]*dv1*gamma+0.23570226039551595*dHdv_surf_C[6]*gamma; 
  out[11] = -(1.6875*H_R[11]*dv1*gamma)-2.6875*H_C[11]*dv1*gamma+1.4072912811497125*H_R[6]*dv1*gamma-1.4072912811497125*H_C[6]*dv1*gamma-0.4082482904638632*dHdv_surf_C[4]*gamma; 
  out[12] = -(1.6875*H_R[12]*dv1*gamma)-2.6875*H_C[12]*dv1*gamma+1.4072912811497125*H_R[8]*dv1*gamma-1.4072912811497125*H_C[8]*dv1*gamma-0.4082482904638632*dHdv_surf_C[5]*gamma; 
  out[13] = -(0.180421959121758*H_R[15]*dv1*gamma)+2.70632938682637*H_C[15]*dv1*gamma+0.1875*H_R[13]*dv1*gamma-0.1875*H_C[13]*dv1*gamma+0.23570226039551595*dHdv_surf_C[7]*gamma; 
  out[14] = -(1.6875*H_R[14]*dv1*gamma)-2.6875*H_C[14]*dv1*gamma+1.4072912811497125*H_R[10]*dv1*gamma-1.4072912811497125*H_C[10]*dv1*gamma-0.4082482904638632*dHdv_surf_C[6]*gamma; 
  out[15] = -(1.6875*H_R[15]*dv1*gamma)-2.6875*H_C[15]*dv1*gamma+1.4072912811497125*H_R[13]*dv1*gamma-1.4072912811497125*H_C[13]*dv1*gamma-0.4082482904638632*dHdv_surf_C[7]*gamma; 
  out[16] = -(0.4034357652299393*H_R[2]*dv1*gamma)-1.6944302139657443*H_C[2]*dv1*gamma+0.4192627457812106*H_R[0]*dv1*gamma-0.4192627457812106*H_C[0]*dv1*gamma+0.5270462766947301*dHdv_surf_C[0]*gamma; 
  out[17] = -(0.4034357652299393*H_R[5]*dv1*gamma)-1.694430213965745*H_C[5]*dv1*gamma+0.4192627457812105*H_R[1]*dv1*gamma-0.4192627457812105*H_C[1]*dv1*gamma+0.5270462766947301*dHdv_surf_C[1]*gamma; 
  out[18] = -(0.4034357652299393*H_R[7]*dv1*gamma)-1.694430213965745*H_C[7]*dv1*gamma+0.4192627457812105*H_R[3]*dv1*gamma-0.4192627457812105*H_C[3]*dv1*gamma+0.5270462766947301*dHdv_surf_C[2]*gamma; 
  out[19] = -(0.4034357652299393*H_R[9]*dv1*gamma)-1.694430213965745*H_C[9]*dv1*gamma+0.4192627457812105*H_R[4]*dv1*gamma-0.4192627457812105*H_C[4]*dv1*gamma+0.5270462766947301*dHdv_surf_C[3]*gamma; 
  out[20] = -(0.4034357652299393*H_R[11]*dv1*gamma)-1.6944302139657443*H_C[11]*dv1*gamma+0.4192627457812106*H_R[6]*dv1*gamma-0.4192627457812106*H_C[6]*dv1*gamma+0.5270462766947301*dHdv_surf_C[4]*gamma; 
  out[21] = -(0.4034357652299393*H_R[12]*dv1*gamma)-1.6944302139657443*H_C[12]*dv1*gamma+0.4192627457812106*H_R[8]*dv1*gamma-0.4192627457812106*H_C[8]*dv1*gamma+0.5270462766947301*dHdv_surf_C[5]*gamma; 
  out[22] = -(0.4034357652299393*H_R[14]*dv1*gamma)-1.6944302139657443*H_C[14]*dv1*gamma+0.4192627457812106*H_R[10]*dv1*gamma-0.4192627457812106*H_C[10]*dv1*gamma+0.5270462766947301*dHdv_surf_C[6]*gamma; 
  out[23] = -(0.4034357652299393*H_R[15]*dv1*gamma)-1.694430213965745*H_C[15]*dv1*gamma+0.4192627457812105*H_R[13]*dv1*gamma-0.4192627457812105*H_C[13]*dv1*gamma+0.5270462766947301*dHdv_surf_C[7]*gamma; 
} 

