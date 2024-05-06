#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH int fpo_drag_coeff_1x3v_vy_ser_p1_upvy(const double *dxv, const double *gamma, 
    const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff, 
    double *drag_coeff_surf, double *sgn_drag_coeff_surf
    ) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 
  // drag_coeff_surf: Surface projection of drag coefficient at lower boundary.
  // sgn_drag_coeff_surf: Sign(drag_coeff_surf) evaluated at quadrature points along lower surface.
  // returns const_sgn_drag_coeff: 1 if sign(drag_coeff_surf) is constant along lower boundary, 0 otherwise. 

  // Use cell-average value for gamma. 
  double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[2]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  
  const double *dHdv_surf_C = &fpo_dhdv_surf[8]; 
  
  double *out = &drag_coeff[40]; 
  double *out_surf = &drag_coeff_surf[8]; 
  double *sgn_alpha_surf = &sgn_drag_coeff_surf[8]; 
  
  out[0] = (-(0.180421959121758*H_L[3])+2.70632938682637*H_C[3]-0.1875*H_L[0]+0.1875*H_C[0])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[0]*gamma_avg; 
  out[1] = (-(0.180421959121758*H_L[6])+2.70632938682637*H_C[6]-0.1875*H_L[1]+0.1875*H_C[1])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[1]*gamma_avg; 
  out[2] = (-(0.180421959121758*H_L[7])+2.70632938682637*H_C[7]-0.1875*H_L[2]+0.1875*H_C[2])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[2]*gamma_avg; 
  out[3] = (1.6875*H_L[3]+2.6875*H_C[3]+1.4072912811497125*H_L[0]-1.4072912811497125*H_C[0])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[0]*gamma_avg; 
  out[4] = (-(0.180421959121758*H_L[10])+2.70632938682637*H_C[10]-0.1875*H_L[4]+0.1875*H_C[4])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[3]*gamma_avg; 
  out[5] = (-(0.180421959121758*H_L[11])+2.70632938682637*H_C[11]-0.1875*H_L[5]+0.1875*H_C[5])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[4]*gamma_avg; 
  out[6] = (1.6875*H_L[6]+2.6875*H_C[6]+1.4072912811497125*H_L[1]-1.4072912811497125*H_C[1])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[1]*gamma_avg; 
  out[7] = (1.6875*H_L[7]+2.6875*H_C[7]+1.4072912811497125*H_L[2]-1.4072912811497125*H_C[2])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[2]*gamma_avg; 
  out[8] = (-(0.180421959121758*H_L[13])+2.70632938682637*H_C[13]-0.1875*H_L[8]+0.1875*H_C[8])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[5]*gamma_avg; 
  out[9] = (-(0.180421959121758*H_L[14])+2.70632938682637*H_C[14]-0.1875*H_L[9]+0.1875*H_C[9])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[6]*gamma_avg; 
  out[10] = (1.6875*H_L[10]+2.6875*H_C[10]+1.4072912811497125*H_L[4]-1.4072912811497125*H_C[4])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[3]*gamma_avg; 
  out[11] = (1.6875*H_L[11]+2.6875*H_C[11]+1.4072912811497125*H_L[5]-1.4072912811497125*H_C[5])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[4]*gamma_avg; 
  out[12] = (-(0.180421959121758*H_L[15])+2.70632938682637*H_C[15]-0.1875*H_L[12]+0.1875*H_C[12])*dv1*gamma_avg+0.23570226039551595*dHdv_surf_C[7]*gamma_avg; 
  out[13] = (1.6875*H_L[13]+2.6875*H_C[13]+1.4072912811497125*H_L[8]-1.4072912811497125*H_C[8])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[5]*gamma_avg; 
  out[14] = (1.6875*H_L[14]+2.6875*H_C[14]+1.4072912811497125*H_L[9]-1.4072912811497125*H_C[9])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[6]*gamma_avg; 
  out[15] = (1.6875*H_L[15]+2.6875*H_C[15]+1.4072912811497125*H_L[12]-1.4072912811497125*H_C[12])*dv1*gamma_avg+0.4082482904638632*dHdv_surf_C[7]*gamma_avg; 
  out[24] = (-(0.4034357652299393*H_L[3])-1.6944302139657443*H_C[3]-0.4192627457812106*H_L[0]+0.4192627457812106*H_C[0])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[0]*gamma_avg; 
  out[25] = (-(0.4034357652299393*H_L[6])-1.694430213965745*H_C[6]-0.4192627457812105*H_L[1]+0.4192627457812105*H_C[1])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[1]*gamma_avg; 
  out[26] = (-(0.4034357652299393*H_L[7])-1.694430213965745*H_C[7]-0.4192627457812105*H_L[2]+0.4192627457812105*H_C[2])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[2]*gamma_avg; 
  out[27] = (-(0.4034357652299393*H_L[10])-1.694430213965745*H_C[10]-0.4192627457812105*H_L[4]+0.4192627457812105*H_C[4])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[3]*gamma_avg; 
  out[28] = (-(0.4034357652299393*H_L[11])-1.6944302139657443*H_C[11]-0.4192627457812106*H_L[5]+0.4192627457812106*H_C[5])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[4]*gamma_avg; 
  out[29] = (-(0.4034357652299393*H_L[13])-1.6944302139657443*H_C[13]-0.4192627457812106*H_L[8]+0.4192627457812106*H_C[8])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[5]*gamma_avg; 
  out[30] = (-(0.4034357652299393*H_L[14])-1.6944302139657443*H_C[14]-0.4192627457812106*H_L[9]+0.4192627457812106*H_C[9])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[6]*gamma_avg; 
  out[31] = (-(0.4034357652299393*H_L[15])-1.694430213965745*H_C[15]-0.4192627457812105*H_L[12]+0.4192627457812105*H_C[12])*dv1*gamma_avg+0.5270462766947301*dHdv_surf_C[7]*gamma_avg; 

  out_surf[0] = -(0.1767766952966368*(8.660254037844386*H_L[3]+8.660254037844386*H_C[3]+9.0*H_L[0]-9.0*H_C[0])*dv1*gamma_avg); 
  out_surf[1] = -(0.1767766952966368*(8.660254037844386*H_L[6]+8.660254037844386*H_C[6]+9.0*H_L[1]-9.0*H_C[1])*dv1*gamma_avg); 
  out_surf[2] = -(0.1767766952966368*(8.660254037844386*H_L[7]+8.660254037844386*H_C[7]+9.0*H_L[2]-9.0*H_C[2])*dv1*gamma_avg); 
  out_surf[3] = -(0.1767766952966368*(8.660254037844386*H_L[10]+8.660254037844386*H_C[10]+9.0*H_L[4]-9.0*H_C[4])*dv1*gamma_avg); 
  out_surf[4] = -(0.1767766952966368*(8.660254037844386*H_L[11]+8.660254037844386*H_C[11]+9.0*H_L[5]-9.0*H_C[5])*dv1*gamma_avg); 
  out_surf[5] = -(0.1767766952966368*(8.660254037844386*H_L[13]+8.660254037844386*H_C[13]+9.0*H_L[8]-9.0*H_C[8])*dv1*gamma_avg); 
  out_surf[6] = -(0.1767766952966368*(8.660254037844386*H_L[14]+8.660254037844386*H_C[14]+9.0*H_L[9]-9.0*H_C[9])*dv1*gamma_avg); 
  out_surf[7] = -(0.1767766952966368*(8.660254037844386*H_L[15]+8.660254037844386*H_C[15]+9.0*H_L[12]-9.0*H_C[12])*dv1*gamma_avg); 

  int const_sgn_alpha_surf = 1;  
  
  if (-(0.3535533905932734*out_surf[7])+0.3535533905932734*(out_surf[6]+out_surf[5]+out_surf[4])-0.3535533905932734*(out_surf[3]+out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[0] = 1.0; 
  else  
    sgn_alpha_surf[0] = -1.0; 
  
  if (0.3535533905932734*out_surf[7]-0.3535533905932734*(out_surf[6]+out_surf[5])+0.3535533905932734*(out_surf[4]+out_surf[3])-0.3535533905932734*(out_surf[2]+out_surf[1])+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[1] = 1.0; 
  else  
    sgn_alpha_surf[1] = -1.0; 
  
  if (sgn_alpha_surf[1] == sgn_alpha_surf[0]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*out_surf[7]-0.3535533905932734*out_surf[6]+0.3535533905932734*out_surf[5]-0.3535533905932734*(out_surf[4]+out_surf[3])+0.3535533905932734*out_surf[2]-0.3535533905932734*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[2] = 1.0; 
  else  
    sgn_alpha_surf[2] = -1.0; 
  
  if (sgn_alpha_surf[2] == sgn_alpha_surf[1]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*out_surf[7])+0.3535533905932734*out_surf[6]-0.3535533905932734*(out_surf[5]+out_surf[4])+0.3535533905932734*(out_surf[3]+out_surf[2])-0.3535533905932734*out_surf[1]+0.3535533905932734*out_surf[0] > 0.) 
    sgn_alpha_surf[3] = 1.0; 
  else  
    sgn_alpha_surf[3] = -1.0; 
  
  if (sgn_alpha_surf[3] == sgn_alpha_surf[2]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(out_surf[7]+out_surf[6])-0.3535533905932734*(out_surf[5]+out_surf[4]+out_surf[3]+out_surf[2])+0.3535533905932734*(out_surf[1]+out_surf[0]) > 0.) 
    sgn_alpha_surf[4] = 1.0; 
  else  
    sgn_alpha_surf[4] = -1.0; 
  
  if (sgn_alpha_surf[4] == sgn_alpha_surf[3]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*(out_surf[7]+out_surf[6]))+0.3535533905932734*out_surf[5]-0.3535533905932734*out_surf[4]+0.3535533905932734*out_surf[3]-0.3535533905932734*out_surf[2]+0.3535533905932734*(out_surf[1]+out_surf[0]) > 0.) 
    sgn_alpha_surf[5] = 1.0; 
  else  
    sgn_alpha_surf[5] = -1.0; 
  
  if (sgn_alpha_surf[5] == sgn_alpha_surf[4]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (-(0.3535533905932734*(out_surf[7]+out_surf[6]+out_surf[5]))+0.3535533905932734*out_surf[4]-0.3535533905932734*out_surf[3]+0.3535533905932734*(out_surf[2]+out_surf[1]+out_surf[0]) > 0.) 
    sgn_alpha_surf[6] = 1.0; 
  else  
    sgn_alpha_surf[6] = -1.0; 
  
  if (sgn_alpha_surf[6] == sgn_alpha_surf[5]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  if (0.3535533905932734*(out_surf[7]+out_surf[6]+out_surf[5]+out_surf[4]+out_surf[3]+out_surf[2]+out_surf[1]+out_surf[0]) > 0.) 
    sgn_alpha_surf[7] = 1.0; 
  else  
    sgn_alpha_surf[7] = -1.0; 
  
  if (sgn_alpha_surf[7] == sgn_alpha_surf[6]) 
    const_sgn_alpha_surf = const_sgn_alpha_surf ? 1 : 0; 
  else  
    const_sgn_alpha_surf = 0; 
  
  return const_sgn_alpha_surf; 
} 

