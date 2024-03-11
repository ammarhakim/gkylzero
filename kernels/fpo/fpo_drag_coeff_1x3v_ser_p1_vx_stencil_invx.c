#include <gkyl_fpo_vlasov_kernels.h> 

GKYL_CU_DH void fpo_drag_coeff_1x3v_vx_ser_p1_invx(const double *dxv, const double *gamma, const double* fpo_h_stencil[3], const double* fpo_dhdv_surf, double *drag_coeff) {
  // dxv[NDIM]: Cell spacing in each direction. 
  // gamma: Scalar factor gamma. 
  // fpo_h_stencil[3]: 3 cell stencil of Rosenbluth potential H. 
  // fpo_dhdv_surf: Surface projection of dH/dv in center cell. 
  // drag_coeff: Output array for drag coefficient. 

  // Use cell-average value for gamma. 
 double gamma_avg = gamma[0]/sqrt(pow(2, 1)); 
  double dv1 = 2.0/dxv[1]; 

  const double* H_L = fpo_h_stencil[0]; 
  const double* H_C = fpo_h_stencil[1]; 
  const double* H_R = fpo_h_stencil[2]; 
  
  const double *fpo_dHdv_surf_C_vx = &fpo_dhdv_surf[0]; 
  const double *fpo_dHdv_surf_C_vy = &fpo_dhdv_surf[8]; 
  const double *fpo_dHdv_surf_C_vz = &fpo_dhdv_surf[16]; 
  
  const double* dHdv_surf_C = fpo_dHdv_surf_C_vx; 
  
  double *drag_coeff_vx = &drag_coeff[0]; 
  double *drag_coeff_vy = &drag_coeff[40]; 
  double *drag_coeff_vz = &drag_coeff[80]; 
  
  double *out = drag_coeff_vx; 
  
  out[0] = -(0.5773502691896257*H_R[2]*dv1*gamma_avg)-0.5773502691896257*H_L[2]*dv1*gamma_avg+1.1547005383792515*H_C[2]*dv1*gamma_avg+0.5*H_R[0]*dv1*gamma_avg-0.5*H_L[0]*dv1*gamma_avg; 
  out[1] = -(0.5773502691896257*H_R[5]*dv1*gamma_avg)-0.5773502691896257*H_L[5]*dv1*gamma_avg+1.1547005383792515*H_C[5]*dv1*gamma_avg+0.5*H_R[1]*dv1*gamma_avg-0.5*H_L[1]*dv1*gamma_avg; 
  out[2] = -(1.0*H_R[2]*dv1*gamma_avg)+1.0*H_L[2]*dv1*gamma_avg+0.8660254037844386*H_R[0]*dv1*gamma_avg+0.8660254037844386*H_L[0]*dv1*gamma_avg-1.7320508075688772*H_C[0]*dv1*gamma_avg; 
  out[3] = -(0.5773502691896257*H_R[7]*dv1*gamma_avg)-0.5773502691896257*H_L[7]*dv1*gamma_avg+1.1547005383792515*H_C[7]*dv1*gamma_avg+0.5*H_R[3]*dv1*gamma_avg-0.5*H_L[3]*dv1*gamma_avg; 
  out[4] = -(0.5773502691896257*H_R[9]*dv1*gamma_avg)-0.5773502691896257*H_L[9]*dv1*gamma_avg+1.1547005383792515*H_C[9]*dv1*gamma_avg+0.5*H_R[4]*dv1*gamma_avg-0.5*H_L[4]*dv1*gamma_avg; 
  out[5] = -(1.0*H_R[5]*dv1*gamma_avg)+1.0*H_L[5]*dv1*gamma_avg+0.8660254037844386*H_R[1]*dv1*gamma_avg+0.8660254037844386*H_L[1]*dv1*gamma_avg-1.7320508075688772*H_C[1]*dv1*gamma_avg; 
  out[6] = -(0.5773502691896257*H_R[11]*dv1*gamma_avg)-0.5773502691896257*H_L[11]*dv1*gamma_avg+1.1547005383792515*H_C[11]*dv1*gamma_avg+0.5*H_R[6]*dv1*gamma_avg-0.5*H_L[6]*dv1*gamma_avg; 
  out[7] = -(1.0*H_R[7]*dv1*gamma_avg)+1.0*H_L[7]*dv1*gamma_avg+0.8660254037844386*H_R[3]*dv1*gamma_avg+0.8660254037844386*H_L[3]*dv1*gamma_avg-1.7320508075688772*H_C[3]*dv1*gamma_avg; 
  out[8] = -(0.5773502691896257*H_R[12]*dv1*gamma_avg)-0.5773502691896257*H_L[12]*dv1*gamma_avg+1.1547005383792515*H_C[12]*dv1*gamma_avg+0.5*H_R[8]*dv1*gamma_avg-0.5*H_L[8]*dv1*gamma_avg; 
  out[9] = -(1.0*H_R[9]*dv1*gamma_avg)+1.0*H_L[9]*dv1*gamma_avg+0.8660254037844386*H_R[4]*dv1*gamma_avg+0.8660254037844386*H_L[4]*dv1*gamma_avg-1.7320508075688772*H_C[4]*dv1*gamma_avg; 
  out[10] = -(0.5773502691896257*H_R[14]*dv1*gamma_avg)-0.5773502691896257*H_L[14]*dv1*gamma_avg+1.1547005383792515*H_C[14]*dv1*gamma_avg+0.5*H_R[10]*dv1*gamma_avg-0.5*H_L[10]*dv1*gamma_avg; 
  out[11] = -(1.0*H_R[11]*dv1*gamma_avg)+1.0*H_L[11]*dv1*gamma_avg+0.8660254037844386*H_R[6]*dv1*gamma_avg+0.8660254037844386*H_L[6]*dv1*gamma_avg-1.7320508075688772*H_C[6]*dv1*gamma_avg; 
  out[12] = -(1.0*H_R[12]*dv1*gamma_avg)+1.0*H_L[12]*dv1*gamma_avg+0.8660254037844386*H_R[8]*dv1*gamma_avg+0.8660254037844386*H_L[8]*dv1*gamma_avg-1.7320508075688772*H_C[8]*dv1*gamma_avg; 
  out[13] = -(0.5773502691896257*H_R[15]*dv1*gamma_avg)-0.5773502691896257*H_L[15]*dv1*gamma_avg+1.1547005383792515*H_C[15]*dv1*gamma_avg+0.5*H_R[13]*dv1*gamma_avg-0.5*H_L[13]*dv1*gamma_avg; 
  out[14] = -(1.0*H_R[14]*dv1*gamma_avg)+1.0*H_L[14]*dv1*gamma_avg+0.8660254037844386*H_R[10]*dv1*gamma_avg+0.8660254037844386*H_L[10]*dv1*gamma_avg-1.7320508075688772*H_C[10]*dv1*gamma_avg; 
  out[15] = -(1.0*H_R[15]*dv1*gamma_avg)+1.0*H_L[15]*dv1*gamma_avg+0.8660254037844386*H_R[13]*dv1*gamma_avg+0.8660254037844386*H_L[13]*dv1*gamma_avg-1.7320508075688772*H_C[13]*dv1*gamma_avg; 
  out[16] = -(1.2909944487358056*H_R[2]*dv1*gamma_avg)-1.2909944487358056*H_L[2]*dv1*gamma_avg-5.163977794943222*H_C[2]*dv1*gamma_avg+1.118033988749895*H_R[0]*dv1*gamma_avg-1.118033988749895*H_L[0]*dv1*gamma_avg; 
  out[17] = -(1.2909944487358056*H_R[5]*dv1*gamma_avg)-1.2909944487358056*H_L[5]*dv1*gamma_avg-5.163977794943222*H_C[5]*dv1*gamma_avg+1.118033988749895*H_R[1]*dv1*gamma_avg-1.118033988749895*H_L[1]*dv1*gamma_avg; 
  out[18] = -(1.2909944487358056*H_R[7]*dv1*gamma_avg)-1.2909944487358056*H_L[7]*dv1*gamma_avg-5.163977794943222*H_C[7]*dv1*gamma_avg+1.118033988749895*H_R[3]*dv1*gamma_avg-1.118033988749895*H_L[3]*dv1*gamma_avg; 
  out[19] = -(1.2909944487358056*H_R[9]*dv1*gamma_avg)-1.2909944487358056*H_L[9]*dv1*gamma_avg-5.163977794943222*H_C[9]*dv1*gamma_avg+1.118033988749895*H_R[4]*dv1*gamma_avg-1.118033988749895*H_L[4]*dv1*gamma_avg; 
  out[20] = -(1.2909944487358056*H_R[11]*dv1*gamma_avg)-1.2909944487358056*H_L[11]*dv1*gamma_avg-5.163977794943222*H_C[11]*dv1*gamma_avg+1.118033988749895*H_R[6]*dv1*gamma_avg-1.118033988749895*H_L[6]*dv1*gamma_avg; 
  out[21] = -(1.2909944487358056*H_R[12]*dv1*gamma_avg)-1.2909944487358056*H_L[12]*dv1*gamma_avg-5.163977794943222*H_C[12]*dv1*gamma_avg+1.118033988749895*H_R[8]*dv1*gamma_avg-1.118033988749895*H_L[8]*dv1*gamma_avg; 
  out[22] = -(1.2909944487358056*H_R[14]*dv1*gamma_avg)-1.2909944487358056*H_L[14]*dv1*gamma_avg-5.163977794943222*H_C[14]*dv1*gamma_avg+1.118033988749895*H_R[10]*dv1*gamma_avg-1.118033988749895*H_L[10]*dv1*gamma_avg; 
  out[23] = -(1.2909944487358056*H_R[15]*dv1*gamma_avg)-1.2909944487358056*H_L[15]*dv1*gamma_avg-5.163977794943222*H_C[15]*dv1*gamma_avg+1.118033988749895*H_R[13]*dv1*gamma_avg-1.118033988749895*H_L[13]*dv1*gamma_avg; 
} 

