#include <gkyl_canonical_pb_kernels.h>
#include <gkyl_basis_tensor_5x_p1_upwind_quad_to_modal.h> 
GKYL_CU_DH double canonical_pb_surfy_2x3v_tensor_p1(const double *w, const double *dxv, const double *hamil, 
  const double *alpha_surf_l, const double *alpha_surf_r, 
  const double *sgn_alpha_surf_l, const double *sgn_alpha_surf_r, 
  const int *const_sgn_alpha_l, const int *const_sgn_alpha_r, 
  const double *fl, const double *fc, const double *fr, double* GKYL_RESTRICT out) 
{ 
  // w[NDIM]: cell-center.
  // dxv[NDIM]: cell length.
  // hamil: hamiltonian.
  // alpha_surf_l: Surface expansion of phase space flux on the left.
  // alpha_surf_r: Surface expansion of phase space flux on the right.
  // sgn_alpha_surf_l: sign(alpha_surf_l) at quadrature points.
  // sgn_alpha_surf_r: sign(alpha_surf_r) at quadrature points.
  // const_sgn_alpha_l: Boolean array true if sign(alpha_surf_l) is only one sign, either +1 or -1.
  // const_sgn_alpha_r: Boolean array true if sign(alpha_surf_r) is only one sign, either +1 or -1.
  // fl,fc,fr: distribution function in left, center and right cells.
  // out: output increment in center cell.

  double rdy2 = 2.0/dxv[1];

  const double *alphaL = &alpha_surf_l[16];
  const double *alphaR = &alpha_surf_r[16];
  const double *sgn_alpha_surfL = &sgn_alpha_surf_l[16];
  const double *sgn_alpha_surfR = &sgn_alpha_surf_r[16];
  const int *const_sgn_alphaL = &const_sgn_alpha_l[1];
  const int *const_sgn_alphaR = &const_sgn_alpha_r[1];

  double fUpL[16] = {0.};
  if (const_sgn_alphaL[0] == 1) {  
    if (sgn_alpha_surfL[0] == 1.0) {  
  fUpL[0] = 1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  fUpL[1] = 1.224744871391589*fl[6]+0.7071067811865475*fl[1]; 
  fUpL[2] = 1.224744871391589*fl[8]+0.7071067811865475*fl[3]; 
  fUpL[3] = 1.224744871391589*fl[10]+0.7071067811865475*fl[4]; 
  fUpL[4] = 1.224744871391589*fl[13]+0.7071067811865475*fl[5]; 
  fUpL[5] = 1.224744871391589*fl[16]+0.7071067811865475*fl[7]; 
  fUpL[6] = 1.224744871391589*fl[17]+0.7071067811865475*fl[9]; 
  fUpL[7] = 1.224744871391589*fl[19]+0.7071067811865475*fl[11]; 
  fUpL[8] = 1.224744871391589*fl[20]+0.7071067811865475*fl[12]; 
  fUpL[9] = 1.224744871391589*fl[22]+0.7071067811865475*fl[14]; 
  fUpL[10] = 1.224744871391589*fl[24]+0.7071067811865475*fl[15]; 
  fUpL[11] = 1.224744871391589*fl[26]+0.7071067811865475*fl[18]; 
  fUpL[12] = 1.224744871391589*fl[27]+0.7071067811865475*fl[21]; 
  fUpL[13] = 1.224744871391589*fl[28]+0.7071067811865475*fl[23]; 
  fUpL[14] = 1.224744871391589*fl[30]+0.7071067811865475*fl[25]; 
  fUpL[15] = 1.224744871391589*fl[31]+0.7071067811865475*fl[29]; 
    } else { 
  fUpL[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[2]; 
  fUpL[1] = 0.7071067811865475*fc[1]-1.224744871391589*fc[6]; 
  fUpL[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[8]; 
  fUpL[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[10]; 
  fUpL[4] = 0.7071067811865475*fc[5]-1.224744871391589*fc[13]; 
  fUpL[5] = 0.7071067811865475*fc[7]-1.224744871391589*fc[16]; 
  fUpL[6] = 0.7071067811865475*fc[9]-1.224744871391589*fc[17]; 
  fUpL[7] = 0.7071067811865475*fc[11]-1.224744871391589*fc[19]; 
  fUpL[8] = 0.7071067811865475*fc[12]-1.224744871391589*fc[20]; 
  fUpL[9] = 0.7071067811865475*fc[14]-1.224744871391589*fc[22]; 
  fUpL[10] = 0.7071067811865475*fc[15]-1.224744871391589*fc[24]; 
  fUpL[11] = 0.7071067811865475*fc[18]-1.224744871391589*fc[26]; 
  fUpL[12] = 0.7071067811865475*fc[21]-1.224744871391589*fc[27]; 
  fUpL[13] = 0.7071067811865475*fc[23]-1.224744871391589*fc[28]; 
  fUpL[14] = 0.7071067811865475*fc[25]-1.224744871391589*fc[30]; 
  fUpL[15] = 0.7071067811865475*fc[29]-1.224744871391589*fc[31]; 
    } 
  } else { 
  double f_lr[16] = {0.};
  double f_cl[16] = {0.};
  double sgn_alphaUpL[16] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(sgn_alpha_surfL, sgn_alphaUpL); 

  f_lr[0] = 1.224744871391589*fl[2]+0.7071067811865475*fl[0]; 
  f_lr[1] = 1.224744871391589*fl[6]+0.7071067811865475*fl[1]; 
  f_lr[2] = 1.224744871391589*fl[8]+0.7071067811865475*fl[3]; 
  f_lr[3] = 1.224744871391589*fl[10]+0.7071067811865475*fl[4]; 
  f_lr[4] = 1.224744871391589*fl[13]+0.7071067811865475*fl[5]; 
  f_lr[5] = 1.224744871391589*fl[16]+0.7071067811865475*fl[7]; 
  f_lr[6] = 1.224744871391589*fl[17]+0.7071067811865475*fl[9]; 
  f_lr[7] = 1.224744871391589*fl[19]+0.7071067811865475*fl[11]; 
  f_lr[8] = 1.224744871391589*fl[20]+0.7071067811865475*fl[12]; 
  f_lr[9] = 1.224744871391589*fl[22]+0.7071067811865475*fl[14]; 
  f_lr[10] = 1.224744871391589*fl[24]+0.7071067811865475*fl[15]; 
  f_lr[11] = 1.224744871391589*fl[26]+0.7071067811865475*fl[18]; 
  f_lr[12] = 1.224744871391589*fl[27]+0.7071067811865475*fl[21]; 
  f_lr[13] = 1.224744871391589*fl[28]+0.7071067811865475*fl[23]; 
  f_lr[14] = 1.224744871391589*fl[30]+0.7071067811865475*fl[25]; 
  f_lr[15] = 1.224744871391589*fl[31]+0.7071067811865475*fl[29]; 

  f_cl[0] = 0.7071067811865475*fc[0]-1.224744871391589*fc[2]; 
  f_cl[1] = 0.7071067811865475*fc[1]-1.224744871391589*fc[6]; 
  f_cl[2] = 0.7071067811865475*fc[3]-1.224744871391589*fc[8]; 
  f_cl[3] = 0.7071067811865475*fc[4]-1.224744871391589*fc[10]; 
  f_cl[4] = 0.7071067811865475*fc[5]-1.224744871391589*fc[13]; 
  f_cl[5] = 0.7071067811865475*fc[7]-1.224744871391589*fc[16]; 
  f_cl[6] = 0.7071067811865475*fc[9]-1.224744871391589*fc[17]; 
  f_cl[7] = 0.7071067811865475*fc[11]-1.224744871391589*fc[19]; 
  f_cl[8] = 0.7071067811865475*fc[12]-1.224744871391589*fc[20]; 
  f_cl[9] = 0.7071067811865475*fc[14]-1.224744871391589*fc[22]; 
  f_cl[10] = 0.7071067811865475*fc[15]-1.224744871391589*fc[24]; 
  f_cl[11] = 0.7071067811865475*fc[18]-1.224744871391589*fc[26]; 
  f_cl[12] = 0.7071067811865475*fc[21]-1.224744871391589*fc[27]; 
  f_cl[13] = 0.7071067811865475*fc[23]-1.224744871391589*fc[28]; 
  f_cl[14] = 0.7071067811865475*fc[25]-1.224744871391589*fc[30]; 
  f_cl[15] = 0.7071067811865475*fc[29]-1.224744871391589*fc[31]; 

  fUpL[0] = (0.125*f_lr[15]-0.125*f_cl[15])*sgn_alphaUpL[15]+(0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[14]+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[13]+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[12]+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[11]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[10]+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[9]+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[8]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[7]+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[6]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[5]+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[4]+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[3]+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[2]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[1]+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[0]+0.5*(f_lr[0]+f_cl[0]); 
  fUpL[1] = (0.125*f_lr[14]-0.125*f_cl[14])*sgn_alphaUpL[15]+sgn_alphaUpL[14]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[13]+sgn_alphaUpL[10]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[12]+sgn_alphaUpL[9]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[11]+sgn_alphaUpL[7]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[8]+sgn_alphaUpL[4]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[6]+sgn_alphaUpL[3]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[5]+sgn_alphaUpL[2]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[1]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[1]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[1]; 
  fUpL[2] = (0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[15]+sgn_alphaUpL[13]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[14]+sgn_alphaUpL[10]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[12]+sgn_alphaUpL[8]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[11]+sgn_alphaUpL[6]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[9]+sgn_alphaUpL[4]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[7]+sgn_alphaUpL[3]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[5]+sgn_alphaUpL[1]*(0.125*f_lr[5]-0.125*f_cl[5])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[2]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[2]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[2]; 
  fUpL[3] = (0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[15]+sgn_alphaUpL[12]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[14]+sgn_alphaUpL[9]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[13]+sgn_alphaUpL[8]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[11]+sgn_alphaUpL[5]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[10]+sgn_alphaUpL[4]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[7]+sgn_alphaUpL[2]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[6]+sgn_alphaUpL[1]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[3]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[3]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[3]; 
  fUpL[4] = (0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[15]+sgn_alphaUpL[11]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[14]+sgn_alphaUpL[7]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[13]+sgn_alphaUpL[6]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[12]+sgn_alphaUpL[5]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[10]+sgn_alphaUpL[3]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[9]+sgn_alphaUpL[2]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[8]+sgn_alphaUpL[1]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[4]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[4]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[4]; 
  fUpL[5] = (0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[15]+sgn_alphaUpL[10]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[13]-0.125*f_cl[13])*sgn_alphaUpL[14]+sgn_alphaUpL[13]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[12]+sgn_alphaUpL[4]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[11]+sgn_alphaUpL[3]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[9]+sgn_alphaUpL[8]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[7]+sgn_alphaUpL[6]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[5]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[5]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[5]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[2]+sgn_alphaUpL[1]*(0.125*f_lr[2]-0.125*f_cl[2]); 
  fUpL[6] = (0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[15]+sgn_alphaUpL[9]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[14]+sgn_alphaUpL[12]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[13]+sgn_alphaUpL[4]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[11]+sgn_alphaUpL[2]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[10]+sgn_alphaUpL[8]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[7]+sgn_alphaUpL[5]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[6]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[6]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[6]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[3]+sgn_alphaUpL[1]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[7] = (0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[15]+sgn_alphaUpL[8]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[14]+sgn_alphaUpL[4]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[12]-0.125*f_cl[12])*sgn_alphaUpL[13]+sgn_alphaUpL[12]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[11]+sgn_alphaUpL[1]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[10]+sgn_alphaUpL[9]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[7]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[7]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[7]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[6]+sgn_alphaUpL[5]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[3]+sgn_alphaUpL[2]*(0.125*f_lr[3]-0.125*f_cl[3]); 
  fUpL[8] = (0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[15]+sgn_alphaUpL[7]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[14]+sgn_alphaUpL[11]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[13]+sgn_alphaUpL[3]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[12]+sgn_alphaUpL[2]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[10]+sgn_alphaUpL[6]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[9]+sgn_alphaUpL[5]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[8]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[8]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[8]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[4]+sgn_alphaUpL[1]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[9] = (0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[15]+sgn_alphaUpL[6]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[14]+sgn_alphaUpL[3]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[13]+sgn_alphaUpL[11]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[12]+sgn_alphaUpL[1]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[10]+sgn_alphaUpL[7]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[9]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[9]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[9]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[8]+sgn_alphaUpL[5]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[4]+sgn_alphaUpL[2]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[10] = (0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[15]+sgn_alphaUpL[5]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[14]+sgn_alphaUpL[2]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[13]+sgn_alphaUpL[1]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[11]-0.125*f_cl[11])*sgn_alphaUpL[12]+sgn_alphaUpL[11]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[10]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[10]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[10]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[9]+sgn_alphaUpL[7]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[8]+sgn_alphaUpL[6]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[4]+sgn_alphaUpL[3]*(0.125*f_lr[4]-0.125*f_cl[4]); 
  fUpL[11] = (0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[15]+sgn_alphaUpL[4]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[14]+sgn_alphaUpL[8]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[13]+sgn_alphaUpL[9]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[12]+sgn_alphaUpL[10]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[11]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[11]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[11]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[7]+sgn_alphaUpL[1]*(0.125*f_lr[7]-0.125*f_cl[7])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[6]+sgn_alphaUpL[2]*(0.125*f_lr[6]-0.125*f_cl[6])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[5]+sgn_alphaUpL[3]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[12] = (0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[15]+sgn_alphaUpL[3]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[14]+sgn_alphaUpL[6]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[13]+sgn_alphaUpL[7]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[12]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[12]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[12]+(0.125*f_lr[10]-0.125*f_cl[10])*sgn_alphaUpL[11]+sgn_alphaUpL[10]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[9]+sgn_alphaUpL[1]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[8]+sgn_alphaUpL[2]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[5]+sgn_alphaUpL[4]*(0.125*f_lr[5]-0.125*f_cl[5]); 
  fUpL[13] = (0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[15]+sgn_alphaUpL[2]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[14]+sgn_alphaUpL[5]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[13]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[13]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[13]+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[12]+sgn_alphaUpL[7]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[9]-0.125*f_cl[9])*sgn_alphaUpL[11]+sgn_alphaUpL[9]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[10]+sgn_alphaUpL[1]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[8]+sgn_alphaUpL[3]*(0.125*f_lr[8]-0.125*f_cl[8])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[6]+sgn_alphaUpL[4]*(0.125*f_lr[6]-0.125*f_cl[6]); 
  fUpL[14] = (0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[15]+sgn_alphaUpL[1]*(0.125*f_lr[15]-0.125*f_cl[15])+(0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[14]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[14]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[14]+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[13]+sgn_alphaUpL[5]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[12]+sgn_alphaUpL[6]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[8]-0.125*f_cl[8])*sgn_alphaUpL[11]+sgn_alphaUpL[8]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[10]+sgn_alphaUpL[2]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[9]+sgn_alphaUpL[3]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[7]+sgn_alphaUpL[4]*(0.125*f_lr[7]-0.125*f_cl[7]); 
  fUpL[15] = (0.125*f_lr[0]-0.125*f_cl[0])*sgn_alphaUpL[15]+(0.125*sgn_alphaUpL[0]+0.5)*f_lr[15]+(0.5-0.125*sgn_alphaUpL[0])*f_cl[15]+(0.125*f_lr[1]-0.125*f_cl[1])*sgn_alphaUpL[14]+sgn_alphaUpL[1]*(0.125*f_lr[14]-0.125*f_cl[14])+(0.125*f_lr[2]-0.125*f_cl[2])*sgn_alphaUpL[13]+sgn_alphaUpL[2]*(0.125*f_lr[13]-0.125*f_cl[13])+(0.125*f_lr[3]-0.125*f_cl[3])*sgn_alphaUpL[12]+sgn_alphaUpL[3]*(0.125*f_lr[12]-0.125*f_cl[12])+(0.125*f_lr[4]-0.125*f_cl[4])*sgn_alphaUpL[11]+sgn_alphaUpL[4]*(0.125*f_lr[11]-0.125*f_cl[11])+(0.125*f_lr[5]-0.125*f_cl[5])*sgn_alphaUpL[10]+sgn_alphaUpL[5]*(0.125*f_lr[10]-0.125*f_cl[10])+(0.125*f_lr[6]-0.125*f_cl[6])*sgn_alphaUpL[9]+sgn_alphaUpL[6]*(0.125*f_lr[9]-0.125*f_cl[9])+(0.125*f_lr[7]-0.125*f_cl[7])*sgn_alphaUpL[8]+sgn_alphaUpL[7]*(0.125*f_lr[8]-0.125*f_cl[8]); 

  } 
  double fUpR[16] = {0.};
  if (const_sgn_alphaR[0] == 1) {  
    if (sgn_alpha_surfR[0] == 1.0) {  
  fUpR[0] = 1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  fUpR[1] = 1.224744871391589*fc[6]+0.7071067811865475*fc[1]; 
  fUpR[2] = 1.224744871391589*fc[8]+0.7071067811865475*fc[3]; 
  fUpR[3] = 1.224744871391589*fc[10]+0.7071067811865475*fc[4]; 
  fUpR[4] = 1.224744871391589*fc[13]+0.7071067811865475*fc[5]; 
  fUpR[5] = 1.224744871391589*fc[16]+0.7071067811865475*fc[7]; 
  fUpR[6] = 1.224744871391589*fc[17]+0.7071067811865475*fc[9]; 
  fUpR[7] = 1.224744871391589*fc[19]+0.7071067811865475*fc[11]; 
  fUpR[8] = 1.224744871391589*fc[20]+0.7071067811865475*fc[12]; 
  fUpR[9] = 1.224744871391589*fc[22]+0.7071067811865475*fc[14]; 
  fUpR[10] = 1.224744871391589*fc[24]+0.7071067811865475*fc[15]; 
  fUpR[11] = 1.224744871391589*fc[26]+0.7071067811865475*fc[18]; 
  fUpR[12] = 1.224744871391589*fc[27]+0.7071067811865475*fc[21]; 
  fUpR[13] = 1.224744871391589*fc[28]+0.7071067811865475*fc[23]; 
  fUpR[14] = 1.224744871391589*fc[30]+0.7071067811865475*fc[25]; 
  fUpR[15] = 1.224744871391589*fc[31]+0.7071067811865475*fc[29]; 
    } else { 
  fUpR[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[2]; 
  fUpR[1] = 0.7071067811865475*fr[1]-1.224744871391589*fr[6]; 
  fUpR[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[8]; 
  fUpR[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[10]; 
  fUpR[4] = 0.7071067811865475*fr[5]-1.224744871391589*fr[13]; 
  fUpR[5] = 0.7071067811865475*fr[7]-1.224744871391589*fr[16]; 
  fUpR[6] = 0.7071067811865475*fr[9]-1.224744871391589*fr[17]; 
  fUpR[7] = 0.7071067811865475*fr[11]-1.224744871391589*fr[19]; 
  fUpR[8] = 0.7071067811865475*fr[12]-1.224744871391589*fr[20]; 
  fUpR[9] = 0.7071067811865475*fr[14]-1.224744871391589*fr[22]; 
  fUpR[10] = 0.7071067811865475*fr[15]-1.224744871391589*fr[24]; 
  fUpR[11] = 0.7071067811865475*fr[18]-1.224744871391589*fr[26]; 
  fUpR[12] = 0.7071067811865475*fr[21]-1.224744871391589*fr[27]; 
  fUpR[13] = 0.7071067811865475*fr[23]-1.224744871391589*fr[28]; 
  fUpR[14] = 0.7071067811865475*fr[25]-1.224744871391589*fr[30]; 
  fUpR[15] = 0.7071067811865475*fr[29]-1.224744871391589*fr[31]; 
    } 
  } else { 
  double f_cr[16] = {0.};
  double f_rl[16] = {0.};
  double sgn_alphaUpR[16] = {0.};
  // Project tensor nodal quadrature basis back onto modal basis. 
  tensor_5x_p1_upwind_quad_to_modal(sgn_alpha_surfR, sgn_alphaUpR); 

  f_cr[0] = 1.224744871391589*fc[2]+0.7071067811865475*fc[0]; 
  f_cr[1] = 1.224744871391589*fc[6]+0.7071067811865475*fc[1]; 
  f_cr[2] = 1.224744871391589*fc[8]+0.7071067811865475*fc[3]; 
  f_cr[3] = 1.224744871391589*fc[10]+0.7071067811865475*fc[4]; 
  f_cr[4] = 1.224744871391589*fc[13]+0.7071067811865475*fc[5]; 
  f_cr[5] = 1.224744871391589*fc[16]+0.7071067811865475*fc[7]; 
  f_cr[6] = 1.224744871391589*fc[17]+0.7071067811865475*fc[9]; 
  f_cr[7] = 1.224744871391589*fc[19]+0.7071067811865475*fc[11]; 
  f_cr[8] = 1.224744871391589*fc[20]+0.7071067811865475*fc[12]; 
  f_cr[9] = 1.224744871391589*fc[22]+0.7071067811865475*fc[14]; 
  f_cr[10] = 1.224744871391589*fc[24]+0.7071067811865475*fc[15]; 
  f_cr[11] = 1.224744871391589*fc[26]+0.7071067811865475*fc[18]; 
  f_cr[12] = 1.224744871391589*fc[27]+0.7071067811865475*fc[21]; 
  f_cr[13] = 1.224744871391589*fc[28]+0.7071067811865475*fc[23]; 
  f_cr[14] = 1.224744871391589*fc[30]+0.7071067811865475*fc[25]; 
  f_cr[15] = 1.224744871391589*fc[31]+0.7071067811865475*fc[29]; 

  f_rl[0] = 0.7071067811865475*fr[0]-1.224744871391589*fr[2]; 
  f_rl[1] = 0.7071067811865475*fr[1]-1.224744871391589*fr[6]; 
  f_rl[2] = 0.7071067811865475*fr[3]-1.224744871391589*fr[8]; 
  f_rl[3] = 0.7071067811865475*fr[4]-1.224744871391589*fr[10]; 
  f_rl[4] = 0.7071067811865475*fr[5]-1.224744871391589*fr[13]; 
  f_rl[5] = 0.7071067811865475*fr[7]-1.224744871391589*fr[16]; 
  f_rl[6] = 0.7071067811865475*fr[9]-1.224744871391589*fr[17]; 
  f_rl[7] = 0.7071067811865475*fr[11]-1.224744871391589*fr[19]; 
  f_rl[8] = 0.7071067811865475*fr[12]-1.224744871391589*fr[20]; 
  f_rl[9] = 0.7071067811865475*fr[14]-1.224744871391589*fr[22]; 
  f_rl[10] = 0.7071067811865475*fr[15]-1.224744871391589*fr[24]; 
  f_rl[11] = 0.7071067811865475*fr[18]-1.224744871391589*fr[26]; 
  f_rl[12] = 0.7071067811865475*fr[21]-1.224744871391589*fr[27]; 
  f_rl[13] = 0.7071067811865475*fr[23]-1.224744871391589*fr[28]; 
  f_rl[14] = 0.7071067811865475*fr[25]-1.224744871391589*fr[30]; 
  f_rl[15] = 0.7071067811865475*fr[29]-1.224744871391589*fr[31]; 

  fUpR[0] = (0.125*f_cr[15]-0.125*f_rl[15])*sgn_alphaUpR[15]+(0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[14]+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[13]+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[12]+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[11]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[10]+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[9]+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[8]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[7]+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[6]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[5]+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[4]+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[3]+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[2]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[1]+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[0]+0.5*(f_rl[0]+f_cr[0]); 
  fUpR[1] = (0.125*f_cr[14]-0.125*f_rl[14])*sgn_alphaUpR[15]+sgn_alphaUpR[14]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[13]+sgn_alphaUpR[10]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[12]+sgn_alphaUpR[9]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[11]+sgn_alphaUpR[7]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[8]+sgn_alphaUpR[4]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[6]+sgn_alphaUpR[3]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[5]+sgn_alphaUpR[2]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[1]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[1]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[1]; 
  fUpR[2] = (0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[15]+sgn_alphaUpR[13]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[14]+sgn_alphaUpR[10]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[12]+sgn_alphaUpR[8]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[11]+sgn_alphaUpR[6]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[9]+sgn_alphaUpR[4]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[7]+sgn_alphaUpR[3]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[5]+sgn_alphaUpR[1]*(0.125*f_cr[5]-0.125*f_rl[5])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[2]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[2]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[2]; 
  fUpR[3] = (0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[15]+sgn_alphaUpR[12]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[14]+sgn_alphaUpR[9]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[13]+sgn_alphaUpR[8]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[11]+sgn_alphaUpR[5]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[10]+sgn_alphaUpR[4]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[7]+sgn_alphaUpR[2]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[6]+sgn_alphaUpR[1]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[3]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[3]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[3]; 
  fUpR[4] = (0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[15]+sgn_alphaUpR[11]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[14]+sgn_alphaUpR[7]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[13]+sgn_alphaUpR[6]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[12]+sgn_alphaUpR[5]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[10]+sgn_alphaUpR[3]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[9]+sgn_alphaUpR[2]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[8]+sgn_alphaUpR[1]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[4]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[4]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[4]; 
  fUpR[5] = (0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[15]+sgn_alphaUpR[10]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[13]-0.125*f_rl[13])*sgn_alphaUpR[14]+sgn_alphaUpR[13]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[12]+sgn_alphaUpR[4]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[11]+sgn_alphaUpR[3]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[9]+sgn_alphaUpR[8]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[7]+sgn_alphaUpR[6]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[5]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[5]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[5]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[2]+sgn_alphaUpR[1]*(0.125*f_cr[2]-0.125*f_rl[2]); 
  fUpR[6] = (0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[15]+sgn_alphaUpR[9]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[14]+sgn_alphaUpR[12]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[13]+sgn_alphaUpR[4]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[11]+sgn_alphaUpR[2]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[10]+sgn_alphaUpR[8]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[7]+sgn_alphaUpR[5]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[6]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[6]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[6]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[3]+sgn_alphaUpR[1]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[7] = (0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[15]+sgn_alphaUpR[8]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[14]+sgn_alphaUpR[4]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[12]-0.125*f_rl[12])*sgn_alphaUpR[13]+sgn_alphaUpR[12]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[11]+sgn_alphaUpR[1]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[10]+sgn_alphaUpR[9]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[7]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[7]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[7]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[6]+sgn_alphaUpR[5]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[3]+sgn_alphaUpR[2]*(0.125*f_cr[3]-0.125*f_rl[3]); 
  fUpR[8] = (0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[15]+sgn_alphaUpR[7]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[14]+sgn_alphaUpR[11]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[13]+sgn_alphaUpR[3]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[12]+sgn_alphaUpR[2]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[10]+sgn_alphaUpR[6]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[9]+sgn_alphaUpR[5]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[8]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[8]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[8]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[4]+sgn_alphaUpR[1]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[9] = (0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[15]+sgn_alphaUpR[6]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[14]+sgn_alphaUpR[3]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[13]+sgn_alphaUpR[11]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[12]+sgn_alphaUpR[1]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[10]+sgn_alphaUpR[7]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[9]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[9]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[9]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[8]+sgn_alphaUpR[5]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[4]+sgn_alphaUpR[2]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[10] = (0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[15]+sgn_alphaUpR[5]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[14]+sgn_alphaUpR[2]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[13]+sgn_alphaUpR[1]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[11]-0.125*f_rl[11])*sgn_alphaUpR[12]+sgn_alphaUpR[11]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[10]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[10]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[10]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[9]+sgn_alphaUpR[7]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[8]+sgn_alphaUpR[6]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[4]+sgn_alphaUpR[3]*(0.125*f_cr[4]-0.125*f_rl[4]); 
  fUpR[11] = (0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[15]+sgn_alphaUpR[4]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[14]+sgn_alphaUpR[8]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[13]+sgn_alphaUpR[9]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[12]+sgn_alphaUpR[10]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[11]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[11]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[11]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[7]+sgn_alphaUpR[1]*(0.125*f_cr[7]-0.125*f_rl[7])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[6]+sgn_alphaUpR[2]*(0.125*f_cr[6]-0.125*f_rl[6])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[5]+sgn_alphaUpR[3]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[12] = (0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[15]+sgn_alphaUpR[3]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[14]+sgn_alphaUpR[6]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[13]+sgn_alphaUpR[7]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[12]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[12]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[12]+(0.125*f_cr[10]-0.125*f_rl[10])*sgn_alphaUpR[11]+sgn_alphaUpR[10]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[9]+sgn_alphaUpR[1]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[8]+sgn_alphaUpR[2]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[5]+sgn_alphaUpR[4]*(0.125*f_cr[5]-0.125*f_rl[5]); 
  fUpR[13] = (0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[15]+sgn_alphaUpR[2]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[14]+sgn_alphaUpR[5]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[13]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[13]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[13]+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[12]+sgn_alphaUpR[7]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[9]-0.125*f_rl[9])*sgn_alphaUpR[11]+sgn_alphaUpR[9]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[10]+sgn_alphaUpR[1]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[8]+sgn_alphaUpR[3]*(0.125*f_cr[8]-0.125*f_rl[8])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[6]+sgn_alphaUpR[4]*(0.125*f_cr[6]-0.125*f_rl[6]); 
  fUpR[14] = (0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[15]+sgn_alphaUpR[1]*(0.125*f_cr[15]-0.125*f_rl[15])+(0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[14]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[14]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[14]+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[13]+sgn_alphaUpR[5]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[12]+sgn_alphaUpR[6]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[8]-0.125*f_rl[8])*sgn_alphaUpR[11]+sgn_alphaUpR[8]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[10]+sgn_alphaUpR[2]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[9]+sgn_alphaUpR[3]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[7]+sgn_alphaUpR[4]*(0.125*f_cr[7]-0.125*f_rl[7]); 
  fUpR[15] = (0.125*f_cr[0]-0.125*f_rl[0])*sgn_alphaUpR[15]+(0.5-0.125*sgn_alphaUpR[0])*f_rl[15]+(0.125*sgn_alphaUpR[0]+0.5)*f_cr[15]+(0.125*f_cr[1]-0.125*f_rl[1])*sgn_alphaUpR[14]+sgn_alphaUpR[1]*(0.125*f_cr[14]-0.125*f_rl[14])+(0.125*f_cr[2]-0.125*f_rl[2])*sgn_alphaUpR[13]+sgn_alphaUpR[2]*(0.125*f_cr[13]-0.125*f_rl[13])+(0.125*f_cr[3]-0.125*f_rl[3])*sgn_alphaUpR[12]+sgn_alphaUpR[3]*(0.125*f_cr[12]-0.125*f_rl[12])+(0.125*f_cr[4]-0.125*f_rl[4])*sgn_alphaUpR[11]+sgn_alphaUpR[4]*(0.125*f_cr[11]-0.125*f_rl[11])+(0.125*f_cr[5]-0.125*f_rl[5])*sgn_alphaUpR[10]+sgn_alphaUpR[5]*(0.125*f_cr[10]-0.125*f_rl[10])+(0.125*f_cr[6]-0.125*f_rl[6])*sgn_alphaUpR[9]+sgn_alphaUpR[6]*(0.125*f_cr[9]-0.125*f_rl[9])+(0.125*f_cr[7]-0.125*f_rl[7])*sgn_alphaUpR[8]+sgn_alphaUpR[7]*(0.125*f_cr[8]-0.125*f_rl[8]); 

  } 
  double GhatL[16] = {0.};
  double GhatR[16] = {0.};
  GhatL[0] = 0.25*(alphaL[12]*fUpL[12]+alphaL[9]*fUpL[9]+alphaL[8]*fUpL[8]+alphaL[5]*fUpL[5]+alphaL[4]*fUpL[4]+alphaL[2]*fUpL[2]+alphaL[1]*fUpL[1]+alphaL[0]*fUpL[0]); 
  GhatL[1] = 0.25*(alphaL[9]*fUpL[12]+fUpL[9]*alphaL[12]+alphaL[4]*fUpL[8]+fUpL[4]*alphaL[8]+alphaL[2]*fUpL[5]+fUpL[2]*alphaL[5]+alphaL[0]*fUpL[1]+fUpL[0]*alphaL[1]); 
  GhatL[2] = 0.25*(alphaL[8]*fUpL[12]+fUpL[8]*alphaL[12]+alphaL[4]*fUpL[9]+fUpL[4]*alphaL[9]+alphaL[1]*fUpL[5]+fUpL[1]*alphaL[5]+alphaL[0]*fUpL[2]+fUpL[0]*alphaL[2]); 
  GhatL[3] = 0.25*(alphaL[12]*fUpL[15]+alphaL[9]*fUpL[14]+alphaL[8]*fUpL[13]+alphaL[5]*fUpL[11]+alphaL[4]*fUpL[10]+alphaL[2]*fUpL[7]+alphaL[1]*fUpL[6]+alphaL[0]*fUpL[3]); 
  GhatL[4] = 0.25*(alphaL[5]*fUpL[12]+fUpL[5]*alphaL[12]+alphaL[2]*fUpL[9]+fUpL[2]*alphaL[9]+alphaL[1]*fUpL[8]+fUpL[1]*alphaL[8]+alphaL[0]*fUpL[4]+fUpL[0]*alphaL[4]); 
  GhatL[5] = 0.25*(alphaL[4]*fUpL[12]+fUpL[4]*alphaL[12]+alphaL[8]*fUpL[9]+fUpL[8]*alphaL[9]+alphaL[0]*fUpL[5]+fUpL[0]*alphaL[5]+alphaL[1]*fUpL[2]+fUpL[1]*alphaL[2]); 
  GhatL[6] = 0.25*(alphaL[9]*fUpL[15]+alphaL[12]*fUpL[14]+alphaL[4]*fUpL[13]+alphaL[2]*fUpL[11]+alphaL[8]*fUpL[10]+alphaL[5]*fUpL[7]+alphaL[0]*fUpL[6]+alphaL[1]*fUpL[3]); 
  GhatL[7] = 0.25*(alphaL[8]*fUpL[15]+alphaL[4]*fUpL[14]+alphaL[12]*fUpL[13]+alphaL[1]*fUpL[11]+alphaL[9]*fUpL[10]+alphaL[0]*fUpL[7]+alphaL[5]*fUpL[6]+alphaL[2]*fUpL[3]); 
  GhatL[8] = 0.25*(alphaL[2]*fUpL[12]+fUpL[2]*alphaL[12]+alphaL[5]*fUpL[9]+fUpL[5]*alphaL[9]+alphaL[0]*fUpL[8]+fUpL[0]*alphaL[8]+alphaL[1]*fUpL[4]+fUpL[1]*alphaL[4]); 
  GhatL[9] = 0.25*(alphaL[1]*fUpL[12]+fUpL[1]*alphaL[12]+alphaL[0]*fUpL[9]+fUpL[0]*alphaL[9]+alphaL[5]*fUpL[8]+fUpL[5]*alphaL[8]+alphaL[2]*fUpL[4]+fUpL[2]*alphaL[4]); 
  GhatL[10] = 0.25*(alphaL[5]*fUpL[15]+alphaL[2]*fUpL[14]+alphaL[1]*fUpL[13]+fUpL[11]*alphaL[12]+alphaL[0]*fUpL[10]+fUpL[7]*alphaL[9]+fUpL[6]*alphaL[8]+fUpL[3]*alphaL[4]); 
  GhatL[11] = 0.25*(alphaL[4]*fUpL[15]+alphaL[8]*fUpL[14]+alphaL[9]*fUpL[13]+fUpL[10]*alphaL[12]+alphaL[0]*fUpL[11]+alphaL[1]*fUpL[7]+alphaL[2]*fUpL[6]+fUpL[3]*alphaL[5]); 
  GhatL[12] = 0.25*(alphaL[0]*fUpL[12]+fUpL[0]*alphaL[12]+alphaL[1]*fUpL[9]+fUpL[1]*alphaL[9]+alphaL[2]*fUpL[8]+fUpL[2]*alphaL[8]+alphaL[4]*fUpL[5]+fUpL[4]*alphaL[5]); 
  GhatL[13] = 0.25*(alphaL[2]*fUpL[15]+alphaL[5]*fUpL[14]+alphaL[0]*fUpL[13]+fUpL[7]*alphaL[12]+alphaL[9]*fUpL[11]+alphaL[1]*fUpL[10]+fUpL[3]*alphaL[8]+alphaL[4]*fUpL[6]); 
  GhatL[14] = 0.25*(alphaL[1]*fUpL[15]+alphaL[0]*fUpL[14]+alphaL[5]*fUpL[13]+fUpL[6]*alphaL[12]+alphaL[8]*fUpL[11]+alphaL[2]*fUpL[10]+fUpL[3]*alphaL[9]+alphaL[4]*fUpL[7]); 
  GhatL[15] = 0.25*(alphaL[0]*fUpL[15]+alphaL[1]*fUpL[14]+alphaL[2]*fUpL[13]+fUpL[3]*alphaL[12]+alphaL[4]*fUpL[11]+alphaL[5]*fUpL[10]+fUpL[6]*alphaL[9]+fUpL[7]*alphaL[8]); 

  GhatR[0] = 0.25*(alphaR[12]*fUpR[12]+alphaR[9]*fUpR[9]+alphaR[8]*fUpR[8]+alphaR[5]*fUpR[5]+alphaR[4]*fUpR[4]+alphaR[2]*fUpR[2]+alphaR[1]*fUpR[1]+alphaR[0]*fUpR[0]); 
  GhatR[1] = 0.25*(alphaR[9]*fUpR[12]+fUpR[9]*alphaR[12]+alphaR[4]*fUpR[8]+fUpR[4]*alphaR[8]+alphaR[2]*fUpR[5]+fUpR[2]*alphaR[5]+alphaR[0]*fUpR[1]+fUpR[0]*alphaR[1]); 
  GhatR[2] = 0.25*(alphaR[8]*fUpR[12]+fUpR[8]*alphaR[12]+alphaR[4]*fUpR[9]+fUpR[4]*alphaR[9]+alphaR[1]*fUpR[5]+fUpR[1]*alphaR[5]+alphaR[0]*fUpR[2]+fUpR[0]*alphaR[2]); 
  GhatR[3] = 0.25*(alphaR[12]*fUpR[15]+alphaR[9]*fUpR[14]+alphaR[8]*fUpR[13]+alphaR[5]*fUpR[11]+alphaR[4]*fUpR[10]+alphaR[2]*fUpR[7]+alphaR[1]*fUpR[6]+alphaR[0]*fUpR[3]); 
  GhatR[4] = 0.25*(alphaR[5]*fUpR[12]+fUpR[5]*alphaR[12]+alphaR[2]*fUpR[9]+fUpR[2]*alphaR[9]+alphaR[1]*fUpR[8]+fUpR[1]*alphaR[8]+alphaR[0]*fUpR[4]+fUpR[0]*alphaR[4]); 
  GhatR[5] = 0.25*(alphaR[4]*fUpR[12]+fUpR[4]*alphaR[12]+alphaR[8]*fUpR[9]+fUpR[8]*alphaR[9]+alphaR[0]*fUpR[5]+fUpR[0]*alphaR[5]+alphaR[1]*fUpR[2]+fUpR[1]*alphaR[2]); 
  GhatR[6] = 0.25*(alphaR[9]*fUpR[15]+alphaR[12]*fUpR[14]+alphaR[4]*fUpR[13]+alphaR[2]*fUpR[11]+alphaR[8]*fUpR[10]+alphaR[5]*fUpR[7]+alphaR[0]*fUpR[6]+alphaR[1]*fUpR[3]); 
  GhatR[7] = 0.25*(alphaR[8]*fUpR[15]+alphaR[4]*fUpR[14]+alphaR[12]*fUpR[13]+alphaR[1]*fUpR[11]+alphaR[9]*fUpR[10]+alphaR[0]*fUpR[7]+alphaR[5]*fUpR[6]+alphaR[2]*fUpR[3]); 
  GhatR[8] = 0.25*(alphaR[2]*fUpR[12]+fUpR[2]*alphaR[12]+alphaR[5]*fUpR[9]+fUpR[5]*alphaR[9]+alphaR[0]*fUpR[8]+fUpR[0]*alphaR[8]+alphaR[1]*fUpR[4]+fUpR[1]*alphaR[4]); 
  GhatR[9] = 0.25*(alphaR[1]*fUpR[12]+fUpR[1]*alphaR[12]+alphaR[0]*fUpR[9]+fUpR[0]*alphaR[9]+alphaR[5]*fUpR[8]+fUpR[5]*alphaR[8]+alphaR[2]*fUpR[4]+fUpR[2]*alphaR[4]); 
  GhatR[10] = 0.25*(alphaR[5]*fUpR[15]+alphaR[2]*fUpR[14]+alphaR[1]*fUpR[13]+fUpR[11]*alphaR[12]+alphaR[0]*fUpR[10]+fUpR[7]*alphaR[9]+fUpR[6]*alphaR[8]+fUpR[3]*alphaR[4]); 
  GhatR[11] = 0.25*(alphaR[4]*fUpR[15]+alphaR[8]*fUpR[14]+alphaR[9]*fUpR[13]+fUpR[10]*alphaR[12]+alphaR[0]*fUpR[11]+alphaR[1]*fUpR[7]+alphaR[2]*fUpR[6]+fUpR[3]*alphaR[5]); 
  GhatR[12] = 0.25*(alphaR[0]*fUpR[12]+fUpR[0]*alphaR[12]+alphaR[1]*fUpR[9]+fUpR[1]*alphaR[9]+alphaR[2]*fUpR[8]+fUpR[2]*alphaR[8]+alphaR[4]*fUpR[5]+fUpR[4]*alphaR[5]); 
  GhatR[13] = 0.25*(alphaR[2]*fUpR[15]+alphaR[5]*fUpR[14]+alphaR[0]*fUpR[13]+fUpR[7]*alphaR[12]+alphaR[9]*fUpR[11]+alphaR[1]*fUpR[10]+fUpR[3]*alphaR[8]+alphaR[4]*fUpR[6]); 
  GhatR[14] = 0.25*(alphaR[1]*fUpR[15]+alphaR[0]*fUpR[14]+alphaR[5]*fUpR[13]+fUpR[6]*alphaR[12]+alphaR[8]*fUpR[11]+alphaR[2]*fUpR[10]+fUpR[3]*alphaR[9]+alphaR[4]*fUpR[7]); 
  GhatR[15] = 0.25*(alphaR[0]*fUpR[15]+alphaR[1]*fUpR[14]+alphaR[2]*fUpR[13]+fUpR[3]*alphaR[12]+alphaR[4]*fUpR[11]+alphaR[5]*fUpR[10]+fUpR[6]*alphaR[9]+fUpR[7]*alphaR[8]); 

  out[0] += (0.7071067811865475*GhatL[0]-0.7071067811865475*GhatR[0])*rdy2; 
  out[1] += (0.7071067811865475*GhatL[1]-0.7071067811865475*GhatR[1])*rdy2; 
  out[2] += ((-1.224744871391589*GhatR[0])-1.224744871391589*GhatL[0])*rdy2; 
  out[3] += (0.7071067811865475*GhatL[2]-0.7071067811865475*GhatR[2])*rdy2; 
  out[4] += (0.7071067811865475*GhatL[3]-0.7071067811865475*GhatR[3])*rdy2; 
  out[5] += (0.7071067811865475*GhatL[4]-0.7071067811865475*GhatR[4])*rdy2; 
  out[6] += ((-1.224744871391589*GhatR[1])-1.224744871391589*GhatL[1])*rdy2; 
  out[7] += (0.7071067811865475*GhatL[5]-0.7071067811865475*GhatR[5])*rdy2; 
  out[8] += ((-1.224744871391589*GhatR[2])-1.224744871391589*GhatL[2])*rdy2; 
  out[9] += (0.7071067811865475*GhatL[6]-0.7071067811865475*GhatR[6])*rdy2; 
  out[10] += ((-1.224744871391589*GhatR[3])-1.224744871391589*GhatL[3])*rdy2; 
  out[11] += (0.7071067811865475*GhatL[7]-0.7071067811865475*GhatR[7])*rdy2; 
  out[12] += (0.7071067811865475*GhatL[8]-0.7071067811865475*GhatR[8])*rdy2; 
  out[13] += ((-1.224744871391589*GhatR[4])-1.224744871391589*GhatL[4])*rdy2; 
  out[14] += (0.7071067811865475*GhatL[9]-0.7071067811865475*GhatR[9])*rdy2; 
  out[15] += (0.7071067811865475*GhatL[10]-0.7071067811865475*GhatR[10])*rdy2; 
  out[16] += ((-1.224744871391589*GhatR[5])-1.224744871391589*GhatL[5])*rdy2; 
  out[17] += ((-1.224744871391589*GhatR[6])-1.224744871391589*GhatL[6])*rdy2; 
  out[18] += (0.7071067811865475*GhatL[11]-0.7071067811865475*GhatR[11])*rdy2; 
  out[19] += ((-1.224744871391589*GhatR[7])-1.224744871391589*GhatL[7])*rdy2; 
  out[20] += ((-1.224744871391589*GhatR[8])-1.224744871391589*GhatL[8])*rdy2; 
  out[21] += (0.7071067811865475*GhatL[12]-0.7071067811865475*GhatR[12])*rdy2; 
  out[22] += ((-1.224744871391589*GhatR[9])-1.224744871391589*GhatL[9])*rdy2; 
  out[23] += (0.7071067811865475*GhatL[13]-0.7071067811865475*GhatR[13])*rdy2; 
  out[24] += ((-1.224744871391589*GhatR[10])-1.224744871391589*GhatL[10])*rdy2; 
  out[25] += (0.7071067811865475*GhatL[14]-0.7071067811865475*GhatR[14])*rdy2; 
  out[26] += ((-1.224744871391589*GhatR[11])-1.224744871391589*GhatL[11])*rdy2; 
  out[27] += ((-1.224744871391589*GhatR[12])-1.224744871391589*GhatL[12])*rdy2; 
  out[28] += ((-1.224744871391589*GhatR[13])-1.224744871391589*GhatL[13])*rdy2; 
  out[29] += (0.7071067811865475*GhatL[15]-0.7071067811865475*GhatR[15])*rdy2; 
  out[30] += ((-1.224744871391589*GhatR[14])-1.224744871391589*GhatL[14])*rdy2; 
  out[31] += ((-1.224744871391589*GhatR[15])-1.224744871391589*GhatL[15])*rdy2; 

  double cflFreq = fmax(fabs(alphaL[0]), fabs(alphaR[0])); 
  return 0.375*rdy2*cflFreq; 

} 
