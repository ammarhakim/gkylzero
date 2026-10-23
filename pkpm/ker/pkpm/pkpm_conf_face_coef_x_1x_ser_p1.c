#include <gkyl_pkpm_kinetic_kernels.h> 
#include <gkyl_pkpm_harm_coeffs.h> 
GKYL_CU_DH void pkpm_conf_face_coef_x_1x_ser_p1(const double *b_l, const double *b_r, const double *u_l, const double *u_r,
  const double *xi_l, const double *xi_r, const double *theta_l, const double *theta_r,
  double* GKYL_RESTRICT fc) 
{ 
  // b_l/b_r:   Magnetic field unit vector [bx, by, bz] in the left/current cell. 
  // u_l/u_r:   Flow velocity [ux, uy, uz] in the left/current cell. 
  // xi_l/xi_r: Thermal speed [xi, xi^2, p] in the left/current cell (the first expansion is used). 
  // theta_l/r: T_perp/T in the left/current cell. 
  // fc:        Output face coefficients of the flux on the lower face of the current cell in 
  //            direction x at the 1 transverse nodes: fc[q*1 + i] for q = PKPM_FC_*. 
  const double *b0_l = &b_l[0], *b0_r = &b_r[0]; 
  const double *b1_l = &b_l[2], *b1_r = &b_r[2]; 
  const double *b2_l = &b_l[4], *b2_r = &b_r[4]; 
  const double *ud_l = &u_l[0], *ud_r = &u_r[0]; 
  const double *xd_l = &xi_l[0], *xd_r = &xi_r[0]; 
 
  // Traces of b, u_x, xi and theta at the transverse nodes: the left cell at its upper 
  // edge, the current cell at its lower edge. 
  double bq0_l[1], bq0_r[1]; 
  bq0_l[0] = 1.224744871391589*b0_l[1]+0.7071067811865475*b0_l[0]; 
  bq0_r[0] = 0.7071067811865475*b0_r[0]-1.224744871391589*b0_r[1]; 
  double bq1_l[1], bq1_r[1]; 
  bq1_l[0] = 1.224744871391589*b1_l[1]+0.7071067811865475*b1_l[0]; 
  bq1_r[0] = 0.7071067811865475*b1_r[0]-1.224744871391589*b1_r[1]; 
  double bq2_l[1], bq2_r[1]; 
  bq2_l[0] = 1.224744871391589*b2_l[1]+0.7071067811865475*b2_l[0]; 
  bq2_r[0] = 0.7071067811865475*b2_r[0]-1.224744871391589*b2_r[1]; 
  double uq_l[1], uq_r[1], xq_l[1], xq_r[1], thq_l[1], thq_r[1]; 
  uq_l[0] = 1.224744871391589*ud_l[1]+0.7071067811865475*ud_l[0]; 
  uq_r[0] = 0.7071067811865475*ud_r[0]-1.224744871391589*ud_r[1]; 
  xq_l[0] = 1.224744871391589*xd_l[1]+0.7071067811865475*xd_l[0]; 
  xq_r[0] = 0.7071067811865475*xd_r[0]-1.224744871391589*xd_r[1]; 
  thq_l[0] = 1.224744871391589*theta_l[1]+0.7071067811865475*theta_l[0]; 
  thq_r[0] = 0.7071067811865475*theta_r[0]-1.224744871391589*theta_r[1]; 
 
  for (int i = 0; i < 1; ++i) { 
    const double bn_l[3] = {bq0_l[i], bq1_l[i], bq2_l[i]}, bn_r[3] = {bq0_r[i], bq1_r[i], bq2_r[i]}; 
    pkpm_harm_face_coef(0, uq_l[i], uq_r[i], xq_l[i], xq_r[i], bn_l, bn_r, thq_l[i], thq_r[i], 1, &fc[i]); 
  } 
} 
