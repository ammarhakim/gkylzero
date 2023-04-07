#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_2x_p1_exp_sq.h> 
#include <gkyl_basis_ser_2x_p1_inv.h> 
#include <gkyl_basis_ser_2x_p1_sqrt_with_sign.h> 
GKYL_CU_DH void em_pkpm_kappa_inv_b_2x_ser_p1(const double *bvar, const double *ExB, double* kappa_inv_b) 
{ 
  // bvar:        Input magnetic field unit vector and unit tensor. 
  // ExB:         Input E x B velocity. 
  // kappa_inv_b: b_i/kappa = B_i/|B| * sqrt(1 - |E x B|^2/(c^2 |B|^4)). 
 
  const double *bx = &bvar[0]; 
  const double *by = &bvar[4]; 
  const double *bz = &bvar[8]; 
  const double *bxbx = &bvar[12]; 
  const double *byby = &bvar[24]; 
  const double *bzbz = &bvar[32]; 
 
  const double *ExB_x = &ExB[0]; 
  const double *ExB_y = &ExB[4]; 
  const double *ExB_z = &ExB[8]; 
 
  double *kappa_inv_bx = &kappa_inv_b[0]; 
  double *kappa_inv_by = &kappa_inv_b[4]; 
  double *kappa_inv_bz = &kappa_inv_b[8]; 
 
  // Calculate ((E x B)/|B|^2)^2. 
  double ExB_x_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(ExB_x, ExB_x_sq); 
 
  double ExB_y_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(ExB_y, ExB_y_sq); 
 
  double ExB_z_sq[4] = {0.0}; 
  ser_2x_p1_exp_sq(ExB_z, ExB_z_sq); 
 
  double kappa_inv_bx_sq[4] = {0.0}; 
  double kappa_inv_by_sq[4] = {0.0}; 
  double kappa_inv_bz_sq[4] = {0.0}; 
 
  kappa_inv_bx_sq[0] = (-0.5*ExB_z_sq[3]*bxbx[3])-0.5*ExB_y_sq[3]*bxbx[3]-0.5*ExB_x_sq[3]*bxbx[3]-0.5*ExB_z_sq[2]*bxbx[2]-0.5*ExB_y_sq[2]*bxbx[2]-0.5*ExB_x_sq[2]*bxbx[2]-0.5*ExB_z_sq[1]*bxbx[1]-0.5*ExB_y_sq[1]*bxbx[1]-0.5*ExB_x_sq[1]*bxbx[1]-0.5*ExB_z_sq[0]*bxbx[0]-0.5*ExB_y_sq[0]*bxbx[0]-0.5*ExB_x_sq[0]*bxbx[0]+bxbx[0]; 
  kappa_inv_bx_sq[1] = (-0.5*ExB_z_sq[2]*bxbx[3])-0.5*ExB_y_sq[2]*bxbx[3]-0.5*ExB_x_sq[2]*bxbx[3]-0.5*bxbx[2]*ExB_z_sq[3]-0.5*bxbx[2]*ExB_y_sq[3]-0.5*bxbx[2]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*bxbx[1]-0.5*ExB_y_sq[0]*bxbx[1]-0.5*ExB_x_sq[0]*bxbx[1]+bxbx[1]-0.5*bxbx[0]*ExB_z_sq[1]-0.5*bxbx[0]*ExB_y_sq[1]-0.5*bxbx[0]*ExB_x_sq[1]; 
  kappa_inv_bx_sq[2] = (-0.5*ExB_z_sq[1]*bxbx[3])-0.5*ExB_y_sq[1]*bxbx[3]-0.5*ExB_x_sq[1]*bxbx[3]-0.5*bxbx[1]*ExB_z_sq[3]-0.5*bxbx[1]*ExB_y_sq[3]-0.5*bxbx[1]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*bxbx[2]-0.5*ExB_y_sq[0]*bxbx[2]-0.5*ExB_x_sq[0]*bxbx[2]+bxbx[2]-0.5*bxbx[0]*ExB_z_sq[2]-0.5*bxbx[0]*ExB_y_sq[2]-0.5*bxbx[0]*ExB_x_sq[2]; 
  kappa_inv_bx_sq[3] = (-0.5*ExB_z_sq[0]*bxbx[3])-0.5*ExB_y_sq[0]*bxbx[3]-0.5*ExB_x_sq[0]*bxbx[3]+bxbx[3]-0.5*bxbx[0]*ExB_z_sq[3]-0.5*bxbx[0]*ExB_y_sq[3]-0.5*bxbx[0]*ExB_x_sq[3]-0.5*ExB_z_sq[1]*bxbx[2]-0.5*ExB_y_sq[1]*bxbx[2]-0.5*ExB_x_sq[1]*bxbx[2]-0.5*bxbx[1]*ExB_z_sq[2]-0.5*bxbx[1]*ExB_y_sq[2]-0.5*bxbx[1]*ExB_x_sq[2]; 

  kappa_inv_by_sq[0] = (-0.5*ExB_z_sq[3]*byby[3])-0.5*ExB_y_sq[3]*byby[3]-0.5*ExB_x_sq[3]*byby[3]-0.5*ExB_z_sq[2]*byby[2]-0.5*ExB_y_sq[2]*byby[2]-0.5*ExB_x_sq[2]*byby[2]-0.5*ExB_z_sq[1]*byby[1]-0.5*ExB_y_sq[1]*byby[1]-0.5*ExB_x_sq[1]*byby[1]-0.5*ExB_z_sq[0]*byby[0]-0.5*ExB_y_sq[0]*byby[0]-0.5*ExB_x_sq[0]*byby[0]+byby[0]; 
  kappa_inv_by_sq[1] = (-0.5*ExB_z_sq[2]*byby[3])-0.5*ExB_y_sq[2]*byby[3]-0.5*ExB_x_sq[2]*byby[3]-0.5*byby[2]*ExB_z_sq[3]-0.5*byby[2]*ExB_y_sq[3]-0.5*byby[2]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*byby[1]-0.5*ExB_y_sq[0]*byby[1]-0.5*ExB_x_sq[0]*byby[1]+byby[1]-0.5*byby[0]*ExB_z_sq[1]-0.5*byby[0]*ExB_y_sq[1]-0.5*byby[0]*ExB_x_sq[1]; 
  kappa_inv_by_sq[2] = (-0.5*ExB_z_sq[1]*byby[3])-0.5*ExB_y_sq[1]*byby[3]-0.5*ExB_x_sq[1]*byby[3]-0.5*byby[1]*ExB_z_sq[3]-0.5*byby[1]*ExB_y_sq[3]-0.5*byby[1]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*byby[2]-0.5*ExB_y_sq[0]*byby[2]-0.5*ExB_x_sq[0]*byby[2]+byby[2]-0.5*byby[0]*ExB_z_sq[2]-0.5*byby[0]*ExB_y_sq[2]-0.5*byby[0]*ExB_x_sq[2]; 
  kappa_inv_by_sq[3] = (-0.5*ExB_z_sq[0]*byby[3])-0.5*ExB_y_sq[0]*byby[3]-0.5*ExB_x_sq[0]*byby[3]+byby[3]-0.5*byby[0]*ExB_z_sq[3]-0.5*byby[0]*ExB_y_sq[3]-0.5*byby[0]*ExB_x_sq[3]-0.5*ExB_z_sq[1]*byby[2]-0.5*ExB_y_sq[1]*byby[2]-0.5*ExB_x_sq[1]*byby[2]-0.5*byby[1]*ExB_z_sq[2]-0.5*byby[1]*ExB_y_sq[2]-0.5*byby[1]*ExB_x_sq[2]; 

  kappa_inv_bz_sq[0] = (-0.5*ExB_z_sq[3]*bzbz[3])-0.5*ExB_y_sq[3]*bzbz[3]-0.5*ExB_x_sq[3]*bzbz[3]-0.5*ExB_z_sq[2]*bzbz[2]-0.5*ExB_y_sq[2]*bzbz[2]-0.5*ExB_x_sq[2]*bzbz[2]-0.5*ExB_z_sq[1]*bzbz[1]-0.5*ExB_y_sq[1]*bzbz[1]-0.5*ExB_x_sq[1]*bzbz[1]-0.5*ExB_z_sq[0]*bzbz[0]-0.5*ExB_y_sq[0]*bzbz[0]-0.5*ExB_x_sq[0]*bzbz[0]+bzbz[0]; 
  kappa_inv_bz_sq[1] = (-0.5*ExB_z_sq[2]*bzbz[3])-0.5*ExB_y_sq[2]*bzbz[3]-0.5*ExB_x_sq[2]*bzbz[3]-0.5*bzbz[2]*ExB_z_sq[3]-0.5*bzbz[2]*ExB_y_sq[3]-0.5*bzbz[2]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*bzbz[1]-0.5*ExB_y_sq[0]*bzbz[1]-0.5*ExB_x_sq[0]*bzbz[1]+bzbz[1]-0.5*bzbz[0]*ExB_z_sq[1]-0.5*bzbz[0]*ExB_y_sq[1]-0.5*bzbz[0]*ExB_x_sq[1]; 
  kappa_inv_bz_sq[2] = (-0.5*ExB_z_sq[1]*bzbz[3])-0.5*ExB_y_sq[1]*bzbz[3]-0.5*ExB_x_sq[1]*bzbz[3]-0.5*bzbz[1]*ExB_z_sq[3]-0.5*bzbz[1]*ExB_y_sq[3]-0.5*bzbz[1]*ExB_x_sq[3]-0.5*ExB_z_sq[0]*bzbz[2]-0.5*ExB_y_sq[0]*bzbz[2]-0.5*ExB_x_sq[0]*bzbz[2]+bzbz[2]-0.5*bzbz[0]*ExB_z_sq[2]-0.5*bzbz[0]*ExB_y_sq[2]-0.5*bzbz[0]*ExB_x_sq[2]; 
  kappa_inv_bz_sq[3] = (-0.5*ExB_z_sq[0]*bzbz[3])-0.5*ExB_y_sq[0]*bzbz[3]-0.5*ExB_x_sq[0]*bzbz[3]+bzbz[3]-0.5*bzbz[0]*ExB_z_sq[3]-0.5*bzbz[0]*ExB_y_sq[3]-0.5*bzbz[0]*ExB_x_sq[3]-0.5*ExB_z_sq[1]*bzbz[2]-0.5*ExB_y_sq[1]*bzbz[2]-0.5*ExB_x_sq[1]*bzbz[2]-0.5*bzbz[1]*ExB_z_sq[2]-0.5*bzbz[1]*ExB_y_sq[2]-0.5*bzbz[1]*ExB_x_sq[2]; 

  // Calculate b_i/kappa = (B_i/|B|)/kappa at quadrature points. 
  // Uses the sign of b_i at quadrature points to get the correct sign of b_i/kappa. 
  ser_2x_p1_sqrt_with_sign(bx, kappa_inv_bx_sq, kappa_inv_bx); 
  ser_2x_p1_sqrt_with_sign(by, kappa_inv_by_sq, kappa_inv_by); 
  ser_2x_p1_sqrt_with_sign(bz, kappa_inv_bz_sq, kappa_inv_bz); 
} 
 
