#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p1_exp_sq.h> 
#include <gkyl_basis_ser_1x_p1_inv.h> 
#include <gkyl_basis_ser_1x_p1_sqrt.h> 
GKYL_CU_DH void em_pkpm_kappa_inv_b_1x_ser_p1(const double *bvar, const double *ExB, double* GKYL_RESTRICT kappa_inv_b) 
{ 
  // bvar:        Input magnetic field unit vector and unit tensor. 
  // ExB:         Input E x B velocity. 
  // kappa_inv_b: b_i/kappa = B_i/|B| * sqrt(1 - |E x B|^2/(c^2 |B|^4)). 
 
  const double *bxbx = &bvar[6]; 
  const double *byby = &bvar[12]; 
  const double *bzbz = &bvar[16]; 
 
  const double *ExB_x = &ExB[0]; 
  const double *ExB_y = &ExB[2]; 
  const double *ExB_z = &ExB[4]; 
 
  double *kappa_inv_bx = &kappa_inv_b[0]; 
  double *kappa_inv_by = &kappa_inv_b[2]; 
  double *kappa_inv_bz = &kappa_inv_b[4]; 
 
  // Calculate ((E x B)/|B|^2)^2. 
  double ExB_x_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(ExB_x, ExB_x_sq); 
 
  double ExB_y_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(ExB_y, ExB_y_sq); 
 
  double ExB_z_sq[2] = {0.0}; 
  ser_1x_p1_exp_sq(ExB_z, ExB_z_sq); 
 
  double kappa_inv_bx_sq[2] = {0.0}; 
  double kappa_inv_by_sq[2] = {0.0}; 
  double kappa_inv_bz_sq[2] = {0.0}; 
 
  kappa_inv_bx_sq[0] = (-0.7071067811865475*ExB_z_sq[1]*bxbx[1])-0.7071067811865475*ExB_y_sq[1]*bxbx[1]-0.7071067811865475*ExB_x_sq[1]*bxbx[1]-0.7071067811865475*ExB_z_sq[0]*bxbx[0]-0.7071067811865475*ExB_y_sq[0]*bxbx[0]-0.7071067811865475*ExB_x_sq[0]*bxbx[0]+bxbx[0]; 
  kappa_inv_bx_sq[1] = (-0.7071067811865475*ExB_z_sq[0]*bxbx[1])-0.7071067811865475*ExB_y_sq[0]*bxbx[1]-0.7071067811865475*ExB_x_sq[0]*bxbx[1]+bxbx[1]-0.7071067811865475*bxbx[0]*ExB_z_sq[1]-0.7071067811865475*bxbx[0]*ExB_y_sq[1]-0.7071067811865475*bxbx[0]*ExB_x_sq[1]; 

  kappa_inv_by_sq[0] = (-0.7071067811865475*ExB_z_sq[1]*byby[1])-0.7071067811865475*ExB_y_sq[1]*byby[1]-0.7071067811865475*ExB_x_sq[1]*byby[1]-0.7071067811865475*ExB_z_sq[0]*byby[0]-0.7071067811865475*ExB_y_sq[0]*byby[0]-0.7071067811865475*ExB_x_sq[0]*byby[0]+byby[0]; 
  kappa_inv_by_sq[1] = (-0.7071067811865475*ExB_z_sq[0]*byby[1])-0.7071067811865475*ExB_y_sq[0]*byby[1]-0.7071067811865475*ExB_x_sq[0]*byby[1]+byby[1]-0.7071067811865475*byby[0]*ExB_z_sq[1]-0.7071067811865475*byby[0]*ExB_y_sq[1]-0.7071067811865475*byby[0]*ExB_x_sq[1]; 

  kappa_inv_bz_sq[0] = (-0.7071067811865475*ExB_z_sq[1]*bzbz[1])-0.7071067811865475*ExB_y_sq[1]*bzbz[1]-0.7071067811865475*ExB_x_sq[1]*bzbz[1]-0.7071067811865475*ExB_z_sq[0]*bzbz[0]-0.7071067811865475*ExB_y_sq[0]*bzbz[0]-0.7071067811865475*ExB_x_sq[0]*bzbz[0]+bzbz[0]; 
  kappa_inv_bz_sq[1] = (-0.7071067811865475*ExB_z_sq[0]*bzbz[1])-0.7071067811865475*ExB_y_sq[0]*bzbz[1]-0.7071067811865475*ExB_x_sq[0]*bzbz[1]+bzbz[1]-0.7071067811865475*bzbz[0]*ExB_z_sq[1]-0.7071067811865475*bzbz[0]*ExB_y_sq[1]-0.7071067811865475*bzbz[0]*ExB_x_sq[1]; 

  ser_1x_p1_sqrt(kappa_inv_bx_sq, kappa_inv_bx); 
  ser_1x_p1_sqrt(kappa_inv_by_sq, kappa_inv_by); 
  ser_1x_p1_sqrt(kappa_inv_bz_sq, kappa_inv_bz); 
} 
 
