#include <gkyl_maxwell_kernels.h> 
#include <gkyl_basis_ser_1x_p2_exp_sq.h> 
#include <gkyl_basis_ser_1x_p2_inv.h> 
#include <gkyl_basis_ser_1x_p2_sqrt_with_sign.h> 
GKYL_CU_DH void em_pkpm_kappa_inv_b_1x_ser_p2(const double *bvar, const double *ExB, double* kappa_inv_b) 
{ 
  // bvar:        Input magnetic field unit vector and unit tensor. 
  // ExB:         Input E x B velocity. 
  // kappa_inv_b: b_i/kappa = B_i/|B| * sqrt(1 - |E x B|^2/(c^2 |B|^4)). 
 
  const double *bxbx = &bvar[9]; 
  const double *byby = &bvar[18]; 
  const double *bzbz = &bvar[24]; 
 
  const double *ExB_x = &ExB[0]; 
  const double *ExB_y = &ExB[3]; 
  const double *ExB_z = &ExB[6]; 
 
  double *kappa_inv_bx = &kappa_inv_b[0]; 
  double *kappa_inv_by = &kappa_inv_b[3]; 
  double *kappa_inv_bz = &kappa_inv_b[6]; 
 
  // Calculate ((E x B)/|B|^2)^2. 
  double ExB_x_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(ExB_x, ExB_x_sq); 
 
  double ExB_y_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(ExB_y, ExB_y_sq); 
 
  double ExB_z_sq[3] = {0.0}; 
  ser_1x_p2_exp_sq(ExB_z, ExB_z_sq); 
 
  double kappa_inv_bx_sq[3] = {0.0}; 
  double kappa_inv_by_sq[3] = {0.0}; 
  double kappa_inv_bz_sq[3] = {0.0}; 
 
  kappa_inv_bx_sq[0] = (-0.7071067811865475*ExB_z_sq[2]*bxbx[2])-0.7071067811865475*ExB_y_sq[2]*bxbx[2]-0.7071067811865475*ExB_x_sq[2]*bxbx[2]-0.7071067811865475*ExB_z_sq[1]*bxbx[1]-0.7071067811865475*ExB_y_sq[1]*bxbx[1]-0.7071067811865475*ExB_x_sq[1]*bxbx[1]-0.7071067811865475*ExB_z_sq[0]*bxbx[0]-0.7071067811865475*ExB_y_sq[0]*bxbx[0]-0.7071067811865475*ExB_x_sq[0]*bxbx[0]+bxbx[0]; 
  kappa_inv_bx_sq[1] = (-0.6324555320336759*ExB_z_sq[1]*bxbx[2])-0.6324555320336759*ExB_y_sq[1]*bxbx[2]-0.6324555320336759*ExB_x_sq[1]*bxbx[2]-0.6324555320336759*bxbx[1]*ExB_z_sq[2]-0.6324555320336759*bxbx[1]*ExB_y_sq[2]-0.6324555320336759*bxbx[1]*ExB_x_sq[2]-0.7071067811865475*ExB_z_sq[0]*bxbx[1]-0.7071067811865475*ExB_y_sq[0]*bxbx[1]-0.7071067811865475*ExB_x_sq[0]*bxbx[1]+bxbx[1]-0.7071067811865475*bxbx[0]*ExB_z_sq[1]-0.7071067811865475*bxbx[0]*ExB_y_sq[1]-0.7071067811865475*bxbx[0]*ExB_x_sq[1]; 
  kappa_inv_bx_sq[2] = (-0.4517539514526256*ExB_z_sq[2]*bxbx[2])-0.4517539514526256*ExB_y_sq[2]*bxbx[2]-0.4517539514526256*ExB_x_sq[2]*bxbx[2]-0.7071067811865475*ExB_z_sq[0]*bxbx[2]-0.7071067811865475*ExB_y_sq[0]*bxbx[2]-0.7071067811865475*ExB_x_sq[0]*bxbx[2]+bxbx[2]-0.7071067811865475*bxbx[0]*ExB_z_sq[2]-0.7071067811865475*bxbx[0]*ExB_y_sq[2]-0.7071067811865475*bxbx[0]*ExB_x_sq[2]-0.6324555320336759*ExB_z_sq[1]*bxbx[1]-0.6324555320336759*ExB_y_sq[1]*bxbx[1]-0.6324555320336759*ExB_x_sq[1]*bxbx[1]; 

  kappa_inv_by_sq[0] = (-0.7071067811865475*ExB_z_sq[2]*byby[2])-0.7071067811865475*ExB_y_sq[2]*byby[2]-0.7071067811865475*ExB_x_sq[2]*byby[2]-0.7071067811865475*ExB_z_sq[1]*byby[1]-0.7071067811865475*ExB_y_sq[1]*byby[1]-0.7071067811865475*ExB_x_sq[1]*byby[1]-0.7071067811865475*ExB_z_sq[0]*byby[0]-0.7071067811865475*ExB_y_sq[0]*byby[0]-0.7071067811865475*ExB_x_sq[0]*byby[0]+byby[0]; 
  kappa_inv_by_sq[1] = (-0.6324555320336759*ExB_z_sq[1]*byby[2])-0.6324555320336759*ExB_y_sq[1]*byby[2]-0.6324555320336759*ExB_x_sq[1]*byby[2]-0.6324555320336759*byby[1]*ExB_z_sq[2]-0.6324555320336759*byby[1]*ExB_y_sq[2]-0.6324555320336759*byby[1]*ExB_x_sq[2]-0.7071067811865475*ExB_z_sq[0]*byby[1]-0.7071067811865475*ExB_y_sq[0]*byby[1]-0.7071067811865475*ExB_x_sq[0]*byby[1]+byby[1]-0.7071067811865475*byby[0]*ExB_z_sq[1]-0.7071067811865475*byby[0]*ExB_y_sq[1]-0.7071067811865475*byby[0]*ExB_x_sq[1]; 
  kappa_inv_by_sq[2] = (-0.4517539514526256*ExB_z_sq[2]*byby[2])-0.4517539514526256*ExB_y_sq[2]*byby[2]-0.4517539514526256*ExB_x_sq[2]*byby[2]-0.7071067811865475*ExB_z_sq[0]*byby[2]-0.7071067811865475*ExB_y_sq[0]*byby[2]-0.7071067811865475*ExB_x_sq[0]*byby[2]+byby[2]-0.7071067811865475*byby[0]*ExB_z_sq[2]-0.7071067811865475*byby[0]*ExB_y_sq[2]-0.7071067811865475*byby[0]*ExB_x_sq[2]-0.6324555320336759*ExB_z_sq[1]*byby[1]-0.6324555320336759*ExB_y_sq[1]*byby[1]-0.6324555320336759*ExB_x_sq[1]*byby[1]; 

  kappa_inv_bz_sq[0] = (-0.7071067811865475*ExB_z_sq[2]*bzbz[2])-0.7071067811865475*ExB_y_sq[2]*bzbz[2]-0.7071067811865475*ExB_x_sq[2]*bzbz[2]-0.7071067811865475*ExB_z_sq[1]*bzbz[1]-0.7071067811865475*ExB_y_sq[1]*bzbz[1]-0.7071067811865475*ExB_x_sq[1]*bzbz[1]-0.7071067811865475*ExB_z_sq[0]*bzbz[0]-0.7071067811865475*ExB_y_sq[0]*bzbz[0]-0.7071067811865475*ExB_x_sq[0]*bzbz[0]+bzbz[0]; 
  kappa_inv_bz_sq[1] = (-0.6324555320336759*ExB_z_sq[1]*bzbz[2])-0.6324555320336759*ExB_y_sq[1]*bzbz[2]-0.6324555320336759*ExB_x_sq[1]*bzbz[2]-0.6324555320336759*bzbz[1]*ExB_z_sq[2]-0.6324555320336759*bzbz[1]*ExB_y_sq[2]-0.6324555320336759*bzbz[1]*ExB_x_sq[2]-0.7071067811865475*ExB_z_sq[0]*bzbz[1]-0.7071067811865475*ExB_y_sq[0]*bzbz[1]-0.7071067811865475*ExB_x_sq[0]*bzbz[1]+bzbz[1]-0.7071067811865475*bzbz[0]*ExB_z_sq[1]-0.7071067811865475*bzbz[0]*ExB_y_sq[1]-0.7071067811865475*bzbz[0]*ExB_x_sq[1]; 
  kappa_inv_bz_sq[2] = (-0.4517539514526256*ExB_z_sq[2]*bzbz[2])-0.4517539514526256*ExB_y_sq[2]*bzbz[2]-0.4517539514526256*ExB_x_sq[2]*bzbz[2]-0.7071067811865475*ExB_z_sq[0]*bzbz[2]-0.7071067811865475*ExB_y_sq[0]*bzbz[2]-0.7071067811865475*ExB_x_sq[0]*bzbz[2]+bzbz[2]-0.7071067811865475*bzbz[0]*ExB_z_sq[2]-0.7071067811865475*bzbz[0]*ExB_y_sq[2]-0.7071067811865475*bzbz[0]*ExB_x_sq[2]-0.6324555320336759*ExB_z_sq[1]*bzbz[1]-0.6324555320336759*ExB_y_sq[1]*bzbz[1]-0.6324555320336759*ExB_x_sq[1]*bzbz[1]; 

  // Calculate b_i/kappa = (B_i/|B|)/kappa at quadrature points. 
  // Uses the sign of b_i at quadrature points to get the correct sign of b_i/kappa. 
  ser_1x_p2_sqrt_with_sign(bx, kappa_inv_bx_sq, kappa_inv_bx); 
  ser_1x_p2_sqrt_with_sign(by, kappa_inv_by_sq, kappa_inv_by); 
  ser_1x_p2_sqrt_with_sign(bz, kappa_inv_bz_sq, kappa_inv_bz); 
} 
 
