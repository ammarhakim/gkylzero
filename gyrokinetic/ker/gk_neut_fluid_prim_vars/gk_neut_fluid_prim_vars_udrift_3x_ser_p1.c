#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
#include <gkyl_binop_mul_ser.h> 
#include <gkyl_basis_ser_3x_p1_inv.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_set_prob_3x_ser_p1(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, 
    const double *moms) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A: preallocated LHS matrix. 
  // rhs: preallocated RHS vector. 
  // moms: moments (rho, rho ux, rho uy, rho uz, totalE).

  // For poly_order = 1, we can analytically invert the matrix and just store the solution 
  struct gkyl_mat rhs_ux = gkyl_nmat_get(rhs, count); 
  struct gkyl_mat rhs_uy = gkyl_nmat_get(rhs, count+1); 
  struct gkyl_mat rhs_uz = gkyl_nmat_get(rhs, count+2); 
  // Clear rhs for each component of flow velocity being solved for 
  gkyl_mat_clear(&rhs_ux, 0.0); 
  gkyl_mat_clear(&rhs_uy, 0.0); 
  gkyl_mat_clear(&rhs_uz, 0.0); 
  const double *rho   = &moms[0]; 
  const double *rhoux = &moms[8]; 
  const double *rhouy = &moms[16]; 
  const double *rhouz = &moms[24]; 

  double rho_inv[8] = {0.0}; 
  ser_3x_p1_inv(rho, rho_inv); 
  // Calculate expansions of flow velocity, which can be calculated free of aliasing errors. 
  double ux[8] = {0.0}; 
  double uy[8] = {0.0}; 
  double uz[8] = {0.0}; 
 
  binop_mul_3d_ser_p1(rho_inv, rhoux, ux); 
  binop_mul_3d_ser_p1(rho_inv, rhouy, uy); 
  binop_mul_3d_ser_p1(rho_inv, rhouz, uz); 
 
  gkyl_mat_set(&rhs_ux,0,0,ux[0]); 
  gkyl_mat_set(&rhs_uy,0,0,uy[0]); 
  gkyl_mat_set(&rhs_uz,0,0,uz[0]); 
  gkyl_mat_set(&rhs_ux,1,0,ux[1]); 
  gkyl_mat_set(&rhs_uy,1,0,uy[1]); 
  gkyl_mat_set(&rhs_uz,1,0,uz[1]); 
  gkyl_mat_set(&rhs_ux,2,0,ux[2]); 
  gkyl_mat_set(&rhs_uy,2,0,uy[2]); 
  gkyl_mat_set(&rhs_uz,2,0,uz[2]); 
  gkyl_mat_set(&rhs_ux,3,0,ux[3]); 
  gkyl_mat_set(&rhs_uy,3,0,uy[3]); 
  gkyl_mat_set(&rhs_uz,3,0,uz[3]); 
  gkyl_mat_set(&rhs_ux,4,0,ux[4]); 
  gkyl_mat_set(&rhs_uy,4,0,uy[4]); 
  gkyl_mat_set(&rhs_uz,4,0,uz[4]); 
  gkyl_mat_set(&rhs_ux,5,0,ux[5]); 
  gkyl_mat_set(&rhs_uy,5,0,uy[5]); 
  gkyl_mat_set(&rhs_uz,5,0,uz[5]); 
  gkyl_mat_set(&rhs_ux,6,0,ux[6]); 
  gkyl_mat_set(&rhs_uy,6,0,uy[6]); 
  gkyl_mat_set(&rhs_uz,6,0,uz[6]); 
  gkyl_mat_set(&rhs_ux,7,0,ux[7]); 
  gkyl_mat_set(&rhs_uy,7,0,uy[7]); 
  gkyl_mat_set(&rhs_uz,7,0,uz[7]); 
 
} 
#include <gkyl_mat.h> 
#include <gkyl_gk_neut_fluid_prim_vars_kernels.h> 
GKYL_CU_DH void gk_neut_fluid_prim_vars_udrift_get_sol_3x_ser_p1(int count, struct gkyl_nmat *xsol, 
    double* GKYL_RESTRICT udrift) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // xsol: Input solution vector. 
  // udrift: Output volume expansion of flow velocity: 
  //         [ux, uy, uz]. 
 
  struct gkyl_mat x_ux = gkyl_nmat_get(xsol, count); 
  struct gkyl_mat x_uy = gkyl_nmat_get(xsol, count+1); 
  struct gkyl_mat x_uz = gkyl_nmat_get(xsol, count+2); 
  double *ux = &udrift[0]; 
  double *uy = &udrift[8]; 
  double *uz = &udrift[16]; 

  ux[0] = gkyl_mat_get(&x_ux,0,0); 
  uy[0] = gkyl_mat_get(&x_uy,0,0); 
  uz[0] = gkyl_mat_get(&x_uz,0,0); 
  ux[1] = gkyl_mat_get(&x_ux,1,0); 
  uy[1] = gkyl_mat_get(&x_uy,1,0); 
  uz[1] = gkyl_mat_get(&x_uz,1,0); 
  ux[2] = gkyl_mat_get(&x_ux,2,0); 
  uy[2] = gkyl_mat_get(&x_uy,2,0); 
  uz[2] = gkyl_mat_get(&x_uz,2,0); 
  ux[3] = gkyl_mat_get(&x_ux,3,0); 
  uy[3] = gkyl_mat_get(&x_uy,3,0); 
  uz[3] = gkyl_mat_get(&x_uz,3,0); 
  ux[4] = gkyl_mat_get(&x_ux,4,0); 
  uy[4] = gkyl_mat_get(&x_uy,4,0); 
  uz[4] = gkyl_mat_get(&x_uz,4,0); 
  ux[5] = gkyl_mat_get(&x_ux,5,0); 
  uy[5] = gkyl_mat_get(&x_uy,5,0); 
  uz[5] = gkyl_mat_get(&x_uz,5,0); 
  ux[6] = gkyl_mat_get(&x_ux,6,0); 
  uy[6] = gkyl_mat_get(&x_uy,6,0); 
  uz[6] = gkyl_mat_get(&x_uz,6,0); 
  ux[7] = gkyl_mat_get(&x_ux,7,0); 
  uy[7] = gkyl_mat_get(&x_uy,7,0); 
  uz[7] = gkyl_mat_get(&x_uz,7,0); 

} 
 
