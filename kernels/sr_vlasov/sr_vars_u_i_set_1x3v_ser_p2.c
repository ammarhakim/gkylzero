#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_u_i_set_1x3v_ser_p2(int count, struct gkyl_nmat *A, struct gkyl_nmat *rhs, const double *M0, const double *M1i, const double *n) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // A:     preallocated LHS matrix. 
  // rhs:   preallocated RHS vector. 
  // M0:    Lab frame density = GammaV*n.
  // M1i:   Lab frame flux = GammaV*n*V_drift.
  // n:     Rest-frame density.

  struct gkyl_mat A0 = gkyl_nmat_get(A, count); 
  struct gkyl_mat rhs0 = gkyl_nmat_get(rhs, count); 
  gkyl_mat_clear(&A0, 0.0); gkyl_mat_clear(&rhs0, 0.0); 
  gkyl_mat_set(&rhs0,0,0,M0[0]); 
  gkyl_mat_set(&rhs0,1,0,M0[1]); 
  gkyl_mat_set(&rhs0,2,0,M0[2]); 
  gkyl_mat_set(&A0,0,0,0.7071067811865475*n[0]); 
  gkyl_mat_set(&A0,0,1,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A0,0,2,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A0,1,0,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A0,1,1,0.6324555320336759*n[2]+0.7071067811865475*n[0]); 
  gkyl_mat_set(&A0,1,2,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A0,2,0,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A0,2,1,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A0,2,2,0.4517539514526256*n[2]+0.7071067811865475*n[0]); 
 
  struct gkyl_mat A1 = gkyl_nmat_get(A, count+1); 
  struct gkyl_mat rhs1 = gkyl_nmat_get(rhs, count+1); 
  gkyl_mat_clear(&A1, 0.0); gkyl_mat_clear(&rhs1, 0.0); 
  const double *M10 = &M1i[0]; 
  gkyl_mat_set(&rhs1,0,0,M10[0]); 
  gkyl_mat_set(&rhs1,1,0,M10[1]); 
  gkyl_mat_set(&rhs1,2,0,M10[2]); 
  gkyl_mat_set(&A1,0,0,0.7071067811865475*n[0]); 
  gkyl_mat_set(&A1,0,1,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A1,0,2,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A1,1,0,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A1,1,1,0.6324555320336759*n[2]+0.7071067811865475*n[0]); 
  gkyl_mat_set(&A1,1,2,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A1,2,0,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A1,2,1,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A1,2,2,0.4517539514526256*n[2]+0.7071067811865475*n[0]); 
 
  struct gkyl_mat A2 = gkyl_nmat_get(A, count+2); 
  struct gkyl_mat rhs2 = gkyl_nmat_get(rhs, count+2); 
  gkyl_mat_clear(&A2, 0.0); gkyl_mat_clear(&rhs2, 0.0); 
  const double *M11 = &M1i[3]; 
  gkyl_mat_set(&rhs2,0,0,M11[0]); 
  gkyl_mat_set(&rhs2,1,0,M11[1]); 
  gkyl_mat_set(&rhs2,2,0,M11[2]); 
  gkyl_mat_set(&A2,0,0,0.7071067811865475*n[0]); 
  gkyl_mat_set(&A2,0,1,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A2,0,2,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A2,1,0,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A2,1,1,0.6324555320336759*n[2]+0.7071067811865475*n[0]); 
  gkyl_mat_set(&A2,1,2,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A2,2,0,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A2,2,1,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A2,2,2,0.4517539514526256*n[2]+0.7071067811865475*n[0]); 
 
  struct gkyl_mat A3 = gkyl_nmat_get(A, count+3); 
  struct gkyl_mat rhs3 = gkyl_nmat_get(rhs, count+3); 
  gkyl_mat_clear(&A3, 0.0); gkyl_mat_clear(&rhs3, 0.0); 
  const double *M12 = &M1i[6]; 
  gkyl_mat_set(&rhs3,0,0,M12[0]); 
  gkyl_mat_set(&rhs3,1,0,M12[1]); 
  gkyl_mat_set(&rhs3,2,0,M12[2]); 
  gkyl_mat_set(&A3,0,0,0.7071067811865475*n[0]); 
  gkyl_mat_set(&A3,0,1,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A3,0,2,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A3,1,0,0.7071067811865475*n[1]); 
  gkyl_mat_set(&A3,1,1,0.6324555320336759*n[2]+0.7071067811865475*n[0]); 
  gkyl_mat_set(&A3,1,2,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A3,2,0,0.7071067811865475*n[2]); 
  gkyl_mat_set(&A3,2,1,0.6324555320336759*n[1]); 
  gkyl_mat_set(&A3,2,2,0.4517539514526256*n[2]+0.7071067811865475*n[0]); 
 
} 
