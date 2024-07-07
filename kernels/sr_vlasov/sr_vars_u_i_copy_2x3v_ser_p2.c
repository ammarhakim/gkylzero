#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_u_i_copy_2x3v_ser_p2(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u_i) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // u_i:   Bulk four-velocity (GammaV, GammaV*V_drift). 
 
  struct gkyl_mat x0 = gkyl_nmat_get(x, count+0); 
  double *u_0 = &u_i[0]; 
  u_0[0] = gkyl_mat_get(&x0,0,0); 
  u_0[1] = gkyl_mat_get(&x0,1,0); 
  u_0[2] = gkyl_mat_get(&x0,2,0); 
  u_0[3] = gkyl_mat_get(&x0,3,0); 
  u_0[4] = gkyl_mat_get(&x0,4,0); 
  u_0[5] = gkyl_mat_get(&x0,5,0); 
  u_0[6] = gkyl_mat_get(&x0,6,0); 
  u_0[7] = gkyl_mat_get(&x0,7,0); 
 
  struct gkyl_mat x1 = gkyl_nmat_get(x, count+1); 
  double *u_1 = &u_i[8]; 
  u_1[0] = gkyl_mat_get(&x1,0,0); 
  u_1[1] = gkyl_mat_get(&x1,1,0); 
  u_1[2] = gkyl_mat_get(&x1,2,0); 
  u_1[3] = gkyl_mat_get(&x1,3,0); 
  u_1[4] = gkyl_mat_get(&x1,4,0); 
  u_1[5] = gkyl_mat_get(&x1,5,0); 
  u_1[6] = gkyl_mat_get(&x1,6,0); 
  u_1[7] = gkyl_mat_get(&x1,7,0); 
 
  struct gkyl_mat x2 = gkyl_nmat_get(x, count+2); 
  double *u_2 = &u_i[16]; 
  u_2[0] = gkyl_mat_get(&x2,0,0); 
  u_2[1] = gkyl_mat_get(&x2,1,0); 
  u_2[2] = gkyl_mat_get(&x2,2,0); 
  u_2[3] = gkyl_mat_get(&x2,3,0); 
  u_2[4] = gkyl_mat_get(&x2,4,0); 
  u_2[5] = gkyl_mat_get(&x2,5,0); 
  u_2[6] = gkyl_mat_get(&x2,6,0); 
  u_2[7] = gkyl_mat_get(&x2,7,0); 
 
  struct gkyl_mat x3 = gkyl_nmat_get(x, count+3); 
  double *u_3 = &u_i[24]; 
  u_3[0] = gkyl_mat_get(&x3,0,0); 
  u_3[1] = gkyl_mat_get(&x3,1,0); 
  u_3[2] = gkyl_mat_get(&x3,2,0); 
  u_3[3] = gkyl_mat_get(&x3,3,0); 
  u_3[4] = gkyl_mat_get(&x3,4,0); 
  u_3[5] = gkyl_mat_get(&x3,5,0); 
  u_3[6] = gkyl_mat_get(&x3,6,0); 
  u_3[7] = gkyl_mat_get(&x3,7,0); 
 
} 
 