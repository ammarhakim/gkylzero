#include <gkyl_mat.h> 
#include <gkyl_sr_Gamma_kernels.h> 
GKYL_CU_DH void sr_vars_u_i_copy_1x2v_ser_p1(int count, struct gkyl_nmat *x, double* GKYL_RESTRICT u_i) 
{ 
  // count: integer to indicate which matrix being fetched. 
  // x:     Input solution vector. 
  // u_i:   Bulk four-velocity (GammaV, GammaV*V_drift). 
 
  struct gkyl_mat x0 = gkyl_nmat_get(x, count+0); 
  double *u_0 = &u_i[0]; 
  u_0[0] = gkyl_mat_get(&x0,0,0); 
  u_0[1] = gkyl_mat_get(&x0,1,0); 
 
  struct gkyl_mat x1 = gkyl_nmat_get(x, count+1); 
  double *u_1 = &u_i[2]; 
  u_1[0] = gkyl_mat_get(&x1,0,0); 
  u_1[1] = gkyl_mat_get(&x1,1,0); 
 
  struct gkyl_mat x2 = gkyl_nmat_get(x, count+2); 
  double *u_2 = &u_i[4]; 
  u_2[0] = gkyl_mat_get(&x2,0,0); 
  u_2[1] = gkyl_mat_get(&x2,1,0); 
 
} 
 
